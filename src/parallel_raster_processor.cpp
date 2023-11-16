

#include <stdexcept>

#include "grid.h"
#include "operation.h"
#include "parallel_raster_processor.h"
#include "processor.h"
#include "stats_registry.h"
#include "operation.h"
#include "raster_sequential_processor.h"
#include "raster_source.h"

#include <map>
#include <memory>
#include <set>

#include <oneapi/tbb.h>

namespace {

void errorHandlerParallel(const char *fmt, ...) {

    char buf[BUFSIZ], *p;
    va_list ap;
    va_start(ap, fmt);
    vsprintf(buf, fmt, ap);
    va_end(ap);
    p = buf + strlen(buf) - 1;
    if(strlen(buf) > 0 && *p == '\n') *p = '\0';

    std::cout << "Error: " << buf << std::endl;
}

}

namespace exactextract {

    struct RasterBlock {
        Operation* _op;
        Grid<bounded_extent> _grid;
        std::unique_ptr<AbstractRaster<double>> _raster;

        RasterBlock() : _grid(Grid<bounded_extent>::make_empty()) {}

        RasterBlock(Operation* o, Grid<bounded_extent> g, std::unique_ptr<AbstractRaster<double>> r)
        : _op(o), _grid(g), _raster(std::move(r)) {}
    };

    typedef std::vector<std::pair<std::string, geom_ptr_r>> GeometryBatch;
    typedef std::shared_ptr<GeometryBatch> GeometryBatchPtr;
    typedef std::pair<GeometryBatchPtr, std::shared_ptr<StatsRegistry>> BatchAndResults;
    const int BatchSize = 1000;

    void ParallelRasterProcessor::read_features() {
        while (m_shp.next()) {
            Feature feature = std::make_pair(
                    m_shp.feature_field(m_shp.id_field()),
                    geos_ptr(m_geos_context, m_shp.feature_geometry(m_geos_context)));
            m_features.push_back(std::move(feature));
        }
    }

    void ParallelRasterProcessor::populate_index() {
        for (const Feature& f : m_features) {
            // TODO compute envelope of dataset, and crop raster by that extent before processing?
            GEOSSTRtree_insert_r(m_geos_context, m_feature_tree.get(), f.second.get(), (void *) &f);
        }
    }

    void ParallelRasterProcessor::process() {

        for (const auto& op : m_operations) {
            m_output.add_operation(*op);
        }
        bool store_values = StatsRegistry::requires_stored_values(m_operations);

        read_features();
        populate_index();
        auto grid = common_grid(m_operations.begin(), m_operations.end());
        auto subdividedGrid = subdivide(grid, m_max_cells_in_memory);

        oneapi::tbb::enumerable_thread_specific<GEOSContextHandle_t> geos_context([] () -> GEOSContextHandle_t {
            return initGEOS_r(errorHandlerParallel, errorHandlerParallel);
        });

        oneapi::tbb::parallel_pipeline(6,
            oneapi::tbb::make_filter<void, RasterBlock>(oneapi::tbb::filter_mode::serial_in_order,
            //read features and setup batch for processing
            [&grid, &subdividedGrid, &geos_context, this] (oneapi::tbb::flow_control& fc) -> RasterBlock {
                //m_shp should only be accessed here
                //iterate through m_shp
                //return name + geom

                if (subdividedGrid.begin() == subdividedGrid.end()) {
                    fc.stop();
                    return RasterBlock();
                }

                auto op = m_operations.begin()->get();
                auto subgrid = *subdividedGrid.begin();

                RasterBlock block{
                    op,
                    subgrid,
                    op->values->read_box(subgrid.extent().intersection(op->values->grid().extent()))
                };

                subdividedGrid.erase(subdividedGrid.begin());
                return block;
            }) &
            //process batch
            oneapi::tbb::make_filter<RasterBlock, StatsRegistry>(oneapi::tbb::filter_mode::parallel,
            [&geos_context, &featureTree = m_feature_tree] (const RasterBlock& block) -> StatsRegistry {
                //intersect with raster
                //return m_reg + name
                auto& threadGeosContext = geos_context.local();
                // bool store_values = StatsRegistry::requires_stored_values(ops);
                bool store_values = true;
                StatsRegistry reg;

                auto query_rect = geos_make_box_polygon(threadGeosContext, block._grid.extent());
                std::vector<const Feature *> hits;
                auto values = block._raster.get();

                GEOSSTRtree_query_r(threadGeosContext,
                    featureTree.get(),
                    query_rect.get(),
                    [](void *hit, void *userdata) {
                        auto feature = static_cast<const Feature *>(hit);
                        auto vec = static_cast<std::vector<const Feature *> *>(userdata);

                        vec->push_back(feature);
                    }, &hits);
                
                for (const auto& f : hits) {
                    std::unique_ptr<Raster<float>> coverage = std::make_unique<Raster<float>>(
                        raster_cell_intersection(block._grid, threadGeosContext, f->second.get()));
                    
                    
                    reg.stats(f->first, *block._op, store_values).process(*coverage, *values);
                }

                return reg;
            }) &
            //write out batch results
            oneapi::tbb::make_filter<StatsRegistry, void>(oneapi::tbb::filter_mode::serial_out_of_order,
            [&output = m_output] (const StatsRegistry& reg) {
                //output
                //m_output should only be accessed here
                output.set_registry(&reg);

                for (auto& feature : reg.m_feature_stats) {
                    output.write(feature.first);
                }

                std::cout << "Wrote " << reg.m_feature_stats.size() << " features." << std::endl;
            }));
    }

}