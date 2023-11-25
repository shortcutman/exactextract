

#include <stdexcept>

#include "feature_source.h"
#include "grid.h"
#include "operation.h"
#include "parallel_raster_processor.h"
#include "processor.h"
#include "stats_registry.h"
#include "operation.h"
#include "raster_source.h"

#include <map>
#include <memory>
#include <set>

#include <oneapi/tbb.h>
#include <unordered_map>

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
        Grid<bounded_extent> _grid;
        std::set<std::pair<RasterSource*, RasterSource*>> processed;
        std::map<RasterSource*, std::unique_ptr<AbstractRaster<double>>> raster_values;
        StatsRegistry _reg;
        std::vector<const Feature *> _hits;

        RasterBlock() : _grid(Grid<bounded_extent>::make_empty()) {}

        RasterBlock(Grid<bounded_extent> g)
        : _grid(g) {}
    };

    void ParallelRasterProcessor::read_features() {
        while (m_shp.next()) {
            const Feature& feature = m_shp.feature();
            MapFeature mf(feature);
            m_features.push_back(std::move(mf));
        }
    }

    void ParallelRasterProcessor::populate_index() {
        for (const auto& f : m_features) {
            // TODO compute envelope of dataset, and crop raster by that extent before processing?
            GEOSSTRtree_insert_r(m_geos_context, m_feature_tree.get(), f.geometry(), (void *) &f);
        }
    }

    void ParallelRasterProcessor::process() {
        read_features();
        populate_index();

        for (const auto& op : m_operations) {
            m_output.add_operation(*op);
        }

        bool store_values = StatsRegistry::requires_stored_values(m_operations);
        auto operationsGridExtent = common_grid(m_operations.begin(), m_operations.end());
        auto subdividedGrid = subdivide(operationsGridExtent, 10000);

        oneapi::tbb::enumerable_thread_specific<GEOSContextHandle_t> geos_context([] () -> GEOSContextHandle_t {
            return initGEOS_r(errorHandlerParallel, errorHandlerParallel);
        });

        oneapi::tbb::parallel_pipeline(6,
            oneapi::tbb::make_filter<void, std::shared_ptr<RasterBlock>>(oneapi::tbb::filter_mode::serial_in_order,
            //read features and setup batch for processing
            [&operationsGridExtent, &subdividedGrid, &geos_context, this] (oneapi::tbb::flow_control& fc) -> std::shared_ptr<RasterBlock> {
                //iterate through subdivisions and do file IO
                if (subdividedGrid.begin() == subdividedGrid.end()) {
                    fc.stop();
                    return std::shared_ptr<RasterBlock>();
                }

                auto subgrid = *subdividedGrid.begin();

                auto rasterBlock = std::make_shared<RasterBlock>(subgrid);

                std::set<std::pair<RasterSource*, RasterSource*>> processed;
                for (const auto &op : m_operations) {
                    auto key = std::make_pair(op->weights, op->values);
                    if (processed.find(key) != processed.end()) {
                        continue;
                    } else {
                        processed.insert(key);
                    }

                    if (!op->values->grid().extent().contains(subgrid.extent())) {
                        continue;
                    }

                    if (op->weighted() && !op->weights->grid().extent().contains(subgrid.extent())) {
                        continue;
                    }

                    auto values = rasterBlock->raster_values[op->values].get();
                    if (values == nullptr) {
                        rasterBlock->raster_values[op->values] = op->values->read_box(subgrid.extent().intersection(op->values->grid().extent()));
                        values = rasterBlock->raster_values[op->values].get();
                    }
                }                

                subdividedGrid.erase(subdividedGrid.begin());
                return rasterBlock;
            }) &
            //process batch
            oneapi::tbb::make_filter<std::shared_ptr<RasterBlock>, std::shared_ptr<RasterBlock>>(oneapi::tbb::filter_mode::parallel,
            [&geos_context, &featureTree = m_feature_tree, this] (std::shared_ptr<RasterBlock> block) -> std::shared_ptr<RasterBlock> {
                //intersect with raster
                //return m_reg + name
                auto& threadGeosContext = geos_context.local();
                // bool store_values = StatsRegistry::requires_stored_values(ops);
                bool store_values = true;

                auto query_rect = geos_make_box_polygon(threadGeosContext, block->_grid.extent());

                GEOSSTRtree_query_r(threadGeosContext,
                    featureTree.get(),
                    query_rect.get(),
                    [](void *hit, void *userdata) {
                        auto feature = static_cast<const Feature *>(hit);
                        auto vec = static_cast<std::vector<const Feature *> *>(userdata);

                        vec->push_back(feature);
                    }, &block->_hits);
                
                for (const auto& f : block->_hits) {
                    std::string fid = f->get_string(m_shp.id_field());
                    std::unique_ptr<Raster<float>> coverage = std::make_unique<Raster<float>>(raster_cell_intersection(block->_grid, threadGeosContext, f->geometry()));
                    
                    for (const auto &op : m_operations) {
                        auto values = block->raster_values[op->values].get();

                        if (op->weighted()) {
                            auto weights = block->raster_values[op->weights].get();
                            block->_reg.stats(fid, *op, store_values).process(*coverage, *values, *weights);
                        } else {
                            block->_reg.stats(fid, *op, store_values).process(*coverage, *values);
                        }
                    }
                }

                return block;
            }) &
            //write out batch results
            oneapi::tbb::make_filter<std::shared_ptr<RasterBlock>, void>(oneapi::tbb::filter_mode::serial_out_of_order,
            [this] (std::shared_ptr<RasterBlock> block) {

                for (const auto& f_in : block->_hits) {
                    auto f_out = m_output.create_feature();
                    std::string fid = f_in->get_string(m_shp.id_field());
                    f_out->set(m_shp.id_field(), fid);
                    for (const auto& col: m_include_cols) {
                        f_out->set(col, *f_in);
                    }
                    for (const auto& op : m_operations) {
                        op->set_result(m_reg, fid, *f_out);
                    }
                    m_output.write(*f_out);
                    m_reg.flush_feature(fid);
                }

                std::cout << "Wrote " << block->_hits.size() << " features." << std::endl;
            }));
    }

}