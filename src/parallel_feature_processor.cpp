

#include <stdexcept>

#include "operation.h"
#include "parallel_feature_processor.h"
#include "processor.h"
#include "stats_registry.h"

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

    std::cerr << buf << std::endl;
}

}

namespace exactextract {

    typedef std::vector<std::pair<std::string, geom_ptr_r>> GeometryBatch;
    typedef std::shared_ptr<GeometryBatch> GeometryBatchPtr;
    typedef std::pair<GeometryBatchPtr, std::shared_ptr<StatsRegistry>> BatchAndResults;
    const int BatchSize = 1000;

    void ParallelFeatureProcessor::process() {
        for (const auto& op : m_operations) {
            m_output.add_operation(*op);
        }
        bool store_values = StatsRegistry::requires_stored_values(m_operations);

        oneapi::tbb::enumerable_thread_specific<GEOSContextHandle_t> geos_context([] () -> GEOSContextHandle_t {
            return initGEOS_r(errorHandlerParallel, errorHandlerParallel);
        });

        oneapi::tbb::parallel_pipeline(4,
            oneapi::tbb::make_filter<void, GeometryBatchPtr>(oneapi::tbb::filter_mode::serial_in_order,
            //read features and setup batch for processing
            [&shp = m_shp, &geos_context] (oneapi::tbb::flow_control& fc) -> GeometryBatchPtr {
                //m_shp should only be accessed here
                //iterate through m_shp
                //return name + geom

                auto batch = std::make_shared<GeometryBatch>();
                auto& threadGeosContext = geos_context.local();

                // while (shp.next()) {
                //     std::string name = shp.feature_field(shp.id_field());
                //     auto geom = geos_ptr(threadGeosContext, shp.feature_geometry(threadGeosContext));
                //     batch->push_back(std::move(std::make_pair(name, std::move(geom))));

                //     if (batch->size() == BatchSize) {
                //         return batch;
                //     }
                // }

                fc.stop();
                return batch;
            }) &
            //process batch
            oneapi::tbb::make_filter<GeometryBatchPtr, BatchAndResults>(oneapi::tbb::filter_mode::parallel,
            [&geos_context, &ops = std::as_const(m_operations), m_max_cells_in_memory = m_max_cells_in_memory] (const GeometryBatchPtr& batch) -> BatchAndResults {
                //intersect with raster
                //return m_reg + name
                auto& threadGeosContext = geos_context.local();
                bool store_values = StatsRegistry::requires_stored_values(ops);
                auto reg = std::make_shared<StatsRegistry>();

                for (auto& feature : *batch) {
                    Box feature_bbox = exactextract::geos_get_box(threadGeosContext, feature.second.get());
                    auto grid = common_grid(ops.begin(), ops.end());

                    if (feature_bbox.intersects(grid.extent())) {
                        // Crop grid to portion overlapping feature
                        auto cropped_grid = grid.crop(feature_bbox);

                        for (const auto &subgrid : subdivide(cropped_grid, m_max_cells_in_memory)) {
                            std::unique_ptr<Raster<float>> coverage;

                            std::set<std::pair<RasterSource*, RasterSource*>> processed;

                            for (const auto &op : ops) {
                                // TODO avoid reading same values/weights multiple times. Just use a map?

                                // Avoid processing same values/weights for different stats
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

                                // Lazy-initialize coverage
                                if (coverage == nullptr) {
                                    coverage = std::make_unique<Raster<float>>(
                                            raster_cell_intersection(subgrid, threadGeosContext, feature.second.get()));
                                }

                                auto values = op->values->read_box(subgrid.extent().intersection(op->values->grid().extent()));

                                if (op->weighted()) {
                                    throw std::runtime_error("not implemented");
                                    auto weights = op->weights->read_box(subgrid.extent().intersection(op->weights->grid().extent()));

                                    reg->stats(feature.first, *op, store_values).process(*coverage, *values, *weights);
                                } else {
                                    reg->stats(feature.first, *op, store_values).process(*coverage, *values);
                                }
                            }
                        }
                    }
                }

                return std::make_pair(batch, reg);
            }) &
            //write out batch results
            oneapi::tbb::make_filter<BatchAndResults, void>(oneapi::tbb::filter_mode::serial_out_of_order,
            [&output = m_output] (const BatchAndResults& results) {
                //output
                //m_output should only be accessed here
                // output.set_registry(results.second.get());

                // for (auto& feature : *(results.first)) {
                //     output.write(feature.first);
                // }

                std::cout << "Wrote " << results.first->size() << " features." << std::endl;
            }));

        throw std::runtime_error("Not implemented.");
    }

}