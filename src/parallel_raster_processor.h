#ifndef EXACTEXTRACT_PARALLEL_RASTER_PROCESSOR_H
#define EXACTEXTRACT_PARALLEL_RASTER_PROCESSOR_H

#include "geos_utils.h"
#include "map_feature.h"
#include "processor.h"

namespace exactextract {

    class ParallelRasterProcessor : public Processor {
    private:
        std::vector<MapFeature> m_features;
        tree_ptr_r m_feature_tree{geos_ptr(m_geos_context, GEOSSTRtree_create_r(m_geos_context, 10))};

    public:
        using Processor::Processor;

        void read_features();
        void populate_index();

        void process() override;
    };
}

#endif