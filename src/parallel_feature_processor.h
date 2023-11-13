#ifndef EXACTEXTRACT_PARALLEL_PROCESSOR_H
#define EXACTEXTRACT_PARALLEL_PROCESSOR_H

#include "processor.h"

namespace exactextract {

    class ParallelFeatureProcessor : public Processor {
    public:
        using Processor::Processor;

        void process() override;
    };
}

#endif