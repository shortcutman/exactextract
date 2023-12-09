// Copyright (c) 2019-2023 ISciences, LLC.
// All rights reserved.
//
// This software is licensed under the Apache License, Version 2.0 (the "License").
// You may not use this file except in compliance with the License. You may
// obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0.
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "stats_registry.h"

#include "operation.h"
#include "raster_stats.h"

namespace exactextract {

RasterStats<double>&
StatsRegistry::stats(const std::string& feature, const Operation& op, bool store_values)
{

    // TODO come up with a better storage method.
    auto& stats_for_feature = m_feature_stats[feature];

    // can't use find because this requires RasterStats to be copy-constructible before C++ 17
    auto exists = stats_for_feature.count(op.key());
    if (!exists) {
        // can't use emplace because this requires RasterStats be copy-constructible before C++17
        RasterStats<double> new_stats(store_values);
        stats_for_feature[op.key()] = std::move(new_stats);
    }

    return stats_for_feature[op.key()];
}

void StatsRegistry::join(StatsRegistry& source) {
    for (auto& feature : source.m_feature_stats) {
        auto statsFeatIt = this->m_feature_stats.find(feature.first);
        if (statsFeatIt == this->m_feature_stats.end()) {
            //move
            auto& thisFeature = this->m_feature_stats[feature.first];
            for (auto& op : feature.second) {
                thisFeature[op.first] = std::move(op.second);
            }
        } else {
            //merge
            auto& thisFeature = this->m_feature_stats[feature.first];
            for (auto& op : feature.second) {
                thisFeature[op.first].merge(op.second);
            }
        }
    }
}

bool
StatsRegistry::contains(const std::string& feature, const Operation& op) const
{
    const auto& m = m_feature_stats;

    auto it = m.find(feature);

    if (it == m.end()) {
        return false;
    }

    const auto& m2 = it->second;

    return m2.find(op.key()) != m2.end();
}

const RasterStats<double>&
StatsRegistry::stats(const std::string& feature, const Operation& op) const
{
    // TODO come up with a better storage method.
    auto stats = m_feature_stats.find(feature);
    if (stats != m_feature_stats.end()) {
        return stats->second.at(op.key());
    } else {
        return RasterStats<double>();
    }
}

}
