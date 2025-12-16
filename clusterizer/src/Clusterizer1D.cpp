//
// Created by dylan on 12/11/25.
//

#include "Clusterizer1D.h"
#include <cmath>

Clusterizer1D::Clusterizer1D(const StripMapping* mapping)
    : mapping_(mapping)
{}

bool Clusterizer1D::isCompatible(const Hit& h, const Cluster& c) const
{
    const StripInfo* s1 = mapping_->getStripInfo(h.channel);
    if (!s1) return false;

    for (const Hit& existing : c.hits) {
        const StripInfo* s2 = mapping_->getStripInfo(existing.channel);
        if (!s2) continue;

        float dx = s1->xGerber - s2->xGerber;
        float dy = s1->yGerber - s2->yGerber;

        float dist_mm = std::sqrt(dx*dx + dy*dy);
        float dt_ns = std::fabs(h.time_ns - existing.time_ns);

        if (dist_mm <= spatialThreshold_ && dt_ns <= timeThreshold_)
            return true;
    }

    return false;
}

std::vector<Cluster> Clusterizer1D::clusterEvent(const std::vector<Hit>& hits)
{
    std::vector<Cluster> clusters;

    // Process each hit
    for (const Hit& h : hits) {
        const StripInfo* si = mapping_->getStripInfo(h.channel);
        if (!si) continue;

        bool added = false;

        // Try existing clusters
        for (Cluster& c : clusters) {
            if (c.axis != si->axis) continue; // only cluster 1D within an axis
            if (isCompatible(h, c)) {
                c.hits.push_back(h);
                added = true;
                break;
            }
        }

        // If not added, create new cluster
        if (!added) {
            Cluster c;
            c.axis = si->axis;
            c.hits.push_back(h);
            clusters.push_back(std::move(c));
        }
    }

    return clusters;
}
