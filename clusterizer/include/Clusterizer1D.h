//
// Created by dylan on 12/11/25.
//

#ifndef CLUSTERIZER1D_H
#define CLUSTERIZER1D_H

#include "StripMapping.h"
#include <vector>

struct Hit {
    int channel;
    float amplitude;
    float time_ns;
};

struct Cluster {
    char axis;  // 'x' or 'y'
    std::vector<Hit> hits;
};

class Clusterizer1D {
public:
    Clusterizer1D(const StripMapping* mapping);

    void setSpatialThreshold(float mm) { spatialThreshold_ = mm; }
    void setTimeThreshold(float ns) { timeThreshold_ = ns; }

    // Run clustering for one event
    std::vector<Cluster> clusterEvent(const std::vector<Hit>& hits);

private:
    const StripMapping* mapping_;

    float spatialThreshold_ = 2.0f; // mm
    float timeThreshold_ = 10.0f;   // ns

    bool isCompatible(const Hit& h, const Cluster& c) const;
};

#endif
