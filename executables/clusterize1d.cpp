//
// Created by dylan on 12/11/25.
//

#include <iostream>
#include <string>
#include <map>

#include "TFile.h"
#include "TTree.h"

#include "StripMapping.h"
#include "Clusterizer1D.h"

int main() {
    // ------------------------------------------------------------------
    // Hard-coded test files
    // ------------------------------------------------------------------
    std::string hits_file = "/media/dylan/data/x17/dream_run/run_69/hits.root";
    std::string mapping_file = "/home/dylan/PycharmProjects/Cosmic_Bench_DAQ_Control/config/detectors/rd5_map.txt";

    // Load strip mapping
    StripMapping mapping;
    if (!mapping.loadFromCSV(mapping_file)) {
        std::cerr << "Failed to load strip mapping.\n";
        return 1;
    }

    // Print strip mapping
    for (int ch = 0; ch < 256; ++ch) {
        const StripInfo* info = mapping.getStripInfo(ch);
        if (info) {
            std::cout << "Channel " << ch
                      << ": Connector " << info->connector
                      << ", ConnCh " << info->connectorChannel
                      << ", Strip " << info->stripNumber
                      << ", Axis " << info->axis
                      << ", Pitch " << info->pitch
                      << ", X " << info->xGerber
                      << ", Y " << info->yGerber
                      << "\n";
        } else {
            std::cout << "Channel " << ch << ": No mapping info.\n";
        }
    }

    // ------------------------------------------------------------------
    // Initialize clusterizer
    // ------------------------------------------------------------------
    Clusterizer1D clusterizer(&mapping);
    clusterizer.setSpatialThreshold(2.0f);   // mm
    clusterizer.setTimeThreshold(10.0f);     // ns

    // ------------------------------------------------------------------
    // Load hits file
    // ------------------------------------------------------------------
    TFile* f = TFile::Open(hits_file.c_str());
    if (!f || f->IsZombie()) {
        std::cerr << "Failed to open hits file.\n";
        return 1;
    }

    TTree* t = nullptr;
    f->GetObject("hits", t);
    if (!t) {
        std::cerr << "Tree 'hits' not found.\n";
        return 1;
    }

    // ------------------------------------------------------------------
    // Set branch addresses
    // ------------------------------------------------------------------
    ULong64_t eventId;
    UShort_t  channel;
    Float_t   amplitude;
    Float_t   time_ns;
    Float_t   sample;

    t->SetBranchAddress("eventId",   &eventId);
    t->SetBranchAddress("channel",   &channel);
    t->SetBranchAddress("amplitude", &amplitude);
    t->SetBranchAddress("time",      &time_ns);
    t->SetBranchAddress("sample",    &sample);

    Long64_t nEntries = t->GetEntries();
    std::cout << "Hits entries = " << nEntries << "\n";

    // ------------------------------------------------------------------
    // Loop over hits and group per event
    // ------------------------------------------------------------------
    ULong64_t current_event = -1;
    std::vector<Hit> eventHits;

    auto flushEvent = [&](ULong64_t evtId) {
        if (eventHits.empty()) return;

        std::cout << "\n==== Event " << evtId << " ====\n";
        std::cout << "Input hits: " << eventHits.size() << "\n";

        // Run clustering
        auto clusters = clusterizer.clusterEvent(eventHits);

        // Print clusters
        int idx = 0;
        for (const auto& c : clusters) {
            std::cout << "-- Cluster " << idx++
                      << " (axis=" << c.axis << ", nhits=" << c.hits.size() << ")\n";

            for (const auto& h : c.hits) {
                std::cout << "    ch=" << h.channel
                          << " amp=" << h.amplitude
                          << " time=" << h.time_ns
                          << "\n";
            }
        }

        eventHits.clear();
    };

    // ------------------------------------------------------------------
    // Loop over all hits
    // ------------------------------------------------------------------
    for (Long64_t i = 0; i < nEntries; ++i) {
        t->GetEntry(i);

        if (eventId != current_event && current_event != (ULong64_t)-1) {
            // Finish previous event
            flushEvent(current_event);
        }

        if (eventId != current_event) {
            current_event = eventId;
        }

        Hit h;
        h.channel   = channel;
        h.amplitude = amplitude;
        h.time_ns   = time_ns;

        eventHits.push_back(h);
    }

    // Flush final event
    flushEvent(current_event);

    std::cout << "Clustering complete.\n";

    return 0;
}
