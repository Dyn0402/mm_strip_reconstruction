//
// Created by dn277127 on 2025-12-05.
//

#ifndef WAVEFORM_ANALYZER_H
#define WAVEFORM_ANALYZER_H
#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"

struct PedestalData {
    double sum = 0;
    double sumsq = 0;
    uint64_t count = 0;
};

struct ChannelPedestal {
    float mean = 0;
    float rms  = 0;
};

class WaveformAnalyzer {
public:
    WaveformAnalyzer(const std::string& inputFileName,
                     const std::string& outputFileName,
                     const std::string& pedestalFileName = "");

    void computePedestals();
    void analyzeWaveforms();
    int run();

    void setAllowMultiplePeaks(bool v) { allowMultiplePeaks = v; }
    void setThresholdSigma(float v)    { thresholdSigma = v; }

private:
    std::string inputFileName;
    std::string outputFileName;
    std::string pedestalFileName;
    bool hasPedestal = false;

    std::unordered_map<int, ChannelPedestal> pedestalMap;

    // configuration
    bool commonNoiseSubtraction = false;  // if true, subtract common noise per event per channel
    bool allowMultiplePeaks = true;  // if false, only the highest peak per channel per event is kept
    float thresholdSigma = 10.0;  // Number of pedestal RMS above which a hit is registered
    int peakMergeDistance = 5;  // number of samples within which peaks are merged
    float timePerSample = 20.0;  // ns per sample. Sampling period
    float timePerFtst = 10.0;  // ns per fine timestamp unit
    float timePerTimestamp = 10.0;  // ns Timestamp is in clock cycles of 10 ns
    float timingPercentMax = 0.3;  // fraction of peak amplitude at which timing is calculated
    std::string timingMethod = "percent_max";  // "percent_max" or "parabola"

    // helpers
    void loadPedestals();
    ChannelPedestal computePedestalForChannel(int ch);
    float subtractPedestal(int ch, float ampl) const;
    void fitParabola(const std::vector<float>& amps,
                     int maxIndex,
                     float& peakAmp, float& peakOffset) const;
    float findxPercentofMax(const std::vector<float>& amps,
                           int maxIndex,
                           float fraction) const;

    std::vector<int> findPeakIndices(const std::vector<float>& wf, float noiseRMS) const;
    std::vector<int> mergePeaks(const std::vector<int>& peaks, const std::vector<float>& wf) const;

};

#endif
