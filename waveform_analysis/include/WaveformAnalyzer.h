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
#include <algorithm>
#include <cmath>
#include <numeric>
#include <limits>
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

// Put this in your WaveformAnalyzer class header (or adjust to your style)
struct PeakInfo {
    int peakIndex;           // integer sample index of the maximum sample (or plateau center)
    float peakAmplitude;     // baseline-subtracted amplitude from parabola fit
    float peakMax;           // maximum sample amplitude (baseline-subtracted)
    float peakSample;        // sub-sample peak position (sample units) from parabola
    float timingSample;      // x% timing sample (sample units)
    int leftCrossIdx;        // integer sample index where waveform fell/rises below threshold on left (exclusive)
    int rightCrossIdx;       // same for right
    float timeOverThreshold; // in sample units (rightCrossIdx - leftCrossIdx)
    float integral;          // baseline-subtracted sum of samples between left+1 .. right-1
    bool saturated;          // whether the peak shows saturation (flat top near ADC max)
    float localBaseline;     // baseline used for this peak (same units as samples)
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
    bool commonNoiseSubtraction = true;  // if true, subtract common noise per event per channel
    bool allowMultiplePeaks = true;  // if false, only the highest peak per channel per event is kept
    bool local_baseline = false;  // if true, use local baseline per peak; if false, use global pedestal mean
    float thresholdSigma = 5.0;  // Number of pedestal RMS above which a hit is registered
    int peakMergeDistance = 5;  // number of samples within which peaks are merged

    int minSamplesForPeak = 3;  // minimum number of samples above threshold to consider a peak
    int minWidthSamples = 2;  // minimum width in samples above threshold to consider a pulse
    int baselineLeftWindow = 8;  // number of samples left of leading edge used to estimate baseline (robust)
    float satFrac = 0.94;  // fraction of max_adc above which samples are considered saturated for peak saturation detection

    float zeroSupressedBaseline = 256.0f; // baseline level for zero-suppressed pedestal-subtracted waveforms

    // float timePerSample = 20.0;  // ns per sample. Sampling period
    float timePerSample = 60.0;  // ns per sample. Sampling period
    float timePerFtst = 10.0;  // ns per fine timestamp unit. Fixed by DREAM clock of 100MHz --> 10 ns. Shift the timestamp by this amount.
    float timePerTimestamp = 10.0;  // ns Timestamp is in clock cycles of 10 ns
    float timingPercentMax = 0.3;  // fraction of peak amplitude at which timing is calculated
    std::string timingMethod = "percent_max";  // "percent_max" or "parabola"
    // std::string timingMethod = "parabola";  // "percent_max" or "parabola"

    int max_adc = 4095;  // maximum ADC value (saturation level)

    // helpers
    ChannelPedestal computePedestalForChannel(int ch);
    float subtractPedestal(int ch, float ampl) const;

    static std::unordered_map<int, std::vector<float>> fillZeroSuppressedSamples(
                                                        const std::unordered_map<int, std::vector<float>>& waves,
                                                        const std::unordered_map<int, std::vector<int>>& samplesByCh);

    std::vector<PeakInfo> analyzeWaveform(
    const std::vector<float>& wf,
    float noiseRMS,
    float adcMax          // ADC saturation value for detection
    ) const;

    void saturatedLinearExtrapolation(
    int satStartIdx,
    int satEndIdx,
    const float* wf,
    int N,
    float baseline,
    float& peakSample,
    float& peakAmpFit
    ) const;

};

#endif
