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

struct PeakInfo {
    int peakIndex;           // integer sample index of the maximum sample (or plateau center)
    float peakAmplitude;     // baseline-subtracted amplitude from parabola fit
    float peakMax;           // maximum sample amplitude (baseline-subtracted)
    float peakSample;        // sub-sample peak position (sample units) from parabola
    float timingSample;      // x% timing sample (sample units)
    float leftCross;         // interpolated threshold-crossing sample on the left (float)
    float rightCross;        // interpolated threshold-crossing sample on the right (float)
    float timeOverThreshold; // rightCross - leftCross, sample units (float, interpolated)
    float integral;          // baseline-subtracted sum over the full pulse span (to baseline return)
    bool saturated;          // whether the peak shows saturation (flat top near ADC max)
    bool truncLeft;          // pulse already above threshold at the first sample (rise not recorded)
    bool truncRight;         // pulse still above threshold at the last sample (tail clipped)
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
    void setTimePerSample(float v)     { timePerSample = v; }
    void setCommonNoiseSubtraction(bool v) { commonNoiseSubtraction = v; }  // toggles CNS on DATA (pedestal RMS is always post-CNS)

private:
    std::string inputFileName;
    std::string outputFileName;
    std::string pedestalFileName;
    bool hasPedestal = false;

    std::unordered_map<int, ChannelPedestal> pedestalMap;

    // configuration
    // TODO: expose these through the yaml ConfigManager (common/) instead of
    // compile-time constants — needed once one binary serves both the bench
    // micro-TPC detectors and P2. Documented, deferred.
    bool commonNoiseSubtraction = true;  // if true, subtract common noise per event per channel (Cosmic Bench default ON)
    bool allowMultiplePeaks = true;  // if false, only the LARGEST peak per channel per event is kept
    bool local_baseline = true;  // if true, use local baseline per peak; if false, use global pedestal mean
    float thresholdSigma = 5.0;  // Number of pedestal RMS above which a hit is registered

    int minSamplesForPeak = 3;  // minimum number of samples above threshold to consider a peak
    int minWidthSamples = 2;  // minimum width in samples above threshold to consider a pulse
    int baselineLeftWindow = 4;  // number of samples used for the median local baseline
    int baselineGapSamples = 2;  // guard gap between the baseline window and the pulse start (skips the sub-threshold rise)
    float satFrac = 0.94;  // fraction of max_adc above which samples are considered saturated
    int minSatSamples = 2;  // consecutive samples above satFrac*adcMax required to flag saturation
    int satSlopeFitSamples = 3;  // samples per side used to fit the edge slopes of a saturated peak

    // --- Pulse-region finding (single unified algorithm) ---
    // Contiguous above-threshold runs are the pulse gate; short sub-threshold
    // gaps are bridged; pile-up inside one run is separated by splitting at
    // valleys whose depth below both flanking maxima is significant vs noise.
    int gapMergeSamples = 2;            // sub-threshold gaps <= this many samples do not end a pulse
    float splitProminenceSigma = 4.0f;  // valley depth (units of noiseRMS) needed to split pile-up

    float zeroSupressedBaseline = 256.0f; // baseline level for zero-suppressed pedestal-subtracted waveforms

    // float timePerSample = 20.0;  // ns per sample. Sampling period
    float timePerSample = 60.0;  // ns per sample. Sampling period -- set from commandline
    float timePerFtst = 10.0;  // ns per fine timestamp unit. Fixed by DREAM clock of 100MHz --> 10 ns. Shift the timestamp by this amount.
    float timePerTimestamp = 10.0;  // ns Timestamp is in clock cycles of 10 ns
    float timingPercentMax = 0.3;  // fraction of peak amplitude at which timing is calculated

    int max_adc = 4095;  // maximum ADC value (saturation level)

    // helpers
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

    // Pulse-region finder. Returns [start, end] threshold cores plus the
    // analysis bounds [loBound, hiBound] a pulse may extend into (limited by
    // its neighbours), one region per detected pulse.
    struct PulseRegion { int start; int end; int loBound; int hiBound; };
    std::vector<PulseRegion> findPulseRegions(
        const std::vector<float>& wf,
        float noiseRMS
    ) const;

    static void splitRegionByValleys(
        const std::vector<float>& wf,
        int start, int end,
        float prominence,
        int minHalfWidth,
        std::vector<std::pair<int,int>>& out);
};

#endif
