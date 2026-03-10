//
// Created by dn277127 on 2025-12-05.
//

#include "WaveformAnalyzer.h"
#include <iostream>
#include <cmath>
#include <limits>
#include <algorithm>

WaveformAnalyzer::WaveformAnalyzer(const std::string& inputFileName,
                                   const std::string& outputFileName,
                                   const std::string& pedestalFileName)
        : inputFileName(inputFileName),
          outputFileName(outputFileName),
          pedestalFileName(pedestalFileName)
{
    hasPedestal = !pedestalFileName.empty();
}

void WaveformAnalyzer::computePedestals() {
    if (pedestalFileName.empty()) {
        std::cout << "No pedestal file provided — assuming zero suppressed and using 256.\n";
        // For zero-suppressed and pedestal subtracted data, the pedestals for each strip are set to 256.
        // If channels is not suppressed, it should be a hit, so set RMS to 1 ADC count
        for (int ch = 0; ch < 512; ch++) {
            pedestalMap[ch] = {zeroSupressedBaseline, 1.0f};
        }
    }
    else  // compute from pedestal file
    {
        std::cout << "Computing pedestals from " << pedestalFileName << "\n";

        TFile f(pedestalFileName.c_str(), "READ");
        if (!f.IsOpen()) {
            std::cerr << "Cannot open pedestal file.\n";
            return;
        }

        TTree* nt = (TTree*)f.Get("nt");
        if (!nt) {
            std::cerr << "Pedestal file has no nt tree.\n";
            return;
        }

        std::vector<UInt_t>* channel = nullptr;
        std::vector<UShort_t>* sample = nullptr;
        std::vector<UShort_t>* amplitude = nullptr;

        nt->SetBranchAddress("channel", &channel);
        nt->SetBranchAddress("sample", &sample);
        nt->SetBranchAddress("amplitude", &amplitude);

        std::unordered_map<int, PedestalData> accum;

        Long64_t nentries = nt->GetEntries();
        for (Long64_t i = 0; i < nentries; i++) {
            nt->GetEntry(i);

            for (size_t j = 0; j < channel->size(); j++) {
                int ch = (*channel)[j];
                float amp = (*amplitude)[j];

                auto& pd = accum[ch];
                pd.sum   += amp;
                pd.sumsq += amp * amp;
                pd.count++;
            }
        }

        // compute mean and RMS for each channel
        pedestalMap.clear();
        for (auto& kv : accum) {
            int ch = kv.first;
            const PedestalData& pd = kv.second;

            float mean = pd.sum / pd.count;
            float rms  = std::sqrt(pd.sumsq / pd.count - mean * mean);

            pedestalMap[ch] = {mean, rms};
        }
        f.Close();
    }

    // Write pedestal tree to output
    TFile fout(outputFileName.c_str(), "UPDATE");
    TTree ped("pedestals", "pedestal values");

    UShort_t ch;
    Float_t mean, rms;

    ped.Branch("channel", &ch, "channel/s");
    ped.Branch("mean",    &mean, "mean/F");
    ped.Branch("rms",     &rms,  "rms/F");

    for (auto& kv : pedestalMap) {
        ch = kv.first;
        mean = kv.second.mean;
        rms  = kv.second.rms;
        ped.Fill();
    }

    ped.Write();

    // Prepare data arrays for TGraph
    int nChannels = pedestalMap.size();
    if (nChannels == 0) {
        std::cout << "No channels found to graph.\n";
        fout.Close();
        return;
    }

    // We must use arrays or vectors of doubles for TGraph
    std::vector<double> chVec, meanVec, rmsVec;
    for (auto const& [channel, data] : pedestalMap) {
        chVec.push_back((double)channel);
        meanVec.push_back(data.mean);
        rmsVec.push_back(data.rms);
    }

    // Create Graphs
    TGraph* gMean = new TGraph(nChannels, chVec.data(), meanVec.data());
    gMean->SetName("g_mean_vs_channel");
    gMean->SetTitle("Pedestal Mean vs. Channel;Channel;Mean Amplitude (ADC)");
    gMean->SetMarkerStyle(20); // Circle marker
    gMean->SetMarkerSize(0.8);
    gMean->SetMarkerColor(kBlue);
    gMean->SetLineColor(kBlue);

    TGraph* gRMS = new TGraph(nChannels, chVec.data(), rmsVec.data());
    gRMS->SetName("g_rms_vs_channel");
    gRMS->SetTitle("Pedestal RMS vs. Channel;Channel;RMS (ADC)");
    gRMS->SetMarkerStyle(20);
    gRMS->SetMarkerSize(0.8);
    gRMS->SetMarkerColor(kRed);
    gRMS->SetLineColor(kRed);

    // Create Canvas and Draw
    TCanvas* c1 = new TCanvas("c_pedestals", "Pedestal Analysis Graphs", 1000, 500);
    c1->Divide(2, 1);

    // Draw Mean Graph
    c1->cd(1);
    gMean->Draw("AP");

    // Draw RMS Graph
    c1->cd(2);
    gRMS->Draw("AP");

    // Write the canvas (which contains both graphs) to the file
    c1->Write("", TObject::kOverwrite);

    // Clean up TObjects (ROOT manages memory, but we delete the canvas)
    delete c1;

    // ----------------------------------------------------
    // END: Graphing Logic
    // ----------------------------------------------------

    fout.Close();

    std::cout << "Pedestals TTree and analysis graphs written to output file.\n";

    fout.Close();

    std::cout << "Pedestals written to output file.\n";
}


float WaveformAnalyzer::subtractPedestal(int ch, float ampl) const {
    auto it = pedestalMap.find(ch);
    if (it == pedestalMap.end()) return ampl;
    return ampl - it->second.mean;
}


// ---------------------------------------------------------
// Common-noise subtraction across 64-channel blocks
// ---------------------------------------------------------
std::unordered_map<int, std::vector<float>>
applyCommonNoiseSubtraction(
        const std::unordered_map<int, std::vector<float>>& waves,
        const std::unordered_map<int, std::vector<int>>& samplesByCh)
{
    // -------- Step 1: build aligned waveforms (dense vectors) --------

    // Get max sample index per channel
    std::unordered_map<int, int> maxSample;
    for (auto &kv : samplesByCh) {
        maxSample[kv.first] =
                *std::max_element(kv.second.begin(), kv.second.end());
    }

    // Build aligned waveforms
    std::unordered_map<int, std::vector<float>> aligned;
    aligned.reserve(waves.size());

    for (auto &kv : waves) {
        int ch = kv.first;
        int maxS = maxSample.at(ch);      // safe because waves and samplesByCh match
        aligned[ch].assign(maxS + 1, 0.0);

        auto &samps = samplesByCh.at(ch);
        auto &amps  = kv.second;

        for (size_t k = 0; k < samps.size(); k++)
            aligned[ch][samps[k]] = amps[k];
    }

    // -------- Step 2: Common-noise subtraction in 64-channel blocks --------

    static constexpr int blockSize = 64;
    static constexpr int maxChannels = 512;   // Max for DREAM FEU

    for (int blockStart = 0; blockStart < maxChannels; blockStart += blockSize)
    {
        int blockEnd = blockStart + blockSize;

        // Find maximum waveform length in this block
        int maxLen = 0;
        for (int ch = blockStart; ch < blockEnd; ch++) {
            if (aligned.count(ch))
                maxLen = std::max(maxLen, (int)aligned[ch].size());
        }

        // For each sample index, compute median
        for (int s = 0; s < maxLen; s++) {
            std::vector<float> vals;
            vals.reserve(blockSize);

            // Collect sample across channels
            for (int ch = blockStart; ch < blockEnd; ch++) {
                if (!aligned.count(ch)) continue;
                auto &wf = aligned[ch];
                if (s < (int)wf.size())
                    vals.push_back(wf[s]);
            }

            if (vals.empty()) continue;

            // Median
            std::nth_element(vals.begin(), vals.begin() + vals.size()/2, vals.end());
            float median = vals[vals.size()/2];

            // Subtract from each waveform
            for (int ch = blockStart; ch < blockEnd; ch++) {
                if (!aligned.count(ch)) continue;
                auto &wf = aligned[ch];
                if (s < (int)wf.size())
                    wf[s] -= median;
            }
        }
    }

    return aligned;   // cleaned waveforms
}


void WaveformAnalyzer::analyzeWaveforms() {
    std::cout << "Analyzing waveforms from " << inputFileName << "\n";

    TFile f(inputFileName.c_str(), "READ");
    TTree* nt = (TTree*)f.Get("nt");

    if (!nt) {
        std::cerr << "Error: input file missing nt tree.\n";
        return;
    }

    // input branches
    ULong64_t eventID, timestamp;
    UShort_t ftst;

    std::vector<UShort_t>* sample = nullptr;
    std::vector<UInt_t>*   channel = nullptr;
    std::vector<UShort_t>* amplitude = nullptr;

    nt->SetBranchAddress("eventId", &eventID);
    nt->SetBranchAddress("timestamp", &timestamp);
    nt->SetBranchAddress("ftst", &ftst);
    nt->SetBranchAddress("sample", &sample);
    nt->SetBranchAddress("channel", &channel);
    nt->SetBranchAddress("amplitude", &amplitude);

    // output file
    TFile fout(outputFileName.c_str(), "UPDATE");

    TTree hitTree("hits", "reconstructed hits");

    ULong64_t out_eventID;
    ULong64_t out_trigger_timestamp_ns;
    UShort_t  out_channel;
    Float_t   out_amp;
    Float_t   out_time_ns;
    Float_t   out_time_of_max_ns;
    Float_t   out_sample;
    Float_t   out_max_sample;
    Float_t   out_local_baseline;
    Float_t   out_local_max;
    Float_t   out_left_sample;
    Float_t   out_right_sample;
    Float_t   out_time_over_threshold;
    Float_t   out_integral;
    Bool_t    out_saturated;

    hitTree.Branch("eventId", &out_eventID, "eventId/l");
    hitTree.Branch("trigger_timestamp_ns", &out_trigger_timestamp_ns, "trigger_timestamp_ns/l");
    hitTree.Branch("channel", &out_channel, "channel/s");
    hitTree.Branch("amplitude", &out_amp, "amplitude/F");
    hitTree.Branch("time", &out_time_ns, "time/F");
    hitTree.Branch("time_of_max", &out_time_of_max_ns, "time_of_max/F");
    hitTree.Branch("sample", &out_sample, "sample/F");
    hitTree.Branch("max_sample", &out_max_sample, "max_sample/F");
    hitTree.Branch("local_baseline", &out_local_baseline, "local_baseline/F");
    hitTree.Branch("local_max", &out_local_max, "local_max/F");
    hitTree.Branch("left_sample", &out_left_sample, "left_sample/F");
    hitTree.Branch("right_sample", &out_right_sample, "right_sample/F");
    hitTree.Branch("time_over_threshold", &out_time_over_threshold, "time_over_threshold/F");
    hitTree.Branch("integral", &out_integral, "integral/F");
    hitTree.Branch("saturated", &out_saturated, "saturated/O");

    Long64_t nentries = nt->GetEntries();

    for (Long64_t i = 0; i < nentries; i++) {
        nt->GetEntry(i);

        // group hits by channel
        std::unordered_map<int, std::vector<float>> waves;
        std::unordered_map<int, std::vector<int>> samplesByCh;

        for (size_t j = 0; j < channel->size(); j++) {
            int ch = (*channel)[j];
            int s  = (*sample)[j];
            float a = subtractPedestal(ch, (*amplitude)[j]);

            waves[ch].push_back(a);
            samplesByCh[ch].push_back(s);
        }
        if(commonNoiseSubtraction) {
            waves = applyCommonNoiseSubtraction(waves, samplesByCh);
        }

        // This line updates the 'waves' map to contain the dense, regular waveforms.
        waves = fillZeroSuppressedSamples(waves, samplesByCh);

        // analyze per channel
        for (auto& kv : waves) {
            int ch = kv.first;
            std::vector<float>& amps = kv.second;

            float noiseRMS = pedestalMap.count(ch) ? pedestalMap[ch].rms : 3.0f;
            float max_adc_ped_sub = max_adc - (pedestalMap.count(ch) ? pedestalMap[ch].mean : 0.0f);
            auto peaks = analyzeWaveform(amps, noiseRMS, max_adc_ped_sub);
            for (auto& peak : peaks) {

                // Correct samples for ftst. ftst in units of clock cycles
                peak.peakSample += static_cast<float>(ftst) * timePerFtst / timePerSample;
                peak.timingSample += static_cast<float>(ftst) * timePerFtst / timePerSample;

                out_eventID        = eventID;
                out_trigger_timestamp_ns = timestamp * static_cast<int>(timePerTimestamp);
                out_channel        = ch;
                out_amp            = peak.peakAmplitude;
                out_local_max      = peak.peakMax;
                out_time_ns        = peak.timingSample * timePerSample;
                out_time_of_max_ns = peak.peakSample * timePerSample;
                out_sample         = peak.timingSample;
                out_max_sample     = peak.peakSample;
                out_local_baseline = peak.localBaseline;
                out_left_sample    = peak.leftCrossIdx + ftst * timePerFtst / timePerSample;
                out_right_sample   = peak.rightCrossIdx + ftst * timePerFtst / timePerSample;
                out_time_over_threshold = peak.timeOverThreshold * timePerSample;
                out_integral       = peak.integral;
                out_saturated      = peak.saturated;

                hitTree.Fill();

                if (!allowMultiplePeaks)
                    break;
            }
        }
    }

    hitTree.Write();
    fout.Close();
}


/**
 * @brief Converts zero-suppressed (sparse) waveform data to regular (dense) data
 * by filling missing samples with a zero amplitude.
 *
 * The length of the final waveform for each channel is determined by the
 * maximum sample index present in the input data + 1.
 *
 * @param waves The input map of channel ID -> vector of amplitudes (ZS data).
 * @param samplesByCh The input map of channel ID -> vector of sample indices (ZS data).
 * @return A new std::unordered_map<int, std::vector<float>> with regular (dense) waveforms.
 */
std::unordered_map<int, std::vector<float>> WaveformAnalyzer::fillZeroSuppressedSamples(
    const std::unordered_map<int, std::vector<float>>& waves,
    const std::unordered_map<int, std::vector<int>>& samplesByCh)
{
    // The new map to hold the dense, regular waveforms
    std::unordered_map<int, std::vector<float>> regularWaves;

    // Iterate over every channel present in the input data
    for (const auto& pair : waves) {
        int channelID = pair.first;
        const std::vector<float>& zsAmplitudes = pair.second;

        // Retrieve corresponding sample indices
        auto samplesIt = samplesByCh.find(channelID);
        if (samplesIt == samplesByCh.end()) {
            std::cerr << "Warning: No sample indices found for channel " << channelID << ". Skipping." << std::endl;
            regularWaves[channelID] = {};
            continue;
        }

        const std::vector<int>& zsIndices = samplesIt->second;

        // Basic validation
        if (zsAmplitudes.size() != zsIndices.size()) {
            std::cerr << "Error: Mismatch in size between amplitudes (" << zsAmplitudes.size()
                      << ") and indices (" << zsIndices.size() << ") for channel " << channelID << ". Skipping." << std::endl;
            regularWaves[channelID] = {};
            continue;
        }

        // --- 1. Determine the maximum sample index to size the new waveform ---
        int maxIndex = -1;
        if (!zsIndices.empty()) {
            // Finds the largest index recorded in the zero-suppressed data
            maxIndex = *std::max_element(zsIndices.begin(), zsIndices.end());
        }

        if (maxIndex < 0) {
            regularWaves[channelID] = {};
            continue;
        }

        // The size of the regular waveform is maxIndex + 1 (0-based indexing)
        size_t regularWaveSize = static_cast<size_t>(maxIndex) + 1;

        // --- 2. Initialize the regular waveform with zeros (0.0f) ---
        // This vector now represents samples [0, 1, ..., maxIndex], all set to 0.0f
        std::vector<float> regularWave(regularWaveSize, 0.0f);

        // --- 3. Fill in the recorded (non-zero-suppressed) values ---
        for (size_t i = 0; i < zsAmplitudes.size(); ++i) {
            int index = zsIndices[i];
            float amplitude = zsAmplitudes[i];

            // Place the recorded amplitude at its correct, absolute index
            // The check below is technically redundant due to maxIndex, but safe.
            if (index < regularWaveSize) {
                regularWave[index] = amplitude;
            }
        }

        // --- 4. Store the new regular waveform ---
        // std::move is used for efficiency to transfer ownership of the vector data.
        regularWaves[channelID] = std::move(regularWave);
    }

    return regularWaves;
}


/// ---------------------------------------------------------------------------
/// findPulseRegions  –  derivative-based, pile-up-aware pulse separator
/// ---------------------------------------------------------------------------
/// Strategy
/// --------
/// 1. Smooth the waveform with a box-car of half-width derivativeSmoothWidth.
/// 2. Differentiate (central differences) to get the instantaneous slope.
/// 3. Locate every local maximum of the derivative that exceeds
///    derivThresholdSigma * noiseRMS / sample ("rising-edge seed").
/// 4. Merge seeds that are closer than derivMergeDistance samples.
/// 5. For each seed, walk left to where the waveform first drops to/below
///    the amplitude threshold (= thresholdSigma * noiseRMS) — this is the
///    pulse start.  Walk right to either:
///      (a) the sample just before the *next* seed's start (pile-up cut), or
///      (b) where the waveform falls back below threshold (isolated pulse).
/// 6. Reject regions narrower than minWidthSamples.
///
/// The resulting [start, end] intervals are non-overlapping and cover every
/// pulse the derivative can resolve, including pulses riding on a tail.
/// ---------------------------------------------------------------------------
std::vector<WaveformAnalyzer::PulseRegion>
WaveformAnalyzer::findPulseRegions(
    const std::vector<float>& wf,
    float noiseRMS) const
{
    std::vector<PulseRegion> regions;
    const int N = (int)wf.size();
    if (N < 3) return regions;

    const float ampThr   = thresholdSigma  * noiseRMS;
    const float derivThr = derivThresholdSigma * noiseRMS; // per-sample slope threshold

    // ------------------------------------------------------------------ //
    // Step 1 – box-car smooth                                             //
    // ------------------------------------------------------------------ //
    const int hw = std::max(0, derivativeSmoothWidth); // half-width
    std::vector<float> smooth(N, 0.0f);
    for (int i = 0; i < N; ++i) {
        int lo = std::max(0, i - hw);
        int hi = std::min(N - 1, i + hw);
        float s = 0.0f;
        for (int k = lo; k <= hi; ++k) s += wf[k];
        smooth[i] = s / float(hi - lo + 1);
    }

    // ------------------------------------------------------------------ //
    // Step 2 – central-difference derivative                              //
    // ------------------------------------------------------------------ //
    std::vector<float> deriv(N, 0.0f);
    for (int i = 1; i < N - 1; ++i)
        deriv[i] = 0.5f * (smooth[i + 1] - smooth[i - 1]);
    deriv[0]     = smooth[1] - smooth[0];
    deriv[N - 1] = smooth[N - 1] - smooth[N - 2];

    // ------------------------------------------------------------------ //
    // Step 3 – find local maxima of deriv above derivThr ("seeds")       //
    // Uses an arm/trigger/reset state machine so a new peak is only      //
    // accepted after the derivative has fallen back below derivResetThr. //
    // ------------------------------------------------------------------ //
    // derivResetThr should be <= derivThr; a good starting value is 0    //
    // (derivative must return to ~flat) or a small positive fraction of  //
    // derivThr (e.g. 0.3 * derivThr) to tolerate a slowly-falling tail. //

    std::vector<int> seeds;
    int    pendingPeak  = -1;     // index of the best candidate seen so far
    float  pendingVal   = -std::numeric_limits<float>::max();
    bool   armed        = true;   // ready to accept a new peak?

    for (int i = 1; i < N - 1; ++i) {
        if (armed) {
            // Track the running maximum while above derivThr
            if (deriv[i] > derivThr) {
                if (deriv[i] > pendingVal) {
                    pendingVal  = deriv[i];
                    pendingPeak = i;
                }
            } else if (pendingPeak >= 0) {
                // Just dropped below derivThr — commit the best peak seen
                seeds.push_back(pendingPeak);
                pendingPeak = -1;
                pendingVal  = -std::numeric_limits<float>::max();
                armed       = false;   // must reset before next peak
            }
        } else {
            // Waiting for derivative to fall back below the reset threshold
            if (deriv[i] <= derivResetThr) {
                armed = true;
            }
        }
    }
    // Flush a pending peak at end-of-waveform
    if (armed && pendingPeak >= 0)
        seeds.push_back(pendingPeak);

    if (seeds.empty()) return regions;

    // ------------------------------------------------------------------ //
    // Step 4 – merge seeds within derivMergeDistance                     //
    // ------------------------------------------------------------------ //
    std::vector<int> merged;
    merged.push_back(seeds[0]);
    for (int k = 1; k < (int)seeds.size(); ++k) {
        if (seeds[k] - merged.back() <= derivMergeDistance) {
            // keep the one with the larger derivative
            if (deriv[seeds[k]] > deriv[merged.back()])
                merged.back() = seeds[k];
        } else {
            merged.push_back(seeds[k]);
        }
    }

    // ------------------------------------------------------------------ //
    // Step 5 – assign pulse regions (two-pass boundary negotiation)       //
    // ------------------------------------------------------------------ //

    // Pass 1: find the natural right boundary for every seed independently.
    // Stop at whichever comes first: waveform drops below ampThr, or we
    // reach the next detected seed (the start of the next pulse's rise).
    std::vector<int> rightEnds(merged.size());
    for (int k = 0; k < (int)merged.size(); ++k) {
        int end        = merged[k];
        int nextSeedAt = (k + 1 < (int)merged.size()) ? merged[k + 1] : N;

        while (end + 1 < N && wf[end + 1] > ampThr && end + 1 < nextSeedAt)
            ++end;

        rightEnds[k] = end;
    }

    // Pass 2: assign left boundaries and resolve overlaps between neighbours.
    // For each seed k:
    //   - Its natural left boundary is "walk left until below threshold".
    //   - If that boundary overlaps with the previous seed's region, the two
    //     regions share a contested zone [prevSeed .. thisSeed].  We place the
    //     split at the local minimum in that valley, giving each pulse the
    //     downslope that belongs to it.
    for (int k = 0; k < (int)merged.size(); ++k) {
        int seed   = merged[k];
        int endIdx = rightEnds[k];

        // --- natural left boundary ---
        // Stop when either:
        //   (a) amplitude drops back below ampThr, or
        //   (b) derivative drops back below derivResetThr (same reset threshold
        //       used in Step 3 — the waveform has stopped rising sharply enough
        //       that we are on the tail of a previous pulse)
        int startIdx = seed;
        while (startIdx > 0
               && wf[startIdx - 1]     > ampThr
               && deriv[startIdx - 1]  > derivResetThr)
            --startIdx;

        // --- negotiate with previous region only if still overlapping ---
        // If either threshold above already stopped the walk inside the previous
        // region, no negotiation is needed. If the walk overshot (neither
        // threshold fired before we hit the previous region), fall back to the
        // valley-minimum split.
        if (k > 0) {
            int prevSeed   = merged[k - 1];
            int prevEndIdx = rightEnds[k - 1];

            if (startIdx <= prevEndIdx) {
                // Neither threshold separated the pulses — find the local minimum
                // in the inter-seed valley as the physical split point.
                int   valleyMin = prevSeed;
                float valleyVal = wf[prevSeed];
                for (int s = prevSeed; s <= seed && s < N; ++s) {
                    if (wf[s] < valleyVal) {
                        valleyVal = wf[s];
                        valleyMin = s;
                    }
                }

                rightEnds[k - 1] = valleyMin - 1;
                startIdx          = valleyMin;
            }
        }

        // --- right pile-up cut with *next* seed (identical to before) ---
        if (k + 1 < (int)merged.size()) {
            int nextSeed  = merged[k + 1];
            int nextStart = nextSeed;
            while (nextStart > 0 && wf[nextStart - 1] > ampThr)
                --nextStart;

            if (nextStart <= endIdx) {
                int   splitAt = seed;
                float minVal  = wf[seed];
                for (int s = seed; s < nextSeed && s < N; ++s) {
                    if (wf[s] < minVal) { minVal = wf[s]; splitAt = s; }
                }
                endIdx = splitAt;
                // Keep rightEnds consistent so the next iteration sees the update.
                rightEnds[k] = endIdx;
            }
        }

        if (endIdx - startIdx + 1 < minWidthSamples) continue;
        regions.push_back({startIdx, endIdx});
    }

    return regions;
}



std::vector<PeakInfo> WaveformAnalyzer::analyzeWaveform(
    const std::vector<float>& wf,
    float noiseRMS,
    float adcMax = 4000.0f          // ADC saturation value for detection
) const
{
    std::vector<PeakInfo> results;
    if (wf.size() < 3) return results;  // Waveform is too short

    const int N = (int)wf.size();
    const float dynThr = thresholdSigma * noiseRMS;

    // ------------------------------------------------------------------
    // Build the list of pulse regions to analyse.
    // In derivative mode we use findPulseRegions() to separate piled-up
    // pulses; otherwise we fall back to the classic threshold scan so
    // behaviour is identical to the original for non-pileup data.
    // ------------------------------------------------------------------
    struct Region { int startIdx; int endIdx; };
    std::vector<Region> pulseRegions;

    if (useDerivativeTrigger) {
        for (auto& r : findPulseRegions(wf, noiseRMS))
            pulseRegions.push_back({r.start, r.end});
    } else {
        // Legacy threshold scan: collect contiguous above-threshold windows
        int i = 0;
        while (i < N) {
            if (wf[i] <= dynThr) { ++i; continue; }
            int startIdx = i;
            int j = i;
            while (j + 1 < N && wf[j + 1] > dynThr) ++j;
            pulseRegions.push_back({startIdx, j});
            i = j + 1;
        }
    }

    // ------------------------------------------------------------------
    // Analyse each region with the existing peak-analysis logic.
    // The only change vs the original code is that startIdx/endIdx now
    // come from the region list above rather than from a threshold scan.
    // All arithmetic below uses absolute sample indices into wf[].
    // ------------------------------------------------------------------
    for (auto& reg : pulseRegions) {
        int startIdx = reg.startIdx;
        int endIdx   = reg.endIdx;

        // Width guard (same as original)
        if (endIdx - startIdx + 1 < minWidthSamples) continue;

        // Within [startIdx..endIdx] find plateau-aware maximum sample (middle of plateau)
        int maxIdx = startIdx;
        float maxVal = wf[startIdx];
        int k = startIdx;
        while (k <= endIdx) {
            // handle plateau runs: detect equal neighbor values
            int runStart = k;
            int runEnd = k;
            while (runEnd + 1 <= endIdx && wf[runEnd + 1] == wf[k]) ++runEnd;
            // choose plateau middle sample as candidate
            int candidate = (runStart + runEnd) / 2;
            if (wf[candidate] > maxVal) {
                maxVal = wf[candidate];
                maxIdx = candidate;
            }
            k = runEnd + 1;
        }

        // Estimate a local baseline robustly from samples left of the leading edge.
        // We'll use the median of up to baselineLeftWindow samples immediately left of startIdx.
        float baseline = 0.0f;
        if (local_baseline)
        {
            // ── Collect the pre-pulse window ─────────────────────────────────────
            // Use samples [startIdx - baselineLeftWindow, maxIdx), clamped to [0, N).
            int leftStart = std::max(0, startIdx - baselineLeftWindow);
            int leftCount = startIdx - leftStart;

            if (leftCount <= 0) {
                baseline = 0.0f;
            } else {
                // Use minimum of these points
                float minVal = wf[leftStart];
                for (int l = leftStart + 1; l < maxIdx; ++l) {
                    if (wf[l] < minVal) minVal = wf[l];
                }
                baseline = minVal;
            }
        }

        // Detect saturation (flat top near adcMax)
        bool saturated = false;
        int satStart = -1, satEnd = -1;
        float satThreshold = satFrac * adcMax;
        {
            for (int t = startIdx; t <= endIdx; ++t) {
                if (wf[t] >= satThreshold) {
                    saturated = true;
                    if (satStart < 0) satStart = t;
                    satEnd = t;
                }
            }
        }

        // Find peak amplitude and time of max
        float peakAmpFit = wf[maxIdx] - baseline;
        float peakSample = maxIdx;
        if (saturated)  // For saturated peaks, do linear extrapolation from edges of saturated region
        {
            saturatedLinearExtrapolation(
                satStart,
                satEnd,
                wf.data(),
                N,
                baseline,
                peakSample,
                peakAmpFit
            );
        }
        else  // Fit parabola around the integer maxIdx if possible to get sub-sample peak
        {
            float peakOffset = 0.0f;
            if (maxIdx > 0 && maxIdx < N-1) {
                float y1 = wf[maxIdx-1] - baseline;
                float y2 = wf[maxIdx]   - baseline;
                float y3 = wf[maxIdx+1] - baseline;
                float denom = (y1 - 2*y2 + y3);
                if (std::abs(denom) > 1e-9f) {
                    peakOffset = 0.5f * (y1 - y3) / denom; // in samples
                    peakAmpFit = y2 - 0.25f * (y1 - y3) * peakOffset;
                } else {
                    peakOffset = 0.0f;
                    peakAmpFit = y2;
                }
            }
            peakSample = maxIdx + peakOffset;
        }

        // find x% timing on leading edge relative to baseline-subtracted peak amplitude
        float fraction = timingPercentMax; // e.g. 0.5
        float timing_amp = std::min(peakAmpFit, adcMax);  // Not toally ideal, want to tweak for baseline and pedestal probably
        float target = fraction * timing_amp + baseline;
        // float target = fraction * peakAmpFit + baseline;
        int leadIdx = maxIdx;
        while (leadIdx > 0 && leadIdx >= startIdx && wf[leadIdx] > target) --leadIdx;
        // linear interp between leadIdx and leadIdx+1
        float timingSample = (float)leadIdx;
        if (leadIdx >= 0 && leadIdx+1 < N) {
            float yL = wf[leadIdx];
            float yH = wf[leadIdx+1];
            if (std::abs(yH - yL) > 1e-9f) {
                timingSample = leadIdx + (target - yL) / (yH - yL);
            } else {
                timingSample = leadIdx + 0.5f;
            }
        } else {
            timingSample = (float)maxIdx;
        }

        // For the left/right threshold crossings we clamp within the region
        // so that pile-up neighbours are not included in TOT / integral.
        // Walk left from maxIdx until waveform falls <= dynamic threshold (or touches baseline) to find left crossing.
        int leftCross = maxIdx;
        while (leftCross > startIdx && wf[leftCross] > dynThr) --leftCross;
        // Walk right from maxIdx until waveform falls <= dynamic threshold, clamped at region end.
        int rightCross = maxIdx;
        while (rightCross + 1 <= endIdx && wf[rightCross] > dynThr) ++rightCross;

        // Compute integral and TOT using baseline-subtracted samples between leftCross+1 .. rightCross-1 (inclusive)
        float integral = 0.0f;
        for (int tt = leftCross + 1; tt <= rightCross - 1; ++tt) {
            if (tt >= 0 && tt < N) integral += (wf[tt] - baseline);
        }

        float tot = float(std::max(0, rightCross - leftCross - 1)); // number of samples above threshold

        // Build PeakInfo
        PeakInfo pi;
        pi.peakIndex = maxIdx;
        pi.peakAmplitude = peakAmpFit;
        pi.peakMax = wf[maxIdx] - baseline;
        pi.peakSample = peakSample;
        pi.timingSample = timingSample;
        pi.leftCrossIdx = leftCross;
        pi.rightCrossIdx = rightCross;
        pi.timeOverThreshold = tot;
        pi.integral = integral;
        pi.saturated = saturated;
        pi.localBaseline = baseline;

        // Only accept if amplitude large enough relative to noise
        if (pi.peakAmplitude >= thresholdSigma * noiseRMS && (endIdx - startIdx + 1) >= minSamplesForPeak) {
            results.push_back(pi);
        }
    }

    return results;
}


// Function to perform linear extrapolation to find peak for saturated pulses
void WaveformAnalyzer::saturatedLinearExtrapolation(
    int satStartIdx,
    int satEndIdx,
    const float* wf,
    int N,
    float baseline,
    float& peakSample,
    float& peakAmpFit
) const {
    // Get linear extrapolations from the edges of the saturated region
    float leftSlope = 0.0f;
    float rightSlope = 0.0f;
    if (satStartIdx > 0) {
        leftSlope = (wf[satStartIdx] - wf[satStartIdx - 1]) ;
    }
    if (satEndIdx + 1 < N) {
        rightSlope = (wf[satEndIdx + 1] - wf[satEndIdx]) ;
    }

    // Ensure we have valid slopes. Left should be positive, right should be negative
    if (leftSlope <= 0 || rightSlope >= 0)
    {
        peakSample = static_cast<float>(satStartIdx);
        peakAmpFit = wf[satStartIdx] - baseline; // fallback to first saturated sample
        return;
    }

    // Estimate peak position by finding intersection of the two lines
    float leftIntercept = wf[satStartIdx] - leftSlope * satStartIdx;
    float rightIntercept = wf[satEndIdx] - rightSlope * satEndIdx;

    peakSample = (rightIntercept - leftIntercept) / (leftSlope - rightSlope);
    peakAmpFit = leftSlope * peakSample + leftIntercept - baseline;
}


int WaveformAnalyzer::run() {
    TFile fout(outputFileName.c_str(), "RECREATE");
    fout.Close();
    computePedestals();  // If no pedestals, this will assume zero-suppressed data and subtract 256
    analyzeWaveforms();
    return 0;
}
