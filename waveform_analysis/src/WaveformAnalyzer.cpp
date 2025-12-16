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


/// New unified analyzer function
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

    int i = 0;
    while (i < N) {
        // Skip until we find a sample above dynamic threshold
        if (wf[i] <= dynThr) { ++i; continue; }

        // Found a sample above threshold: treat this as a candidate pulse start
        int startIdx = i;

        // Walk forward until waveform goes below threshold to locate the pulse window.
        int j = i;
        while (j + 1 < N && wf[j + 1] > dynThr) ++j;
        int endIdx = j;

        // Ensure width is meaningfully wide, otherwise skip past
        if (endIdx - startIdx + 1 < minWidthSamples) {
            i = endIdx + 1;
            continue;
        }

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
            int leftStart = std::max(0, startIdx - baselineLeftWindow);
            int leftCount = startIdx - leftStart;
            if (leftCount <= 0) {
                // fallback: use very local min before peak or zero
                baseline = 0.0f;
            } else {
                // copy into temp and take median (robust against previous pulses)
                std::vector<float> tmp;
                tmp.reserve(leftCount);
                for (int tt = leftStart; tt < startIdx; ++tt) tmp.push_back(wf[tt]);
                std::nth_element(tmp.begin(), tmp.begin() + tmp.size()/2, tmp.end());
                baseline = tmp[tmp.size()/2];
                if (tmp.size() > 1 && tmp.size()%2 == 0) {
                    // average with lower median for even count (not necessary but fine)
                    auto it = std::max_element(tmp.begin(), tmp.begin() + tmp.size()/2);
                    baseline = 0.5f * (baseline + *it);
                }
            }
        }

        // Subtract baseline when computing peak amplitude and integrals
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
            std::cout << "Max ADC: " << adcMax << ", satThreshold: " << satThreshold << "\n";
            saturatedLinearExtrapolation(
                satStart,
                satEnd,
                wf.data(),
                N,
                baseline,
                peakSample,
                peakAmpFit
            );

            std::cout << "Saturated peak detected at index " << maxIdx << " with sat range [" << satStart << ", " << satEnd << "]\n";
            std::cout << "Peak sample is " << peakSample << " samples, fitted amplitude = " << peakAmpFit << "\n";
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
        float target = fraction * peakAmpFit + baseline; // careful: function findxPercentofMax expects raw samples
        // We'll implement a local leading-edge crossing search to be robust with baseline-subtraction.
        int leadIdx = maxIdx;
        while (leadIdx > 0 && wf[leadIdx] > target) --leadIdx;
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

        // Walk left from maxIdx until waveform falls <= dynamic threshold (or touches baseline) to find left crossing.
        int leftCross = maxIdx;
        while (leftCross > 0 && wf[leftCross] > dynThr) --leftCross;
        // Walk right from maxIdx until waveform falls <= dynamic threshold
        int rightCross = maxIdx;
        while (rightCross + 1 < N && wf[rightCross] > dynThr) ++rightCross;

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

        // Advance i to after the rightCross to avoid merging
        i = rightCross + 1;
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

    // If this exact waveform is found, pause for user input [  -8.    8.  708. 3720. 3712. 3700. 3692. 3696.  939. -256. -256.]
    std::vector<float> exampleWaveform = {  -8.0f,    8.0f,  708.0f, 3720.0f, 3712.0f, 3700.0f, 3692.0f, 3696.0f,  939.0f, -256.0f, -256.0f };
    bool match = false;
    if (N == (int)exampleWaveform.size())
    {
        match = true;
        for (int t = 0; t < N; ++t)
        {
            if (std::abs(wf[t] - exampleWaveform[t]) > 1e-3f)
            {
                match = false;
                break;
            }
        }
    }


    if (match)
    {
        std::cout << "Waveform found:" << std::endl;

        // Cout startIdx and endIdx samples
        std::cout << "Saturated region from index " << satStartIdx << " to " << satEndIdx << "\n";
        // Cout all the data points from startIdx - 1 to endIdx + 1 for debugging
        std::cout << "All samples:\n";
        for (int t = 0; t<N; ++t) {
            if (t >= 0 && t < N) {
                std::cout << "  Sample " << t << ": " << wf[t] << "\n";
            }
        }

        // Cout everything and then pause until user hits enter
        std::cout << "Saturated peak extrapolation:\n";
        std::cout << "  Left slope: " << leftSlope << "\n";
        std::cout << "  Right slope: " << rightSlope << "\n";
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

    if (match)
    {
        std::cout << "Saturated peak extrapolation:\n";
        std::cout << "Saturated region from index " << satStartIdx << " to " << satEndIdx << "\n";
        std::cout << "  Left slope: " << leftSlope << ", intercept: " << leftIntercept << "\n";
        std::cout << "  Right slope: " << rightSlope << ", intercept: " << rightIntercept << "\n";
        std::cout << "  Estimated peak sample: " << peakSample << "\n";
        std::cout << "  Estimated peak amplitude: " << peakAmpFit << "\n";

        std::cin.get();  // wait for user input
    }



}


int WaveformAnalyzer::run() {
    TFile fout(outputFileName.c_str(), "RECREATE");
    fout.Close();
    computePedestals();  // If no pedestals, this will assume zero-suppressed data and subtract 256
    analyzeWaveforms();
    return 0;
}
