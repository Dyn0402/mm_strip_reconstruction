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
        std::cout << "No pedestal file provided — skipping pedestal calculation.\n";
        return;
    }

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


void WaveformAnalyzer::loadPedestals() {
    // pedestalMap already filled by computePedestals()
    if (pedestalMap.empty()) {
        std::cout << "Warning: No pedestal data loaded.\n";
    }
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



// simple 3-point parabolic interpolation
// maxIndex is the bin of the max sample
void WaveformAnalyzer::fitParabola(const std::vector<float>& amps,
                                   int maxIndex,
                                   float& peakAmp, float& peakOffset) const
{
    if (maxIndex <= 0 || maxIndex >= (int)amps.size() - 1) {
        peakAmp = amps[maxIndex];
        peakOffset = 0.0;
        return;
    }

    float y1 = amps[maxIndex - 1];
    float y2 = amps[maxIndex];
    float y3 = amps[maxIndex + 1];

    float denom = (y1 - 2*y2 + y3);
    if (std::abs(denom) < 1e-9) {
        peakAmp = y2;
        peakOffset = 0;
        return;
    }

    peakOffset = 0.5f * (y1 - y3) / denom; // in units of sample spacing
    peakAmp = y2 - 0.25f * (y1 - y3) * peakOffset;
}


// Find the time (in sample units) at which the waveform rises above a fraction of the max
float WaveformAnalyzer::findxPercentofMax(const std::vector<float>& amps,
                                          int maxIndex,
                                          float fraction) const
{
    if (amps.empty() || maxIndex <= 0 || maxIndex >= (int)amps.size())
        return maxIndex; // fallback

    float maxVal = amps[maxIndex];
    float target = fraction * maxVal;

    // Walk backwards from the peak until waveform drops below target
    int i = maxIndex;
    while (i > 0 && amps[i] > target) {
        i--;
    }

    // i is now the last sample below target, i+1 is above target
    float y1 = amps[i];
    float y2 = amps[i + 1];

    // Linear interpolation: t_cross = i + (target - y1)/(y2 - y1)
    if (std::abs(y2 - y1) < 1e-9f)
        return i + 0.5f; // avoid divide by zero, take mid-sample

    float i_cross = i + (target - y1) / (y2 - y1);
    return i_cross;
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
    UShort_t  out_channel;
    Float_t   out_amp;
    Float_t   out_time_ns;
    Float_t   out_sample;

    hitTree.Branch("eventId", &out_eventID, "eventId/l");
    hitTree.Branch("channel", &out_channel, "channel/s");
    hitTree.Branch("amplitude", &out_amp, "amplitude/F");
    hitTree.Branch("time", &out_time_ns, "time/F");
    hitTree.Branch("sample", &out_sample, "sample/F");

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

        // analyze per channel
        for (auto& kv : waves) {
            int ch = kv.first;
            std::vector<float>& amps = kv.second;

            float noiseRMS = pedestalMap.count(ch) ? pedestalMap[ch].rms : 3.0f;
            auto peakIndices = findPeakIndices(amps, noiseRMS);
            if (!peakIndices.empty()) {
                peakIndices = mergePeaks(peakIndices, amps);
            }

            // Fallback: no peaks found → continue
            if (peakIndices.empty()) continue;

            // If multiple peaks not allowed, peakIndices will have size=1.
            // Otherwise we loop over them.
            for (int idx : peakIndices) {

                int sampleOfMax = samplesByCh[ch][idx];
                float peakAmp, peakSample;
                if (timingMethod == "parabola") {
                    float peakOffset;
                    fitParabola(amps, idx, peakAmp, peakOffset);
                    peakSample = sampleOfMax + peakOffset;
                }
                else if (timingMethod == "percent_max")
                {
                    peakAmp = amps[idx];
                    peakSample = findxPercentofMax(amps, idx, timingPercentMax);
                } else {
                    std::cerr << "Unknown timing method: " << timingMethod << "\n";
                    peakAmp = amps[idx];
                    peakSample = sampleOfMax;
                }

                float hitTime_ns =
                        timePerTimestamp * timestamp +
                        timePerFtst * ftst +
                        timePerSample * peakSample;

                out_eventID = eventID;
                out_channel = ch;
                out_amp     = peakAmp;
                out_time_ns = hitTime_ns;
                out_sample  = peakSample;

                hitTree.Fill();

                if (!allowMultiplePeaks)
                    break;
            }
        }
    }

    hitTree.Write();
    fout.Close();
}

std::vector<int> WaveformAnalyzer::findPeakIndices(const std::vector<float>& wf,
                                                   float noiseRMS) const
{
    std::vector<int> peaks;
    if (wf.size() < 3) return peaks;

    // dynamic threshold
    float thr = thresholdSigma * noiseRMS;

    int bestIdx = -1;
    float bestVal = thr;

    int i = 1; // start from second sample
    bool peakActive = false; // flag to hold until waveform drops below threshold

    while (i < (int)wf.size() - 1) {
        float y0 = wf[i-1];
        float y1 = wf[i];
        float y2 = wf[i+1];

        if (y1 > y0 && y1 >= y2 && y1 > thr) {
            if (!peakActive) {
                // find plateau end
                int start = i;
                int end = i;
                while (end + 1 < (int)wf.size() && wf[end + 1] == y1) end++;

                int peakIdx = (start + end) / 2; // middle of plateau

                if (allowMultiplePeaks) {
                    peaks.push_back(peakIdx);
                } else {
                    if (y1 > bestVal) {
                        bestVal = y1;
                        bestIdx = peakIdx;
                    }
                }

                peakActive = true;    // now hold
                i = end + 1;          // skip plateau
            } else {
                i++; // peak still active, skip
            }
        } else {
            if (y1 < thr) peakActive = false; // waveform dropped below threshold
            i++;
        }
    }

    if (!allowMultiplePeaks && bestIdx != -1) {
        peaks.push_back(bestIdx);
    }

    return peaks;
}




std::vector<int> WaveformAnalyzer::mergePeaks(const std::vector<int>& peaks,
                                              const std::vector<float>& wf) const
{
    std::vector<int> result;
    if (peaks.empty()) return result;

    // copy and ensure peaks are sorted ascending
    std::vector<int> p = peaks;
    std::sort(p.begin(), p.end());

    auto valueAt = [&](int idx) -> float {
        if (idx < 0 || idx >= (int)wf.size()) return -std::numeric_limits<float>::infinity();
        return wf[idx];
    };

    // Start first group with the first peak
    int bestIdx = p[0];
    float bestVal = valueAt(bestIdx);

    for (size_t i = 1; i < p.size(); ++i) {
        // if within merge distance from previous peak, treat as same group
        if (p[i] - p[i - 1] <= peakMergeDistance) {
            float val = valueAt(p[i]);
            if (val > bestVal) {
                bestVal = val;
                bestIdx = p[i];
            }
        } else {
            // group ended, keep best of that group
            result.push_back(bestIdx);
            bestIdx = p[i];
            bestVal = valueAt(bestIdx);
        }
    }

    // push last group's best
    result.push_back(bestIdx);
    return result;
}


int WaveformAnalyzer::run() {
    TFile fout(outputFileName.c_str(), "RECREATE");
    fout.Close();
    if (hasPedestal) computePedestals();
    analyzeWaveforms();
    return 0;
}
