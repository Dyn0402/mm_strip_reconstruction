//
// Created by dn277127 on 2025-12-05.
//

#include <iostream>
// #include "../waveform_analysis/include/WaveformAnalyzer.h"
#include "WaveformAnalyzer.h"

int main(int argc, char **argv) {
//    WaveformAnalyzer wf("decoded.root", "hits.root", "pedestal.root");
    // WaveformAnalyzer wf("/local/home/dn277127/x17/decoder_test/ftest.root", "hits.root", "/local/home/dn277127/x17/dream_run/ped_thresh_1_12_25_18_30/Mx17_ped_pedthr_251201_18H27_000_05.root");
    // WaveformAnalyzer wf("Mx17_run_datrun_251204_18H17_000_05.root", "hits.root", "/home/dylan/CLionProjects/mm_strip_reconstruction/test/ped_thresh_1_12_25_18_30/Mx17_ped_pedthr_251201_18H27_000_05.root");
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " <input.root> [output.root] [pedestal.root]" << std::endl;
        return 1;
    }
    std::string inputFile = argv[1];
    std::string outputFile = "hits.root";
    std::string pedestalFile = "";
    if (argc >= 3) outputFile = argv[2];
    if (argc >= 4) pedestalFile = argv[3];

    WaveformAnalyzer wf(inputFile, outputFile, pedestalFile);
    wf.run();
    return 0;
}
