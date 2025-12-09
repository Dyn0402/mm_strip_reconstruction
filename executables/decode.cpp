//
// Created by dn277127 on 2025-12-02.
//

#include "DreamDecoder.h"
#include <iostream>

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input.fdf> [output.root]" << std::endl;
        return 1;
    }
    std::string in = argv[1];
    std::string out = "ftest.root";
    if (argc >= 3) out = argv[2];

    try {
        DreamDecoder dec(in, out);
        int n = dec.run();
        std::cout << "Decoder finished, events: " << n << std::endl;
    } catch (const std::exception &e) {
        std::cerr << "Decoder error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

