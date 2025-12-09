//
// Created by dn277127 on 2025-12-02.
//

#ifndef DREAM_DECODER_H
#define DREAM_DECODER_H

#include <string>

class DreamDecoder {
public:
    // Construct with input filename and optional output ROOT filename
    // Throws std::runtime_error if file can't be opened or is empty.
    DreamDecoder(const std::string &input_filename, const std::string &output_root = "ftest.root");

    // Run the decode. Returns number of events processed.
    // Throws std::runtime_error on fatal IO errors.
    int run();

private:
    // Note: read16 replicates legacy behavior: it reads 2 bytes, ntohs them,
    // and returns is.eof() (true if EOF). This matches your old read16.
    bool read16(std::ifstream &is, uint16_t &data);

    void print_data(uint16_t data);

    std::string input_filename_;
    std::string output_root_;
};

#endif // DREAM_DECODER_H

