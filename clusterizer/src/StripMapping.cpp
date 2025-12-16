//
// Created by dylan on 12/11/25.
//

#include "StripMapping.h"
#include <fstream>
#include <sstream>
#include <iostream>

std::vector<std::string> StripMapping::split(const std::string& s, char delim)
{
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim))
        out.push_back(item);
    return out;
}

bool StripMapping::loadFromCSV(const std::string& filename)
{
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "ERROR: Cannot open mapping file " << filename << "\n";
        return false;
    }

    std::string line;
    bool isHeader = true;

    while (std::getline(infile, line)) {
        if (line.empty()) continue;

        // Skip header row
        if (isHeader) {
            isHeader = false;
            continue;
        }

        auto cols = split(line, ',');
        if (cols.size() < 9) {
            std::cerr << "WARNING: Bad line, skipping: " << line << "\n";
            continue;
        }

        StripInfo info;
        info.connector        = std::stoi(cols[0]);
        info.connectorChannel = std::stoi(cols[1]);
        info.stripNumber      = std::stoi(cols[2]);
        info.axis             = cols[3][0];
        info.pitch            = std::stof(cols[4]);
        info.interpitch       = std::stof(cols[5]);
        info.xGerber          = std::stof(cols[7]);
        info.yGerber          = std::stof(cols[8]);

        // Parse neighbours — colon-separated list
        for (auto& n : split(cols[6], ':')) {
            if (!n.empty())
                info.neighbours.push_back(std::stoi(n));
        }

        // Compute global channel = connector*64 + connectorChannel
        int globalChannel = info.connector * 64 + info.connectorChannel;

        mapping_[globalChannel] = info;

        if (info.axis == 'x') xChannels_.push_back(globalChannel);
        else if (info.axis == 'y') yChannels_.push_back(globalChannel);
    }

    return true;
}

const StripInfo* StripMapping::getStripInfo(int channel) const
{
    auto it = mapping_.find(channel);
    return (it == mapping_.end() ? nullptr : &it->second);
}

const std::vector<int>& StripMapping::getChannelsForAxis(char axis) const
{
    return (axis == 'x' ? xChannels_ : yChannels_);
}
