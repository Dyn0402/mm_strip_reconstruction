//
// Created by dylan on 12/11/25.
//

#ifndef STRIPMAPPING_H
#define STRIPMAPPING_H

#include <string>
#include <unordered_map>
#include <vector>

struct StripInfo {
    int connector = -1;
    int connectorChannel = -1;
    int stripNumber = -1;

    char axis = 'x';      // 'x' or 'y'
    float pitch = 0.0f;
    float interpitch = 0.0f;

    std::vector<int> neighbours;

    float xGerber = 0.0f;
    float yGerber = 0.0f;
};

class StripMapping {
public:
    bool loadFromCSV(const std::string& filename);

    // Lookup by global channel index
    const StripInfo* getStripInfo(int channel) const;

    // Convenience access by axis
    const std::vector<int>& getChannelsForAxis(char axis) const;

private:
    std::unordered_map<int, StripInfo> mapping_; // key = global channel
    std::vector<int> xChannels_;
    std::vector<int> yChannels_;

    // Helper
    static std::vector<std::string> split(const std::string& s, char delim);
};

#endif
