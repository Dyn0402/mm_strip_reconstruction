//
// Created by dn277127 on 2025-12-02.
//

// include/ConfigManager.h
#pragma once
#include <string>
#include <yaml-cpp/yaml.h>

class ConfigManager {
public:
    explicit ConfigManager(const std::string &fname);
    int getSampleRate() const;
    std::string getOutputDir() const;

private:
    YAML::Node cfg;
};
