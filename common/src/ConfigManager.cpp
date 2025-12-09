//
// Created by dn277127 on 2025-12-02.
//

// src/ConfigManager.cpp
#include "../include/ConfigManager.h"

ConfigManager::ConfigManager(const std::string &fname) {
    cfg = YAML::LoadFile(fname);
}

int ConfigManager::getSampleRate() const {
    return cfg["sample_rate"].as<int>();
}
std::string ConfigManager::getOutputDir() const {
    return cfg["output_dir"].as<std::string>();
}
