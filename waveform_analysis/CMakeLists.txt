add_library(waveform_analysis STATIC
        src/WaveformAnalyzer.cpp
)

target_include_directories(common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(waveform_analysis
        PUBLIC common
        PUBLIC ${ROOT_LIBRARIES}
)

# Executable
add_executable(analyze_waveforms ../executables/analyze_waveforms.cpp)

target_link_libraries(analyze_waveforms
        PRIVATE waveform_analysis
        PRIVATE ${ROOT_LIBRARIES}
)
