add_library( strip_mapping STATIC
    src/StripMapping.cpp)

target_include_directories(strip_mapping PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(strip_mapping
    PUBLIC common
    PUBLIC ${ROOT_LIBRARIES}
)

add_library(clusterizer STATIC
        src/Clusterizer1D.cpp)

target_include_directories(clusterizer PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(clusterizer
        PUBLIC common
        PUBLIC strip_mapping
        PUBLIC ${ROOT_LIBRARIES}
)

# ----------------------------------------------------------------------
# Executable
# ----------------------------------------------------------------------
add_executable(clusterize1d ../executables/clusterize1d.cpp)

# Link ROOT
target_link_libraries(clusterize1d
        PRIVATE strip_mapping
        PRIVATE clusterizer
        PRIVATE ${ROOT_LIBRARIES}
)
