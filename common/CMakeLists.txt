add_library(common STATIC
        src/ConfigManager.cpp
)

target_include_directories(common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(common
        PRIVATE ${ROOT_LIBRARIES}
)
