add_library(decoder STATIC
        src/DreamDecoder.cpp
)

target_include_directories(decoder PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(decoder
        PUBLIC common
        PUBLIC ${ROOT_LIBRARIES}
)

# Executable
add_executable(decode ../executables/decode.cpp)

target_link_libraries(decode
        PRIVATE decoder
        PRIVATE ${ROOT_LIBRARIES}
)
