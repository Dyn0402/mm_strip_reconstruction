add_executable(combine_feus_hits src/combine_feus_hits.cpp)

target_link_libraries(combine_feus_hits
        PUBLIC common
        PUBLIC ${ROOT_LIBRARIES}
)

target_include_directories(combine_feus_hits
        PRIVATE ${ROOT_INCLUDE_DIRS}
)
