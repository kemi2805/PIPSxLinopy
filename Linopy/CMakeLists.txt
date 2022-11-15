############## callbackExample ##############
add_executable(pipsipmLinopyCallback LinopyCallback.cpp)

add_library(Linopy
            Linopy_Source/LinopyModel.cpp
           )

target_include_directories(pipsipmLinopyCallback PRIVATE ${includes_for_interface})
target_compile_options(pipsipmLinopyCallback PRIVATE ${OpenMP_CXX_FLAGS})
target_link_libraries(pipsipmLinopyCallback
        PRIVATE pips-ipmpp
        PRIVATE OpenMP::OpenMP_CXX
        PRIVATE MPI::MPI_CXX
        PRIVATE Linopy
        )

           
