pybind11_add_module(pyCGMF pyCGMF.cpp)

target_link_libraries(pyCGMF PRIVATE cgmf)
target_compile_options(pyCGMF PRIVATE -fsanitize=address)
