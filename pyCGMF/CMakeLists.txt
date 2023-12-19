pybind11_add_module(pyCGMF pyCGMF.cpp)

target_link_libraries(pyCGMF PRIVATE cgmf nlohmann_json::nlohmann_json osiris_lib)
target_compile_definitions(pyCGMF PRIVATE VERSION_INFO=${PYCGMF_VERSION_INFO})
