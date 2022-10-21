#ifndef __GKDN__
#define __GKDN__

#ifdef __cplusplus
#include <string>
using std::string;
#endif

#include "nlohmann/json.hpp"
using nlohmann::json;

#include "potential/params.hpp"
constexpr auto n = omplib::Proj::neutron;

using OMPFile = omplib::OMParams<n>;

#endif
