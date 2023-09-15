#ifndef PARAM_H
#define PARAM_H

#include <vector>

struct Params
{
    std::vector<std::string> parameter_names;
    std::vector<std::string> differentiable_parameters;
};

#endif