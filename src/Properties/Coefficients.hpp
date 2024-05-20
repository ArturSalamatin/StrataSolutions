#pragma once

#include "../Grid/Defines.h"

namespace EqSolver
{
    namespace Properties
    {
        struct Capacity
        {   
            float_t operator()(float_t x, float_t y)
            {
                return x+y+1.0;
            }
        };
        
        struct Conductivity
        {   
            float_t operator()(float_t x, float_t y)
            {
                return 1.0;
            }
        };
    } // Coefficients
} // EqSolver