#pragma once


namespace EqSolver
{
    using float_t = double;
    
    struct Box
    {
        Box(float_t x_a, float_t x_b, float_t y_a, float_t y_b)
        : x_a{x_a}, x_b{x_b}, y_a{y_a}, y_b{y_b}
        {
            assert(x_a < x_b);
            assert(y_a < y_b);
        };

        float_t x_a, x_b, y_a, y_b;
    };

    
    struct Steps
    {
        Steps(float_t step_x, float_t step_y)
        : step_x{step_x}, step_y{step_y}
        {
        }

        float_t step_x, step_y;
    };
} // EqSolver