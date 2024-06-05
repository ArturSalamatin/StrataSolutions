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

    namespace BoundaryConditions
    {
        struct BoundaryCondition
        {
            enum BCType
            {
                first,
                second,
                third
            };
            BoundaryCondition(BCType type) : type{type} {}

            BCType type;
        };

        struct BCSouth : public BoundaryCondition
        {
            BCSouth(float_t fixed_x, BCType type = BCType::first)
                : BoundaryCondition{type}, fixed_x{fixed_x}
            {
            }

            float_t fixed_x;
        };

        struct BCNorth : public BoundaryCondition
        {
            BCNorth(float_t fixed_x, BCType type = BCType::first)
                : BoundaryCondition{type}, fixed_x{fixed_x}
            {
            }

            float_t fixed_x;
        };

        struct BCEast : public BoundaryCondition
        {
            BCEast(float_t fixed_y, BCType type = BCType::first)
                : BoundaryCondition{type}, fixed_y{fixed_y}
            {
            }

            float_t fixed_y;
        };

        struct BCWest : public BoundaryCondition
        {
            BCWest(float_t fixed_y, BCType type = BCType::first)
                : BoundaryCondition{type}, fixed_y{fixed_y}
            {
            }

            float_t fixed_y;
        };

        struct BCFunctorBase
        {
            virtual float_t operator()(float_t x, float_t y, float_t t) const = 0;
        };
    } // BoundaryConditions
} // EqSolver