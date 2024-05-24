#pragma once

#include "../Grid/Defines.h"

namespace EqSolver
{
    namespace Problem
    {
        struct BoundaryCondition
        {
            enum BCType
            {
                first,
                second,
                third
            };
            BoundaryCondition(BCType type) : type{type}{}

            BCType type;
        };

        struct BCSouth : public BoundaryCondition
        {
            BCSouth(BCType type, float_t fixed_x)
            : BoundaryCondition{type}, fixed_x{fixed_x}
            {}

            float_t fixed_x;
        };

        struct BCNorth : public BoundaryCondition
        {
            BCNorth(BCType type, float_t fixed_x)
            : BoundaryCondition{type}, fixed_x{fixed_x}
            {}
            
            float_t fixed_x;
        };

        struct BCEast : public BoundaryCondition
        {
            BCEast(BCType type, float_t fixed_y)
            : BoundaryCondition{type}, fixed_y{fixed_y}
            {}
            
            float_t fixed_y;
        };

        struct BCWest : public BoundaryCondition
        {
            BCWest(BCType type, float_t fixed_y)
            : BoundaryCondition{type}, fixed_y{fixed_y}
            {}
            
            float_t fixed_y;
        };

        struct BoundaryConditions
        {
            template<typename Grid_t>
            BoundaryConditions(const Grid_t& grid)
                : south{BoundaryCondition::BCType::first, grid.X_nodes.front()},
                  east{BoundaryCondition::BCType::first, grid.Y_nodes.back()},
                  north{BoundaryCondition::BCType::first, grid.X_nodes.back()},
                  west{BoundaryCondition::BCType::first, grid.Y_nodes.front()}
            {
            }

            BCSouth south;
            BCEast east;
            BCNorth north;
            BCWest west;
        };
    } // ConcreteProblem
} // EqSolver