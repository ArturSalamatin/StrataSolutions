#pragma once

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

            BCType type;
        };

        struct BCSouth : public BoundaryCondition
        {
        };

        struct BCNorth : public BoundaryCondition
        {
        };

        struct BCEast : public BoundaryCondition
        {
        };

        struct BCWest : public BoundaryCondition
        {
        };

        struct BCEastWest
        {
            BCEast east;
            BCWest west;
        };
        
        struct BCNorthSouth
        {
            BCNorth north;
            BCSouth south;
        };

        struct BoundaryConditions
        {
            BoundaryConditions()
                : south{BoundaryCondition::BCType::first},
                  east{BoundaryCondition::BCType::first},
                  north{BoundaryCondition::BCType::first},
                  west{BoundaryCondition::BCType::first}
            {
            }

            BCSouth south;
            BCEast east;
            BCNorth north;
            BCWest west;
        };
    } // ConcreteProblem
} // EqSolver