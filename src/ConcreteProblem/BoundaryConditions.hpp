#pragma once

#include "../Grid/Defines.h"
#include "../Grid/UniformGrid1D.h"
#include "../Grid/UniformGrid2D.h"

namespace EqSolver
{
    namespace Problem
    {
        struct BCSouthNorth
        {
            BCSouthNorth(
                const BCSouth &south, 
                const BCNorth &north,
                 const Grid::UniformGrid1D& grid)
                : south{south}, north{north}, 
                grid{grid}
            {
            }

            BCSouth south;
            BCNorth north;
            Grid::UniformGrid1D grid;
        };

        struct BCEastWest
        {
            BCEastWest(const BCEast &east, 
            const BCWest &west, 
            const Grid::UniformGrid1D& grid)
                : east{east}, west{west}, grid{grid}
            {
            }

            BCEast east;
            BCWest west;
            Grid::UniformGrid1D grid;
        };

        template <
        typename Functor>
        struct BoundaryConditions
        {
            BoundaryConditions(
                const Grid::UniformGrid2D &grid,
                const Functor& functor)
                : south_north{
                    BCSouth{grid.X_nodes.front()},
                    BCNorth{grid.X_nodes.back()},
                    grid.Y_nodes},
                  east_west{
                    BCEast{grid.Y_nodes.back()},
                    BCWest{grid.Y_nodes.front()},
                    grid.X_nodes},
                  functor{functor}
            {
            }

            BoundaryConditions(const BoundaryConditions&) = default;

            BCEastWest east_west;
            BCSouthNorth south_north;
            Functor functor;
        };
    } // ConcreteProblem
} // EqSolver