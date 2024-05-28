#pragma once

#include "../Grid/Defines.h"
#include "../Grid/UniformGrid1D.h"

namespace EqSolver
{
    namespace Problem
    {
        template <typename Grid_t>
        struct BCSouthNorth
        {
            BCSouthNorth(const BCSouth &south, const BCNorth &north, const Grid_t gtid)
                : south{south}, north{north}, grid{grid}
            {
            }

            BCSouth south;
            BCNorth north;
            Grid_t grid;
        };

        template <typename Grid_t>
        struct BCEastWest
        {
            BCEastWest(const BCEast &east, const BCWest &west, const Grid_t gtid)
                : east{east}, west{west}, grid{grid}
            {
            }

            BCEast east;
            BCWest west;
            Grid_t grid;
        };

        template <typename Grid_1D_t, typename Functor>
        struct BoundaryConditions
        {
            template <typename Grid_2D_t>
            BoundaryConditions(
                const Grid_2D_t &grid,
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

            BCSouthNorth<Grid_1D_t> south_north;
            BCEastWest<Grid_1D_t> east_west;
            Functor functor;
        };
    } // ConcreteProblem
} // EqSolver