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
            template <typename Functor>
            BCSouthNorth(
                const BCSouth &south,
                const BCNorth &north,
                const Grid::UniformGrid1D &grid,
                const Functor &functor)
                : south{south}, north{north},
                  grid{grid},
                  south_vals(grid.size(), 0.0),
                  north_vals(grid.size(), 0.0)
            {
                for (size_t i = 0; i < grid.size(); ++i)
                {
                    south_vals[i] =
                        functor(south.fixed_x, grid[i]);
                    north_vals[i] =
                        functor(north.fixed_x, grid[i]);
                }
            }

            BCSouth south;
            BCNorth north;
            Grid::UniformGrid1D grid;

            std::vector<float_t> south_vals, north_vals;
        };

        struct BCEastWest
        {
            template <
                typename Functor>
            BCEastWest(const BCEast &east,
                       const BCWest &west,
                       const Grid::UniformGrid1D &grid,
                       const Functor &functor)
                : east{east}, west{west}, grid{grid}, east_vals(grid.size()), west_vals(grid.size())
            {
                for (size_t i = 0; i < grid.size(); ++i)
                {
                    east_vals[i] = functor(grid[i], east.fixed_y);
                    west_vals[i] = functor(grid[i], west.fixed_y);
                }
            }

            BCEast east;
            BCWest west;
            Grid::UniformGrid1D grid;

            std::vector<float_t> east_vals, west_vals;
        };

        struct BoundaryConditions
        {
            template <
                typename Functor>
            BoundaryConditions(
                const Grid::UniformGrid2D &grid,
                const Functor &functor)
                : south_north{
                      BCSouth{grid.X_nodes.front()},
                      BCNorth{grid.X_nodes.back()},
                      grid.Y_nodes,
                      functor},
                  east_west{BCEast{grid.Y_nodes.back()}, BCWest{grid.Y_nodes.front()}, grid.X_nodes, functor}
            {
            }

            BoundaryConditions(const BoundaryConditions &) = default;

            BCEastWest east_west;
            BCSouthNorth south_north;
        };
    } // ConcreteProblem
} // EqSolver