#pragma once

#include "../Grid/Defines.h"
#include "../Grid/UniformGrid1D.h"
#include "../Grid/UniformGrid2D.h"

namespace EqSolver
{
    namespace BoundaryConditions
    {
        struct BCSouthNorth
        {
            BCSouthNorth(
                const BCSouth &south,
                const BCNorth &north,
                const Grid::UniformGrid1D &grid,
                std::shared_ptr<BCFunctorBase> functor,
                float_t t = 0.0)
                : south{south}, north{north},
                  grid{grid},
                  south_vals(grid.size(), 0.0),
                  north_vals(grid.size(), 0.0),
                  functor{functor}
            {
                set_vals(t);
            }

            BCSouth south;
            BCNorth north;
            Grid::UniformGrid1D grid;

            std::vector<float_t> south_vals, north_vals;
            std::shared_ptr<BCFunctorBase> functor;

            void set_vals(float_t t)
            {
                for (size_t i = 0; i < grid.size(); ++i)
                {
                    south_vals[i] =
                        (*functor)(south.fixed_x, grid[i], t);
                    north_vals[i] =
                        (*functor)(north.fixed_x, grid[i], t);
                }
            }
        };

        struct BCEastWest
        {
            BCEastWest(const BCEast &east,
                       const BCWest &west,
                       const Grid::UniformGrid1D &grid,
                       std::shared_ptr<BCFunctorBase> functor,
                       float_t t = 0.0)
                : east{east}, west{west}, grid{grid}, east_vals(grid.size()), west_vals(grid.size()),
                  functor{functor}
            {
                set_vals(t);
            }

            BCEast east;
            BCWest west;
            Grid::UniformGrid1D grid;

            std::shared_ptr<BCFunctorBase> functor;

            std::vector<float_t> east_vals, west_vals;

            void set_vals(float_t t)
            {
                for (size_t i = 0; i < grid.size(); ++i)
                {
                    east_vals[i] = (*functor)(grid[i], east.fixed_y, t);
                    west_vals[i] = (*functor)(grid[i], west.fixed_y, t);
                }
            }
        };

        struct BoundaryConditions
        {
            BoundaryConditions(
                const Grid::UniformGrid2D &grid,
                std::shared_ptr<BCFunctorBase> functor)
                : south_north{
                      BCSouth{grid.X_nodes.front()},
                      BCNorth{grid.X_nodes.back()},
                      grid.Y_nodes,
                      functor},
                  east_west{BCEast{grid.Y_nodes.back()}, BCWest{grid.Y_nodes.front()}, grid.X_nodes, functor}
            {
            }

            BoundaryConditions(const BoundaryConditions &) = default;

            float_t east_vals(size_t i) const
            {
                return east_west.east_vals[i];
            }
            float_t west_vals(size_t i) const
            {
                return east_west.west_vals[i];
            }

            float_t south_vals(size_t i) const
            {
                return south_north.south_vals[i];
            }
            float_t north_vals(size_t i) const
            {
                return south_north.north_vals[i];
            }

            void set_vals(float_t t)
            {
                east_west.set_vals(t);
                south_north.set_vals(t);
            }

        protected:
            BCEastWest east_west;
            BCSouthNorth south_north;
        };
    } // BoundaryConditions
} // EqSolver