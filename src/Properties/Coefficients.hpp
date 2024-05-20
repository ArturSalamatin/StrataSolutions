#pragma once

#include <omp.h>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "../Grid/Defines.h"
#include "../Grid/UniformGrid2D.h"

namespace EqSolver
{
    namespace Properties
    {
        struct Fields
        {
            using Volume_f = Eigen::MatrixX<float_t>;
            using Capacity_f = Eigen::MatrixX<float_t>;
            using Conductivity_f = Eigen::MatrixX<float_t>;

            Volume_f volume_f;
            Capacity_f capacity_f;
            Conductivity_f conductivity_f;

            template <typename Grid_t, typename Capacity_t, typename Conductivity_t>
            Fields(const Grid_t &grid,
                   const Capacity_t &capacity,
                   const Conductivity_t &conductivity)
                : volume_f{set_volumes(grid)},
                  capacity_f{set_functor(grid, capacity)},
                  conductivity_f{set_functor(grid, conductivity)}
            {
            }

        protected:
            template <typename Grid_t>
            static auto set_volumes(const Grid_t &grid)
            {
                Volume_f volumes{
                    grid.X_nodes.size(),
                    grid.Y_nodes.size()};

#pragma omp parallel for
                for (ptrdiff_t j = 0; j < volumes.cols(); ++j)
                    for (ptrdiff_t i = 0; i < volumes.rows(); ++i)
                        volumes(i, j) =
                            grid.X_nodes.step(i) * grid.Y_nodes.step(j);

                // account for boundary cells
                for (ptrdiff_t j = 0; j < volumes.cols(); ++j)
                {
                    volumes(0, j) /= 2.0;
                    volumes(volumes.rows()-1, j) /= 2.0;
                }
                for (ptrdiff_t i = 0; i < volumes.rows(); ++i)
                {
                    volumes(i, 0) /= 2.0;
                    volumes(i, volumes.cols()-1) /= 2.0;
                }
                return volumes;
            }

            template <typename Grid_t, typename Functor>
            static auto set_functor(
                const Grid_t &grid,
                const Functor &functor)
            {
                Eigen::MatrixX<float_t> functor_nodes{
                    grid.X_nodes.size(),
                    grid.Y_nodes.size()};

#pragma omp parallel for
                for (ptrdiff_t j = 0; j < functor_nodes.cols(); ++j)
                    for (ptrdiff_t i = 0; i < functor_nodes.rows(); ++i)
                    {
                        functor_nodes(i, j) =
                            functor(
                                grid.X_nodes[i],
                                grid.Y_nodes[j]);
                    }

                return functor_nodes;
            }
        };

    } // Coefficients
} // EqSolver