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
            using Volume_f = Eigen::Matrix<float_t>;
            using Capacity_f = Eigen::Matrix<float_t>;
            using Conductivity_f = Eigen::Matrix<float_t>;

            Volume_f volume_f;
            Capacity_f capacity_f;
            Conductivity_f conductivity_f;

            template <typename Grid_t, typename Capacity_t, typename Conductivity_t>
            Fields(const Grid_t &grid,
                   const Capacity_t &capacity,
                   const Conductivity_t &conductivity)
                : volumes{set_volumes(grid)},
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
                for (size_t j = 0; j < volumes.cols(); ++j)
                    for (size_t i = 0; i < volumes.rows(); ++i)
                        volumes(i, j) =
                            grid.X_nodes.step(i) * grid.Y_nodes.step(j);

                // account for boundary cells
                for (size_t j = 0; j < volumes.cols(); ++j)
                {
                    volumes(0, j) /= 2.0;
                    volumes(end, j) /= 2.0;
                }
                for (size_t i = 0; i < volumes.rows(); ++i)
                {
                    volumes(i, 0) /= 2.0;
                    volumes(i, end) /= 2.0;
                }
                return volumes;
            }

            template <typename Grid_t, typename Functor>
            static auto set_functor(
                const Grid_t &grid,
                const Functor &functor)
            {
                Eigen::Matrix<float_t> functor_nodes{
                    grid.X_nodes.size(),
                    grid.Y_nodes.size()};

#pragma omp parallel for
                for (size_t j = 0; j < functor_nodes.cols(); ++j)
                    for (size_t i = 0; i < functor_nodes.rows(); ++i)
                    {
                        functor_nodes(i, j) =
                            properties.conductivity(
                                grid.X_nodes[i],
                                grid.Y_nodes[j]);
                    }

                return functor_nodes;
            }
        };

    } // Coefficients
} // EqSolver