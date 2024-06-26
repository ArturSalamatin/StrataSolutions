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
            using Capacity_f = Eigen::ArrayXX<float_t>;
            using Source_f = Eigen::ArrayXX<float_t>;
            using Conductivity_f = Eigen::ArrayXX<float_t>;

            Capacity_f capacity_vol_f;
            Source_f source_vol_f;
            Conductivity_f conductivity_f;

            template <
                typename Grid_t,
                typename Factory_t>
            Fields(const Grid_t &grid,
                   const Factory_t &factory)
                : Fields{grid,
                         Factory_t::Capacity{},
                         Factory_t::Conductivity{},
                         Factory_t::Source{}}
            {
            }

            template <
                typename Grid_t,
                typename Capacity_t,
                typename Conductivity_t,
                typename Source_t>
            Fields(const Grid_t &grid,
                   const Capacity_t &capacity,
                   const Conductivity_t &conductivity,
                   const Source_t &source)
                : capacity_vol_f{set_functor(grid, capacity) * grid.cell_volume},
                  conductivity_f{set_functor(grid, conductivity)},
                  source_vol_f{set_functor(grid, source) * grid.cell_volume}
            {
            }

        protected:
            template <typename Grid_t, typename Functor>
            static auto set_functor(
                const Grid_t &grid,
                const Functor &functor)
            {
                Eigen::ArrayXX<float_t> functor_nodes{
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