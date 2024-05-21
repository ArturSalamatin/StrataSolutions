#pragma once

#include <vector>
#include <memory>

#include <omp.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

#include "../../Grid/Defines.h"
#include "../../Properties/Coefficients.hpp"

#include "BaseSplit.hpp"

namespace EqSolver
{
    namespace SplittingMethod
    {

        struct SplitY
            : public BaseSplit
        {
            template <typename Grid_t>
            SplitY(
                std::shared_ptr<Properties::Fields> properties,
                const Grid_t &grid)
                : BaseSplit{properties, grid.Y_nodes.size(), grid.X_nodes.size()}
            {
                FillLaplaceTerm(grid);
            }

        protected:
            template <typename Grid_t>
            void FillLaplaceTerm(const Grid_t &grid)
            {
                Conductivity_f conductivity_x_bounds{
                    set_conductivity_x_bounds(properties->conductivity_f, grid)};

                const auto &gr_x = grid.X_nodes;
                const auto &y_vol = grid.Y_nodes.cell_volume();

#pragma omp parallel for
                for (ptrdiff_t m_id = 0; m_id < ptrdiff_t(LaplaceTerm.size()); ++m_id)
                {
                    auto &matrix = LaplaceTerm[m_id];
                    std::vector<Eigen::Triplet<float_t, ptrdiff_t>> tripletList;
                    tripletList.reserve(matrix.rows() * 3ull - 2ull);

                    tripletList.emplace_back(
                        0, 0,
                        conductivity_x_bounds(0, m_id) / gr_x.step(0) * y_vol(0));
                    tripletList.emplace_back(
                        0, 1,
                        conductivity_x_bounds(1, m_id) / gr_x.step(0) * y_vol(0));
                    for (ptrdiff_t row = 1; row < ptrdiff_t(matrix.rows()) - 1; ++row)
                    {
                        tripletList.emplace_back(
                            row, row - 1,
                            conductivity_x_bounds(row, m_id) / gr_x.step(row - 1) * y_vol(row));
                        tripletList.emplace_back(
                            row, row,
                            (conductivity_x_bounds(row, m_id) / gr_x.step(row - 1) +
                             conductivity_x_bounds(row + 1, m_id) / gr_x.step(row)) *
                                y_vol(row));
                        tripletList.emplace_back(
                            row, row + 1,
                            conductivity_x_bounds(row + 1, m_id) / gr_x.step(row) * y_vol(row));
                    }
                    ptrdiff_t end = ptrdiff_t(matrix.rows()) - 1;
                    tripletList.emplace_back(
                        end, end - 1,
                        conductivity_x_bounds(end, m_id) / gr_x.step(end) * y_vol(end));
                    tripletList.emplace_back(
                        end, end,
                        conductivity_x_bounds(end + 1, m_id) / gr_x.step(end) * y_vol(end));

                    matrix.setFromTriplets(tripletList.begin(), tripletList.end());
                }
            }

            /**
             * @brief Average conductivity in y-direction, to get values on y-boundaries
             *
             * @tparam Grid_t
             * @param conductivity_nodes
             * @param grid
             * @return auto
             */
            template <typename Grid_t>
            Conductivity_f set_conductivity_x_bounds(
                const Conductivity_f &conductivity_nodes,
                const Grid_t &grid) const
            {
                Conductivity_f conductivity_x_bounds{
                    grid.X_nodes.size() + 1,
                    grid.Y_nodes.size()};

#pragma omp parallel for
                for (ptrdiff_t j = 0; j < ptrdiff_t(conductivity_x_bounds.cols()); ++j)
                {
                    conductivity_x_bounds(0, j) =
                        conductivity_nodes(0, j);
                    for (ptrdiff_t i = 1; i < ptrdiff_t(conductivity_x_bounds.rows()) - 1; ++i)
                    {
                        conductivity_x_bounds(i, j) =
                            (conductivity_nodes(i - 1, j) + conductivity_nodes(i, j)) / 2.0;
                    }
                    ptrdiff_t end = conductivity_x_bounds.rows() - 1;
                    conductivity_x_bounds(end, j) =
                        conductivity_nodes(end - 1, j);
                }

                return conductivity_x_bounds;
            }
        };

    } // Splitting method
} // EqSolver