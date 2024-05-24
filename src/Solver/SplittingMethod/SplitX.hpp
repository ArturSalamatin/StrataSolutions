#pragma once

#include <vector>
#include <memory>

#include <omp.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

#include "../../Grid/Defines.h"
#include "../../Properties/Coefficients.hpp"
#include "../State2D.hpp"

#include "BaseSplit.hpp"

namespace EqSolver
{
    namespace SplittingMethod
    {
        struct SplitX
            : public BaseSplit
        {
            template <typename Grid_t>
            SplitX(
                std::shared_ptr<Properties::Fields> properties,
                const Grid_t &grid)
                : BaseSplit{
                      properties,
                      grid.X_nodes.size(),
                      grid.Y_nodes.size()}
            {
                FillLaplaceTerm(grid);
            }

        protected:
            template <typename Grid_t>
            void FillLaplaceTerm(const Grid_t &grid)
            {
                Conductivity_f conductivity_y_bounds{
                    set_conductivity_y_bounds(properties->conductivity_f, grid)};

                const auto &gr_y = grid.Y_nodes;
                const auto &x_vol = grid.X_nodes.cell_volume();

#pragma omp parallel for
                for (ptrdiff_t m_id = 0; m_id < ptrdiff_t(its_LaplaceTerm.size()); ++m_id)
                {
                    auto &matrix = its_LaplaceTerm[m_id];
                    std::vector<Eigen::Triplet<float_t, ptrdiff_t>> tripletList;
                    tripletList.reserve(matrix.rows() * 3ull - 2ull);

                    tripletList.emplace_back(
                        0, 0,
                        conductivity_y_bounds(m_id, 0) / gr_y.step(0) * x_vol(0));
                    tripletList.emplace_back(
                        0, 1,
                        conductivity_y_bounds(m_id, 0) / gr_y.step(0) * x_vol(0));
                    for (ptrdiff_t row = 1; row < ptrdiff_t(matrix.rows()) - 1; ++row)
                    {
                        tripletList.emplace_back(
                            row, row - 1,
                            conductivity_y_bounds(m_id, row) / gr_y.step(row - 1) * x_vol(row));
                        tripletList.emplace_back(
                            row, row,
                            (conductivity_y_bounds(m_id, row) / gr_y.step(row - 1) +
                             conductivity_y_bounds(m_id, row + 1) / gr_y.step(row)) *
                                x_vol(row));
                        tripletList.emplace_back(
                            row, row + 1,
                            conductivity_y_bounds(m_id, row + 1) / gr_y.step(row) * x_vol(row));
                    }
                    ptrdiff_t end = ptrdiff_t(matrix.rows()) - 1;
                    tripletList.emplace_back(
                        end, end - 1,
                        conductivity_y_bounds(m_id, end + 1) / gr_y.step(end) * x_vol(end));
                    tripletList.emplace_back(
                        end, end,
                        conductivity_y_bounds(m_id, end + 1) / gr_y.step(end) * x_vol(end));

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
            auto set_conductivity_y_bounds(
                const Conductivity_f &conductivity_nodes,
                const Grid_t &grid) const
            {
                Conductivity_f conductivity_y_bounds{
                    grid.X_nodes.size(),
                    grid.Y_nodes.size() + 1};

#pragma omp parallel for
                for (ptrdiff_t i = 0; i < conductivity_y_bounds.rows(); ++i)
                {
                    conductivity_y_bounds(i, 0) =
                        conductivity_nodes(i, 0);
                    for (ptrdiff_t j = 1; j < conductivity_y_bounds.cols() - 1; ++j)
                    {
                        conductivity_y_bounds(i, j) =
                            (conductivity_nodes(i, j - 1) + conductivity_nodes(i, j)) / 2.0;
                    }
                    ptrdiff_t end = conductivity_y_bounds.cols() - 1;
                    conductivity_y_bounds(i, end) =
                        conductivity_nodes(i, end - 1);
                }
                return conductivity_y_bounds;
            }
        };
    } // SplittingMethod
} // EqSolver