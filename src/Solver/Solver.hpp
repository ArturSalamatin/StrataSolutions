#pragma once

#include <omp.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

#include "../Grid/Defines.h"
#include "../Grid/UniformGrid2D.h"
#include "../Properties/Coefficients.hpp"

namespace EqSolver
{
    struct State2D
    {
        static State2D FillWithZeros(const Grid::UniformGrid2D &grid)
        {
            Eigen::MatrixX<float_t> cur_state{grid.X_nodes.size(), grid.Y_nodes.size()};
            cur_state.fill(0.0);

            return {cur_state};
        }

        State2D(const State2D &) noexcept = default;
        State2D(State2D &&) noexcept = default;

    protected:
        Eigen::MatrixX<float_t> cur_state;

    protected:
        State2D(const Eigen::MatrixX<float_t> &cur_state)
            : cur_state{cur_state}
        {
        }
    };

    namespace SplittingMethod
    {
        using SpMatrix = Eigen::SparseMatrix<float_t>;
        using VectSpMatrix = std::vector<Eigen::SparseMatrix<float_t>>;

        struct BaseSplit
        {
            BaseSplit(
                std::shared_ptr<Properties::Fields> properties,
                size_t serial_nodes,
                size_t matrix_size)
                : properties{properties},
                  LaplaceTerm(
                      serial_nodes,                // number of matricies
                      Eigen::SparseMatrix<float_t>{// ctor per matrix
                                                   ptrdiff_t(matrix_size),
                                                   ptrdiff_t(matrix_size)})
            {
                for (auto &m : LaplaceTerm)
                {
                    m.reserve(matrix_size * 3ull - 2ull);
                }
            }

        public:
            VectSpMatrix LaplaceTerm;
            std::shared_ptr<Properties::Fields> properties;

            using Conductivity_f = Properties::Fields::Conductivity_f;
        };

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

#pragma omp parallel for
                for (ptrdiff_t m_id = 0; m_id < ptrdiff_t(LaplaceTerm.size()); ++m_id)
                {
                    auto &matrix = LaplaceTerm[m_id];
                    std::vector<Eigen::Triplet<float_t, ptrdiff_t>> tripletList;
                    tripletList.reserve(matrix.rows() * 3ull - 2ull);

                    tripletList.emplace_back(0, 0, conductivity_y_bounds(m_id, 0) / gr_y.step(0));
                    tripletList.emplace_back(0, 1, conductivity_y_bounds(m_id, 0) / gr_y.step(0));
                    for (ptrdiff_t row = 1; row < ptrdiff_t(matrix.rows()) - 1; ++row)
                    {
                        tripletList.emplace_back(row, row - 1,
                                                 conductivity_y_bounds(m_id, row) / gr_y.step(row - 1));
                        tripletList.emplace_back(row, row,
                                                 conductivity_y_bounds(m_id, row) / gr_y.step(row - 1) +
                                                     conductivity_y_bounds(m_id, row + 1) / gr_y.step(row));
                        tripletList.emplace_back(row, row + 1,
                                                 conductivity_y_bounds(m_id, row + 1) / gr_y.step(row));
                    }
                    ptrdiff_t end = ptrdiff_t(matrix.rows()) - 1;
                    tripletList.emplace_back(end, end - 1,
                                             conductivity_y_bounds(m_id, end + 1) / gr_y.step(end));
                    tripletList.emplace_back(end, end,
                                             conductivity_y_bounds(m_id, end + 1) / gr_y.step(end));

                    matrix.setFromTriplets(tripletList.begin(), tripletList.end());
                }
            }

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
                        conductivity_nodes(i, end-1);
                }
                return conductivity_y_bounds;
            }
        };

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

#pragma omp parallel for
                for (ptrdiff_t m_id = 0; m_id < ptrdiff_t(LaplaceTerm.size()); ++m_id)
                {
                    auto &matrix = LaplaceTerm[m_id];
                    std::vector<Eigen::Triplet<float_t, ptrdiff_t>> tripletList;
                    tripletList.reserve(matrix.rows() * 3ull - 2ull);

                    tripletList.emplace_back(0, 0,
                                             conductivity_x_bounds(0, m_id) / gr_x.step(0));
                    tripletList.emplace_back(0, 1,
                                             conductivity_x_bounds(1, m_id) / gr_x.step(0));
                    for (ptrdiff_t row = 1; row < ptrdiff_t(matrix.rows()) - 1; ++row)
                    {
                        tripletList.emplace_back(row, row - 1,
                                                 conductivity_x_bounds(row, m_id) / gr_x.step(row - 1));
                        tripletList.emplace_back(row, row,
                                                 conductivity_x_bounds(row, m_id) / gr_x.step(row - 1) +
                                                     conductivity_x_bounds(row + 1, m_id) / gr_x.step(row));
                        tripletList.emplace_back(row, row + 1,
                                                 conductivity_x_bounds(row + 1, m_id) / gr_x.step(row));
                    }
                    ptrdiff_t end = ptrdiff_t(matrix.rows()) - 1;
                    tripletList.emplace_back(end, end - 1,
                                             conductivity_x_bounds(end, m_id) / gr_x.step(end));
                    tripletList.emplace_back(end, end,
                                             conductivity_x_bounds(end + 1, m_id) / gr_x.step(end));

                    matrix.setFromTriplets(tripletList.begin(), tripletList.end());
                }
            }

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
                        conductivity_nodes(end-1, j);
                }

                return conductivity_x_bounds;
            }
        };

    } // Splitting method
} // EqSolver