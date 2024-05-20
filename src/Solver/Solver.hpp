#pragma once

#include <omp.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

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
                std::shared_ptr<Fields> properties,
                size_t serial_nodes,
                size_t matrix_size)
                : properties{properties},
                  LaplaceTerm(
                      serial_nodes, // number of matricies
                      Eigen::SparseMatrix<float_t>{ // ctor per matrix
                          matrix_size,
                          matrix_size}),
                  grid{grid}
            {
                for (auto &m : LaplaceTerm)
                {
                    m.reserve(matrix_size * 3ull - 2ull);
                }
            }

        protected:
            VectSpMatrix LaplaceTerm;
            std::shared_ptr<Fields> properties;
        };

        struct SplitX
            : public BaseSplit
        {
            SplitX(
                std::shared_ptr<Fields> properties,
                const Grid::UniformGrid2D &grid)
                : BaseSplit{properties, grid.X_nodes.size(), grid.Y_nodes.size()}
            {

                FillLaplaceTerm();
            }

        protected:
            void FillLaplaceTerm()
            {
                Eigen::Matrix<float_t> conductivity_y_bounds{
                    set_conductivity_y_bounds(properties->conductivity_f)};
            }

            auto set_conductivity_y_bounds(const auto &conductivity_nodes) const
            {
                Eigen::Matrix<float_t> conductivity_y_bounds{
                    grid.X_nodes.size(),
                    grid.Y_nodes.size() + 1};

#pragma omp parallel for
                for (size_t i = 0; i < conductivity_y_bounds.rows(); ++i)
                {
                    conductivity_y_bounds(i, 0) =
                        conductivity_nodes(i, 0);
                    for (size_t j = 1; j < conductivity_y_bounds.cols() - 1; ++j)
                    {
                        conductivity_y_bounds(i, j) =
                            (conductivity_nodes(i, j - 1) + conductivity_nodes(i, j)) / 2.0;
                    }
                    conductivity_y_bounds(i, end) =
                        conductivity_nodes(i, end);
                }
                return conductivity_y_bounds;
            }
        };

        struct SplitY
            : public BaseSplit
        {
            SplitY(
                std::shared_ptr<Fields> properties,
                const Grid::UniformGrid2D &grid)
                : properties{properties},
                  LaplaceTerm(
                      grid.X_nodes.size(),
                      Eigen::SparseMatrix<float_t>{
                          grid.X_nodes.size(),
                          grid.X_nodes.size()}),
                  grid{grid}
            {
                FillLaplaceTerm();
            }

        protected:

            void FillLaplaceTerm()
            {
                Eigen::Matrix<float_t> conductivity_x_bounds{
                    set_conductivity_x_bounds(properties->conductivity_f)};
            }

            auto set_conductivity_x_bounds(const auto &conductivity_nodes) const
            {
                Eigen::Matrix<float_t> conductivity_x_bounds{
                    grid.X_nodes.size() + 1,
                    grid.Y_nodes.size()};

#pragma omp parallel for
                for (size_t j = 0; j < conductivity_x_bounds.cols(); ++j)
                {
                    conductivity_x_bounds(0, j) =
                        conductivity_nodes(0, j);
                    for (size_t i = 1; i < conductivity_x_bounds.rows() - 1; ++i)
                    {
                        conductivity_x_bounds(i, j) =
                            (conductivity_nodes(i - 1, j) + conductivity_nodes(i, j)) / 2.0;
                    }
                    conductivity_x_bounds(end, j) =
                        conductivity_nodes(end, j);
                }

                return conductivity_x_bounds;
            }
        };

    } // EqSolver