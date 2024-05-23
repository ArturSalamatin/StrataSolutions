#pragma once

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SparseCholesky>

#include "../../Grid/Defines.h"

#include "SplittingMethod/SplitX.hpp"
#include "SplittingMethod/SplitY.hpp"

namespace EqSolver
{
    namespace SplittingMethod
    {
        template <typename Grid_t>
        struct Solver
        {
            using Map1D = Eigen::Map<Eigen::ArrayX<float_t>>;
            using Map2D = Eigen::Map<Eigen::ArrayXX<float_t>>;

            Solver(
                std::shared_ptr<Properties::Fields> properties,
                const Grid_t &grid,
                const State2D &state)
                : splitX{properties, grid}, 
                splitY{properties, grid}, factor{properties, grid}, grid{grid}, 
                state{state}
            {
            }

            void advance(float_t tau)
            {
                // tau_factor multiplies Delta_u at different time moments,
                // t and t+tau
                Eigen::ArrayXX<float_t> tau_factor{factor.DivideByTemporalStep(tau)};
                // solve a set of 1D problems in y-direction, for various x-coords
                solve_split_x(tau, tau_factor.data());
                // solve a set of 1D problems in x-direction, for various y-coords
                solve_split_y(tau, tau_factor.data());
            }

protected:
            void solve_split_x(float_t tau, float_t *tau_factor)
            {
                // nmbr of nodes in the 1D problem
                ptrdiff_t chunk_size = grid.Y_nodes.size();
                // stride in a 1D layout of 2D unknown temperature values
                ptrdiff_t stride_size = grid.X_nodes.size();
                // to be proided to Eigen::Map
                Eigen::InnerStride<> stride{stride_size};

#pragma omp parallel for
                // take every line along x-direction. A line per y-node
                for (ptrdiff_t i = 0; i < grid.X_nodes.size(); ++i)
                {
                    // memory chunk in capacity-container, corresponding to x-line
                    auto time_factor =
                        Map1D{
                            tau_factor + i,
                            chunk_size, stride};

                    // memory chunk in temperature, corresponding to x-line
                    auto data =
                        Map1D{
                            state.cur_state.data() + i,
                            chunk_size, stride};

                    // right handside of Au = b problem
                    Eigen::ArrayX<float_t> rhs = data * time_factor;

                    SpMatrix A{splitX.LaplaceTerm[i]};
                    A = A + time_factor;

                    data = solve_3Diag_problem(A, rhs);
                }
            }

            void solve_split_y(float_t tau, float_t *tau_factor)
            {
#pragma omp parallel for
                // take every line along x-direction. A line per y-node
                for (ptrdiff_t j = 0; j < grid.Y_nodes.size(); ++j)
                {
                    // memory chunk in capacity-container, corresponding to x-line
                    auto time_factor =
                        Map1D{
                            tau_factor + grid.X_nodes.size() * j,
                            grid.X_nodes.size()};

                    // memory chunk in temperature, corresponding to x-line
                    auto data =
                        Map1D{
                            state.cur_state.data() + grid.X_nodes.size() * j,
                            grid.X_nodes.size()};

                    // right handside of Au = b problem
                    Eigen::ArrayX<float_t> rhs = data * time_factor;

                    SpMatrix A{splitY.LaplaceTerm[j]};
                    A = A + time_factor;

                    data = solve_3Diag_problem(A, rhs);
                }
            }

            SplitX splitX;
            SplitY splitY;
            TemporalTerm factor;
            Grid_t &grid;
            State2D state;

        protected:
            auto solve_linear_problem(auto &A, const auto &b)
            {
                A.makeCompressed();

                Eigen::SimplicialCholesky<SpMatrix> chol(A);

                return chol.solve(b);
            }
        };

    } // SplittingMethod

} // EqSolver