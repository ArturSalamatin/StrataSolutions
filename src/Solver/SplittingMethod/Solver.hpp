#pragma once

#include <vector>
#include <tuple>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

#include "../../Grid/Defines.h"
#include "../BoundaryConditions.hpp"

#include "SplitX.hpp"
#include "SplitY.hpp"

namespace EqSolver
{
    namespace SplittingMethod
    {
        template <typename Grid_t>
        struct Solver
        {
            using Map1D =
                Eigen::Map<
                    Eigen::ArrayX<float_t>>;

            using Map1D_Stride =
                Eigen::Map<
                    Eigen::VectorX<float_t>,
                    0,
                    Eigen::OuterStride<Eigen::Dynamic>>;

            using RHS_t = Eigen::VectorX<float_t>;

            template <typename ProblemFactory_t>
            Solver(
                const Grid_t &grid,
                std::shared_ptr<Properties::Fields>
                    properties,
                const ProblemFactory_t &factory,
                float_t initial_moment = 0.0)
                : splitX{properties, grid},
                  splitY{properties, grid},
                  factor{properties, grid},
                  grid{grid},
                  properties{properties},
                  state{factory.zero_state},
                  bc{factory.bc},
                  states(),
                  time_moments()
            {
                time_moments.push_back(initial_moment);
                states.emplace_back(state);
            }

            void advance(float_t tau)
            {
                // tau_factor multiplies Delta_u at different time moments,
                // t and t+tau
                Eigen::ArrayXX<float_t> tau_factor{
                    factor.DivideByTemporalStep(tau)};
                // solve a set of 1D problems in y-direction, for various x-coords
                solve_split_x(tau, tau_factor.data());
                // solve a set of 1D problems in x-direction, for various y-coords
                solve_split_y(tau, tau_factor.data());

                time_moments.push_back(time_moments.back() + tau);
                states.emplace_back(state);
            }

            struct Solution
            {
                Solution(
                    const std::vector<float_t> &times,
                    const std::vector<State::State2D> &states)
                    : times{times}, states{states}
                {
                }

                auto back() const{
                    return std::pair{times.back(), states.back()};
                }

                const std::vector<float_t> &times;
                const std::vector<State::State2D> &states;
            };

            Solution solution() const
            {
                return {time_moments, states};
            }

        protected:
            void solve_split_x(float_t tau, float_t *tau_factor)
            {
                // nmbr of nodes in the 1D problem
                ptrdiff_t chunk_size = grid.Y_nodes.size();
                // stride in a 1D layout of 2D unknown temperature values
                ptrdiff_t stride_size = grid.X_nodes.size();
                // to be provided to Eigen::Map
                Eigen::OuterStride<Eigen::Dynamic> stride{stride_size};

#pragma omp parallel for
                //  take every line along x-direction. A line per y-node
                for (ptrdiff_t i = 0; i < (ptrdiff_t)grid.X_nodes.size(); ++i)
                {
                    // memory chunk in capacity-container, corresponding to x-line
                    const auto time_factor =
                        Map1D_Stride{
                            tau_factor + i,
                            chunk_size, stride};

                    // memory chunk in temperature, corresponding to x-line
                    auto data =
                        Map1D_Stride{
                            state.cur_state.data() + i,
                            chunk_size, stride};

                    // right handside of Au = b problem
                    RHS_t rhs =
                        (data.array() * time_factor.array()).matrix();

                    SpMatrix A{splitX.LaplaceTerm(i)};
                    A.diagonal() = A.diagonal() + time_factor;
                    applyBC_split_x(A, rhs, i);
                    data = solve_linear_problem(A, rhs);
                }
            }

            void solve_split_y(float_t tau, float_t *tau_factor)
            {
#pragma omp parallel for
                // take every line along x-direction. A line per y-node
                for (ptrdiff_t j = 0; j < (ptrdiff_t)grid.Y_nodes.size(); ++j)
                {
                    // memory chunk in capacity-container, corresponding to x-line
                    const auto time_factor =
                        Map1D{
                            tau_factor + grid.X_nodes.size() * j,
                            (ptrdiff_t)grid.X_nodes.size()};

                    // memory chunk in temperature, corresponding to x-line
                    auto data =
                        Map1D{
                            state.cur_state.data() + grid.X_nodes.size() * j,
                            (ptrdiff_t)grid.X_nodes.size()};

                    auto source =
                        Map1D{
                            properties->source_vol_f.data() + grid.X_nodes.size() * j,
                            (ptrdiff_t)grid.X_nodes.size()};

                    // right handside of Au = b problem
                    RHS_t rhs = (data * time_factor + source).matrix();

                    SpMatrix A{splitY.LaplaceTerm(j)};
                    A.diagonal() = A.diagonal() + time_factor.matrix();
                    applyBC_split_y(A, rhs, j);
                    data = solve_linear_problem(A, rhs);
                }
            }

            SplitX splitX;
            SplitY splitY;
            Problem::BoundaryConditions bc;
            std::shared_ptr<Properties::Fields> properties;
            TemporalTerm factor;
            const Grid_t &grid;
            State::State2D state;
            std::vector<State::State2D> states;
            std::vector<float_t> time_moments;

        protected:
            Eigen::VectorXd solve_linear_problem(
                const SpMatrix &A,
                const Eigen::VectorX<float_t> &b)
            {
                // A.makeCompressed();
                Eigen::SparseLU<SpMatrix> lu;
                lu.analyzePattern(A); // this is common for every matrix A. Can be optimized
                lu.factorize(A);
                return lu.solve(b);
            }

            void applyBC_split_x(SpMatrix &A, RHS_t &b, ptrdiff_t i)
            {
                A.coeffRef(0, 0) = 1;
                A.coeffRef(0, 1) = 0;
                b(0) = bc.east_west.west_vals[i];
                ptrdiff_t n = A.outerSize() - 1;
                A.coeffRef(n, n) = 1.0;
                A.coeffRef(n, n - 1) = 0.0;
                b(n) = bc.east_west.east_vals[i];
            }

            void applyBC_split_y(SpMatrix &A, RHS_t &b, ptrdiff_t j)
            {
                A.coeffRef(0, 0) = 1;
                A.coeffRef(0, 1) = 0;
                b(0) = bc.south_north.south_vals[j];
                ptrdiff_t n = A.outerSize() - 1;
                A.coeffRef(n, n) = 1.0;
                A.coeffRef(n, n - 1) = 0.0;
                b(n) = bc.south_north.north_vals[j];
            }
        };
    } // SplittingMethod
} // EqSolver