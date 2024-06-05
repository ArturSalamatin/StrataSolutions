#pragma once

#include "../Grid/Defines.h"
#include "../Grid/UniformGridFactory.hpp"
#include "../Solver/State2D.hpp"

#include "../Solver/InitialCondition.hpp"
#include "../Solver/BoundaryConditions.hpp"

namespace EqSolver
{
    namespace ConcreteProblem
    {
        /**
         * @brief Quite simple problem for tests
         * du/dt - d^2 u/dx^2 - d^2 u/dy^2 = 0
         *
         * solution: u = 1/(t-tau)*exp(-(r-R)^2/4*(t-tau))
         */

        struct BaseFactory
        {
            struct Point
            {
                float_t tau, X, Y;
            };

            struct Solution
            {
                float_t operator()(
                    float_t x, float_t y, float_t t) const
                {
                    float_t r2 = (x - X) * (x - X) + (y - Y) * (y - Y);

                    return (1 / (t - tau)) * exp(-r2 / (4 * (t - tau)));
                }
                Solution(const Point &p)
                    : X{p.X}, Y{p.Y}, tau{p.tau}
                {
                }

            protected:
                float_t tau, X, Y;
            };

            Point P;
            Solution solution;

            BaseFactory(const Point &P)
                : P{P}, solution{P}
            {
            }
        };

        struct SelfSimilarProblem : public BaseFactory
        {
            struct Capacity
            {
                float_t operator()(
                    float_t x, float_t y) const
                {
                    return 1.0;
                }
            };

            struct Conductivity
            {
                float_t operator()(
                    float_t x, float_t y) const
                {
                    return 1.0;
                }
            };

            struct Source
            {
                float_t operator()(
                    float_t, float_t) const
                {
                    return 0.0;
                }
            };

            struct BCFunctor : public BoundaryConditions::BCFunctorBase
            {
                float_t operator()(
                    float_t x, float_t y, float_t t) const override
                {
                    return solution(x, y, t);
                }
                BCFunctor(const Solution &s)
                    : solution{s}
                {
                }

            protected:
                const Solution &solution;
            };

            Problem::InitialCondition
                zero_state;
            BoundaryConditions::BoundaryConditions
                bc;

            SelfSimilarProblem(
                const Box &box, const Steps &steps, const Point &p)
                : SelfSimilarProblem{
                      GridFactory::CreateGridFromStep(box, steps), p}
            {
            }

            template <typename Grid_t>
            SelfSimilarProblem(const Grid_t &grid, const Point &p)
                : BaseFactory{p},
                  zero_state{initial_state(grid, solution, p.tau)},
                  bc{grid, std::make_shared<BCFunctor>(solution)}
            {
            }

        protected:
            template <typename Grid_t>
            static Problem::InitialCondition initial_state(
                const Grid_t &grid, const Solution &solution, float_t initial_moment)
            {
                return {State::State2D::FillWithFunctor(grid, solution, initial_moment)};
            }
        };
    } // ConcreteProblem
} // EqSolver