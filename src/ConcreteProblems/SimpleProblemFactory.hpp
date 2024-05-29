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
         * u(t=0) = 0
         * u(boundary) = 0
         */
        struct SimpleProblemFactory
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

            struct BCFunctor
            {
                float_t operator()(
                    float_t x, float_t y) const
                {
                    return 0.0;
                }
            };

            Problem::InitialCondition
                zero_state;
            Problem::BoundaryConditions
                bc;

            SimpleProblemFactory(
                const Box &box, const Steps &steps)
                : SimpleProblemFactory{
                      GridFactory::CreateGridFromStep(box, steps)}
            {
            }

            template <typename Grid_t>
            SimpleProblemFactory(const Grid_t &grid)
                : zero_state{
                      State::State2D::FillWithZeros(
                          grid)},
                  bc{grid, BCFunctor{}}
            {
            }
        };
    } // ConcreteProblem
} // EqSolver