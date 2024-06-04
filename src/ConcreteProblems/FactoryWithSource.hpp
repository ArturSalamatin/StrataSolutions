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
         * du/dt - d^2 u/dx^2 - d^2 u/dy^2 = 1
         *
         * u(t=0) = 0
         * u(boundary) = 0
         * 
         * solution: u = t
         */
        struct FactoryWithSource
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
                    return 1.0;
                }
            };

            struct BCFunctor : public BoundaryConditions::BCFunctorBase
            {
                float_t operator()(
                    float_t, float_t, float_t t) const override
                {
                    return t;
                }
            };

            Problem::InitialCondition
                zero_state;
            BoundaryConditions::BoundaryConditions
                bc;

            FactoryWithSource(
                const Box &box, const Steps &steps)
                : FactoryWithSource{
                      GridFactory::CreateGridFromStep(box, steps)}
            {
            }

            template <typename Grid_t>
            FactoryWithSource(const Grid_t &grid)
                : zero_state{
                      State::State2D::FillWithZeros(
                          grid)},
                  bc{grid, std::make_shared<BCFunctor>()}
            {
            }
        };
    } // ConcreteProblem
} // EqSolver