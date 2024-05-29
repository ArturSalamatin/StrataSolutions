#pragma once

#include "../Grid/Defines.h"
#include "../Solver/State2D.hpp"

#include "InitialCondition.hpp"
#include "BoundaryConditions.hpp"

namespace EqSolver
{
    namespace ConcreteProblem
    {
        /**
         * @brief Quite simple problem for tests
         * All coefficients are unities
         * Source term is zero
         * 
         * Inital state is zero
         * 
         * BCs are of first kind: zero values 
         * 
         */
        struct SimpleProblemFactory
        {
            struct Capacity
            {
                float_t operator()(float_t x, float_t y) const
                {
                    return 1.0;
                }
            };

            struct Conductivity
            {
                float_t operator()(float_t x, float_t y) const
                {
                    return 1.0;
                }
            };

            struct Source
            {
                float_t operator()(float_t, float_t) const
                {
                    return 0.0;
                }
            };

            struct BCFunctor
            {
                float_t operator()(float_t x, float_t y) const
                {
                    return 1+x+y;
                }
            };
            
            Problem::InitialCondition zero_state;
            Problem::BoundaryConditions<BCFunctor> bc;

            SimpleProblemFactory(
                const Box &box, const Steps &steps)
                : SimpleProblemFactory{Grid::UniformGrid2D{
                      Grid::UniformGrid1D::CreateFromStep(
                          box.x_a, box.x_b, steps.step_x),
                      Grid::UniformGrid1D::CreateFromStep(
                          box.y_a, box.y_b, steps.step_y)}}
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