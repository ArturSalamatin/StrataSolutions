#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>

#include "../Grid/Defines.h"
#include "../Grid/UniformGrid2D.h"

namespace EqSolver
{
    namespace State
    {
        struct State2D
        {
            using State_Container =
                Eigen::ArrayXX<float_t>;

            State2D(
                const State_Container &cur_state)
                : cur_state{cur_state}
            {
            }

            static State2D FillWithZeros(
                const Grid::UniformGrid2D &grid)
            {
                State_Container cur_state{
                    grid.X_nodes.size(),
                    grid.Y_nodes.size()};
                cur_state.fill(0.0);

                return {cur_state};
            }
            
            template<typename Functor>
            static State2D FillWithFunctor(
                const Grid::UniformGrid2D &grid, const Functor& f, float_t initial_moment)
            {
                State_Container cur_state{
                    grid.X_nodes.size(),
                    grid.Y_nodes.size()};

                for(ptrdiff_t j = 0; j < cur_state.outerSize(); ++j)
                for(ptrdiff_t i = 0; i < cur_state.innerSize(); ++i)
                {
                    cur_state(i,j) = f(grid.X_nodes[i], grid.Y_nodes[j], initial_moment);
                }

                return {cur_state};
            }


            State2D(const State2D &) noexcept = default;
            State2D(State2D &&) noexcept = default;

            float_t operator()(ptrdiff_t i, ptrdiff_t j) const
            {
                return cur_state(i, j);
            }
            
            float_t& operator()(ptrdiff_t i, ptrdiff_t j) 
            {
                return cur_state(i, j);
            }

        public:
            State_Container cur_state;
        };

    } // Coefficients
} // EqSolver