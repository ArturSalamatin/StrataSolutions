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

            State2D(const State2D &) noexcept = default;
            State2D(State2D &&) noexcept = default;

        public:
            State_Container cur_state;
        };

    } // Coefficients
} // EqSolver