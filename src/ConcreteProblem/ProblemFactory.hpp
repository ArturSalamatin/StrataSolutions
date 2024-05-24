#pragma once

#include "../Grid/Defines.h"
#include "../Solver/State2D.hpp"

namespace EqSolver
{
    namespace ConcreteProblem
    {
        struct InitialCondition : public State::State2D
        {
            InitialCondition(const State::State2D &state)
                : State::State2D{state}
            {
            }
        };

        struct BoundaryCondition
        {
            enum BCType
            {
                first,
                second,
                third
            };

            BCType type;
        };

        struct BCSouth : public BoundaryCondition
        {
        };

        struct BCNorth : public BoundaryCondition
        {
        };

        struct BCEast : public BoundaryCondition
        {
        };

        struct BCWest : public BoundaryCondition
        {
        };

        struct BoundaryConditions
        {
            BoundaryConditions()
                : south{BoundaryCondition::BCType::first},
                  east{BoundaryCondition::BCType::first},
                  north{BoundaryCondition::BCType::first},
                  west{BoundaryCondition::BCType::first}
            {
            }

            BCSouth south;
            BCEast east;
            BCNorth north;
            BCWest west;
        };

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

        struct ProblemFactory
        {
            InitialCondition zero_state;
            BoundaryConditions bc;
            Capacity capacity;
            Conductivity conductivity;

            ProblemFactory(const Box &box, const Steps &steps)
                : ProblemFactory{Grid::UniformGrid2D{
                      Grid::UniformGrid1D::CreateFromStep(
                          box.x_a, box.x_b, steps.step_x),
                      Grid::UniformGrid1D::CreateFromStep(
                          box.y_a, box.y_b, steps.step_y)}}
            {
            }

            template <typename Grid_t>
            ProblemFactory(const Grid_t &grid)
                : zero_state{
                      State::State2D::FillWithZeros(
                          grid)},
                  bc{}, capacity{}, conductivity{}
            {
            }
        };
    } // ConcreteProblem
} // EqSolver