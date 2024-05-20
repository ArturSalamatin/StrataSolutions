#pragma once

#include "Grid/Defines.h"
#include "Grid/StructuredGrid2D.h"

namespace EqSolver
{
    struct GridFactory
    {
        static Grid::StructuredGrid2D CreateGridFromStep(const Box& box, const Steps& steps)
        {
            return {
                Grid::StructuredGrid1D::CreateFromStep(
                    box.x_a, box.x_b, steps.step_x), 
                Grid::StructuredGrid1D::CreateFromStep(
                    box.y_a, box.y_b, steps.step_y)
            };
        }
    };
} // EqSolver
