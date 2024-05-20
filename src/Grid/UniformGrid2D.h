#pragma once

#include "UniformGrid1D.h"

namespace EqSolver
{
    namespace Grid
    {
        struct UniformGrid2D
        {
            public:

            UniformGrid2D(const UniformGrid1D& X_nodes, const UniformGrid1D& Y_nodes)
            : X_nodes{X_nodes}, Y_nodes{Y_nodes}
            {}

            UniformGrid2D(const UniformGrid2D&) noexcept = default;
            UniformGrid2D(UniformGrid2D&&) noexcept = default;

            public:
                UniformGrid1D X_nodes, Y_nodes;
        };
    } // Grid
} // EqSolver