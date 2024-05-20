#pragma once

#include "StructuredGrid1D.h"

namespace EqSolver
{
    namespace Grid
    {
        struct StructuredGrid2D
        {
            public:

            StructuredGrid2D(const UniformGrid1D& X_nodes, const UniformGrid1D& Y_nodes)
            : X_nodes{X_nodes}, Y_nodes{Y_nodes}
            {}

            StructuredGrid2D(const StructuredGrid2D&) noexcept = default;
            StructuredGrid2D(StructuredGrid2D&&) noexcept = default;

            public:
                UniformGrid1D X_nodes, Y_nodes;
        };
    } // Grid
} // EqSolver