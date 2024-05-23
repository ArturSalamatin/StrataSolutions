#pragma once

#include "UniformGrid1D.h"

namespace EqSolver
{
    namespace Grid
    {
        struct UniformGrid2D
        {
            using Cell_Volume = Eigen::ArrayXX<float_t>;
            public:

            UniformGrid2D(const UniformGrid1D& X_nodes, const UniformGrid1D& Y_nodes)
            : X_nodes{X_nodes}, Y_nodes{Y_nodes}
            {
                cell_volume = 
                X_nodes.cell_volume().matrix() * 
                Y_nodes.cell_volume().transpose().matrix();
            }

            UniformGrid2D(const UniformGrid2D&) noexcept = default;
            UniformGrid2D(UniformGrid2D&&) noexcept = default;

            public:
                UniformGrid1D X_nodes, Y_nodes;
                Cell_Volume cell_volume;
        };
    } // Grid
} // EqSolver