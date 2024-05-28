#pragma once

#include "UniformGrid1D.h"

namespace EqSolver
{
    namespace Grid
    {
        struct UniformGrid2D
        {
            struct Point
            {
                float_t x,y;
            };
            using Cell_Volume = Eigen::ArrayXX<float_t>;

            UniformGrid2D(const UniformGrid1D &X_nodes, const UniformGrid1D &Y_nodes)
                : X_nodes{X_nodes}, Y_nodes{Y_nodes}
            {
                cell_volume =
                    X_nodes.cell_volume().matrix() *
                    Y_nodes.cell_volume().transpose().matrix();
            }

            UniformGrid2D(const UniformGrid2D &) noexcept = default;
            UniformGrid2D(UniformGrid2D &&) noexcept = default;

            Point operator()(size_t i, size_t j)
            {
                return {X_nodes[i], Y_nodes[j]};
            }

            UniformGrid1D X_nodes, Y_nodes;
            Cell_Volume cell_volume;
        };
    } // Grid
} // EqSolver