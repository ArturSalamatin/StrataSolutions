#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "../Grid/Defines.h"
#include "../Grid/StructuredGrid2D.h"
#include "../Properties/Coefficients.hpp"

namespace EqSolver
{
    struct State2D
    {
        static State2D FillWithZeros(const Grid::StructuredGrid2D& grid)
        {
            Eigen::MatrixX<float_t> cur_state{grid.X_nodes.size(), grid.Y_nodes.size()};
            cur_state.fill(0.0);

            return {cur_state};

        }

        State2D(const State2D&) noexcept = default;
        State2D(State2D&&) noexcept = default;

        protected:
            Eigen::MatrixX<float_t> cur_state;

        protected:
            State2D(const Eigen::MatrixX<float_t>& cur_state) 
                : cur_state{cur_state}
            { }
    };

    namespace SplittingMethod
    {
        using SpMatrix = Eigen::SparseMatrix<float_t>;
        using VectSpMatrix = std::vector<Eigen::SparseMatrix<float_t>>;

        template<typename Properties_t>
        struct SplitX
        {
            SplitX(
                const Properties_t& properties,
                const Grid::StructuredGrid& grid)
            : LaplaceTerm(
                grid.X_nodes.size(), 
                Eigen::SparseMatrix<float_t>{
                    grid.Y_nodes.size(), 
                    grid.Y_nodes.size()
                }
            ),
            grid{grid}
            {
                for(auto& m : matrcies)
                {
                    m.reserve(grid.Y_nodes.size()*3);
                }
            }

            protected:
            VectSpMatrix LaplaceTerm;
            Grid::StructuredGrid grid;
            Grid::StructuredGrid grid;
            Properties_t properties;

            void FillLaplaceTerm()
            {
                Eigen::Matrix<float_t> conductivity_nodes{
                    grid.X_nodes.size(), 
                    grid.Y_nodes.size()};

#pragma omp parallel for
                for(size_t j = 0; j < conductivity_nodes.cols(); ++j)
                    for(size_t i = 0; i < conductivity_nodes.rows(); ++i)
                    {
                        conductivity_nodes(i,j) = 
                            properties(
                                grid.X_nodes[i], 
                                grid.Y_nodes[j]);
                    }

                
                Eigen::Matrix<float_t> conductivity_shift_y{
                    grid.X_nodes.size(), 
                    grid.Y_nodes.size()+1};
                    
                for(size_t i = 0; i < conductivity_shift_y.rows(); ++i)
                {
                    conductivity_shift_y(i,0) = 
                        conductivity_nodes(i, 0);
                    for(size_t j = 1; j < conductivity_shift_y.cols()-1; ++j)
                    {
                        conductivity_shift_y(i,j) = 
                            (conductivity_nodes(i, j-1)+conductivity_nodes(i,j))/2.0;
                    }
                    conductivity_shift_y(i,end) = 
                        conductivity_nodes(i, end);
                }
        };
    }

} // EqSolver