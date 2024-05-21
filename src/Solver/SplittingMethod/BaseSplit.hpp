#pragma once

#include <vector>
#include <memory>

// #include <Eigen/Core>
// #include <Eigen/Dense>
#include <Eigen/SparseCore>

#include "../../Grid/Defines.h"
#include "../../Properties/Coefficients.hpp"

namespace EqSolver
{
    namespace SplittingMethod
    {
        struct TemporalTerm
        {
            using Factor = Properties::Fields::Conductivity_f;

            template <typename Grid_t>
            TemporalTerm(
                std::shared_ptr<Properties::Fields> properties,
                const Grid_t &grid)
                : factor{properties->capacity_vol_f} // volumes are taken into account
            {
            }

            auto UseTemporalStep(float_t tau) const
            {
                return (1.0/tau)*factor;
            }

protected:
            Factor factor;
        };

        using SpMatrix = Eigen::SparseMatrix<float_t>;
        using VectSpMatrix = std::vector<Eigen::SparseMatrix<float_t>>;

        struct BaseSplit
        {
            BaseSplit(
                std::shared_ptr<Properties::Fields> properties,
                size_t serial_nodes,
                size_t matrix_size)
                : properties{properties},
                  LaplaceTerm(
                      serial_nodes,                // number of matricies
                      Eigen::SparseMatrix<float_t>{// ctor per matrix
                                                   ptrdiff_t(matrix_size),
                                                   ptrdiff_t(matrix_size)})
            {
                for (auto &m : LaplaceTerm)
                {
                    m.reserve(matrix_size * 3ull - 2ull);
                }
            }

        public:
            VectSpMatrix LaplaceTerm;
            std::shared_ptr<Properties::Fields> properties;

            using Conductivity_f = Properties::Fields::Conductivity_f;
        };
    } // SplittingMethod
} // EqSolver