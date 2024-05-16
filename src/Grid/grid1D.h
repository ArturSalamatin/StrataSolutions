#pragma once

#include <vector>
#include <algorithm>
#include <cassert>

namespace EqSolver
{
    namespace Grid
    {
        using float_t = double;

        struct StructuredGrid1D
        {
            public:
            static StructuredGrid1D CreateFromStep(float_t a, float_t b, float_t step)
            {
                assert(b > a);
                assert(step > 0.0);
                assert(step < b-a);

                size_t segm_nmbr = (b-a)/step + 1;
                size_t nodes_nmbr{segm_nmbr + 1ull};
                step = (b-a)/segm_nmbr;
                std::vector<float_t> nodes(nodes_nmbr, a);

                 std::generate(nodes.begin(), nodes.end(), [n = 0, &step, &a]() mutable { return a + n++ * step; });

                 return {nodes};
            }

            static StructuredGrid1D CreateFromNodes(float_t a, float_t b, size_t nodes_nmbr)
            {
                assert(b > a);
                assert(nodes_nmbr > 1ull);

                float_t step = (b-a)/(nodes_nmbr-1ull);
                std::vector<float_t> nodes(nodes_nmbr, a);

                 std::generate(nodes.begin(), nodes.end(), [n = 0, &step, &a]() mutable { return a+ n++ * step; });

                 return {nodes};
            }

            StructuredGrid1D(const StructuredGrid1D&) noexcept = default;
            StructuredGrid1D(StructuredGrid1D&&) noexcept = default;


            public:
                std::vector<float_t> data;

            protected:
                StructuredGrid1D(const std::vector<float_t>& nodes) 
                    : data{nodes}
                {}
        };
    } // Grid
} // EqSolver