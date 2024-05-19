#pragma once

#include <vector>
#include <algorithm>
#include <cassert>

#include "Defines.h"

namespace EqSolver
{
    namespace Grid
    {
        struct StructuredGrid1D : private std::vector<float_t>
        {
            public:
            static StructuredGrid1D CreateFromStep(float_t a, float_t b, float_t step)
            {
                assert(b > a);
                assert(step > 0.0);
                assert(step < b-a);

                size_t segm_nmbr = static_cast<size_t>((b-a)/step) + 1;
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
            //     std::vector<float_t> data;

            using std::vector<float_t>::operator[];
            using std::vector<float_t>::data;
            using std::vector<float_t>::begin;
            using std::vector<float_t>::end;
            using std::vector<float_t>::front;
            using std::vector<float_t>::back;

            protected:
                StructuredGrid1D(const std::vector<float_t>& nodes) 
                    : std::vector<float_t>{nodes}
                {}
        };
    } // Grid
} // EqSolver