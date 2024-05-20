#pragma once

#include <vector>
#include <algorithm>
#include <cassert>

#include "Defines.h"

namespace EqSolver
{
    namespace Grid
    {
        struct UniformGrid1D : private std::vector<float_t>
        {
            public:
            static UniformGrid1D CreateFromStep(float_t a, float_t b, float_t step)
            {
                assert(b > a);
                assert(step > 0.0);
                assert(step < b-a);

                size_t segm_nmbr = static_cast<size_t>((b-a)/step) + 1;

                return CreateFromNodes(a, b, segm_nmbr+1);
            }

            static UniformGrid1D CreateFromNodes(float_t a, float_t b, size_t nodes_nmbr)
            {
                assert(b > a);
                assert(nodes_nmbr > 1ull);

                float_t step = (b-a)/(nodes_nmbr-1ull);
                std::vector<float_t> nodes(nodes_nmbr, a);

                std::generate(nodes.begin(), nodes.end(), [n = 0, &step, &a]() mutable { return a + n++ * step; });

                 return {nodes};
            }

            UniformGrid1D(const UniformGrid1D&) noexcept = default;
            UniformGrid1D(UniformGrid1D&&) noexcept = default;

            float_t step(size_t) const
            {
                return its_step;
            }

            public:
            using std::vector<float_t>::operator[];
            using std::vector<float_t>::data;
            using std::vector<float_t>::size;
            using std::vector<float_t>::begin;
            using std::vector<float_t>::end;
            using std::vector<float_t>::front;
            using std::vector<float_t>::back;

            protected:
                UniformGrid1D(const std::vector<float_t>& nodes) 
                    : std::vector<float_t>{nodes}
                {
                    assert(nodes.size() > 1);
                    its_step = nodes[1]-nodes[0];
                }

            float_t its_step;
        };
    } // Grid
} // EqSolver