#include <gtest/gtest.h>

#include <iostream>
#include <src/Grid/StructuredGrid1D.h>

using namespace EqSolver;
using namespace EqSolver::Grid;

TEST(Grid, ctor) 
{
    StructuredGrid1D grid_x{StructuredGrid1D::CreateFromNodes(0.0, 1.0, 11)};
    for(size_t i = 0; i < grid_x.size(); ++i)
    {
      ASSERT_NEAR(grid_x[i], i / 10.0, 1E-10);
    }
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}