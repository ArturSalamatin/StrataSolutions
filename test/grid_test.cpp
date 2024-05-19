#include <gtest/gtest.h>

#include <iostream>
#include <src/Grid/StructuredGrid1D.h>

using namespace EqSolver;
using namespace EqSolver::Grid;


TEST(Grid, ctor) {

StructuredGrid1D grid_x{StructuredGrid1D::CreateFromNodes(0.0, 1.0, 11)};
    ASSERT_EQ(grid_x[0], 0.0);
    ASSERT_EQ(grid_x.back(), 1.0);
    ASSERT_EQ(grid_x[5], 0.5);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}