#include <gtest/gtest.h>

#include <iostream>
#include <src/Grid/UniformGrid1D.h>

using namespace EqSolver;
using namespace EqSolver::Grid;

TEST(Grid1D, ctor)
{
  const double tol = 1E-10;

  UniformGrid1D grid_x{UniformGrid1D::CreateFromNodes(0.0, 1.0, 11)};
  for (size_t i = 0; i < grid_x.size(); ++i)
    ASSERT_NEAR(grid_x[i], i / 10.0, tol);
}

TEST(Grid1D, ctor2)
{
  const double tol = 1E-10;

  UniformGrid1D grid_x{UniformGrid1D::CreateFromNodes(0.0, 1.0, 17)};
  for (size_t i = 0; i < grid_x.size(); ++i)
    ASSERT_NEAR(grid_x[i], i / 16.0, tol);
}

int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}