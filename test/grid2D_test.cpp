#include <gtest/gtest.h>

#include <iostream>
#include <src/Grid/UniformGrid2D.h>
#include <src/GridFactory.hpp>

using namespace EqSolver;
using namespace EqSolver::Grid;

TEST(Grid2D, ctor)
{
  const double tol = 1E-8;

  UniformGrid2D grid{
      GridFactory::CreateGridFromStep(
          Box{0.0, 1.0, 0.0, 1.0},
          Steps{0.065, 0.11})};

  // test steps of 1D grids
  for (size_t i = 0; i < grid.X_nodes.size(); ++i)
    ASSERT_NEAR(grid.X_nodes.step(i), 0.0625, tol);
  for (size_t j = 0; j < grid.Y_nodes.size(); ++j)
    ASSERT_NEAR(grid.Y_nodes.step(j), 0.1, tol);

  // test volumes of internal cells
  for (size_t i = 1; i < grid.X_nodes.size() - 1; ++i)
    for (size_t j = 1; j < grid.Y_nodes.size() - 1; ++j)
    {
      ASSERT_NEAR(
          grid.cell_volume(i, j),
          grid.X_nodes.step(i) * grid.Y_nodes.step(j),
          tol);
    }

  // test volumes of boundary cells except corners
  for (size_t i = 0; i < grid.X_nodes.size(); i += grid.X_nodes.size() - 1)
    for (size_t j = 1; j < grid.Y_nodes.size() - 1; ++j)
    {
      ASSERT_NEAR(
          grid.cell_volume(i, j),
          grid.X_nodes.step(i) * grid.Y_nodes.step(j) / 2.0,
          tol);
    }
  for (size_t j = 0; j < grid.Y_nodes.size(); j += grid.Y_nodes.size() - 1)
    for (size_t i = 1; i < grid.X_nodes.size() - 1; ++i)
    {
      ASSERT_NEAR(
          grid.cell_volume(i, j),
          grid.X_nodes.step(i) * grid.Y_nodes.step(j) / 2.0,
          tol);
    }

  // test volumes of corner cells
  for (size_t j = 0; j < grid.Y_nodes.size(); j += grid.Y_nodes.size() - 1)
    for (size_t i = 0; i < grid.X_nodes.size(); i += grid.X_nodes.size() - 1)
    {
      ASSERT_NEAR(
          grid.cell_volume(i, j),
          grid.X_nodes.step(i) * grid.Y_nodes.step(j) / 4.0,
          tol);
    }
}

int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}