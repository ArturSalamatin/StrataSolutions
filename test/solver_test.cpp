#include <gtest/gtest.h>

#include <iostream>
#include <memory>

#include <src/Grid/Defines.h>
#include <src/Solver/Solver.hpp>

#include <src/GridFactory.hpp>

#include <src/ConcreteProblem/ProblemFactory.hpp>

using namespace EqSolver;
using namespace EqSolver::Grid;
using namespace EqSolver::Properties;
using namespace EqSolver::SplittingMethod;

TEST(Solver, splitX)
{
  const double tol = 1E-8;

  UniformGrid2D grid{
      GridFactory::CreateGridFromStep(
          Box{0.0, 1.0, 0.0, 1.0},
          Steps{0.11, 0.11})};

  for (size_t i = 0; i < grid.X_nodes.size(); ++i)
    ASSERT_NEAR(grid.X_nodes.step(i), 0.1, tol);
  for (size_t j = 0; j < grid.Y_nodes.size(); ++j)
    ASSERT_NEAR(grid.Y_nodes.step(j), 0.1, tol);

  for (size_t i = 1; i < grid.X_nodes.size() - 1; ++i)
    for (size_t j = 1; j < grid.Y_nodes.size() - 1; ++j)
    {
      ASSERT_NEAR(
          grid.cell_volume(i, j),
          grid.X_nodes.step(i) * grid.Y_nodes.step(j),
          tol);
    }

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
    
  for (size_t j = 0; j < grid.Y_nodes.size(); j += grid.Y_nodes.size() - 1)
    for (size_t i = 0; i < grid.X_nodes.size(); i += grid.X_nodes.size() - 1)
    {
      ASSERT_NEAR(
          grid.cell_volume(i, j),
          grid.X_nodes.step(i) * grid.Y_nodes.step(j) / 4.0,
          tol);
    }


  std::shared_ptr<Properties::Fields>
      properties{
          std::make_shared<Properties::Fields>(
              grid,
              Capacity{},
              Conductivity{})};

  SplitX splitx{properties, grid};
}

int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}