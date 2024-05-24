#include <gtest/gtest.h>

#include <iostream>
#include <memory>

#include <src/Grid/Defines.h>
#include <src/Solver/SplittingMethod/SplitX.hpp>

#include <src/GridFactory.hpp>

// #include <src/ConcreteProblem/ProblemFactory.hpp>

using namespace EqSolver;
using namespace EqSolver::Grid;
using namespace EqSolver::Properties;
// using namespace EqSolver::ConcreteProblem;
using namespace EqSolver::SplittingMethod;

struct Capacity
{
  EqSolver::float_t operator()(EqSolver::float_t x, EqSolver::float_t y) const
  {
    return 1.0;
  }
};

struct Conductivity
{
  EqSolver::float_t operator()(EqSolver::float_t x, EqSolver::float_t y) const
  {
    return 1.0;
  }
};

struct Source
{
  EqSolver::float_t operator()(EqSolver::float_t x, EqSolver::float_t y) const
  {
    return 1.0;
  }
};

TEST(Solver, splitX)
{
  const double tol = 1E-8;

  UniformGrid2D grid{
      GridFactory::CreateGridFromStep(
          Box{0.0, 1.0, 0.0, 1.0},
          Steps{0.11, 0.11})};

  std::shared_ptr<Properties::Fields>
      properties{
          std::make_shared<Properties::Fields>(
              grid,
              Capacity{},
              Conductivity{},
              Source{})};

  SplitX splitx{properties, grid};
}

int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}