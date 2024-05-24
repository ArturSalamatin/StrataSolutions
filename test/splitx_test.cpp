#include <gtest/gtest.h>

#include <iostream>
#include <memory>

#include <src/Grid/Defines.h>
#include <src/Solver/SplittingMethod/SplitX.hpp>

#include <src/GridFactory.hpp>

#include <src/ConcreteProblem/ProblemFactory.hpp>

using namespace EqSolver;
using namespace EqSolver::Grid;
using namespace EqSolver::Properties;
using namespace EqSolver::ConcreteProblem;
using namespace EqSolver::SplittingMethod;

TEST(Solver, splitX)
{
  const double tol = 1E-8;

  UniformGrid2D grid{
      GridFactory::CreateGridFromStep(
          Box{0.0, 1.0, 0.0, 1.0},
          Steps{0.11, 0.11})};

//  ConcreteProblem::ProblemFactory factory{grid};

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