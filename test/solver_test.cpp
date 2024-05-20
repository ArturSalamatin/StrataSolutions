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

TEST(Solver, splitX)
{
  UniformGrid2D grid{
      GridFactory::CreateGridFromStep(
          Box{0.0, 1.0, 0.0, 1.0},
          Steps{0.1, 0.1})};

   std::shared_ptr<Properties::Fields>
     properties{std::make_shared<Properties::Fields>(grid, Capacity{}, Conductivity{})};
}

int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}