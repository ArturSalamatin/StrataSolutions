#include <gtest/gtest.h>

#include <iostream>

#include <src/Grid/Defines.h>
#include <src/GridFactory.hpp>

using namespace EqSolver;
using namespace EqSolver::Grid;

TEST(Grid, ctor)
{
  UniformGrid2D grid{
      GridFactory::CreateGridFromStep(
          Box{0.0, 1.0, 0.0, 1.0},
          Steps{0.1, 0.1})};
}

int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}