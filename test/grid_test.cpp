#include <gtest/gtest.h>

#include <iostream>

using namespace OptLib;


TEST(Grid, ctor) {

    // testing Point ctors
    auto rp1{Point<3>{1.0, 2.0, 3.0}};
    auto rp2 = rp1;
    auto rp3 = Point<3>{1.0, 2.0, 3.0};
    auto rp4{std::move(rp2)};
    auto rp5 = std::move(rp4);
    Point<3> rp6{2.0, 1.0, 3.0};

    ASSERT_EQ(0, 0);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}