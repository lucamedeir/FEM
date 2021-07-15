#include <gtest/gtest.h>

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    auto test = RUN_ALL_TESTS();
   
    int result = 0;


    return test;
}