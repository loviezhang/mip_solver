#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;

#if 0
TEST(Eigen, resize) {
    MatrixXi m(3, 3);
    m << 1, 2, 3,
         4, 5, 6,
         7, 8, 9;
    cout << "m:" << endl << m << endl;

    MatrixXi m2(m.rows(), m.cols()-1);
    int col = 2;
    cout << "m.block(0, 0, 3, col):" << endl << m.block(0, 0, m.rows(), 1) << endl;
    cout << "m.block(0, col, 3, m.cols() - col - 1):" << endl << m.block(0, col+1, m.rows(), m.cols() - col - 1) << endl;
    m2 << m.block(0, 0, m.rows(), col), m.block(0, col+1, m.rows(), m.cols() - col - 1);
    cout << "m2:" << endl << m2 << endl;

    ASSERT_TRUE(true);
}
#endif
