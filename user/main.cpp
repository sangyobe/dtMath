#include <dtMath/dtMath.h>

int main()
{
    dt::Math::Matrix<3, 3, double> m1;
    dt::Math::Matrix<3, 3, double> m2;

    m1 << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;
    m2 = m1.Transpose();

    m1.Print();
    m2.Print();

    return 0;
}