#include "testVector6.h"
#include "testPrint.h"

#include <dtMath/dtMath.h>

using dtMath::dtMatrix;
using dtMath::dtVector;
using dtMath::dtVector3;
using dtMath::dtVector6;

void Test_Vector6()
{
    PrintTitle("Test Vector6");
    Vec6Init();
    Vec6MemberFunc();
    Vec6MemberAccessOperator();
    Vec6ArithmeticSum();
    Vec6ArithmeticProduct();
    Vec6ComparisonOperator();
}

void Vec6Init()
{
    PrintHeading("Initialize Vector6 ");
    Printf("CommaInit, vec << 1, 2, 3, 4, 5, 6;\n");
    Printf("vec = \n");

    dtVector6<> vec;
    vec << 1, 2, 3, 4, 5, 6;

    vec.Print();
    Println;
}

void Vec6MemberFunc()
{
    PrintHeading("Vector6 Member Functions ");

    float v[6] = {1, 2, 3, 4, 5, 6};

    Printf("/* Class Create: () */\n");
    dtVector6<> v1;
    v1.Print();
    Println;

    Printf("/* Class Create: with array */\n");
    dtVector6<> v2(v, sizeof(v));
    Printf("array[6] = {1,2,3,4,5,6}\n");
    Printf("vec =\n");
    v2.Print();
    Println;

    Printf("/* Class Create: with element */\n");
    dtVector6<> v3(1, 2, 3, 4, 5, 6);
    Printf("elements are 1, 2, 3, 4, 5, 6\n");
    Printf("vec =\n");
    v3.Print();
    Println;

    Printf("/* Class Create: with dtVector6 */\n");
    dtVector6<> v4(v3);
    Printf("vec = [1,2,3,4,5,6]T\n");
    Printf("vec =\n");
    v4.Print();
    Println;

    Printf("/* Class Create: with dtVector3 and dtVector3 */\n");
    dtVector3<> vec1, vec2;
    vec1.SetElement(1, 2, 3);
    vec2.SetElement(4, 5, 6);
    dtVector6<> v5(vec1, vec2);
    Printf("vec1 = [1,2,3]T, vec2 = [4,5,6]T\n");
    Printf("vec =\n");
    v5.Print();
    Println;

    Printf("/* Class Create: with dtVector */\n");
    dtVector<6> v61(v, sizeof(v));
    dtVector6<> v6(v61);
    Printf("vec = [1,2,3,4,5,6]T\n");
    Printf("vec =\n");
    v6.Print();
    Println;

    Printf("/* Class Create: with dtMatrix */\n");
    dtMatrix<6, 1> m61(v, sizeof(v));
    dtVector6<> v7(m61);
    Printf("mat61 = [1,2,3,4,5,6]T\n");
    Printf("vec =\n");
    v7.Print();
    Println;

    Printf("/* Function: SetZero() */\n");
    Printf("vec =\n");
    v2.SetZero();
    v2.Print();
    Println;

    Printf("/* Function: SetFill(value) */\n");
    Printf("value is 1\n");
    Printf("vec =\n");
    v2.SetFill(1);
    v2.Print();
    Println;

    Printf("/* Function: SetElement(* element) */\n");
    Printf("element[3] = {1,2,3,4,5,6}\n");
    Printf("vec =\n");
    v2.SetElement(v, sizeof(v));
    v2.Print();
    Println;

    Printf("/* Function: SetElement(px, py, pz, ox, oy, oz) */\n");
    Printf("px=1, py=2, pz=3, ox=4, oy=5, oz=6\n");
    Printf("vec =\n");
    v2.SetElement(1, 2, 3, 4, 5, 6);
    v2.Print();
    Println;

    Printf("/* Function: SetElement(dtVector6) */\n");
    Printf("vec = [1,2,3,4,5,6]T\n");
    Printf("vec =\n");
    v2.SetElement(v3);
    v2.Print();
    Println;

    Printf("/* Function: SetElement(dtVector3, dtVector3) */\n");
    Printf("vec1 = [1,2,3]T, vec2 = [4,5,6]T\n");
    Printf("vec =\n");
    v2.SetElement(vec1, vec2);
    v2.Print();
    Println;

    Printf("/* Function: SetElement(dtVector) */\n");
    Printf("vec = [1,2,3,4,5,6]T\n");
    Printf("vec =\n");
    v2.SetElement(v61);
    v2.Print();
    Println;

    Printf("/* Function: SetElement(dtMatrix) */\n");
    Printf("mat61 = [1,2,3,4,5,6]T\n");
    Printf("vec =\n");
    v2.SetElement(m61);
    v2.Print();
    Println;

    Printf("/* Function: SetSwap(i, j) */\n");
    Printf("vec =\n");
    v2.Print();
    Printf("SetSwap(0,2);\n");
    Printf("vec =\n");
    v2.SetSwap(0, 2);
    v2.Print();
    Println;

    Printf("/* Function: SetNormalize() */\n");
    v2.SetFill(1);
    Printf("vec =\n");
    v2.Print();
    v2.SetNormalize();
    Printf("Normalized vec =\n");
    v2.Print();
    Println;

    Printf("/* Function: GetPos() */\n");
    v2.SetElement(1, 2, 3, 4, 5, 6);
    Printf("vec =\n");
    v2.Print();
    Printf("GetPos = \n");
    v2.GetPos().Print();
    Println;

    Printf("/* Function: GetOri() */\n");
    v2.SetElement(1, 2, 3, 4, 5, 6);
    Printf("vec =\n");
    v2.Print();
    Printf("GetOri = \n");
    v2.GetOri().Print();
    Println;

    Printf("/* Function: GetNorm() */\n");
    v2.SetFill(1);
    Printf("vec =\n");
    v2.Print();
    Printf("Norm is %f\n", v2.GetNorm());
    Println;

    Printf("/* Function: GetSqNorm() */\n");
    v2.SetFill(1);
    Printf("vec =\n");
    v2.Print();
    Printf("Squared Norm is %f\n", v2.GetSqNorm());
    Println;

    Printf("/* Function: GetSum() */\n");
    v2.SetFill(1);
    Printf("vec =\n");
    v2.Print();
    Printf("Sum of elements is %f\n", v2.GetSum());
    Println;

    Printf("/* Function: GetNormalized() */\n");
    v2.SetFill(1);
    Printf("vec =\n");
    v2.Print();
    Printf("Normalized vec =\n");
    v2.GetNormalized().Print();
    Println;

    Printf("/* Function: Transpose() */\n");
    v2.SetElement(1, 2, 3, 4, 5, 6);
    Printf("vec =\n");
    v2.Print();
    Printf("Transposed vec =\n");
    v2.Transpose().Print();
    Println;
}

void Vec6MemberAccessOperator()
{
    dtVector6<> v1;

    PrintHeading("Vector6 Access operator ");
    Printf("/* Operator: () */\n");
    Printf("vec(i) = i; Printf("
           "%%f"
           ", vec(i));\n");

    for (uint16_t i = 0; i < 6; i++)
    {
        v1(i) = i;
        Printf("%7.3f ", v1(i));
    }

    Println;
    Println;
}

void Vec6ArithmeticSum()
{
    PrintHeading("Vector6 Arithmetic operators - Sum ");

    float v[6] = {1, 2, 3, 4, 5, 6};
    dtVector6<> v1(1, 2, 3, 4, 5, 6);
    dtVector<6> v61(v, sizeof(v));
    dtMatrix<6, 1> m61(v, sizeof(v));

    Printf("vec =\n");
    v1.Print();
    Println;

    Printf("/* Operator: -() */\n");
    Printf("-vec =\n");
    (-v1).Print();
    Println;

    Printf("/* Operator: +(dtVector6) */\n");
    Printf("vec + (-vec) =\n");
    (v1 + (-v1)).Print();
    Println;

    Printf("/* Operator: -(dtVector6) */\n");
    Printf("vec - vec =\n");
    (v1 - v1).Print();
    Println;

    Printf("/* Operator: +(dtVector) */\n");
    Printf("vec = [1,2,3,4,5,6]T, vec61 = [1,2,3,4,5,6]T\n");
    Printf("vec + (-vec61) =\n");
    (v1 + (-v61)).Print();
    Println;

    Printf("/* Operator: -(dtVector) */\n");
    Printf("vec = [1,2,3,4,5,6]T, vec61 = [1,2,3,4,5,6]T\n");
    Printf("vec - vec61 =\n");
    (v1 - v61).Print();
    Println;

    Printf("/* Operator: +(dtMatrix) */\n");
    Printf("vec = [1,2,3,4,5,6]T, mat61 = [1,2,3,4,5,6]T\n");
    Printf("vec + (-mat61) =\n");
    (v1 + (-m61)).Print();
    Println;

    Printf("/* Operator: -(dtMatrix) */\n");
    Printf("vec = [1,2,3,4,5,6]T, mat61 = [1,2,3,4,5,6]T\n");
    Printf("vec - mat61 =\n");
    (v1 - m61).Print();
    Println;
}

void Vec6ArithmeticProduct()
{
    PrintHeading("Vector6 Arithmetic operators - Product ");

    float v[6] = {1, 2, 3, 4, 5, 6};
    dtVector6<> v1(1, 2, 3, 4, 5, 6);
    dtVector6<> v2(v, sizeof(v));
    dtVector<6> v61(v, sizeof(v));
    dtMatrix<6, 1> m61(v, sizeof(v));
    dtMatrix<1, 6> m16(v, sizeof(v));

    Printf("/* Operator: *(scalar) */\n");
    Printf("vec = [1,1,1,1,1,1]T\n");
    Printf("vec * 0.1 =\n");
    v1.SetFill(1);
    (v1 * 0.1f).Print();
    Println;

    Printf("/* Operator: /(scalar) */\n");
    Printf("vec = [1,1,1,1,1,1]T\n");
    Printf("vec / 10 =\n");
    v1.SetFill(1);
    (v1 / 10).Print();
    Println;

    Printf("/* Operator: *(scalar, dtVector) */\n");
    Printf("vec = [1,1,1,1,1,1]T\n");
    Printf("0.1 * vec =\n");
    v1.SetFill(1);
    (0.1f * v1).Print();
    Println;

    Printf("/* Operator: *(dtMatrix<1, col>) */\n");
    Printf("vec = [1,1,1,1,1,1]T, mat16 = [1,1,1,1,1,1]\n");
    Printf("vec * mat16 =\n");
    v1.SetFill(1);
    m16.SetFill(1);
    (v1 * m16).Print();
    Println;

    Printf("/* Operator: dot(dtVector6) */\n");
    Printf("vec1 = [1,1,1,1,1,1]T, vec2 = [1,1,1,1,1,1]T\n");
    v1.SetFill(1);
    v2.SetFill(1);
    Printf("vec1 dot vec2 = %f\n", v1.dot(v2));
    Println;

    Printf("/* Operator: dot(dtVector) */\n");
    Printf("vec1 = [1,1,1,1,1,1]T, vec61 = [1,1,1,1,1,1]T\n");
    v1.SetFill(1);
    v61.SetFill(1);
    Printf("vec1 dot vec61 = %f\n", v1.dot(v61));
    Println;

    Printf("/* Operator: dot(dtMatrix) */\n");
    Printf("vec1 = [1,1,1,1,1,1]T, mat61 = [1,1,1,1,1,1]T\n");
    v1.SetFill(1);
    m61.SetFill(1);
    Printf("vec1 dot mat61 = %f\n", v1.dot(m61));
    Println;
}

void Vec6ComparisonOperator()
{
    PrintHeading("Vector6 Comparison operators ");

    float v[6] = {1, 2, 3, 4, 5, 6};
    dtVector6<> v1(1, 2, 3, 4, 5, 6);
    dtVector<6> v61(v, sizeof(v));
    dtMatrix<6, 1> m61(v, sizeof(v));

    Printf("/* Operator: ==(dtVector6 &v) and !=(dtVector6 &v) */\n");
    if (v1 == v1) Printf("true");
    else Printf("false");
    Println;

    if (v1 != v1) Printf("false");
    else Printf("true");
    Println;

    Printf("/* Operator: ==(dtVector &v) and !=(dtVector &v) */\n");
    if (v1 == v61) Printf("true");
    else Printf("false");
    Println;

    if (v1 != v61) Printf("false");
    else Printf("true");
    Println;

    Printf("/* Operator: ==(dtMatrix &v) and !=(dtMatrix &v) */\n");
    if (v1 == m61) Printf("true");
    else Printf("false");
    Println;

    if (v1 != m61) Printf("false");
    else Printf("true");
    Println;
    Println;
}