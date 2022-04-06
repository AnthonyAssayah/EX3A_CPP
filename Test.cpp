#include <iostream>
#include <vector>
#include "Matrix.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include "doctest.h"

using namespace std;
using namespace zich;


TEST_CASE("Constructor takingt vector and matrix"){

    vector<double> vect0 = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    vector<double> vect1 = {1, 2, 3, 1, 2, 3, 1, 2, 3};
    vector<double> vect2 = {2.5, 4.1, 3.0, 8.2};
    vector<double> vect3 = {3, 11, -5, 4.5, 3.2, 0.5, -10, 12};
    vector<double> vect4 = {-5.0, 0.5, -3.5, 2.2};
    vector<double> vect5 = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};

    CHECK_NOTHROW(Matrix mat1(vect0,3,3));
    CHECK_NOTHROW(Matrix mat2(vect2,2,2));
    CHECK_NOTHROW(Matrix mat3(vect5,3,4));
    CHECK_NOTHROW(Matrix mat4(vect1,3,3));
    CHECK_NOTHROW(Matrix mat5(vect3,2,4));
    CHECK_NOTHROW(Matrix mat6(vect4,2,2));

} 

TEST_CASE("Constructors with illegal arguments"){

    vector<double> vect0 = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    vector<double> vect1 = {1, 2, 3, 1, 2, 3, 1, 2, 3};
    vector<double> vect2 = {2.5, 4.1, 3.0, 8.2};
    vector<double> vect3 = {3, 11, -5, 4.5, 3.2, 0.5, -10, 12};
    vector<double> vect4 = {-5.0, 0.5, -3.5, 2.2};
    vector<double> vect5 = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};

    CHECK_THROWS(Matrix mat4(vect1,0,3));       // row 0
    CHECK_THROWS(Matrix mat5(vect3,2,0));       // col 0
    CHECK_THROWS(Matrix mat6(vect4,0,0));       // row and col 0

}

TEST_CASE( "Addition & Soustraction operators ") {

    vector<double> vect0 = {1, 0, 0, 0, 1, 0, 0, 0, 1};  
    vector<double> vect1 = {1, 2, 3, 1, 2, 3, 1, 2, 3};
    vector<double> vect2 = {2.5, 4.1, 3.0, 8.2};   
    vector<double> vect3 = {3, 11, -5, 4.5, 3.2, 0.5, -10, 12, 9.2};
    vector<double> vect4 = {-5.0, 0.5, -3.5, 2.2}; 
    vector<double> vect5 = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};

    // Construct matrix
    Matrix mat1(vect0,3,3);
    Matrix mat2(vect2,2,2);
    Matrix mat3(vect5,3,4);
    Matrix mat4(vect1,3,3);
    Matrix mat5(vect3,3,3);
    Matrix mat6(vect4,2,2);
   

    // new vectors after addition between two matrix
    vector<double> vect_plus1 = {2, 2, 3, 1, 3, 3, 1, 2, 4};
    vector<double> vect_plus2 = {-2,5, 4.6, -0.5, 10.4}; 
    vector<double> vect_plus3 = {2, 0, 3, 1, -6, 3, 10, 2, 4};
    vector<double> vect_plus4 = {2,5, 5.5, -6.5, -8.4};
    

    // Two matrix plus operator (a + b)

    CHECK((mat1+mat4) == Matrix(vect_plus1,3,3));           // vect0 + vect1  ==> mat1 + mat4 = vect_add1
    CHECK_NE(mat1+mat4, Matrix(vect_plus3,3,3));            // vect0 + vect1  ==> mat1 + mat4 != vect_add3
    CHECK_NE(mat2+mat6, Matrix(vect_plus4,2,2));            // vect2 + vect4  ==> mat2 + mat6 != vect_add4
    
    // new vectors after soustraction between two matrix
    vector<double> vect_minus1 = {0, -2, -3, -1, -1, -3, -1, -2, -2};
    vector<double> vect_minus2 = {7.5, 3.6, 6.5, 6};  
    vector<double> vect_minus3 = {2.2, -6, 1.2, 0.3, -6.3, 30, 11, 2.0, 4.1};
    vector<double> vect_minus4 = {2,5, 5.5, -6.5, -8.4};

    // Two matrix minus operator (a - b)
    CHECK_EQ(mat1-mat4, Matrix(vect_minus1,3,3));     // vect0 - vect1  ==> mat1 - mat4 = vect_minus1
    CHECK_NE(mat1-mat4, Matrix(vect_minus3,3,3));     // vect0 - vect1  ==> mat1 - mat4 != vect_minus3
    CHECK_NE(mat2-mat6, Matrix(vect_minus4,2,2));     // vect2 - vect4  ==> mat2 - mat6 != vect_minus4

   
     // Minus Unary operator (-)
    CHECK_EQ(-mat1, Matrix({-1, 0, 0, 0, -1, 0, 0, 0, -1},3,3));
    CHECK_EQ(-mat4, Matrix({-1, -2, -3, -1, -2, -3, -1, -2, -3},3,3));
    CHECK_EQ(-mat2, Matrix({-2.5, -4.1, -3.0, -8.2},2,2));

     // Plus Unary operator (+)
    CHECK_EQ(+mat5, Matrix({3, 11, -5, 4.5, 3.2, 0.5, -10, 12, 9.2},3,3));  
    CHECK_EQ(+mat6, Matrix({-5.0, 0.5, -3.5, 2.2},2,2));
    CHECK_EQ(+mat3, Matrix({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12},3,4));

   
     // Plus equal operator (+=) 

    Matrix mat7({-1, 0, 0, 0, -1, 0, 0, 0, -1},3,3);
    Matrix mat8({1, 0, 0, 0, 1, 0, 0, 0, 1},3,3);
    

    CHECK((mat1+=mat4) == Matrix({2, 2, 3, 1, 3, 3, 1, 2, 4},3,3));     
    CHECK((mat8+=mat7) == Matrix({0, 0, 0, 0, 0, 0, 0, 0, 0},3,3));    

    // Minus equal operator (-=)


    Matrix mat9({1, 0, 0, 0, 1, 0, 0, 0, 1},3,3);
    Matrix mat10({1, 2, 3, 1, 2, 3, 1, 2, 3},3,3);
    Matrix mat11({1, 0, 0, 0, 1, 0, 0, 0, 1},3,3);
 
    CHECK((mat9-=mat7) == Matrix({2, 0, 0, 0, 2, 0, 0, 0, 2},3,3)); 
    CHECK((mat10-=mat11) == Matrix({0, 2, 3, 1, 1, 3, 1, 2, 2},3,3));
}


 TEST_CASE( "Increment & Decrement operators ") {

    Matrix mat1({1, 1, 1, 1, 1, 1, 1, 1, 1},3,3);
    Matrix mat2({2.5, 4.1, 3.0, 8.2},2,2);
    Matrix mat3({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12},3,4);

    //prefix increment operator (++)
    CHECK_EQ(++mat1, Matrix({2, 2, 2, 2, 2, 2, 2, 2, 2},3,3));   //CHECK_NOTHROW
    CHECK_EQ(++mat2, Matrix({3.5, 5.1, 4.0, 9.2},2,2));
    CHECK_EQ(++mat3, Matrix({2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13},3,4));

    //postfix incremement operator (++)
    mat1++;
    CHECK_EQ(mat1, Matrix({3, 3, 3, 3, 3, 3, 3, 3, 3},3,3)); 
    CHECK_EQ(mat1++, Matrix({3, 3, 3, 3, 3, 3, 3, 3, 3},3,3)); 
    CHECK_EQ(mat1,Matrix({4, 4, 4, 4, 4, 4, 4, 4, 4},3,3)); 
  
    //prefix decrement operator (--)
    CHECK_EQ(--mat1, Matrix({3, 3, 3, 3, 3, 3, 3, 3, 3},3,3));   //CHECK_NOTHROW  {2.5, 4.1, 3.0, 8.2}; 
    CHECK_EQ(--mat2, Matrix({2.5, 4.1, 3.0, 8.2},2,2));
    CHECK_EQ(--mat3, Matrix({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12},3,4));

    //postfix decrement operator (--)
    mat1--;
    CHECK_EQ(mat1, Matrix({2, 2, 2, 2, 2, 2, 2, 2, 2},3,3)); 
    CHECK_EQ(mat1--, Matrix({2, 2, 2, 2, 2, 2, 2, 2, 2},3,3)); 
    CHECK_EQ(mat1, Matrix({1, 1, 1, 1, 1, 1, 1, 1, 1},3,3));  
 }

 TEST_CASE( "Comparaison operators ") {

    // Construct matrix
    Matrix mat1({1, 0, 0, 0, 1, 0, 0, 0, 1},3,3);
    Matrix mat2({2.5, 4.1, 3.0, 8.2},2,2);
    Matrix mat3({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12},3,4);
    Matrix mat4({1, 2, 3, 1, 2, 3, 1, 2, 3},3,3);
    Matrix mat5({3, 11, -5, 4.5, 3.2, 0.5, -10, 12},2,4);
    Matrix mat6({-5.0, 0.5, -3.5, 2.2},2,2);
    Matrix mat7({1, 0, 0, 0, 1, 0, 0, 0, 1},3,3);
    Matrix mat8({2.5, 4.1, 3.0, 8.2},2,2);

    // == operator
    CHECK(mat1.operator==(mat7));
    CHECK_NOTHROW(mat2.operator==(mat8));
    CHECK_EQ(mat1,mat7);
    CHECK_EQ(mat2,mat8);


    // != operator
    CHECK_NOTHROW(mat1.operator!=(mat4));
    CHECK_NOTHROW(mat2.operator!=(mat6));
    CHECK_NE(mat1,mat4);
    CHECK_NE(mat2,mat6);

    // < and <= operator
    CHECK_NOTHROW(mat1.operator<(mat4));
    CHECK_NOTHROW(mat2.operator<=(mat8));
    CHECK_LT(mat1,mat4);
    CHECK_LE(mat2,mat8);

    // > and >= operator
    CHECK_NOTHROW(mat8.operator>(mat6));
    CHECK_NOTHROW(mat2.operator>=(mat8));
    CHECK_GT(mat8,mat6);
    CHECK_GE(mat2,mat8);
    
  
}

TEST_CASE( "Multiplication operators ") {
  
    Matrix mat1({1, 1, 1, 1, 1, 1, 1, 1, 1},3,3);
    Matrix mat2({2.5, 4.1, 3.0, 8.2},2,2);
    Matrix mat3({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12},3,4);
    Matrix mat4({2.5, 4.5, 6.5, 8.5},2,2);

    
    CHECK((mat4*=1.5) ==  Matrix({3.75, 6.75, 9.75, 12.75},2,2));
    CHECK((mat3*=0) == Matrix({0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},3,4));
}

TEST_CASE("Output operators") {

    ostringstream out;
    Matrix mat1({1, 1, 1, 1, 1, 1, 1, 1, 1},3,3);
    out << mat1;
    CHECK_EQ(out.str(), "[1 1 1]\n[1 1 1]\n[1 1 1]\n");

    ostringstream out1;
    Matrix mat2({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12},3,4);
    out1 << mat2;
    CHECK_EQ(out1.str(), "[1 2 3 4]\n[5 6 7 8]\n[9 10 11 12]\n");


}
