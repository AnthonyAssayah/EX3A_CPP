#include <iostream>
#include <stdexcept>
#include "Matrix.hpp"
#include <vector>
using namespace std;


namespace zich {
    // Constructor
    Matrix::Matrix(vector<double> mat, const size_t &my_rows,const size_t &my_cols) {

        if ( my_rows <= 0 || my_cols <= 0) {
        throw invalid_argument("Error, row and col must be positive numbers");
        }

        this->rows = my_rows;
        this->cols = my_cols;
        
        this->mat=new double* [my_rows];

        for(int i =0; i<my_rows; i++) {
            this->mat[i]=new double[my_cols]; 
        }
        size_t k =0;
        for(int i=0;i<my_rows;i++) {
            for(int j=0;j<my_cols;j++) {
                this->mat[i][j] = mat[k++];
            }
        }
    }

    // Constructor
    Matrix::Matrix(const size_t &my_rows,const size_t &my_cols) {

        if ( my_rows <= 0 || my_cols <= 0) {
            throw invalid_argument("Error, row and col must be positive numbers");
        }

        this->rows = my_rows;
        this->cols = my_cols;

        this->mat=new double* [my_rows];
        
        for(int i =0;i<my_rows;i++)
        {
            this->mat[i]=new double[my_cols]; 
        
        }

        size_t k =0;
        for(size_t i=0;i<my_rows;i++) {
            for(size_t j=0;j<my_cols;j++) {
                this->mat[i][j] +=1;
            }
        }
    }

    // Copy constructor
    Matrix::Matrix(const Matrix &other){ 
        rows = other.rows;
        cols = other.cols;

        mat = new double* [other.rows];
        
        for(int i =0;i<other.rows;i++)
        {
        mat[i]=new double[other.cols];
        
        }
        for(int i=0;i<other.rows;i++) {
        for(int j=0;j<other.cols;j++) {
                mat[i][j] = other.mat[i][j];
            }
        }

    }

    // Desconstructor
    Matrix::~Matrix(){

        for(int i =0;i<rows;i++) {
            delete [] mat[i]; 
        }
        delete [] mat;
    }


    ///////// 1 - Addition / Soustraction operations /////////////

    // Addition of two matrix  - Return a new matrix
    Matrix Matrix::operator+(const Matrix& my_mat) {

        if ( this->rows != my_mat.rows || this->cols != my_mat.cols) {
            throw invalid_argument("Error, the two matrix must have the same size");
        }

        // Create a new matrix for addition 
        Matrix plus_mat(my_mat.rows, my_mat.cols);

        for(int i=0; i<plus_mat.cols ;i++) {
            for(int j=0; j<plus_mat.cols ;j++) {
                plus_mat.mat[i][j] = this->mat[i][j] + my_mat.mat[i][j];
            }
        }

        return plus_mat;
    }

    // Soustraction of two matrix - Return a new matrix
    Matrix Matrix::operator-(const Matrix& my_mat) {

        if ( this->rows != my_mat.rows || this->cols != my_mat.cols) {
            throw invalid_argument("Error, the two matrix must have the same size");
        }

        // Create a new matrix for soustraction 
        Matrix minus_mat(my_mat.rows, my_mat.cols);
        
        for(int i=0; i<minus_mat.cols ;i++) {
            for(int j=0; j<minus_mat.cols ;j++) {
                minus_mat.mat[i][j] = this->mat[i][j] - my_mat.mat[i][j];
            }
        }
        return minus_mat;
    }

    // Plus/equal operator
    Matrix& Matrix::operator+=(const Matrix &my_mat) {                         ///maybe Matrix& my_mat

        if ( this->rows != my_mat.rows || this->cols != my_mat.cols) {
            throw invalid_argument("Error, the two matrix must have the same size");
        }

        for (int i = 0; i < this->rows; i++) {
            for ( int j = 0; j < this->cols; j++) {
                this->mat[i][j] += my_mat.mat[i][j];
            }
        }
        return *this;
    }

    // Minus/equal operator
    Matrix &Matrix::operator-=(const Matrix &my_mat) {

        if ( this->rows != my_mat.rows || this->cols != my_mat.cols) {
            throw invalid_argument("Error, the two matrix must have the same size");
        }

    
        for (int i = 0; i < this->rows; i++) {
            for ( int j = 0; j < this->cols; j++) {
                this->mat[i][j] -= my_mat.mat[i][j];
            }
        }
        return *this;
    }

    // Plus unary operator
    Matrix Matrix::operator+() {

        Matrix plusUnary(*this);
        return plusUnary;

    }

    // Minus unary operator
    Matrix Matrix::operator-(){

        
        // Create a new matrix for soustraction 
        Matrix temp_mat(this->rows, this->cols);

        for (int i =0; i < temp_mat.rows; i++) {
            for (int j =0; j < temp_mat.cols; j++) {
                if (this->mat[i][j] == 0) {
                    temp_mat.mat[i][j] = 0;
                }
                else {
                    temp_mat.mat[i][j] = this->mat[i][j]*(-1);           /// check for negative value
                }
            }
        }
        return temp_mat;
    }

 

    //////////  2 - Compare operation   /////////////

    bool Matrix::operator==(const Matrix& my_mat) const {     /// const ??

        if ( this->rows != my_mat.rows || this->cols != my_mat.cols) {
            throw invalid_argument("Error, the two matrix must have the same size");
        }

        for (int i = 0; i < this->rows; i++) {
            for (int j = 0; j < this->cols; j++) {
                if ( this->mat[i][j] != my_mat.mat[i][j]) {
                    return false;     
                    }
                }
            }
            return true;
    }

    bool Matrix::operator!=(const Matrix& my_mat) const{

        return !this->operator==(my_mat);
    }

    bool Matrix::operator>(const Matrix& my_mat) const {

        if ( this->rows != my_mat.rows || this->cols != my_mat.cols) {
            throw invalid_argument("Error, the two matrix must have the same size");
        }

        double sum1 =0;
        double sum2 =0;

        for (int i = 0; i < this->rows; i++) {
            for (int j = 0; j < this->cols; j++) {
                sum1 += this->mat[i][j];
            }
        }

        for (int i = 0; i < my_mat.rows; i++) {
            for (int j = 0; j < my_mat.cols; j++) {
                sum2 += my_mat.mat[i][j];
            }
        }

        if (sum1 > sum2) {
            return true;
        }
        return false;

    }

    bool Matrix::operator>=(const Matrix& my_mat) const {

        if ( this->rows != my_mat.rows || this->cols != my_mat.cols) {
            throw invalid_argument("Error, the two matrix must have the same size");
        }

        double sum1 =0;
        double sum2 =0;

        for (int i = 0; i < this->rows; i++) {
            for (int j = 0; j < this->cols; j++) {
                sum1 += this->mat[i][j];
            }
        }

        for (int i = 0; i < my_mat.rows; i++) {
            for (int j = 0; j < my_mat.cols; j++) {
                sum2 += my_mat.mat[i][j];
            }
        }

        if (sum1 >= sum2) {
            return true;
        }
        return false;

    }

    bool Matrix::operator<(const Matrix& my_mat) const {

        if ( this->rows != my_mat.rows || this->cols != my_mat.cols) {
            throw invalid_argument("Error, the two must matrix have the same size");
        }

        double sum1 =0;
        double sum2 =0;

        for (int i = 0; i < this->rows; i++) {
            for (int j = 0; j < this->cols; j++) {
                sum1 += this->mat[i][j];
            }
        }

        for (int i = 0; i < my_mat.rows; i++) {
            for (int j = 0; j < my_mat.cols; j++) {
                sum2 += my_mat.mat[i][j];
            }
        }

        if (sum1 < sum2) {
            return true;
        }
        return false;

    }

    bool Matrix::operator<=(const Matrix& my_mat) const{

        if ( this->rows != my_mat.rows || this->cols != my_mat.cols) {
            throw invalid_argument("Error, the two matrix must have the same size");
        }

        double sum1 =0;
        double sum2 =0;

        for (int i = 0; i < this->rows; i++) {
            for (int j = 0; j < this->cols; j++) {
                sum1 += this->mat[i][j];
            }
        }

        for (int i = 0; i < my_mat.rows; i++) {
            for (int j = 0; j < my_mat.cols; j++) {
                sum2 += my_mat.mat[i][j];
            }
        }

        if (sum1 <= sum2) {
            return true;
        }
        return false;

    }


    ///////////  3 - Increment / Decrement operations /////////////

    // Prefix increment
    Matrix &Matrix::operator++() {

        for (int i = 0; i < this->rows; i++) {
            for ( int j = 0; j < this->cols; j++) {
                this->mat[i][j]++;
            }
        }
        return *this;
    }

    // Postfix increment
    Matrix Matrix::operator++(int) {
        Matrix mat_copy(*this);
        for (int i = 0; i < this->rows; i++) {
            for ( int j = 0; j < this->cols; j++) {
                this->mat[i][j]++;
            }
        }
        return mat_copy;
    }

    // Prefix decrement
    Matrix &Matrix::operator--() {

        for (int i = 0; i < this->rows; i++) {
            for ( int j = 0; j < this->cols; j++) {
                (this->mat[i][j])--;
            }
        }
        return *this;
    }

    // Postfix decrement
    Matrix Matrix::operator--(int) {
        Matrix mat_copy(*this);
        for (int i = 0; i < this->rows; i++) {
            for ( int j = 0; j < this->cols; j++) {
                this->mat[i][j]--;
            }
        }
        return mat_copy;
    }


    ///////////  4 - Multiplication operations between matrix and double /////////


    Matrix &Matrix::operator*=(const double num) {

        for ( int i = 0; i < this->rows; i++) {
            for ( int j = 0; j < this->cols; j++) {
                this->mat[i][j] *= num;
            }
        }
      return *this;
    }

    Matrix operator*(const double num, const Matrix &matrix) {

        Matrix mul_mat(matrix.rows, matrix.cols);

        for ( int i = 0; i < matrix.rows; i++) {
            for ( int j = 0; j < matrix.cols; j++) {
                mul_mat.mat[i][j] = num*matrix.mat[i][j];
            }
        }
        return mul_mat;
    }

     Matrix operator*(const Matrix &matrix, const double num) {

        Matrix mul_mat(matrix.rows, matrix.cols);

        for ( int i = 0; i < matrix.rows; i++) {
            for ( int j = 0; j < matrix.cols; j++) {
                mul_mat.mat[i][j] *= num;
            }
        }
        return mul_mat;
    }


    ///////////  5 -  Multiplication operation between two matrix /////////


    Matrix operator*(const Matrix &mat_1, const Matrix &mat_2) {

        if ( mat_1.cols != mat_2.rows) {
            throw invalid_argument(" Error, you cannot multiply those two matrix");
        }

        Matrix mul_mat = Matrix(mat_1.rows, mat_2.cols);

        for(int i  = 0; i < mat_1.rows; i++) {
            for(int j = 0; j < mat_2.cols; j++) {
                for(int k = 0; k < mat_1.cols; k++) {
                    mul_mat.mat[i][j] += mat_1.mat[i][k] * mat_2.mat[k][j];
                }
            }
        }
        return mul_mat;
    }

    Matrix Matrix::operator*=(const Matrix &my_mat) {

        if (this->cols != my_mat.rows) {
            throw invalid_argument(" Error, you cannot multiply those two matrix");
        }

        for(int i  = 0; i < this->rows; i++) {
            for(int j = 0; j < my_mat.cols; j++) {
                for(int k = 0; k < this->cols; k++) {
                    this->mat[i][j] = this->mat[i][k] * my_mat.mat[k][j];
                }
            }
        }
        return *this;    //// return 
    }


    /////////// 6 - Input / Output operations /////////

    ostream& operator<<(ostream& out, const Matrix &matrix) {

        for ( int i = 0; i < matrix.rows; i++) {
            out << "["; 
            for ( int j = 0; j < matrix.cols; j++) {
                if ( j != matrix.cols - 1) {
                    out << matrix.mat[i][j] << ' ';
                }
                else {
                    out << matrix.mat[i][j]; 
                }
            }
            out << ']';
            out << '\n';
            
        }
        return out;
    }


    istream& operator>>(istream& in, Matrix &matrix) {

    //     string input ="";
    //     for ( int i = 0; i < matrix.rows; i++) {
    //         in >>'[';
    //         for (int j =0; j <matrix.cols; j++) {

    //         }
    // }

        return in;
    }
}