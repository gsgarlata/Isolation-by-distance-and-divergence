#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix MatrixMultiplication(NumericMatrix A_matrix, NumericMatrix B_matrix){

  int rows = B_matrix.nrow();
  int cols = B_matrix.ncol();

NumericMatrix product(rows,cols);

    for (unsigned int i = 0; i < rows; i++) {
        for (unsigned int j = 0; j < rows; j++) {
            product(i,j) = 0;
              for (unsigned int k = 0; k < rows; k++) {
                  product(i,j) += A_matrix(i,k) * B_matrix(k,j);
              }
        }
    }
return product;
}

// [[Rcpp::export]]
NumericMatrix IdentityMatrix(NumericMatrix Matrix){

const unsigned int rows = Matrix.nrow();
const unsigned int cols = Matrix.ncol();

NumericMatrix I_matrix(rows,cols);

for(unsigned int t = 0; t < rows; t++)
    I_matrix(t,t) = 1;

return  I_matrix;
}

// [[Rcpp::export]]
NumericMatrix PowerMatrix(NumericMatrix Matrix, int time_bit){

int num_rows =  Matrix.nrow();
int num_cols =  Matrix.ncol();

NumericMatrix R_matrix(num_rows,num_cols);


if (time_bit & 1){
        R_matrix = Matrix;
        time_bit &= ~1;
}else{
        R_matrix = IdentityMatrix(Matrix);
      }
        int i=1;
    NumericMatrix B_matrix = Matrix;                //B will always be M^i, where i is a power of 2
      while(time_bit){
          i *= 2;        //Advance i to the next power of 2
          B_matrix = MatrixMultiplication(B_matrix, B_matrix);       //B was M^(i/2) and is now M^i

          if (time_bit & i){       //i is of the form 2^j. Is the j-th bit of t set?
                R_matrix = MatrixMultiplication(R_matrix,B_matrix);      //Add B=A^i to the result
                time_bit = time_bit & ~i;   //Clear the j-th bit of t
          }

 }
          return R_matrix;
}
