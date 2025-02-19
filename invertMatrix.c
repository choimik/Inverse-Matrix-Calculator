#include <stdio.h> //C standard library for standard input and output like printf
#include <stdlib.h> //C standard library gives access to srand() to generate random numbers
#include <cs50.h> //Harvard's library with basic input functions such as get_double, get_int and more
#include <time.h> //C standard library which gives access to time (used in this program to generate random numbers based off of time)
#include <math.h> //C standard library which gives access to function like pow

//below are prototypes for C functions (C requires these functions to all be listed in the beginning)

void copyMatrix(int size, double matrix[size][size], double copy[size][size]);
void fillMatrixValues(int row, int column, int matrix[row][column], int choice);
void printMatrix(int row, int column, int matrix[row][column]);
void printDoubleMatrix(int rows, int columns, double matrix[rows][columns], int choice);
double calculateDeterminant(int matrixSize, double matrix[matrixSize][matrixSize], double determinantArray[matrixSize]);
void swapRows(int matrixSize, int original_row, int replacement_row, double matrix[matrixSize][matrixSize]);
void addRows(int matrixSize, int row_one, int row_two, double matrix[matrixSize][matrixSize]);
void scalarMultiplyRow(int matrixSize, int row, double scalar, double matrix[matrixSize][matrixSize]);
void fillMinorValues(int matrixSize, int minor_row, int minor_col, double minorMatrix[matrixSize-1][matrixSize-1], double matrix[matrixSize][matrixSize]);
void fillCofactorValues(int determinant, int matrixSize, double adjointMatrix[matrixSize][matrixSize], double matrix[matrixSize][matrixSize]);
void transposeCofactorValues(int matrixSize, double adjointMatrix[matrixSize][matrixSize]);
void divideAdjointMatrixByDeterminant(int matrixSize, int determinant, double adjointMatrix[matrixSize][matrixSize]);

//main function

int main(void){

    //setting up random number generator using seed (seed or value which random number is based off of in this case is just the time currently)
    long t;
    srand((unsigned) time(&t));


    printf("\nWelcome to my program that solves for an inverse matrix for an inputted square matrix. This program is intended to both be used as a calculator to solve for an inverse matrix if a solution is possible and also teach the basics of matrices and their inverses to students.");

    //use get_int to get 0 or other integer so that there can be a calculator teaching/explaining mode and without teaching mode
    int option = get_int("\n\nIf you would like a more in-depth guide about the solving process and to learn about the basics of matrices and their inverses... press '0' otherwise press any other number to simply get the answer: ");

    //use conditional to determine if in calculator teaching mode or not (code in if statement will just print out basic information on matrices)
    if(option == 0){

        //create four example matrices to explain when in teaching mode and fill them with values

        int ex_matrix_one[2][3];
        fillMatrixValues(2, 3, ex_matrix_one, 0);

        int ex_matrix_two[2][2];
        fillMatrixValues(2, 2, ex_matrix_two, 0);

        int ex_identity_matrix_one[3][3];
        fillMatrixValues(3, 3, ex_identity_matrix_one, 1);

        int ex_identity_matrix_two[5][5];
        fillMatrixValues(5, 5, ex_identity_matrix_two, 1);

        //explain what matrices are and teach some basic concepts

        printf("\n\n\nYou pressed '0' for the educational option.");
        printf("\n\nA matrix is a two dimensional array with rows and columns. A matrix can have as many rows and columns as long as the amount of rows or columns are greater than or equal to one.");

        printf("\n\nFor example a 2x3 matrix with two rows and three columns would look like this... \n");

        printMatrix(2,3, ex_matrix_one);

        printf("\nNotice a matrix has individual values corresponding to every row and column pair. These values are called the elements of a matrix and they can be negative or positive rational numbers (including decimals and fractions).\n");

        printf("\n\nOn the other hand a 2x2 matrix with two rows columns would look like this... \n");

        printMatrix(2,2, ex_matrix_two);
        printf("\nNotice this matrix has the same number of rows and columns. These types of matrices are called square matrices. A specific property that only square matrices share is that these matrices have a special matrix linked to them known as the inverse matrix. When a square matrix is multiplied by its corresponding inverse matrix it results in the identity matrix.");

        printf("\n\nHere are two examples of the identity matrix of different sizes. ");
        printf("Notice the identity matrix has a diagonal of '1's through the middle.\n");
        printMatrix(3,3, ex_identity_matrix_one);
        printMatrix(5,5, ex_identity_matrix_two);

        printf("\nThis calculator attempts to find the inverse matrix, the unique corresponding matrix that when multiplied to a square matrix will yield the inverse matrix\n");
    }
    else{
        printf("\nYou selected no explanation and just the calculation of the inverse matrix.\n");
    }

    //get user input for matrix size and element values for the matrix to be created

    printf("\n\nThe program takes in two main parameters: a matrix size and a matrix with values filled in for every element.\n");
    int matrixSize = get_int("\nPlease Enter Matrix Size: ");
    while(matrixSize <= 0){
        matrixSize = get_int("Please Enter Matrix Size (input a number larger than 0): ");
    }
    printf("\n");

    double matrix[matrixSize][matrixSize];
    double temp;

    printf("Please Enter Element Values for Each Row and Column Position.\n\n");

    //repeatedly ask for element values for every row and column pair

    for(int i = 0; i < matrixSize; i++){
        for(int j = 0; j < matrixSize; j++){
            temp = get_double("Row %i and Column %i: ", i+1, j+1);
            matrix[i][j] = temp;
        }
        printf("\n");
    }

    printf("\n\nThe Inverse Matrix is defined as (1/determinant) * adjoint-matrix\n\n");
    printf("Original Square Matrix:\n");
    printDoubleMatrix(matrixSize, matrixSize, matrix, 1);

    // create matrix and variables needed to run the calculate determinant function
    double determinantArray[matrixSize];
    double determinant = calculateDeterminant(matrixSize, matrix, determinantArray);
    if (determinant == 0){
        printf("\nNon-invertible Matrix. Determinant equals zero and no inverse matrix exists.");
        return 1;
    }
    printf("\n\nDeterminant: %.5f\n\n", determinant);
    if(option == 0){
        printf("\nThe determinant is a scalar (where a scalar is just a one dimensional value) or number associated to a square matrix and its inverse.\n\n\n");
    }

    double adjointMatrix[matrixSize][matrixSize];
    if(option == 0){
        printf("\nThe Adjoint Matrix is a square matrix strongly linked with the inverse matrix. It consists of the transposition (reversal of the rows and columns) of its cofactors (the determinants of smaller matrices within the larger original square matrix) of a square matrix.\n");
    }
    printf("\nAdjoint Matrix filled with Cofactor Values (no transposition yet):\n");
    //call function to fill Adjoint Matrix with cofactor values
    fillCofactorValues(determinant, matrixSize, adjointMatrix, matrix);
    printDoubleMatrix(matrixSize, matrixSize, adjointMatrix, 1);
    printf("\n");
    if(option == 0){
        printf("\n\nWithin this step, the calculator fills the Adjoint Matrix with its smaller matrix cofactors. It will calculate the determinant of every smaller submatrix within the larger matrix. If a matrix is of a size nxn then its minors or submatrices are of a size of (n-1)x(n-1).These minors or submatrices and their determinants are used to find the cofactor values of the matrix.\n\n");
    }

    printf("\nComplete Adjoint Matrix After Transposition:\n");
    //call function to transpose/switch order of Adjoint Matrix
    transposeCofactorValues(matrixSize, adjointMatrix);
    printDoubleMatrix(matrixSize, matrixSize, adjointMatrix, 1);
    printf("\n");
    if(option == 0){
        printf("\nWithin this step, the calculator 'shuffles' or transposes the order of the cofactors. It will turn every row of the adjoint matrix into a column and every column into a row. This is the final step of finding the adjoint matrix.\n");
    }

    double inverseMatrix[matrixSize][matrixSize];
    //copies
    copyMatrix(matrixSize, adjointMatrix, inverseMatrix);

    printf("\nInverse Matrix (The Result of Dividing Adjoint Matrix by Determinant):\n");
    divideAdjointMatrixByDeterminant(matrixSize, determinant, inverseMatrix);
    printDoubleMatrix(matrixSize, matrixSize, inverseMatrix, 0);
    if(option == 0){
        printf("\nHaving determined the adjoint and determinant of the square matrix, the final step to find the inverse matrix is to divide every element within the adjoint matrix by the determinant. The calculator above loops through every row and column of the matrix and then divides every element by the determinant to finally get the inverse matrix.\n");
    }
    printf("\n");

    return 0;
}

void printMatrix(int rows, int columns, int matrix[rows][columns]){
    // loops through a function and prints a matrix with a type of integer

    printf("\n");
    for (int i = 0; i < rows; i++){
        printf("|");
        for (int j = 0; j < columns; j++){
            printf("%4i ", matrix[i][j]);
        }
        printf("  |");
        printf("\n");
    }
}

void printDoubleMatrix(int rows, int columns, double matrix[rows][columns], int choice){
    // loops through a matrix of data type double float and prints it out onto terminal
    //creates two choices with choice parameter... choice = 0 will print more decimal values which is especially useful when printing the inverse matrix
    //choice = 1 will print less decimal values and create more padding between each of the matrix' elements
    if (choice == 0){
        printf("\n");

        for (int i = 0; i < rows; i++){
            printf("|");
            for (int j = 0; j < columns; j++){
                printf("  %4.9f  ", matrix[i][j]);
            }
            printf("  |");
            printf("\n");
        }
    }
    else if (choice == 1){
        printf("\n");
        for (int i = 0; i < rows; i++){
            printf("|");
            for (int j = 0; j < columns; j++){
                printf("%10.2f ", matrix[i][j]);
            }
            printf("     |");
            printf("\n");
        }
    }
}

void fillMatrixValues(int row, int column, int matrix[row][column], int choice){

    if(choice == 0){
        for(int i = 0; i < row; i++){
            for (int j = 0; j < column; j++){
                matrix[i][j] = (rand() % 9)-4;
            }
        }
    }

    if(choice == 1){
        for(int i = 0; i < row; i++){
            for (int j = 0; j < column; j++){
                if(i == j){
                    matrix[i][j] = 1;
                }
                else{
                    matrix[i][j] = 0;
                }
            }
        }
    }
}

void copyMatrix(int size, double matrix[size][size], double copy[size][size]){
    //function above just copies the contents of one matrix into another
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            copy[i][j] = matrix[i][j];
        }
    }
}

void swapRows(int matrixSize, int original_row, int replacement_row, double matrix[matrixSize][matrixSize]){
    // matrix operation which swaps two rows positions
    double copy[matrixSize][matrixSize];
    copyMatrix(matrixSize, matrix, copy);

    for(int col = 0; col < matrixSize; col++){
        matrix[original_row][col] = copy[replacement_row][col];
        matrix[replacement_row][col] = copy[original_row][col];
    }
}

void addRows(int matrixSize, int row_one, int row_two, double matrix[matrixSize][matrixSize]){
    // adds two rows within a matrix specifically row_one to row_two

    for(int j = 0; j < matrixSize; j++){
        matrix[row_two][j] = matrix[row_one][j] + matrix[row_two][j];
    }

}

void scalarMultiplyRow(int matrixSize, int row, double scalar, double matrix[matrixSize][matrixSize]){
    //multiply row with scalar value (scalar is just a value that can be represented by a single number like height)

    for(int j = 0; j < matrixSize; j++){
        matrix[row][j] = matrix[row][j]*scalar;
    }

}

double calculateDeterminant(int matrixSize, double matrix[matrixSize][matrixSize], double determinantArray[matrixSize]){
    // calculates the determinant of a given matrix through row echelon form (matrix operations like add, swap and multiply)

    //creates a copy of the original matrix so that there is a copy of the original matrix without the row transformations that will occur later
    double copy_matrix[matrixSize][matrixSize];
    copyMatrix(matrixSize, matrix, copy_matrix);

    int index = 0;
    double scalar;
    int sign = 1;

    //this program uses row echelon form along with row operations to find the determinant
    //the determinant is the product of all the diagonal elements when the matrix is in row echelon form

    // | 3 8 9 |.  the matrix on the left is in row echelon form so as long as the bottom left triangle is made up of all zeros
    // | 0 2 6 |.  and the diagonal consists of nonzero numbers the determinant will be the product of the diagonal which in this case
    // | 0 0 4 |.  is 3*2*4*(any scalar values factored out during solving)

    //iterate through every element of 2d array / matrix
    for(int col = 0; col < matrixSize; col++){

        if(copy_matrix[col][col] == 0){
        //search for row that does not have zero in the middle diagonal and then swap rows of the matrix if there is a zero there (if there is a zero in diagonal cannot work)
        for(int i = col+1; i < matrixSize; i++){
            if(copy_matrix[i][col] != 0){
                swapRows(matrixSize, col, i, copy_matrix);
                sign = (-1)*sign;
                //changes the sign of the determinant
                break;
                }
            }
        }
        //set the first value of the determinant array to the current diagonal value
        //(eventually at the end of iterating through the whole matrix there will be an array with all the diagonal values within the matrix that can just be multiplied together at the end)
        determinantArray[index] = copy_matrix[col][col];
        //multiply the row so that the diagonal value is set to 1 after having saved the data into an array
        scalarMultiplyRow(matrixSize, col, 1.0/determinantArray[index], copy_matrix);
        index++;
        for(int row = 1+col; row < matrixSize; row++){
            //in this part of the loop the algorithm will attempt to create the bottom triangle of zeros needed to be in row echelon form
            //it will do this by subtracting rows from each other so that the element that it is currently positioned on becomes zero
            //also an important fact is that it will only iterate through the bottom triangle of zeros as row is set to 1+col

            //set the scalar to the value that must be subtracted out through row addition
            scalar = copy_matrix[row][col];

            //since scalar could be zero set a condition that will skip the element if the value is already zero
            if(scalar == 0 && row != col){
                continue;
            }
            //if the scalar is zero and also along the diagonal the product of the diagonal will always be zero so just return zero as the determinant
            else if(scalar == 0 && row == col){
                printDoubleMatrix(matrixSize, matrixSize, copy_matrix, 1);
                return 0;
            }
            //multiply row by -scalar value (prior to this it is already 1) so when rows are added the specific element will be turned into zero
            scalarMultiplyRow(matrixSize, col, -scalar, copy_matrix);
            //add rows together leading to the one of the bottom triangle values turning into zero
            addRows(matrixSize, col, row, copy_matrix);
            //revert the original row which was multiplied by -scalar back to its original values so that it can be used again to subtract a different value if there is another elements that needs to be turned into a zero
            scalarMultiplyRow(matrixSize, col, (-1.0/scalar), copy_matrix);

        }
    }

    double determinant = 1;

    //multiply all the diagonal values together and return the product
    for (int i = 0; i < matrixSize; i++){
        determinant = determinant*determinantArray[i];
    }
    determinant = determinant*sign;
    return determinant;
}

void fillMinorValues(int matrixSize, int minor_row, int minor_col, double minorMatrix[matrixSize-1][matrixSize-1], double matrix[matrixSize][matrixSize]){
    // identifies a small sub matrix from a larger matrix (matrix) and then fills those values into the minorMatrix

    //a minor is a sub matrix that is exactly one size smaller than the original matrix and is denoted with two letters i and j which are the row and column taken away to create the submatrix
    //
    // this is an example of a matrix... the minor (1, 1) is the submatrix created when the first row and column are removed
    // in this example the minor (1, 1) is a 2x2 matrix seen to the right

    // | 2 5 8 |.
    // | 3 7 9 |.    | 7 9 |
    // | 1 0 3 |.    | 0 3 |

    //This function will simply fill a submatrix denoted by minorMatrix with its minor values so that it can be used later in a different function

    int row_tracker = 0;
    int col_tracker = 0;

    //this will iterate through the adjointMatrix without having to create another loop using the the variables of row_tracker and col_tracker

    for(int row = 0; row < matrixSize; row++){
        if(row == minor_row){
            continue;
        }
        for(int col = 0; col < matrixSize; col++){
            if(col == minor_col){
                continue;
            }
            else{
                minorMatrix[row_tracker][col_tracker] = matrix[row][col];
                //the conditional below will reset the columns to 0 and add 1 to the row once it reaches the minorSize
                if(col_tracker == matrixSize-2){
                    col_tracker = 0;
                    row_tracker++;
                }
                else{
                    col_tracker++;
                }
            }
        }
    }
}

void fillCofactorValues(int determinant, int matrixSize, double adjointMatrix[matrixSize][matrixSize], double matrix[matrixSize][matrixSize]){
    // calculates the determinant of minor (smaller sub matrix) and adds a positive or negative sign

    //a minor is a sub matrix that is exactly one size smaller than the original matrix and is denoted with two letters i and j which are the row and column taken away to create the submatrix
    //
    // this is an example of a matrix... the minor (1, 1) is the submatrix created when the first row and column are removed
    // in this example the minor (1, 1) is a 2x2 matrix seen to the right

    // | 2 5 8 |.
    // | 3 7 9 |.    | 7 9 |
    // | 1 0 3 |.    | 0 3 |

    //the determinant of this submatrix multiplied by a sign (-1 or 1) is known as the cofactor which is necessary to find the adjoint matrix

    //create variables necessary to call the calculateDeterminant function
    int minorSize = matrixSize - 1;
    double minorDeterminantArray[matrixSize];
    double minorDeterminant;

    double minorMatrix[minorSize][minorSize];

    int row_tracker = 0;
    int col_tracker = 0;

    //find cofactor values by calling fillMinorValues function and then calling the calculateDeterminant function on the submatrix to find the cofactor
    //this will iterate through the entire matrix and find the corresponding cofactor values for all of the possible minors and elements within the matrix
    for(int row = 0; row < matrixSize; row++){
        for(int col = 0; col < matrixSize; col++){
            //find new submatrix for every element
            fillMinorValues(matrixSize, row, col, minorMatrix, matrix);
            minorDeterminant = calculateDeterminant(minorSize, minorMatrix, minorDeterminantArray);
            //multiply by sign (-1 or 1) depending on position of row and col
            adjointMatrix[row_tracker][col_tracker] = pow(-1, (row+1+col+1)) * minorDeterminant;
            //this will iterate through the adjointMatrix without having to create another loop using the the variables of row_tracker and col_tracker
            //the conditional below will reset the columns to 0 and add 1 to the row once it reaches the minorSize
            if(col_tracker == matrixSize - 1){
                col_tracker = 0;
                row_tracker++;
            }
            else{
                col_tracker++;
            }
        }
    }


}

void transposeCofactorValues(int matrixSize, double adjointMatrix[matrixSize][matrixSize]){
    //changes order of rows and columns by making all columns into rows and all rows into columns

    //first creates a copy of the adjointMatrix array so that order can be reversed later
    double copy[matrixSize][matrixSize];
    for(int row = 0; row < matrixSize; row++){
        for(int col = 0; col < matrixSize; col++){
            copy[row][col] = adjointMatrix[row][col];
        }
    }

    //switches the rows for columns and columns for rows in the matrix
    for(int row = 0; row < matrixSize; row++){
        for(int col = 0; col < matrixSize; col++){
            adjointMatrix[row][col] = copy[col][row];
        }
    }
}

void divideAdjointMatrixByDeterminant(int matrixSize, int determinant, double adjointMatrix[matrixSize][matrixSize]){
    // divides all elements of a matrix by a number (specifically adjoint matrix by determinant)

    for(int row = 0; row < matrixSize; row++){
        for(int col = 0; col < matrixSize; col++){
            adjointMatrix[row][col] = adjointMatrix[row][col]/determinant;
        }
    }
}
