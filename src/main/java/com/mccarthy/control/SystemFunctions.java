package com.mccarthy.control;

import org.ejml.simple.SimpleMatrix;

class SystemFunctions {
    private SystemFunctions() {
    }

    protected static void isSystemValid(SimpleMatrix A, SimpleMatrix B, SimpleMatrix Q, SimpleMatrix R) {
        if (A.numCols() != Q.numCols() || A.numRows() != Q.numRows()) {
            throw new InvalidSystemException("Q and A must be of the same dimentions.");
        } else if (R.numCols() != R.numRows()) {
            throw new InvalidSystemException("R must be a square matrix.");
        } else if (B.numCols() != R.numRows()) {
            throw new InvalidSystemException("R must have the same number of rows as B's cols.");
        }
    }

    protected static SimpleMatrix HamiltonianMatrix(SimpleMatrix A, SimpleMatrix B, SimpleMatrix Q, SimpleMatrix R) {
        SimpleMatrix H12 = B.mult(R.invert().mult(B.transpose())).negative();
        SimpleMatrix H1 = A.combine(0, A.numCols(), H12);
        SimpleMatrix H2 = Q.negative().combine(0, Q.numCols(), A.transpose().negative());
        SimpleMatrix H = H1.combine(H1.numRows(), 0, H2);
        return H;
    }

    protected static SimpleMatrix SymplecticMatrix(SimpleMatrix A, SimpleMatrix B, SimpleMatrix Q, SimpleMatrix R) {
        SimpleMatrix S11 = A.plus(B.mult(R.invert().mult(B.transpose().mult(A.invert().transpose().mult(Q)))));
        SimpleMatrix S12 = B.mult(R.invert().mult(B.transpose().mult(A.invert().transpose()))).negative();
        SimpleMatrix S21 = A.invert().transpose().mult(Q).negative();
        SimpleMatrix S22 = A.invert().transpose();
        SimpleMatrix S = S11.combine(0, S11.numCols(), S12);
        S = S.combine(S11.numRows(), 0, S21);
        S = S.combine(S11.numRows(), S11.numCols(), S22);
        return S;
    }
}