package com.mccarthy.control;

import java.util.Comparator;
import java.util.List;

import org.ejml.data.Complex_F64;
import org.ejml.simple.SimpleEVD;
import org.ejml.simple.SimpleMatrix;

public class Control {
    private Control() {
    }

    public static Complex_F64[] pole(ss sys){
        return (Complex_F64[]) sys._A.eig().getEigenvalues().toArray();
    }

    public static StateSolution lqr(ss sys, SimpleMatrix Q, SimpleMatrix R) throws UnableToEvaluateStateSolution {
        return care(sys._A, sys._B, Q, R);
    }

    /**
     * Find the state space solution K, S and E
     * 
     * @param A - A matrix
     * @param B - B matrix
     * @param Q - Q matrix
     * @param R - R matrix
     * @return StateSolution class containing K gain matrix, S riccati solution and
     *         E eiganvalues of control system
     */
    public static StateSolution lqr(SimpleMatrix A, SimpleMatrix B, SimpleMatrix Q, SimpleMatrix R)
            throws UnableToEvaluateStateSolution {
        return care(A, B, Q, R);
    }

    public static StateSolution dare(SimpleMatrix A, SimpleMatrix B, SimpleMatrix Q, SimpleMatrix R)
            throws UnableToEvaluateStateSolution {
        isSystemValid(A, B, Q, R);
        SimpleMatrix S = SymplecticMatrix(A, B, Q, R);
        SimpleMatrix U = null;
        SimpleEVD<SimpleMatrix> sEvd = S.eig();
        List<Complex_F64> EiganValues = sEvd.getEigenvalues();
        EiganValues.sort(new Comparator<Complex_F64>() {

            @Override
            public int compare(Complex_F64 o1, Complex_F64 o2) {
                if (o1.real < o2.real)
                    return 1;
                else if (o1.real > o2.real)
                    return -1;
                return 0;
            }
        });
        for (int i = 0; i < EiganValues.size() / 2; i++) {
            double tmp = EiganValues.get(i).real;
            if (tmp < 0) {
                tmp = tmp * -1;
            }
            if (!(tmp < 1)) {
                throw new UnableToEvaluateStateSolution("The system has eiganvalues on the unit circle.");
            }
        }
        for (int i = 0; i < sEvd.getNumberOfEigenvalues(); i++) {
            int place = 0;
            for (Complex_F64 c : sEvd.getEigenvalues()) {
                if (c.real == EiganValues.get(i).real) {
                    if (U == null)
                        U = sEvd.getEigenVector(place);
                    else {
                        U = U.combine(0, i, sEvd.getEigenVector(place));
                    }
                    break;
                }
                place++;
            }
        }
        int row = U.numRows();
        int col = U.numCols();
        SimpleMatrix U11 = U.extractMatrix(0, row / 2, 0, col / 2);
        SimpleMatrix U21 = U.extractMatrix(row / 2, row, 0, col / 2);
        double[][] tmp = new double[1][sEvd.getNumberOfEigenvalues() / 2];
        for (int i = 0; i < sEvd.getNumberOfEigenvalues() / 2; i++) {
            tmp[0][i] = EiganValues.get(i).real;
        }
        SimpleMatrix P = U21.mult(U11.invert());
        SimpleMatrix E = new SimpleMatrix(tmp);
        SimpleMatrix K = R.plus(B.transpose().mult(P.mult(B))).invert().mult(B.transpose().mult(P.mult(A)));
        return new StateSolution(K, P, E);
    }

    public static StateSolution care(SimpleMatrix A, SimpleMatrix B, SimpleMatrix Q, SimpleMatrix R)
            throws UnableToEvaluateStateSolution {
        isSystemValid(A, B, Q, R);
        SimpleMatrix H = HamiltonianMatrix(A, B, Q, R);
        SimpleMatrix U = null;
        SimpleEVD<SimpleMatrix> hEvd = H.eig();
        List<Complex_F64> EiganValues = hEvd.getEigenvalues();
        for (Complex_F64 c : EiganValues) {
            if (c.imaginary != 0) {
                throw new UnableToEvaluateStateSolution("The system has eigenvalues on the imaginary axis.");
            }
        }
        EiganValues.sort(new Comparator<Complex_F64>() {

            @Override
            public int compare(Complex_F64 o1, Complex_F64 o2) {
                if (o1.real > o2.real)
                    return 1;
                else if (o1.real < o2.real)
                    return -1;
                return 0;
            }
        });
        for (int i = 0; i < hEvd.getNumberOfEigenvalues(); i++) {
            int place = 0;
            for (Complex_F64 c : hEvd.getEigenvalues()) {
                if (c.real == EiganValues.get(i).real) {
                    if (U == null)
                        U = hEvd.getEigenVector(place);
                    else {
                        U = U.combine(0, i, hEvd.getEigenVector(place));
                    }
                    break;
                }
                place++;
            }
        }
        int row = U.numRows();
        int col = U.numCols();
        SimpleMatrix U11 = U.extractMatrix(0, row / 2, 0, col / 2);
        SimpleMatrix U21 = U.extractMatrix(row / 2, row, 0, col / 2);
        double[][] tmp = new double[1][hEvd.getNumberOfEigenvalues() / 2];
        for (int i = 0; i < hEvd.getNumberOfEigenvalues() / 2; i++) {
            tmp[0][i] = EiganValues.get(i).real;
        }
        SimpleMatrix P = U21.mult(U11.invert());
        SimpleMatrix E = new SimpleMatrix(tmp);
        SimpleMatrix K = R.invert().mult(B.transpose().mult(P));
        return new StateSolution(K, P, E);
    }

    public static SimpleMatrix eye(int size) {
        return SimpleMatrix.identity(size);
    }

    private static void isSystemValid(SimpleMatrix A, SimpleMatrix B, SimpleMatrix Q, SimpleMatrix R) {
        if (A.numCols() != Q.numCols() || A.numRows() != Q.numRows()) {
            throw new InvalidSystemException("Q and A must be of the same dimentions.");
        } else if (R.numCols() != R.numRows()) {
            throw new InvalidSystemException("R must be a square matrix.");
        } else if (B.numCols() != R.numRows()) {
            throw new InvalidSystemException("R must have the same number of rows as B's cols.");
        }
    }

    private static SimpleMatrix HamiltonianMatrix(SimpleMatrix A, SimpleMatrix B, SimpleMatrix Q, SimpleMatrix R) {
        SimpleMatrix H12 = B.mult(R.invert().mult(B.transpose())).negative();
        SimpleMatrix H1 = A.combine(0, A.numCols(), H12);
        SimpleMatrix H2 = Q.negative().combine(0, Q.numCols(), A.transpose().negative());
        SimpleMatrix H = H1.combine(H1.numRows(), 0, H2);
        return H;
    }

    private static SimpleMatrix SymplecticMatrix(SimpleMatrix A, SimpleMatrix B, SimpleMatrix Q, SimpleMatrix R) {
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