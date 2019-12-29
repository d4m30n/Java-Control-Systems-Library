package com.mccarthy.control;

import java.util.Comparator;
import java.util.List;

import org.ejml.data.Complex_F64;
import org.ejml.simple.SimpleEVD;
import org.ejml.simple.SimpleMatrix;

public class Care implements StateSolution {

    private final SimpleMatrix _K;
    private final SimpleMatrix _P;
    private final Complex_F64[] _E;

    public Care(SS sys, SimpleMatrix Q, SimpleMatrix R) throws UnableToEvaluateStateSolution {
        SimpleMatrix A = sys._A;
        SimpleMatrix B = sys._B;
        SystemFunctions.isSystemValid(A, B, Q, R);
        SimpleMatrix H = SystemFunctions.HamiltonianMatrix(A, B, Q, R);
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
        Complex_F64[] E = new Complex_F64[hEvd.getNumberOfEigenvalues() / 2];
        for (int i = 0; i < hEvd.getNumberOfEigenvalues() / 2; i++) {
            E[i] = EiganValues.get(i);
        }
        SimpleMatrix P = U21.mult(U11.invert());
        SimpleMatrix K = R.invert().mult(B.transpose().mult(P));
        this._K = K;
        this._P = P;
        this._E = E;
    }

    @Override
    public SimpleMatrix getK() {
        return _K.copy();
    }

    @Override
    public SimpleMatrix getS() {
        return _P.copy();
    }

    @Override
    public Complex_F64[] getE() {
        return _E.clone();
    }
}