package com.mccarthy.control;

import java.util.Comparator;
import java.util.List;

import org.ejml.data.Complex_F64;
import org.ejml.simple.SimpleEVD;
import org.ejml.simple.SimpleMatrix;

public class Dare implements StateSolution {

    private final SimpleMatrix _K;
    private final SimpleMatrix _P;
    private final Complex_F64[] _E;

    public Dare(SS sys, SimpleMatrix Q, SimpleMatrix R) throws UnableToEvaluateStateSolution {
        SimpleMatrix A = sys._A;
        SimpleMatrix B = sys._B;
        SystemFunctions.isSystemValid(A, B, Q, R);
        SimpleMatrix S = SystemFunctions.SymplecticMatrix(A, B, Q, R);
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
        Complex_F64[] E = new Complex_F64[sEvd.getNumberOfEigenvalues() / 2];
        for (int i = 0; i < sEvd.getNumberOfEigenvalues() / 2; i++) {
            E[i] = EiganValues.get(i);
        }
        SimpleMatrix P = U21.mult(U11.invert());
        SimpleMatrix K = R.plus(B.transpose().mult(P.mult(B))).invert().mult(B.transpose().mult(P.mult(A)));
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