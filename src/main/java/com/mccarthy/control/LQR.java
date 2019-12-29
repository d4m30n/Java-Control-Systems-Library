package com.mccarthy.control;

import org.ejml.data.Complex_F64;
import org.ejml.simple.SimpleMatrix;

public class LQR implements StateSolution {

    private final Care _care;

    public LQR(SS sys, SimpleMatrix Q, SimpleMatrix R) throws UnableToEvaluateStateSolution {
        this._care = new Care(sys, Q, R);
    }

    @Override
    public SimpleMatrix getK() {
        return _care.getK();
    }

    @Override
    public SimpleMatrix getS() {
        return _care.getS();
    }

    @Override
    public Complex_F64[] getE() {
        return _care.getE();
    }
}