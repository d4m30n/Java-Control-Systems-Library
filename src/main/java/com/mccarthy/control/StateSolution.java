package com.mccarthy.control;

import org.ejml.simple.SimpleMatrix;
import org.ejml.data.Complex_F64;

public class StateSolution {
    private SimpleMatrix _K;
    private SimpleMatrix _S;
    private Complex_F64[] _E;

    public StateSolution(SimpleMatrix K, SimpleMatrix S, Complex_F64[] E) {
        _K = K;
        _S = S;
        _E = E;
    }

    public SimpleMatrix GetK() {
        return _K.copy();
    }

    public SimpleMatrix GetS() {
        return _S.copy();
    }

    public Complex_F64[] GetE() {
        return _E.clone();
    }
}