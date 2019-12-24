package com.mccarthy.control;

import org.ejml.simple.SimpleMatrix;

public class StateSolution {
    private SimpleMatrix _K;
    private SimpleMatrix _S;
    private SimpleMatrix _E;

    public StateSolution(SimpleMatrix K, SimpleMatrix S, SimpleMatrix E) {
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

    public SimpleMatrix GetE() {
        return _E.copy();
    }
}