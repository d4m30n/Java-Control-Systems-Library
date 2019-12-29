package com.mccarthy.control;

import org.ejml.simple.SimpleMatrix;
import org.ejml.data.Complex_F64;

public class StateSolution {
    private SimpleMatrix _K;
    private SimpleMatrix _S;
    private Complex_F64[] _E;

    /**
     * Takes and store the state solution
     * 
     * @param K - The gain matrix K
     * @param S - The solution to the algebraic riccati equation
     * @param E - The eiganvalues of the system
     */
    public StateSolution(SimpleMatrix K, SimpleMatrix S, Complex_F64[] E) {
        _K = K;
        _S = S;
        _E = E;
    }

    /**
     * Creates a copy (preventing change to the original) and returns the gain
     * matrix K
     * 
     * @return - The gain matrix K
     */
    public SimpleMatrix getK() {
        return _K.copy();
    }

    /**
     * Creates a copy (preventing change to the original) and returns the solution
     * to the algebraic riccati equation
     * 
     * @return - The algebraic riccati equation
     */
    public SimpleMatrix getS() {
        return _S.copy();
    }

    /**
     * Creates a copy (preventing change to the original) and returns the eiganvalus
     * of the system
     * 
     * @return - The eiganvalues of the system
     */
    public Complex_F64[] getE() {
        return _E.clone();
    }
}