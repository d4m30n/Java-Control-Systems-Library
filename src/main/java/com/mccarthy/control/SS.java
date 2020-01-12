package com.mccarthy.control;

import org.ejml.data.Complex_F64;
import org.ejml.simple.SimpleMatrix;

public class SS {
    protected final SimpleMatrix _A;
    protected final SimpleMatrix _B;
    protected final SimpleMatrix _C;
    protected final SimpleMatrix _D;

    protected double _dt = -1; // Discrete Time state space

    /**
     * Create a continus time state space system
     * 
     * @param A - The system matrix
     * @param B - The input matrix
     * @param C - The output matrix
     * @param D - The feedthrough (or feedforward) matrix
     */
    public SS(SimpleMatrix A, SimpleMatrix B, SimpleMatrix C, SimpleMatrix D) {
        isSystemValid(A, B, C, D);
        _A = A;
        _B = B;
        _C = C;
        _D = D;
    }

    /**
     * Creates a descret time state space system
     * 
     * @param A  - The system matrix
     * @param B  - The input matrix
     * @param C  - The output matrix
     * @param D  - The feedthrough (or feedforward) matrix
     * @param dt - The descret time
     */
    public SS(SimpleMatrix A, SimpleMatrix B, SimpleMatrix C, SimpleMatrix D, double dt) {
        this(A, B, C, D);
        _dt = dt;
    }

    private static void isSystemValid(SimpleMatrix A, SimpleMatrix B, SimpleMatrix C, SimpleMatrix D) {
        if (A.numCols() != A.numRows()) {
            throw new InvalidSystemException("A must be a square matrix.");
        } else if (A.numCols() != B.numRows()) {
            throw new InvalidSystemException("B must have the same number of rows as A's cols.");
        } else if (C.numCols() != A.numCols() || C.numRows() != A.numRows()) {
            throw new InvalidSystemException("C must have the same number of rows and cols as A.");
        } else if (D.numCols() != B.numCols() || D.numRows() != B.numRows()) {
            throw new InvalidSystemException("D must have the same number of rows and cols as B.");
        }
    }

    private void isSystemValid(SimpleMatrix x, SimpleMatrix u) {
        if (x.numCols() > 1)
            throw new InvalidSystemException("x must only have 1 col.");
        if (u.numCols() > 1)
            throw new InvalidSystemException("u must only have 1 col.");
        if (x.numRows() != _A.numCols())
            throw new InvalidSystemException("x must have the same number of rows as A has cols.");
        if (u.numRows() != _B.numCols())
            throw new InvalidSystemException("u must have the same number of rows as B has cols.");
    }

    /**
     * Evaluate the current system given x(k) and u(k)
     * 
     * @param x - The state vector
     * @param u - The input (or control) vector
     * @return The state vector x(k+1)
     */
    public SimpleMatrix stepSystem(SimpleMatrix x, SimpleMatrix u) {
        isSystemValid(x, u);
        return _A.mult(x).plus(_B.mult(u));
    }

    /**
     * Evaluates the output vector of the system given x(k) and u(k)
     * 
     * @param x - The state vector
     * @param u - The input (or control) vector
     * @return The output vector y(k)
     */
    public SimpleMatrix getOutputVector(SimpleMatrix x, SimpleMatrix u) {
        isSystemValid(x, u);
        return _C.mult(x).plus(_D.mult(u));
    }

    /**
     * Given a system the eiganvalues of that system are returned
     * 
     * @param sys - The system
     * @return - The eiganvalues of the system
     */
    public Complex_F64[] pole() {
        return (Complex_F64[]) _A.eig().getEigenvalues().toArray();
    }

    public SimpleMatrix copyA() {
        return _A.copy();
    }

    public SimpleMatrix copyB() {
        return _B.copy();
    }

    public SimpleMatrix copyC() {
        return _C.copy();
    }

    public SimpleMatrix copyD() {
        return _D.copy();
    }

    /**
     * Get the descret time of the state space
     * 
     * @return The descret time
     * @throws DtNotSetException If the state space system is not a descret time
     *                           system
     */
    public double getDt() throws DtNotSetException {
        if (_dt == -1) {
            throw new DtNotSetException("The dt is not set for this system");
        }
        return _dt;
    }
}
