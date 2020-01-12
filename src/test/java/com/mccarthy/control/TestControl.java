package com.mccarthy.control;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import org.ejml.simple.SimpleMatrix;
import org.ejml.data.Complex_F64;
import org.junit.Test;

public class TestControl {

    SimpleMatrix A = null;
    SimpleMatrix B = null;
    SimpleMatrix C = null;
    SimpleMatrix D = null;
    SimpleMatrix R = null;
    SimpleMatrix Q = null;

    public void createSystem() {
        if (A != null)
            return;
        double[][] a = { { -0.5, 0 }, { 0, -0.1 } };
        double[][] b = { { 0.1, 0.2, 0.3 }, { 0.4, 0.5, 0.6 } };
        double[][] c = { { 1, 0 }, { 0, 1 } };
        double[][] d = { { 0, 0, 0 }, { 0, 0, 0 } };
        double[][] r = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
        double[][] q = { { 1, 0 }, { 0, 1 } };
        A = new SimpleMatrix(a);
        B = new SimpleMatrix(b);
        C = new SimpleMatrix(c);
        D = new SimpleMatrix(d);
        R = new SimpleMatrix(r);
        Q = new SimpleMatrix(q);
    }

    @Test
    public void testCare() {
        createSystem();

        testCareOne();
        testCareTwo();
    }

    private void testCareOne() {
        double[][] a = { { -0.5, 0 }, { 0, -0.1 } };
        double[][] b = { { 0.1, 0.2, 0.3 }, { 0.4, 0.5, 0.6 } };
        double[][] c = { { 1, 0 }, { 0, 1 } };
        double[][] d = { { 0, 0, 0 }, { 0, 0, 0 } };
        double[][] r = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
        double[][] q = { { 1, 0 }, { 0, 1 } };
        SimpleMatrix A = new SimpleMatrix(a);
        SimpleMatrix B = new SimpleMatrix(b);
        SimpleMatrix C = new SimpleMatrix(c);
        SimpleMatrix D = new SimpleMatrix(d);
        SimpleMatrix R = new SimpleMatrix(r);
        SimpleMatrix Q = new SimpleMatrix(q);
        SS sys = new SS(A, B, C, D);

        double[][] ssolution = { { 0.96930845, -0.22590614 }, { -0.22590614, 1.10011044 } };
        double[][] ksolution = { { 0.00656839, 0.41745356 }, { 0.08090862, 0.50487399 }, { 0.15524885, 0.59229442 } };
        double[] esolution = { -0.97965974, -0.45854855 };
        SimpleMatrix S = new SimpleMatrix(ssolution);
        SimpleMatrix K = new SimpleMatrix(ksolution);
        Complex_F64[] E = { new Complex_F64(esolution[0], 0d), new Complex_F64(esolution[1], 0d) };

        try {
            Care care = new Care(sys, Q, R);
            validateSolutions(care, K, S, E);
            LQR lqr = new LQR(sys, Q, R);
            validateSolutions(lqr, K, S, E);

        } catch (UnableToEvaluateStateSolution utess) {
            fail("Unable to evaluate the state solution");
        }
    }

    private void testCareTwo() {
        double[][] a = { { -0.5, 0 }, { 0, -0.1 } };
        double[][] b = { { 0.7, 0.8, 0.9 }, { 0.4, 0.5, 0.6 } };
        double[][] c = { { 1, 0 }, { 0, 1 } };
        double[][] d = { { 0, 0, 0 }, { 0, 0, 0 } };
        double[][] r = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 4 } };
        double[][] q = { { 5, 0 }, { 0, 3 } };
        SimpleMatrix A = new SimpleMatrix(a);
        SimpleMatrix B = new SimpleMatrix(b);
        SimpleMatrix C = new SimpleMatrix(c);
        SimpleMatrix D = new SimpleMatrix(d);
        SimpleMatrix R = new SimpleMatrix(r);
        SimpleMatrix Q = new SimpleMatrix(q);
        SS sys = new SS(A, B, C, D);

        double[][] ssolution = { { 3.01005001, -2.92989202 }, { -2.92989202, 6.59071229 } };
        double[][] ksolution = { { 0.9350782, 0.5853605 }, { 0.943094, 0.95144253 }, { 0.23777745, 0.32938114 } };
        double[] esolution = { -2.892366, -0.23815784 };
        SimpleMatrix S = new SimpleMatrix(ssolution);
        SimpleMatrix K = new SimpleMatrix(ksolution);
        Complex_F64[] E = { new Complex_F64(esolution[0], 0d), new Complex_F64(esolution[1], 0d) };

        try {
            Care care = new Care(sys, Q, R);
            validateSolutions(care, K, S, E);
            LQR lqr = new LQR(sys, Q, R);
            validateSolutions(lqr, K, S, E);

        } catch (UnableToEvaluateStateSolution utess) {
            fail("Unable to evaluate the state solution");
        }
    }

    private void validateSolutions(StateSolution ss, SimpleMatrix K, SimpleMatrix S, Complex_F64[] E) {
        assertEquals(K.numCols(), ss.getK().numCols());
        assertEquals(K.numRows(), ss.getK().numRows());
        assertEquals(S.numCols(), ss.getS().numCols());
        assertEquals(S.numRows(), ss.getS().numRows());
        assertEquals(E.length, ss.getE().length);
        for (int i = 0; i < K.numCols(); i++) {
            for (int j = 0; j < K.numRows(); j++) {
                assertEquals(K.get(j, i), ss.getK().get(j, i), 0.00000001);
            }
        }
        for (int i = 0; i < S.numCols(); i++) {
            for (int j = 0; j < S.numRows(); j++) {
                assertEquals(S.get(j, i), ss.getS().get(j, i), 0.00000001);
            }
        }
        for (int i = 0; i < E.length; i++) {
            assertEquals(E[i].real, ss.getE()[i].real, 0.0000001);
        }
    }

    @Test
    public void testDare() {
        createSystem();
        double[][] ksolution = { { -0.01498429, -0.02211717 }, { -0.06273946, -0.02600792 },
                { -0.11049463, -0.02989867 } };
        double[][] ssolution = { { 1.293394, -0.01092363 }, { -0.01092363, 1.00607547 } };
        double[] esolution = { -0.0559309, -0.45708437 };
        SimpleMatrix K = new SimpleMatrix(ksolution);
        SimpleMatrix S = new SimpleMatrix(ssolution);
        Complex_F64[] E = new Complex_F64[2];
        E[0] = new Complex_F64(esolution[0], 0d);
        E[1] = new Complex_F64(esolution[1], 0d);
        try {
            SS sys = new SS(A, B, C, D);
            Dare dare = new Dare(sys, Q, R);
            validateSolutions(dare, K, S, E);

        } catch (UnableToEvaluateStateSolution e) {
            fail("Unable to evaluate the system");
        }
    }
}