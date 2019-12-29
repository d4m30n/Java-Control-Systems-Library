package com.mccarthy.control;

import org.ejml.simple.SimpleMatrix;
import org.ejml.data.Complex_F64;

public interface StateSolution {
    SimpleMatrix getK();

    SimpleMatrix getS();

    Complex_F64[] getE();
}