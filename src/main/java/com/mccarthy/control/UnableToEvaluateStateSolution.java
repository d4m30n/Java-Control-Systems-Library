package com.mccarthy.control;

public class UnableToEvaluateStateSolution extends Exception {

    /**
     *
     */
    private static final long serialVersionUID = 1L;

    public UnableToEvaluateStateSolution(String message) {
        super("UnableToEvalueateStateSolution - " + message);
    }
}