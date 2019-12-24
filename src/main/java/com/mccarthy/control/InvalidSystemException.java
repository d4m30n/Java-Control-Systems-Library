package com.mccarthy.control;

public class InvalidSystemException extends RuntimeException {

    /**
     *
     */
    private static final long serialVersionUID = 1L;

    public InvalidSystemException(String message) {
        super("InvalidSystemException - " + message);
    }

}