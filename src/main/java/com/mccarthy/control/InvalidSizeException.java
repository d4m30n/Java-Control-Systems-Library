package com.mccarthy.control;

public class InvalidSizeException extends RuntimeException{

    /**
     *
     */
    private static final long serialVersionUID = 1L;
    
    public InvalidSizeException(String message){
        super("InvalidSizeException - "+message);
    }
}