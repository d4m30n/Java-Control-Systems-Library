package com.mccarthy.control;

public class DtNotSetException extends Exception {
	/** 
	 *
	 */
	private static final long serialVersionUID = 7807210433754119748L;

	public DtNotSetException(String meassage) {
		super("DtNotSetException - " + meassage);
	}
}