# Java Linear Control Library

This library was designed as a project that required the ability to implement a Linear Quadratic Regulator in Java, however, to my knowledge there is no existing library that can evaluate an LQR controller in Java. The main focused for this library is to implement basic functions that are used in similar libraries such as the Python Control Systems Library by implementing the basic functions of implementing and validating a system along with the evaluation of the gain matrix *K*. Currently this library is a rough draft and may be missing many other functions that are required by a control library such a pole placement and a number of system stability functions.

## Library Usage

This library to keep things simple utilises the Efficient Java Matrix Library (EJML) as a way of managing the matrix equations and layouts, a lot of the results are also given using the EJML’s `SimpleMatrix` class. 

### Creating and Evaluating a System

The main class that is used in this library is the `SS` class which allows for the implementation and checking of an LQR controller made up of the *A*, *B*, *C*, and *D* matrices, this class is implemented using the example code below.

```
SS system = new SS(A, B, C, D);
```

The code above shows a simple implementation of the SS class where the matrices *A*, *B*, *C*, and *D* are passed through as a `SimpleMatrix` class from the EJML. Once a system is created there are a few class’s that can be used to evaluate the stability of the system through the eigenvalues and implement a simple step of the system given the vectors *x* and *u*.

```
system.pole();
system.stepSystem(x, u);
system.getOutputVector(x, u);
```

The three examples above allow for simple operation to be run on the given system where the first gets the eigenvalues of the given system, the second allows the system to be evaluated given the vectors *x* and *u*, and third allows for evaluating the *y* vector in the system also given *x* and *u*. Finally, the last important class is the `LQR` class which allows for evaluating the gain matrix *K* given the system above below shows the main operations of the `LQR` class and its usage.

```
LQR lqr = new LQR(system, Q, R);
lqr.getK();
lqr.getS();
lqr.getE();
```

Given the `SS` class of a system and the two matrices *Q* and *R* the `LQR` class allows for evaluating the gain matrix *K*, and can be retrieved using the `getK()` method witch returns the calculated *K* matrix. The other two methods allow for retrieving the intermediary calculated matrices and vectors where `getS()` returns the calculated solution to the riccati equations and `getE()` returns the eigenvectors of the closed system.

## Notes

This system was designed to serve the need of a specific project but done in a way that can allow for the use within other projects that need to evaluate LQR controllers. However, because it was built for a specific project you may find that somethings that you require are missing or implemented in a way that dose not work for your project so feel free to take this as a base and modify it in anyway that you feel is necessary.
