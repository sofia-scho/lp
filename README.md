# CSC 445 Operations Research: Linear Programming

## Running the Code
The code is written in the Python interpretted language. I used Python3. I also use three libraries in my code: numpy, sys, and math.

You can run the code by using the command python3 simplexMethod.py < inputFile.txt

## About the Program
I implemented the linear algebraic method. In the main function it checks to see if the linear program is primal feasible. If yes, it runs the simplex function (which is just the primal simplex method). Otherwise it checks if the linear program is dual feasible and runs the dualSimplex function titled "dual". If neither way is initially feasible main calls dual with c=0 and returns the value of B and N to then use in the primal simplex function. This method will not pivot because the rule I chose to implement for pivotting is Bland's Rule which does not cycle.

The sigFigs function I got from https://www.delftstack.com/howto/python/round-to-significant-digits-python/ I learned how to use the Python math module to round a number to specified significant figures.

## Extra Features
1. Primal-Dual Methods
2. Linear Algebraic Simplex Method
    Though I still used Matrix inverses.