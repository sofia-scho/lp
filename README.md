# CSC 445 Operations Research: Linear Programming
### Sofia Scholefield
### V00196009


## Running the Code
The code is written in the Python interpretted language. I used Python3. I also use two modules in my code: numpy and sys.

You can run the code by using the command python3 simplexMethod.py < inputFile.txt

## About the Program
I implemented the linear algebraic method. In the main function it checks to see if the linear program is primal feasible. If yes, it runs the simplex function (which is just the primal simplex method). Otherwise it checks if the linear program is dual feasible and runs the dualSimplex function titled "dual". If neither way is initially feasible main calls dual with c=0 and returns the value of B and N to then use in the primal simplex function. This method will not pivot because the rule I chose to implement for pivotting is Bland's Rule which does not cycle.

I took the logic/flow of the program from The Revised Simplex Slides (slide 70 and 94) courtesy of Bill Bird 2022.


## Extra Features
1. Primal-Dual Methods
2. Linear Algebraic Simplex Method
    Though I still used Matrix inverses.

## Future Improvements
I would like to implement something to prevent error accumulation when working with floatin point arithmetic. I would also like to implement a different primary pivotting rule that is "better" or usually faster than Bland's and instead only use Bland's when I detect the program is cycling ie. after three pivots and the objective has not increased use Bland's instead of Steepest Edge.