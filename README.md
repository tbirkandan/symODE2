# symODE2: Symbolic analysis of second-order ordinary differential equations with polynomial coefficients

by Tolga Birkandan
E-mail: birkandant@itu.edu.tr
Department of Physics, Istanbul Technical University, 34469 Istanbul, Turkey.

An open-source package for symbolic analysis of second-order ordinary differential equations with polynomial coefficients is proposed. The approach is mainly based on the singularity structure of the equation and the routines are written under the open-source computer algebra system SageMath. The code is able to obtain the singularity structure, indices and recurrence relations associated with the regular singular points, and symbolic solutions of the hypergeometric equation, Heun equation, and their confluent forms.

The symODE2 package is written under SageMath 9.1 using a laptop computer with Intel(R) Core(TM) i7-6500U CPU @ 2.50GHz and 8 GB memory. The operating system is Windows 10 Enterprise ver.1909.

The package consists of two main parts: ode2analyzer.sage for the general analysis and hypergeometric_heun.sage for the symbolic solutions of the equations. hypergeometric_heun.sage calls the routines defined in ode2analyzer.sage when needed. We suggest the user put these two files in the same directory. 

A sample SageMath session is also given.

## General analysis: ode2analyzer

The ﬁrst part, ode2analyzer contains the routines,

* find_singularities(diffeqn,y,z): Finds the singularity structure of the input ODE. The output is an array that involves the locations of the singularities, indices of the regular singularities and the ranks of the irregular singularities.

* find_indices_recurrence(diffeqn,y,z,point,index,operation): Finds the indices and/or the recurrence relation with respect to a regular singular point using the theta-operator method.

* ode_change_of_variable(diffeqn,y,z,transformation): Performs a change of variables.

* normal_form_ode2(diffeqn,y,z): Finds the normal form of a second order ODE. 

Here, the input arguments are

− diffeqn: A 2nd order ODE with polynomial coeﬃcients.

− y: The dependent variable of the ODE.

− z: The independent variable of the ODE.

− point: The regular singular point for which the indices and/or the recurrence relation will be found.

− index: The name of the parameter that will be used for the characteristic exponent.

− operation: The result of the routine (indices, recurrence, full (both)).

− transformation: The transformation function that deﬁnes the change of variables.

## Symbolic solutions of special ODEs: hypergeometric heun

The second part, hypergeometric_heun contains the routines,

* find_special_ode(diffeqn,y,z): Finds the type of the ODE using its singularity structure and solves it using the routines deﬁned below.

* ode_finder_bruteforce(diffeqn,y,z): Uses a change of variables list in order to bring the input ODE into a special form that is recognized in the package.

* find_2F1(diffeqn,y,z): Solves a hypergeometric equation.

* find_1F1(diffeqn,y,z): Solves a conﬂuent hypergeometric equation.

* find_HG(diffeqn,y,z): Solves a general Heun equation.

* find_HC(diffeqn,y,z): Solves a (singly) conﬂuent Heun equation.

* find_HD(diffeqn,y,z): Solves a double conﬂuent Heun equation.

* find_HB(diffeqn,y,z): Solves a biconﬂuent Heun equation.

* find_HT(diffeqn,y,z): Solves a triconﬂuent Heun equation.

- The input arguments (diffeqn,y,z) are the same as deﬁned in the ﬁrst part.

We change the indices of the regular singular points and move the locations of the singular points in order to obtain a standard form. After reaching the standard singularity structure of an equation that can be recognized by the code, the parameters are read either from the characteristic exponents or by matching the ﬁnal form of the input equation with the standard equation in their normal forms.

The polynomial coefficients of the normal forms are matched in the conﬂuent cases. The parameters of these ODEs can be found by solving single equations, i.e. the code ﬁnds some parameters by solving algebraic equations that depend only on one parameter. The rest of the parameters are found by substitution. For the Fuschian ODEs, the parameters are read from the characteristic exponents, e.g. the non-zero exponent of the singular point at zero in the hypergeometric equation yields 1 − c, etc.

We ﬁnd the parameters with this method and use the Maple or Mathematica forms of the solutions to substitute these parameters. The hypergeometric equation has three pairs of Frobenius solutions around its three regular singular points and these solutions can be transformed into other solutions via speciﬁc transformations. The number of all solutions of the hypergeometric equation is 24. The number of total solutions is 192 for the general Heun equation. The user of our code may need to use some transformations or function identities in order to obtain the desired form of the solution.

The results of the hypergeometric and conﬂuent hypergeometric equations are numerically usable as these functions are deﬁned in SageMath. However, the numerical solutions of the Heun-type functions are not deﬁned. The numerical solutions of the general Heun and (singly) conﬂuent Heun functions are deﬁned by Motygin for GNU Octave/MATLAB. GNU Octave/MATLAB commands can be run in a SageMath session. However, this procedure is not straightforward and it is beyond the scope of this work.

**Regarding our code, although the results of our limited number of tests seem proper, it should be emphasized that the solutions should be taken with great care as such routines should be tested heavily before using them for new problems. The user should approach these routines as an experimental eﬀort and follow the updates from the
GitHub page regularly.**
