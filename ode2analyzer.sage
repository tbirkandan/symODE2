####################################################################
def simplify_fullfull(theinput):
    #This routine is used at some points in the functions 
    return theinput.expand().canonicalize_radical().full_simplify()

####################################################################
def normal_form_ode2(diffeqn,y,z):
    # This function brings the ODE into its normal form
    # by removing its first derivative term.
    # THe output form is y''+derless*y
    # diffeqn: The differential equation to be analyzed
    # y: The dependent function
    # z: The argument of y
    
    f1 = function('f1')(z) #f1*y''+f2*y'+f3*y=0
    f2 = function('f2')(z) #f1*y''+f2*y'+f3*y=0
    f3 = function('f3')(z) #f1*y''+f2*y'+f3*y=0
    p = function('p')(z) #y''+p*y'+q*y=0
    q = function('q')(z) #y''+p*y'+q*y=0
    derless=function('derless')(z) #f''+derless*f=0
    
    f1(z)=diffeqn.coefficient(diff(y(z),z,z))
    f2(z)=diffeqn.coefficient(diff(y(z),z))
    f3(z)=diffeqn.coefficient(y(z))
    p(z)=simplify_fullfull((f2(z)/f1(z)))
    q(z)=simplify_fullfull((f3(z)/f1(z)).full_simplify())
    derless=simplify_fullfull((q(z)-((p(z)^2)/4)-(diff(p(z),z)/2)))
    return diff(y(z),z,z)+derless*y(z)

####################################################################
def ode_change_of_variable(diffeqn,y,z,transformation):
    # This function transforms the 2nd order ode as z->z'(z)
    # diffeqn: The differential equation to be analyzed
    # y: The dependent function
    # z: The argument of y
    # transformation: transformation formula

    f1 = function('f1')(z) #f1*y''+f2*y'+f3*y=0
    f2 = function('f2')(z) #f1*y''+f2*y'+f3*y=0
    f3 = function('f3')(z) #f1*y''+f2*y'+f3*y=0
    f1(z)=diffeqn.coefficient(diff(y(z),z,z))
    f2(z)=diffeqn.coefficient(diff(y(z),z))
    f3(z)=diffeqn.coefficient(y(z))
    
    #Apply the transformation to the equation
    f1(z)=f1(z).subs(z==transformation)
    f2(z)=f2(z).subs(z==transformation)
    f3(z)=simplify_fullfull(f3(z).subs(z==transformation))
    f2(z)=simplify_fullfull((f2(z)*(diff(transformation, z))^2-f1(z)*(diff(transformation, z, z)))/(diff(transformation, z))^3)
    f1(z)=simplify_fullfull(f1(z)/(diff(transformation, z))^2)
                            
    return f1(z)*diff(y(z),z,z)+f2(z)*diff(y(z),z)+f3(z)*y(z)
            
####################################################################
def find_indices_recurrence(diffeqn,y,z,point,index,operation):
    # This function finds the indices and/or recurrence relation 
    # for a point for a 2nd order ODE with polynomial coefficients.
    # diffeqn: The differential equation to be analyzed
    # y: The dependent function
    # z: The independent parameter
    # point: The resular singular point (z=z_0)
    # index: The variable to assign the indices
    # operation: "indices", "recurrence" or "full"
    ## indices: returns the indices
    ## recurrence: returns the recurrence relation
    ## full: returns the indices and the recurrence relation
    
    # For the $\theta$-operator method, please see:
    # Gabriel Allen, "Some Efficient Methods for Obtaining...",
    # NASA Technical Report, NASA TR R-390 (1972).
    
    var("theta,n,thegcd")
    f1 = function('f1')(z) #f1*y''+f2*y'+f3*y=0
    f2 = function('f2')(z) #f1*y''+f2*y'+f3*y=0
    f3 = function('f3')(z) #f1*y''+f2*y'+f3*y=0
    C = function('C')(n) # Recurrence elements
    
    # Make transformation according to z=point
    if point==infinity:
        f1(z)=(z^4)*(diffeqn.coefficient(diff(y(z),z,z))).subs(z==1/z).expand()
        f2(z)=2*(z^3)*(diffeqn.coefficient(diff(y(z),z,z))).subs(z==1/z)-(z^2)*(diffeqn.coefficient(diff(y(z),z))).subs(z==1/z).expand()
        f3(z)=(diffeqn.coefficient(y(z))).subs(z==1/z).expand()
    else:
        f1(z)=(diffeqn.coefficient(diff(y(z),z,z))).subs(z==z+point).expand()
        f2(z)=(diffeqn.coefficient(diff(y(z),z))).subs(z==z+point).expand()
        f3(z)=(diffeqn.coefficient(y(z))).subs(z==z+point).expand()
        
    #Get rid of the common factors coming from some transformations outside this func.
    f2(z)=simplify_fullfull(f2(z)/f1(z))
    f3(z)=simplify_fullfull(f3(z)/f1(z))
    f1(z)=1
    
    while (f3(z).denominator()).has(z) or (f2(z).denominator()).has(z) or (f1(z).denominator()).has(z):
        if (f3(z).denominator()).has(z):
            f1(z)=(f1(z)*f3(z).denominator()).full_simplify().expand()
            f2(z)=(f2(z)*f3(z).denominator()).full_simplify()
            f3(z)=f3.numerator()       
        if (f2(z).denominator()).has(z):
            f1(z)=(f1(z)*f2(z).denominator()).full_simplify().expand()
            f3(z)=(f3(z)*f2(z).denominator()).full_simplify()
            f2(z)=f2.numerator()
        if (f1(z).denominator()).has(z):
            f2(z)=(f2(z)*f1(z).denominator()).full_simplify().expand()
            f3(z)=(f3(z)*f1(z).denominator()).full_simplify()
            f1(z)=f1.numerator()        
    
    # Make the lowest degree of f1 = 2 
    while f1(z).low_degree(z)>2:
        f1(z)=(f1(z)/z).expand()
        f2(z)=(f2(z)/z).expand()
        f3(z)=(f3(z)/z).expand()
    if f2(z)==0:
        while f1(z).low_degree(z)<2 or f3(z).low_degree(z)<0:
            f1(z)=(f1(z)*z)
            f3(z)=(f3(z)*z)
    else:
        while f1(z).low_degree(z)<2 or f2(z).low_degree(z)<1 or f3(z).low_degree(z)<0:
            f1(z)=(f1(z)*z)
            f2(z)=(f2(z)*z)
            f3(z)=(f3(z)*z)
    
    # Apply the $\theta$-operator method to get A[i]
    myeqn2=(((f1(z)/(z*z)).full_simplify())*theta*(theta-1)+((f2(z)/z).full_simplify())*theta+f3(z)).expand()
    
    A=[]
    for i in range(myeqn2.degree(z)+1):
        A.append(myeqn2.coefficient(z,i))

    # Find the indicial equation and solve it
    indices=solve(A[0].subs(theta==index)==0,index)
    for ii in range(len(indices)):
        indices[ii]=indices[ii].rhs().expand().full_simplify().canonicalize_radical()
    if operation=="indices":
        return (indices)
    # Find the recurrence relation
    else:
        recurrence_relation=0
        C(n)=C
        for ii in range(len(A)):
            recurrence_relation+=(C(n-ii)*(A[ii].subs(theta=n-ii+index)).full_simplify()).collect(C(n-ii))
        if operation=="recurrence":
            return (recurrence_relation)
        elif operation=="full":
            return (indices,[recurrence_relation])    
        
####################################################################
def find_singularities(diffeqn,y,z):
    # This function finds the singularity structure of 
    # a 2nd order ODE with polynomial coefficients.
    # diffeqn: The differential equation to be analyzed
    # y: The dependent function
    # z: The argument of y
    
    # The initial form of the eqn is f1*y''+f2*y'+f3*y=0
    # Then we will have the form y''+p*y'+q*y=0
    # We will define y(z)=g(z)*f(z) to get rid of the first derivative
    # The final form will be f''+derless*f=0
    
    # The result will be returned as an array:
    # ([regular singular point(s)],
    #  [indices of regular singular point(s) (respectively)]
    #  [irregular singular point(s)], 
    #  [Rank(s) of irregular singular point(s) (respectively)]) 
    
    var('t,myindices') #t=1/z, myindices will be the index argument of the regular singularities
    f1 = function('f1')(z) #f1*y''+f2*y'+f3*y=0
    f2 = function('f2')(z) #f1*y''+f2*y'+f3*y=0
    f3 = function('f3')(z) #f1*y''+f2*y'+f3*y=0
    p = function('p')(z) #y''+p*y'+q*y=0
    q = function('q')(z) #y''+p*y'+q*y=0
    derless=function('derless')(z) #f''+derless*f=0
    p2 = function('p2')(t) # For analyzing the infinity
    q2 = function('q2')(t) # For analyzing the infinity
    derless2 = function('derless2')(t) # For analyzing the infinity
    
    regularsingularities=[] # The array of the regular singular points
    indicesofregularsingularities=[] # The array of the indices of the regular sing.s
    irregularsingularities=[] # The array of the irregular singular points
    ranksofirregularsingularities=[] # The array of the ranks of the irrreg. sing.s
    
    f1(z)=diffeqn.coefficient(diff(y(z),z,z))
    f2(z)=diffeqn.coefficient(diff(y(z),z))
    f3(z)=diffeqn.coefficient(y(z))
    p(z)=simplify_fullfull((f2(z)/f1(z)))
    q(z)=simplify_fullfull((f3(z)/f1(z)).full_simplify())
    derless=simplify_fullfull((q(z)-((p(z)^2)/4)-(diff(p(z),z)/2)))

    myroots=solve(derless.denominator()==0,z) # Singular points
    
    for theroot in myroots:
        testterim=simplify_fullfull((((z-theroot.rhs())^2)*derless))
        test=simplify_fullfull((testterim.taylor(z,theroot.rhs(),0)))
        #if maxima.freeof(z,test):
        if test.has(z)==False:
            regularsingularities.append(theroot.rhs())
            indicesofregularsingularities.append(find_indices_recurrence(diffeqn=diffeqn,y=y,z=z,point=theroot.rhs(),index=myindices,operation="indices"))
        else:
            irregularsingularities.append(theroot.rhs())
            #Find the rank (see F.W.J. Olver, "Asymptotics and Special Functions, 1997"):
            kk=2
            while (True):
                testterim=simplify_fullfull((((z-theroot.rhs())^(2*kk))*derless))
                test=simplify_fullfull((testterim.taylor(z,theroot.rhs(),0)))
                #if maxima.freeof(z,test):
                if test.has(z)==False:
                    ranksofirregularsingularities.append(kk-1)
                    break
                else:
                    kk=kk+1           
            
    # Let us test z = infinity
    p2(t)=p(z).subs(z==1/t)
    q2(t)=q(z).subs(z==1/t)
    derless2=simplify_fullfull((1/(2*t^2))*(diff(p2(t),t))-(1/(4*t^4))*p2(t)^2+(1/(t^4))*q2(t))
    myroots=solve(derless2.denominator()==0,t)
    testterim=simplify_fullfull(((t^2)*derless2))
    test=simplify_fullfull((testterim.taylor(t,0,0)))
    #if maxima.freeof(t,test):
    if 0 in [ii.rhs() for ii in myroots]:
        if test.find(t)==[]:
            regularsingularities.append(oo)
            indicesofregularsingularities.append(find_indices_recurrence(diffeqn=diffeqn,y=y,z=z,point=oo,index=myindices,operation="indices"))
        else:
            irregularsingularities.append(oo)
            kk=2
            while (True):
                testterim=simplify_fullfull((((t)^(2*kk))*derless2))
                test=simplify_fullfull((testterim.taylor(t,0,0)))
                #if maxima.freeof(t,test):
                if test.find(t)==[]:
                    ranksofirregularsingularities.append(kk-1)
                    break
                else:
                    kk=kk+1

    return (regularsingularities,indicesofregularsingularities,irregularsingularities,ranksofirregularsingularities)