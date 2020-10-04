# The following functions use the package "ode2analyzer.sage"
# Thus it is loaded initially
load('ode2analyzer.sage')

####################################################################
def change_indices(newvariable,mystructure):
    # This function finds the multiplier factor that changes the indices.
    # We try to make one of the indices zero to obtain
    # the standard form of the singularity structure of the ODE.
    
    multiplier=1
    if oo in mystructure[0] or oo in mystructure[2]:
        for ii in range(len(mystructure[0])):
            if len(mystructure[1][ii])==1:
                mystructure[1][ii].append(mystructure[1][ii][0])
            if mystructure[0][ii] == oo:
                pass
            elif mystructure[1][ii][0] != 0 and mystructure[1][ii][1] != 0:
                multiplier=multiplier*((newvariable-mystructure[0][ii])^(mystructure[1][ii][1]))   
    else:
        for ii in range(len(mystructure[0])-1):
            if len(mystructure[1][ii])==1:
                mystructure[1][ii].append(mystructure[1][ii][0])
            if mystructure[1][ii][0] != 0 and mystructure[1][ii][1] != 0:
                termnumer=newvariable-mystructure[0][ii]
                termdenom=newvariable-mystructure[0][-1]
                multiplier=multiplier*((termnumer/termdenom)^(mystructure[1][ii][1]))
    return multiplier

####################################################################
def rootselector(roots):
    # This function selects the positive root if two roots are the same
    # or the smallest one in size.
    var('myroot')
    try:
        if len(str(roots[0].rhs()))<len(str(roots[1].rhs())):
            myroot=roots[0].rhs()
        elif len(str(roots[1].rhs()))<len(str(roots[0].rhs())):
            myroot=roots[1].rhs()
        elif roots[0].rhs()==-roots[1].rhs():
            myroot=roots[1].rhs()
        else:
            myroot=roots[0].rhs()
    except:
        myroot=roots[0].rhs()
    return myroot

####################################################################
def find_special_ode(diffeqn,y,z):
    # This function matches the singularity structure with 
    # some previously defined ODEs.
    
    # diffeqn: The differential equation to be analyzed
    # y: The dependent function
    # z: The argument of y
    # The function shows and returns the "oderesult"
    
    # Find the singularity structure
    sing_struct=find_singularities(diffeqn,y,z)
    # Numbers of regular and irregular singularities
    numregular=len(sing_struct[0])
    numirregular=len(sing_struct[2])
    
    # Matching using the singularity structure 
    # (If a special form is found, the solution 
    # will be attempted.)
    if numregular==3 and numirregular==0:
        oderesult=find_2F1(diffeqn,y,z)
        show(oderesult)
    elif numregular==1 and numirregular==1 and sing_struct[3][0]==1:
        oderesult=find_1F1(diffeqn,y,z)
        show(oderesult)
    elif numregular==4 and numirregular==0:
        oderesult=find_HG(diffeqn,y,z)
        show(oderesult)
    elif numregular==2 and numirregular==1 and sing_struct[3][0]==1:
        oderesult=find_HC(diffeqn,y,z)
        show(oderesult)
    elif numregular==1 and numirregular==1 and sing_struct[3][0]==2:
        oderesult=find_HB(diffeqn,y,z)
        show(oderesult)
    elif numregular==0 and numirregular==1 and sing_struct[3][0]==3:
        oderesult=find_HT(diffeqn,y,z)
        show(oderesult)
    elif numregular==0 and numirregular==2 and sing_struct[3][0]==1 and sing_struct[3][1]==1:
        oderesult=find_HD(diffeqn,y,z)
        show(oderesult)
    else:
        show("No special form is found")
        oderesult=NaN
    return oderesult

####################################################################
def ode_finder_bruteforce(diffeqn,y,z):
    # This function tries to change the singularity structure
    # of a 2nd order ODE with polynomial coefficients using
    # some transformations of the form z->z^n
    
    # diffeqn: The differential equation to be analyzed
    # y: The dependent function
    # z: The argument of y
    # This function returns nothing. 
    
    var('counter') # will count the transformations
    counter=1
    exponents=[1,-1/2,1/2,2,-1/4,1/4,4] #n in z^n
    show("The number of transformations that will be tried:  ",len(exponents))
    for ii in exponents:
        show("Try (",counter,"):")
        show("Trying for z->z**(",ii,")")
        diffeqn2=ode_change_of_variable(diffeqn=diffeqn,y=y,z=z,transformation=z^ii)
        try:
            find_special_ode(diffeqn=diffeqn2,y=y,z=z)
        except:
            show("Unable to find a special form for this transformation.")
            show("*** The transformation may have caused some errors. ***")
        counter+=1
    show("The process is finished.")

####################################################################
# The following functions solve some 2nd order ODEs with 
# polynomial coefficients in terms of the associated special functions.
#
# Please note the following:
# i. The forms of the equations are given in the functions in the form
#    f1(z)*diff(u(z),z,z)+f2(z)*diff(u(z),z)+f3(z)*u(z)
# ii. The forms of the equtions and the solutions can be different
#     from other computer algebra systems (Maple, Mathematica, etc.).
# iii. *** The functions are not heavily tested. They should be used with caution.***
# iv. Please report any issues to Tolga Birkandan (birkandant@itu.edu.tr).
#
# Input arguments:
# diffeqn: The differential equation to be analyzed
# y: The dependent function
# z: The argument of y
# Output: 
# The solution of diffeqn in terms of the associated special function
#
# The method:
# 1. Find the singularity structure.
# 2. If the indices are not in the standard form,
#    i.e. one of the indices in non-zero,
#    redefine y(z) as
#    y(z)->y(z)*[(z-z[i])^(non-zero index)].
# 3. If the locations of the singular points are not
#    in the standard form, apply a transformation
#    to relocate them.
# 4. Use the indices or the normal form of the ODE 
#    to read the parameters.  
# 5. Form the solution using the parameters and
#    the redefinition/transformation defined above.
# 
# *** The functions are NOT heavily tested. They should be used with caution!***
####################################################################
def find_HT(diffeqn,y,z):
    # The standard form of the equation is
    # f1(z)*diff(u(z),z,z)+f2(z)*diff(u(z),z)+f3(z)*u(z)
    # where
    # var('alpha,gamma,delta,epsilon,q');
    # f1(z)=1;f2(z)=gamma+delta*z+epsilon*z^2;f3(z)=(alpha*z-q);
    
    mystructure=find_singularities(diffeqn,y,z)
    
    f1 = function('f1')(z) 
    f2 = function('f2')(z) 
    f3 = function('f3')(z)
    newvariable=z

    f1(z)=diffeqn.coefficient(diff(y(z),z,z))
    f2(z)=((diffeqn.coefficient(diff(y(z),z)))/f1(z)).full_simplify()
    f3(z)=((diffeqn.coefficient(y(z)))/f1(z)).full_simplify()
    diffeqn=diff(y(z),z,z)+f2(z)*diff(y(z),z)+f3(z)*y(z)
    
    var('x11,x22,bb,C1,C2')
    if not(oo in mystructure[2]):
        mytransformation=(x11*z)/(z+1)
        mytransformation=mytransformation.limit(x11=mystructure[2][0])
        newvariable=solve(mytransformation==bb,z)[0].rhs().subs(bb==z)
        diffeqn=ode_change_of_variable(diffeqn,y,z,mytransformation)
    
    mystructure=find_singularities(diffeqn,y,z)   
    
    f1(z)=diffeqn.coefficient(diff(y(z),z,z))
    f2(z)=((diffeqn.coefficient(diff(y(z),z)))/f1(z)).full_simplify()
    f3(z)=((diffeqn.coefficient(y(z)))/f1(z)).full_simplify()
    
    mult1=exp((-1/2)*integrate(f2(z),z))
    
    diffeqn=normal_form_ode2(diffeqn,y,z)

    coef4=diffeqn.coefficient(y(z)).coefficient(z,4)
    coef3=diffeqn.coefficient(y(z)).coefficient(z,3)
    coef2=diffeqn.coefficient(y(z)).coefficient(z,2)
    coef1=diffeqn.coefficient(y(z)).coefficient(z,1)
    coef0=diffeqn.coefficient(y(z)).coefficient(z,0)

    var("myq,myalpha,mygamma,mydelta,myepsilon,myz")
    var("myq2,myalpha2,mygamma2,mydelta2,myepsilon2")
    
    myepsilon=rootselector(solve(-4*coef4==myepsilon2^2,myepsilon2))
    mydelta=-2*coef3/myepsilon
    mygamma=(-4*coef2-mydelta^2)/(2*myepsilon)
    myalpha=coef1+(1/2)*mydelta*mygamma+myepsilon
    myq=-coef0-(1/4)*(mygamma^2)-(1/2)*mydelta
    
    originalmult=exp((-1/2)*integrate(mygamma+mydelta*z+myepsilon*z^2,z))
    mult=(mult1/originalmult).expand().canonicalize_radical().full_simplify()
    mult=mult.subs(z==newvariable)
    
    HT = function('HT')(myq2,myalpha2,mygamma2,mydelta2,myepsilon2,myz)
    myHTresult1=HT(myq2=myq,myalpha2=myalpha,mygamma2=mygamma,\
                   mydelta2=mydelta,myepsilon2=myepsilon,myz=newvariable)
    myHTresult2=HT(myq2=myq+mydelta,myalpha2=myalpha-2*myepsilon,mygamma2=-mygamma,\
                   mydelta2=-mydelta,myepsilon2=-myepsilon,myz=newvariable)
    myHTresult=C1*myHTresult1+\
    C2*(exp((-mygamma*newvariable)-(mydelta/2)*newvariable^2-(myepsilon/3)*newvariable^3))*myHTresult2
    myHTresult=(myHTresult*mult).expand()
    return (myHTresult)

####################################################################
def find_HB(diffeqn,y,z):
    # The standard form of the equation is
    # f1(z)*diff(u(z),z,z)+f2(z)*diff(u(z),z)+f3(z)*u(z)
    # where
    # var('alpha,gamma,delta,epsilon,q');
    # f1(z)=1;f2(z)=((gamma/z)+delta+epsilon*z);f3(z)=(alpha*z-q)/(z);
    
    mystructure=find_singularities(diffeqn,y,z)
    
    f1 = function('f1')(z) 
    f2 = function('f2')(z) 
    f3 = function('f3')(z)
    newvariable=z
    
    f1(z)=diffeqn.coefficient(diff(y(z),z,z))
    f2(z)=diffeqn.coefficient(diff(y(z),z))
    f3(z)=diffeqn.coefficient(y(z))
    
    multiplier=change_indices(newvariable,mystructure)
    
    diffeqn=(f1(z)*diff(multiplier*y(z),z,z)+f2(z)*diff(multiplier*y(z),z)+f3(z)*multiplier*y(z))\
    .normalize().numerator().full_simplify()

    f1(z)=diffeqn.coefficient(diff(y(z),z,z))
    f2(z)=((diffeqn.coefficient(diff(y(z),z)))/f1(z)).full_simplify()
    f3(z)=((diffeqn.coefficient(y(z)))/f1(z)).full_simplify()
    diffeqn=diff(y(z),z,z)+f2(z)*diff(y(z),z)+f3(z)*y(z)
    
    var('x11,x22,bb,C1,C2')
    if not(0 in mystructure[0] and oo in mystructure[2]):
        mytransformation=(x11*z+x22)/(z+1)
        mytransformation=mytransformation.limit(x11=mystructure[2][0])
        mytransformation=mytransformation.limit(x22=mystructure[0][0])
        newvariable=solve(mytransformation==bb,z)[0].rhs().subs(bb==z)
        diffeqn=ode_change_of_variable(diffeqn,y,z,mytransformation)
    
    mystructure=find_singularities(diffeqn,y,z) 
    
    mult1=exp((-1/2)*integrate(f2(z),z))

    diffeqn=normal_form_ode2(diffeqn,y,z)
    
    var("myq,myalpha,mygamma,mydelta,myepsilon,myz")
    var("myqsol,myalphasol,mygammasol,mydeltasol,myepsilonsol,myz")
    var("myq2,myalpha2,mygamma2,mydelta2,myepsilon2")
    var('coef2,coef1,coef0,coefm1,coefm2')
    
    coef2=diffeqn.expand().coefficient(y(z)).coefficient(z,2)
    coef1=diffeqn.expand().coefficient(y(z)).coefficient(z,1)
    coef0=diffeqn.expand().coefficient(y(z)).coefficient(z,0)
    coefm1=diffeqn.expand().coefficient(y(z)).coefficient(z,-1)
    coefm2=diffeqn.expand().coefficient(y(z)).coefficient(z,-2)

    myepsilonsol=solve(-4*coef2==myepsilon2^2,myepsilon2)
    try:
        if myepsilonsol[0].rhs()==-2:
            myepsilon=myepsilonsol[0].rhs()
        elif myepsilonsol[1].rhs()==-2:
            myepsilon=myepsilonsol[1].rhs()
        else:
            myepsilon=rootselector(myepsilonsol)
    except:
        myepsilon=rootselector(myepsilonsol)
        
    mydelta=-2*coef1/myepsilon
    mygamma=rootselector(solve(-4*coefm2==mygamma2^2-2*mygamma2,mygamma2))
    myq=-(1/2)*mydelta*mygamma-coefm1
    myalpha=coef0+(1/4)*mydelta^2+(1/2)*myepsilon*mygamma+myepsilon/2
    
    originalmult=exp((-1/2)*integrate(((mygamma/z)+mydelta+myepsilon*z),z))
    mult=(mult1/originalmult).expand().canonicalize_radical().full_simplify()
    mult=mult.subs(z==newvariable)
    
    HB = function('HB')(myq2,myalpha2,mygamma2,mydelta2,myepsilon2,myz)
    myHBresult1=HB(myq2=myq,myalpha2=myalpha,mygamma2=mygamma,\
                   mydelta2=mydelta,myepsilon2=myepsilon,myz=newvariable)
    myHBresult2=HB(myq2=myq-(1-mygamma)*mydelta,myalpha2=myalpha+(1-mygamma)*myepsilon,mygamma2=2-mygamma,\
                   mydelta2=mydelta,myepsilon2=myepsilon,myz=newvariable)
    myHBresult=C1*myHBresult1+C2*(newvariable^(1-mygamma))*myHBresult2
    myHBresult=(myHBresult*mult).expand()
    return (myHBresult)

####################################################################
def find_HD(diffeqn,y,z):
    # The standard form of the equation is
    # f1(z)*diff(u(z),z,z)+f2(z)*diff(u(z),z)+f3(z)*u(z)
    # where
    # var('alpha,gamma,delta,epsilon,q');
    # f1(z)=1;f2(z)=(delta/z)+(gamma/(z^2))+epsilon;f3(z)=(alpha*z-q)/(z^2);
    
    mystructure=find_singularities(diffeqn,y,z)
    
    f1 = function('f1')(z) 
    f2 = function('f2')(z) 
    f3 = function('f3')(z)
    newvariable=z

    f1(z)=diffeqn.coefficient(diff(y(z),z,z))
    f2(z)=((diffeqn.coefficient(diff(y(z),z)))/f1(z)).full_simplify()
    f3(z)=((diffeqn.coefficient(y(z)))/f1(z)).full_simplify()
    diffeqn=diff(y(z),z,z)+f2(z)*diff(y(z),z)+f3(z)*y(z)
    
    var('x11,x22,bb,C1,C2')
    if not(0 in mystructure[2] and oo in mystructure[2]):
        mytransformation=(x11*z+x22)/(z+1)
        mytransformation=mytransformation.limit(x11=mystructure[2][1])
        mytransformation=mytransformation.limit(x22=mystructure[2][0])
        newvariable=solve(mytransformation==bb,z)[0].rhs().subs(bb==z)
        diffeqn=ode_change_of_variable(diffeqn,y,z,mytransformation)
    
    mystructure=find_singularities(diffeqn,y,z)   
    
    f1(z)=diffeqn.coefficient(diff(y(z),z,z))
    f2(z)=((z^2)*(diffeqn.coefficient(diff(y(z),z)))/f1(z)).full_simplify()
    f3(z)=((z^2)*(diffeqn.coefficient(y(z)))/f1(z)).full_simplify()
    
    ##################################
    mult1=exp((-1/2)*integrate(f2(z),z))

    diffeqn=(normal_form_ode2(diffeqn,y,z)*(z^4)).full_simplify()
    
    var("myq,myalpha,mygamma,mydelta,myepsilon,myz")
    var("myqsol,myalphasol,mygammasol,mydeltasol,myepsilonsol,myz")
    var("myq2,myalpha2,mygamma2,mydelta2,myepsilon2")
    var('coef2,coef1,coef0,coefm1,coefm2')
    
    coef0=diffeqn.expand().coefficient(y(z)).coefficient(z,0)
    coef1=diffeqn.expand().coefficient(y(z)).coefficient(z,1)
    coef2=diffeqn.expand().coefficient(y(z)).coefficient(z,2)
    coef3=diffeqn.expand().coefficient(y(z)).coefficient(z,3)
    coef4=diffeqn.expand().coefficient(y(z)).coefficient(z,4)
    
    myepsilon=((-2*I)*sqrt(coef4)).canonicalize_radical().full_simplify()
    mydelta=(2-(I*coef1)/sqrt(coef0)).canonicalize_radical().full_simplify()
    mygamma=((-2*I)*sqrt(coef0)).canonicalize_radical().full_simplify()
    myalpha=(coef3 + (-2*I - coef1/sqrt(coef0))*sqrt(coef4)).canonicalize_radical().full_simplify()
    myq=(((2*I)*sqrt(coef0)*coef1 + coef1^2 - 4*coef0*coef2 + 8*coef0^(3/2)*sqrt(coef4))/(4*coef0)).canonicalize_radical().full_simplify()
        
    originalmult=exp((-1/2)*integrate(((mydelta/z)+(mygamma/(z^2))+myepsilon)*z^2,z))
    mult=(mult1/originalmult).expand().canonicalize_radical().full_simplify()
    mult=mult.subs(z==newvariable)
                    
    HD = function('HD')(myq2,myalpha2,mygamma2,mydelta2,myepsilon2,myz)
    myHDresult1=HD(myq2=myq,myalpha2=myalpha,mygamma2=mygamma,\
                   mydelta2=mydelta,myepsilon2=myepsilon,myz=newvariable)
    myHDresult2=HD(myq2=myq+mydelta-2,myalpha2=myalpha-2*myepsilon,mygamma2=-mygamma,\
                   mydelta2=4-mydelta,myepsilon2=-myepsilon,myz=newvariable)
    myHDresult=C1*myHDresult1+C2*(exp((mygamma/newvariable)-myepsilon*newvariable))*(newvariable^(2-mydelta))*myHDresult2
    myHDresult=(myHDresult*mult).expand()
    return (myHDresult)

####################################################################
def find_HC(diffeqn,y,z):
    # The standard form of the equation is
    # f1(z)*diff(u(z),z,z)+f2(z)*diff(u(z),z)+f3(z)*u(z)
    # where
    # var('alpha,beta,gamma,mu,nu');
    # f1(z)=1;f2(z)=((beta+1)/z)+((gamma+1)/(z-1))+alpha;f3(z)=(mu/z)+(nu/(z-1))
    
    mystructure=find_singularities(diffeqn,y,z)
    
    f1 = function('f1')(z) 
    f2 = function('f2')(z) 
    f3 = function('f3')(z)
    newvariable=z
    
    f1(z)=diffeqn.coefficient(diff(y(z),z,z))
    f2(z)=diffeqn.coefficient(diff(y(z),z))
    f3(z)=diffeqn.coefficient(y(z))
    
    multiplier=change_indices(newvariable,mystructure)
    
    diffeqn=(f1(z)*diff(multiplier*y(z),z,z)+f2(z)*diff(multiplier*y(z),z)+f3(z)*multiplier*y(z))\
    .normalize().numerator().full_simplify()
    
    f1(z)=diffeqn.coefficient(diff(y(z),z,z))
    f2(z)=((diffeqn.coefficient(diff(y(z),z)))/f1(z)).full_simplify()
    f3(z)=((diffeqn.coefficient(y(z)))/f1(z)).full_simplify()
    diffeqn=diff(y(z),z,z)+f2(z)*diff(y(z),z)+f3(z)*y(z)
    
    var('x11,x22,x33,bb,C1,C2')
    if not(0 in mystructure[0] and 1 in mystructure[0] and oo in mystructure[2]):
        mytransformation=((x22+(z-1)*x33)*x11-x22*x33*z)/(x11*z+(-z+1)*x22-x33)
        mytransformation=mytransformation.limit(x33=mystructure[2][0])
        mytransformation=mytransformation.limit(x11=mystructure[0][1])
        mytransformation=mytransformation.limit(x22=mystructure[0][0])
        newvariable=solve(mytransformation==bb,z)[0].rhs().subs(bb==z)
        diffeqn=ode_change_of_variable(diffeqn,y,z,mytransformation)
    
    mystructure=find_singularities(diffeqn,y,z)   
    
    f1(z)=diffeqn.coefficient(diff(y(z),z,z))
    f2(z)=((diffeqn.coefficient(diff(y(z),z)))/f1(z)).full_simplify()
    f3(z)=((diffeqn.coefficient(y(z)))/f1(z)).full_simplify()
    
    mult1=exp((-1/2)*integrate(f2(z),z))

    diffeqn=(normal_form_ode2(diffeqn,y,z)*(z^2)*(z-1)^2).full_simplify()

    var("myalpha,mybeta,mygamma,mydelta,myeta,mymu,mynu,myz")
    var("myalpha2,mybeta2,mygamma2,mydelta2,myeta2")
    
    coef0=diffeqn.expand().coefficient(y(z)).coefficient(z,0)
    coef1=diffeqn.expand().coefficient(y(z)).coefficient(z,1)
    coef2=diffeqn.expand().coefficient(y(z)).coefficient(z,2)
    coef3=diffeqn.expand().coefficient(y(z)).coefficient(z,3)
    coef4=diffeqn.expand().coefficient(y(z)).coefficient(z,4)
        
    myalpha=rootselector(solve(-4*coef4==myalpha2^2,myalpha2))
    mybeta=rootselector(solve(1-4*coef0==mybeta2^2,mybeta2))
    mygamma=rootselector(solve(myalpha^2+mybeta^2-4*(coef1+coef2+coef3)==mygamma2^2,mygamma2))
    mymu=coef1+((mybeta*(myalpha-1))/2)-((mybeta^2)/2)-((mygamma*(mybeta+1))/2)+((myalpha)/2)
    mynu=coef3-((myalpha^2)/2)+((myalpha*mybeta)/2)+((myalpha*mygamma)/2)+myalpha-mymu
    
    originalmult=exp((-1/2)*integrate((((mybeta+1)/z)+((mygamma+1)/(z-1))+myalpha),z))
    mult=(mult1/originalmult).expand().canonicalize_radical().full_simplify()
    mult=mult.subs(z==newvariable)
    
    mydelta=(mymu+mynu-myalpha*((mybeta+mygamma+2)/2)).full_simplify()
    myeta=(((myalpha*(mybeta+1))/2)-mymu-((mybeta+mygamma+mybeta*mygamma)/2)).full_simplify()
                    
    HC = function('HC')(myalpha2,mybeta2,mygamma2,mydelta2,myeta2,myz)
    myHCresult1=HC(myalpha2=myalpha,mybeta2=mybeta,mygamma2=mygamma,\
                   mydelta2=mydelta,myeta2=myeta,myz=newvariable)
    myHCresult2=HC(myalpha2=myalpha,mybeta2=-mybeta,mygamma2=mygamma,\
                   mydelta2=mydelta,myeta2=myeta,myz=newvariable)
    myHCresult=C1*multiplier*myHCresult1+C2*multiplier*(newvariable^(-mybeta))*myHCresult2
    myHCresult=(myHCresult*mult).expand()
    return (myHCresult)

####################################################################
def find_HG(diffeqn,y,z):
    # The standard form of the equation is
    # f1(z)*diff(u(z),z,z)+f2(z)*diff(u(z),z)+f3(z)*u(z)
    # where
    # var('alpha,beta,gamma,delta,epsilon,q,a1');
    # epsilon=alpha+beta+1-gamma-delta;
    # f1(z)=1;f2(z)=(gamma/z)+(delta/(z-1))+(epsilon/(z-a1));
    # f3(z)=(alpha*beta*z-q)/(z*(z-1)*(z-a1))
    
    mystructure=find_singularities(diffeqn,y,z)
    
    f1 = function('f1')(z) 
    f2 = function('f2')(z) 
    f3 = function('f3')(z)
    newvariable=z
    
    f1(z)=diffeqn.coefficient(diff(y(z),z,z))
    f2(z)=diffeqn.coefficient(diff(y(z),z))
    f3(z)=diffeqn.coefficient(y(z))
    
    multiplier=change_indices(newvariable,mystructure)
    
    diffeqn=(f1(z)*diff(multiplier*y(z),z,z)+f2(z)*diff(multiplier*y(z),z)+f3(z)*multiplier*y(z))\
    .normalize().numerator().full_simplify()
    
    f1(z)=diffeqn.coefficient(diff(y(z),z,z))
    f2(z)=((diffeqn.coefficient(diff(y(z),z)))/f1(z)).full_simplify()
    f3(z)=((diffeqn.coefficient(y(z)))/f1(z)).full_simplify()
    diffeqn=diff(y(z),z,z)+f2(z)*diff(y(z),z)+f3(z)*y(z)
    
    del(f1,f2,f3)
    
    var('x11,x22,x33,bb,C1,C2')
    if not(0 in mystructure[0] and 1 in mystructure[0] and oo in mystructure[0]):
        mytransformation=((x22+(z-1)*x33)*x11-x22*x33*z)/(x11*z+(-z+1)*x22-x33)
        mytransformation=mytransformation.limit(x33=mystructure[0][3])
        mytransformation=mytransformation.limit(x11=mystructure[0][1])
        mytransformation=mytransformation.limit(x22=mystructure[0][0])
        #show(mytransformation)
        newvariable=solve(mytransformation==bb,z)[0].rhs().subs(bb==z)
        #show(newvariable)
        diffeqn=ode_change_of_variable(diffeqn,y,z,mytransformation)
    
    mystructure=find_singularities(diffeqn,y,z)   

    var("mythirdsing,myq,myalpha,mybeta,mygamma,mydelta,myepsilon,myz")
    var("mythirdsing2,myq2,myalpha2,mybeta2,mygamma2,mydelta2,myepsilon2")
    regularsings=mystructure[0]
    indices=mystructure[1]
    
    for ii in range(len(regularsings)):
        if len(indices[ii])==1:
            indices[ii].append(indices[ii][0])
    
    myq=diffeqn.coefficient(y(z))
    for ii in range(len(regularsings)):
        if regularsings[ii]==0:
            myq=myq*(z-regularsings[ii])
            if indices[ii][0]==0:
                mygamma=simplify_fullfull(-indices[ii][1]+1)
            else:
                mygamma=simplify_fullfull(-indices[ii][0]+1)
        if regularsings[ii]==1:
            myq=myq*(z-regularsings[ii])
            if indices[ii][0]==0:
                mydelta=simplify_fullfull(-indices[ii][1]+1)
            else:
                mydelta=simplify_fullfull(-indices[ii][0]+1)
        if not(regularsings[ii] in [0,1,oo]):
            myq=myq*(z-regularsings[ii])
            mythirdsing=regularsings[ii]
            if indices[ii][0]==0:
                myepsilon=simplify_fullfull(-indices[ii][1]+1)
            else:
                myepsilon=simplify_fullfull(-indices[ii][0]+1)
        if regularsings[ii]==infinity:
            mybeta=simplify_fullfull(indices[ii][0])
            myalpha=simplify_fullfull(indices[ii][1])
    
    myq=-simplify_fullfull(myq).coefficient(z,0)
    
    myfuschiancondition=simplify_fullfull(myepsilon-(myalpha+mybeta+1-mygamma-mydelta))
    if myfuschiancondition==0:
        HG = function('HG')(mythirdsing2,myq2,myalpha2,mybeta2,mygamma2,mydelta2,myz)
        myHGresult1=HG(mythirdsing2=mythirdsing,myq2=myq,myalpha2=myalpha,\
                       mybeta2=mybeta,mygamma2=mygamma,mydelta2=mydelta,myz=newvariable)
        myHGresult2=HG(mythirdsing2=mythirdsing,\
                       myq2=(mythirdsing*mydelta+myepsilon)*(1-mygamma)+myq,\
                       myalpha2=myalpha+1-mygamma,mybeta2=mybeta+1-mygamma,mygamma2=2-mygamma,\
                       mydelta2=mydelta,myz=newvariable)
        myHGresult=C1*multiplier*myHGresult1+C2*multiplier*(newvariable^(1-mygamma))*myHGresult2
        return (myHGresult)
    else:
        return (Nan)

####################################################################
def find_1F1(diffeqn,y,z):
    # The standard form of the equation is
    # f1(z)*diff(u(z),z,z)+f2(z)*diff(u(z),z)+f3(z)*u(z)
    # where
    # var('a,b');f1(z)=z;f2(z)=b-z;f3(z)=-a;
    
    mystructure=find_singularities(diffeqn,y,z)
    
    f1 = function('f1')(z) 
    f2 = function('f2')(z) 
    f3 = function('f3')(z)
    newvariable=z
    
    multiplier=change_indices(newvariable,mystructure)
    diffeqn=diffeqn.substitute_function(y(z),multiplier*u(x)).normalize().numerator().full_simplify()
   
    var('x11,x22,bb,C1,C2')
    if not(0 in mystructure[0] and oo in mystructure[2]):
        mytransformation=(x11*z+x22)/(z+1)
        mytransformation=mytransformation.limit(x11=mystructure[2][0])
        mytransformation=mytransformation.limit(x22=mystructure[0][0])
        newvariable=solve(mytransformation==bb,z)[0].rhs().subs(bb==z)
        diffeqn=ode_change_of_variable(diffeqn,y,z,mytransformation)
    
    f1(z)=diffeqn.coefficient(diff(y(z),z,z))
    f2(z)=((diffeqn.coefficient(diff(y(z),z)))/f1(z)).full_simplify()
    f3(z)=((diffeqn.coefficient(y(z)))/f1(z)).full_simplify()
    diffeqn=diff(y(z),z,z)+f2(z)*diff(y(z),z)+f3(z)*y(z)
    
    mult1=exp((-1/2)*integrate(f2(z),z))
    
    diffeqn=normal_form_ode2(diffeqn,y,z)
    
    var("mya,myb,mya2,myb2")
    var('coef0,coefm1,coefm2')
    
    coef0=diffeqn.expand().coefficient(y(z)).coefficient(z,0)
    newvariable=newvariable.subs(z==(2*sqrt(coef0)*z)/I)
    
    diffeqn=ode_change_of_variable(diffeqn,y,z,(I*z)/(2*sqrt(coef0)))    
    diffeqn=(diffeqn/diffeqn.coefficient(diff(y(z),z,z))).expand()
    
    coef0=diffeqn.expand().coefficient(y(z)).coefficient(z,0)
    coefm1=diffeqn.expand().coefficient(y(z)).coefficient(z,-1)
    coefm2=diffeqn.expand().coefficient(y(z)).coefficient(z,-2)
   
    myb=rootselector(solve(-4*coefm2==myb2^2-2*myb2,myb2))
    mya=(myb/2)-coefm1
    
    originalmult=exp((-1/2)*integrate(f2(z),z))
    mult=(mult1/originalmult).expand().canonicalize_radical().full_simplify()
    mult=mult.subs(z==newvariable)
        
    my1F1result=simplify_fullfull(C1*multiplier*hypergeometric_M(mya,myb,newvariable)\
                +C2*multiplier*hypergeometric_U(mya,myb,newvariable))
    my1F1result=(my1F1result*mult).expand()
    return (my1F1result)
####################################################################
def find_2F1(diffeqn,y,z):
    # The standard form of the equation is
    # f1(z)*diff(u(z),z,z)+f2(z)*diff(u(z),z)+f3(z)*u(z)
    # where
    # var('a,b,c');f1(z)=z*(1-z);f2(z)=(c-(a+b+1)*z);f3(z)=-a*b;
    
    mystructure=find_singularities(diffeqn,y,z)
    var('x11,x22,x33,bb,C1,C2')
    f1 = function('f1')(z) 
    f2 = function('f2')(z) 
    f3 = function('f3')(z)
    newvariable=z
    
    f1(z)=diffeqn.coefficient(diff(y(z),z,z))
    f2(z)=diffeqn.coefficient(diff(y(z),z))
    f3(z)=diffeqn.coefficient(y(z))
    
    multiplier=change_indices(newvariable,mystructure)

    if not(0 in mystructure[0] and 1 in mystructure[0] and oo in mystructure[0]):
        mytransformation=((x22+(z-1)*x33)*x11-x22*x33*z)/(x11*z+(-z+1)*x22-x33)
        mytransformation=mytransformation.limit(x33=mystructure[0][2])
        mytransformation=mytransformation.limit(x11=mystructure[0][1])
        mytransformation=mytransformation.limit(x22=mystructure[0][0])
        newvariable=solve(mytransformation==bb,z)[0].rhs().subs(bb==z)
        diffeqn=ode_change_of_variable(diffeqn,y,z,mytransformation)
    
    mystructure=find_singularities(diffeqn,y,z)
        
    diffeqn=(f1(z)*diff(multiplier*y(z),z,z)+f2(z)*diff(multiplier*y(z),z)+f3(z)*multiplier*y(z))\
    .normalize().numerator().full_simplify()
    
    f1(z)=diffeqn.coefficient(diff(y(z),z,z))
    f2(z)=((diffeqn.coefficient(diff(y(z),z)))/f1(z)).full_simplify()
    f3(z)=((diffeqn.coefficient(y(z)))/f1(z)).full_simplify()
    diffeqn=diff(y(z),z,z)+f2(z)*diff(y(z),z)+f3(z)*y(z)
    
    var("mya,myb,myc")
    regularsings=mystructure[0]
    indices=mystructure[1]
    
    for ii in range(len(regularsings)):
        if len(indices[ii])==1:
            indices[ii].append(indices[ii][0])
    
    for ii in range(len(regularsings)):
        if regularsings[ii]==0:
            if indices[ii][0]==0:
                myc=simplify_fullfull(-indices[ii][1]+1)
            else:
                myc=simplify_fullfull(-indices[ii][0]+1)
        if regularsings[ii]==infinity:
            myb=simplify_fullfull(indices[ii][0])
            mya=simplify_fullfull(indices[ii][1])

    my2F1result=simplify_fullfull(C1*multiplier*hypergeometric([mya,myb],[myc],newvariable)\
                +C2*(newvariable^(-myc+1))*multiplier*hypergeometric([mya-myc+1,myb-myc+1],[2-myc],newvariable))
    return (my2F1result)