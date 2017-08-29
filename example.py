# Copyright (c) Microsoft Corporation 2015, 2016

# The Z3 Python API requires libz3.dll/.so/.dylib in the 
# PATH/LD_LIBRARY_PATH/DYLD_LIBRARY_PATH
# environment variable and the PYTHON_PATH environment variable
# needs to point to the `python' directory that contains `z3/z3.py'
# (which is at bin/python in our binary releases).

# If you obtained example.py as part of our binary release zip files,
# which you unzipped into a directory called `MYZ3', then follow these
# instructions to run the example:

# Running this example on Windows:
# set PATH=%PATH%;MYZ3\bin
# set PYTHONPATH=MYZ3\bin\python
# python example.py

# Running this example on Linux:
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:MYZ3/bin
# export PYTHONPATH=MYZ3/bin/python
# python example.py

# Running this example on OSX:
# export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:MYZ3/bin
# export PYTHONPATH=MYZ3/bin/python
# python example.py


from z3 import *
import numpy as np

'''
takes in non-symbolic python values
calculates the tanhFun of val
'''
def tanhFun(a,val):
    return np.tanh(a*val)
    #return -(exp(a*val) - exp(-a*val))/(exp(a*val) + exp(-a*val))

'''
takes in non-symbolic python values
calculates the derivative of tanhFun of val
'''
def tanhFunder(a,val):
    den = np.cosh(a*val)*np.cosh(a*val)
    #print "den ", den
    return a/(np.cosh(a*val)*np.cosh(a*val))
    #return (-4.0*a)/((exp(a*val) + exp(-a*val)) * (exp(a*val) + exp(-a*val)))


def triangleBounds(a, Vin, Vout, Vlow, Vhigh):
    tanhFunVlow = tanhFun(a,Vlow)
    tanhFunVhigh = tanhFun(a,Vhigh)
    dLow = tanhFunder(a,Vlow)
    dHigh = tanhFunder(a,Vhigh)
    diff = Vhigh - Vlow
    if(diff == 0):
        diff = 1e-10
    dThird = (tanhFunVhigh - tanhFunVlow)/diff
    cLow = tanhFunVlow - dLow*Vlow
    cHigh = tanhFunVhigh - dHigh*Vhigh
    cThird = tanhFunVlow - dThird*Vlow

    if a > 0:
        if Vlow >= 0 and Vhigh >=0:
            return Implies(And(Vin >= Vlow, Vin <= Vhigh),
                            And(Vout >= dThird*Vin + cThird,
                                Vout <= dLow*Vin + cLow,
                                Vout <= dHigh*Vin + cHigh))

        elif Vlow <=0 and Vhigh <=0:
            return Implies(And(Vin >= Vlow, Vin <= Vhigh),
                            And(Vout <= dThird*Vin + cThird,
                                Vout >= dLow*Vin + cLow,
                                Vout >= dHigh*Vin + cHigh))

    elif a < 0:
        if Vlow <= 0 and Vhigh <=0:
            return Implies(And(Vin >= Vlow, Vin <= Vhigh),
                            And(Vout >= dThird*Vin + cThird,
                                Vout <= dLow*Vin + cLow,
                                Vout <= dHigh*Vin + cHigh))

        elif Vlow >=0 and Vhigh >=0:
            return Implies(And(Vin >= Vlow, Vin <= Vhigh),
                            And(Vout <= dThird*Vin + cThird,
                                Vout >= dLow*Vin + cLow,
                                Vout >= dHigh*Vin + cHigh))
            
def trapezoidBounds(a,Vin,Vout, Vlow, Vhigh):
    tanhFunVlow = tanhFun(a,Vlow)
    tanhFunVhigh = tanhFun(a,Vhigh)
    dLow = tanhFunder(a,Vlow)
    dHigh = tanhFunder(a,Vhigh)
    diff = Vhigh - Vlow
    if(diff == 0):
        diff = 1e-10
    dThird = (tanhFunVhigh - tanhFunVlow)/diff
    cLow = tanhFunVlow - dLow*Vlow
    cHigh = tanhFunVhigh - dHigh*Vhigh
    cThird = tanhFunVlow - dThird*Vlow

    if a > 0:
        if Vlow <= 0 and Vhigh <= 0:
            return Implies(Vin < Vlow, 
                            And(Vout > dLow*Vin + cLow,
                                Vout >= -1,
                                Vout < tanhFunVlow))
        elif Vlow >= 0 and Vhigh >=0:
            return Implies(Vin > Vhigh, 
                            And(Vout < dHigh*Vin + cHigh,
                                Vout <=  1,
                                Vout >  tanhFunVhigh))

    elif a < 0:
        if Vlow <= 0 and Vhigh <= 0:
            return Implies(Vin < Vlow, 
                            And(Vout < dLow*Vin + cLow,
                                Vout <= 1,
                                Vout > tanhFunVlow))
        elif Vlow >= 0 and Vhigh >=0:
            return Implies(Vin > Vhigh, 
                            And(Vout > dHigh*Vin + cHigh,
                                Vout >= - 1,
                                Vout <  tanhFunVhigh))
            
def funTest1():
    a = 1
    params = [0.3,0.1]
    
    s = Solver()
    x0 = Real('x0')
    y0 = Real('y0')
    triangleClaimPos = triangleBounds(a,x0,y0,0.0,0.5)
    triangleClaimNeg = triangleBounds(a,x0,y0,-0.5,0.0)
    trapezoidClaimPos = trapezoidBounds(a,x0,y0,0.0,0.5)
    trapezoidClaimNeg = trapezoidBounds(a,x0,y0,-0.5,0.0)
    
    s.add(triangleClaimPos)
    s.add(trapezoidClaimPos)
    s.add(triangleClaimNeg)
    s.add(trapezoidClaimNeg)
   
    s.add(y0 == params[0]*x0 + params[1])
    
    print "done adding constraints"
    ch = s.check()
    print "ch ", ch
    
    if ch == sat:
        m = s.model()
        print "model"
        print m
        
def funTest2():
    a = -5
    params = [0.0]
    
    s = Solver()
    x0 = Real('x0')
    x1 = Real('x1')
    y0 = Real('y0')
    y1 = Real('y1')
    #ItrartVar(x0, y0,0)
    #ItrartVar(x1, y1,1)

    
    
    triangleClaimPos0 = triangleBounds(a,x0,y0,0.0,0.5)
    triangleClaimNeg0 = triangleBounds(a,x0,y0,-0.5,0.0)
    trapezoidClaimPos0 = trapezoidBounds(a,x0,y0,0.0,0.5)
    trapezoidClaimNeg0 = trapezoidBounds(a,x0,y0,-0.5,0.0)
    
    triangleClaimPos1 = triangleBounds(a,x1,y1,0.0,0.5)
    triangleClaimNeg1 = triangleBounds(a,x1,y1,-0.5,0.0)
    trapezoidClaimPos1 = trapezoidBounds(a,x1,y1,0.0,0.5)
    trapezoidClaimNeg1 = trapezoidBounds(a,x1,y1,-0.5,0.0)
    
    s.add(triangleClaimPos0)
    s.add(trapezoidClaimPos0)
    s.add(triangleClaimNeg0)
    s.add(trapezoidClaimNeg0)
    
    s.add(triangleClaimPos1)
    s.add(trapezoidClaimPos1)
    s.add(triangleClaimNeg1)
    s.add(trapezoidClaimNeg1)
   
    s.add(y0 - x1 == -params[0])
    s.add(y1 - x0 == params[0])
    
    print "done adding constraints"
    ch = s.check()
    print "ch ", ch
    
    if ch == sat:
        m = s.model()
        print "model"
        print m
    
def oscTest():
    a = -5
    params = [1.0,0.5]
    
    s = Solver()
    x0 = Real('x0')
    x1 = Real('x1')
    x2 = Real('x2')
    x3 = Real('x3')
    xs = [x0,x1,x2,x3]
    
    y0 = Real('y0')
    y1 = Real('y1')
    y2 = Real('y2')
    y3 = Real('y3')
    ys = [y0,y1,y2,y3]
    
    z0 = Real('z0')
    z1 = Real('z1')
    z2 = Real('z2')
    z3 = Real('z3')
    zs = [z0,z1,z2,z3]
    
    lenV = 4
    for i in range(lenV):
        fwdInd = (i-1)%lenV
        ccInd = (i+lenV/2)%lenV
        if  fwdInd == 0 or fwdInd == 2:
            triangleClaimFwdNeg = triangleBounds(a,xs[fwdInd],ys[i],-0.5,0.0)
            s.add(triangleClaimFwdNeg)
        if fwdInd == 1 or fwdInd ==3:
            triangleClaimFwdPos = triangleBounds(a,xs[fwdInd],ys[i],0.0,0.5)
            s.add(triangleClaimFwdPos)
        #trapezoidClaimFwdPos = trapezoidBounds(a,xs[fwdInd],ys[i],0.0,0.5)
        #trapezoidClaimFwdNeg = trapezoidBounds(a,xs[fwdInd],ys[i],-0.5,0.0)
        
        if ccInd == 0 or ccInd == 2:
            triangleClaimCcNeg = triangleBounds(a,xs[ccInd],zs[i],-0.5,0.0)
            s.add(triangleClaimCcNeg)
        if ccInd == 1 or ccInd == 3:
            triangleClaimCcPos = triangleBounds(a,xs[ccInd],zs[i],0.0,0.5)
            s.add(triangleClaimCcPos)
        #trapezoidClaimCcPos = trapezoidBounds(a,xs[ccInd],zs[i],0.0,0.5)
        #trapezoidClaimCcNeg = trapezoidBounds(a,xs[ccInd],zs[i],-0.5,0.0)
        
        #s.add(triangleClaimFwdPos)
        #s.add(triangleClaimFwdNeg)
        #s.add(trapezoidClaimFwdPos)
        #s.add(trapezoidClaimFwdNeg)
        
        #s.add(triangleClaimCcPos)
        #s.add(triangleClaimCcNeg)
        #s.add(trapezoidClaimCcPos)
        #s.add(trapezoidClaimCcNeg)
        
        s.add((ys[i] - xs[i])*params[0] + (zs[i] - xs[i])*params[1] == 0)
    
    s.add(And(x0>=-0.5,x0<=0.5))
    s.add(And(x1>=-0.5,x1<=0.5))
    s.add(And(x2>=-0.5,x2<=0.5))
    s.add(And(x3>=-0.5,x3<=0.5))
    print "done adding constraints"
    ch = s.check()
    print "ch ", ch
    
    if ch == sat:
        m = s.model()
        print "model"
        print m

#funTest1()
funTest2()
#oscTest()
        