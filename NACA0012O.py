import pandas as pd
import datetime
import numpy as np
import matplotlib.pyplot as plt
def BoundrayCondition(m,n):
    X = np.ones((m,n))
    Y = np.ones((m,n))
    for i in range((m-1)//2):
        X[i][0] = 0.5*(1+np.cos(np.pi/((m-1)//2)*i))
        Y[i][0] = -0.6*(-0.1015*X[i][0]**4+0.2843*X[i][0]**3-0.3576*X[i][0]**2-0.1221*X[i][0]+0.2969*X[i][0]**(0.5))
        X[i][n-1] = 0.5+15*np.cos(np.pi/((m-1)//2)*i)
        Y[i][n-1] = -15*np.sin(np.pi/((m-1)//2)*i)  
    for j in range(((m-1)//2),m):
        X[j][0] = 0.5*(1+np.cos(np.pi/((m-1)//2)*j))
        Y[j][0] = 0.6*(-0.1015*X[j][0]**4+0.2843*X[j][0]**3-0.3576*X[j][0]**2-0.1221*X[j][0]+0.2969*X[j][0]**(0.5))
        X[j][n-1] = 0.5+15*np.cos(np.pi/((m-1)//2)*j)
        Y[j][n-1] = -15*np.sin(np.pi/((m-1)//2)*j)
    for i in range(m-1):
        for j in range(1,n-1):
            X[i][j] = X[i][0] + ((X[i][n-1]-X[i][0])/(n-1)) * j
            Y[i][j] = Y[i][0] + ((Y[i][n-1]-Y[i][0])/(n-1)) * j
    for k in range(n):
        X[m-1][k] = X[0][k]
        Y[m-1][k] = Y[0][k]
    return X,Y,m,n

def iteration(X,Y,m,n):
    X0 = np.ones((m,n))
    Y0 = np.ones((m,n))
    XZeta = np.ones((m,n));YZeta = np.ones((m,n));
    XEta = np.ones((m,n));YEta = np.ones((m,n));
    Alphf = np.ones((m,n));Beta = np.ones((m,n));gamma = np.ones((m,n))   
    Zeta = 0.025;Eta = 0.025
    k = 0
    while True:
        X0 = X.copy(); Y0 = Y.copy()
        for j in range(1,n-1):
            XZeta[m-1][j] = XZeta[0][j] = XZeta = (X[1][j]-X[m-2][j])/(2 * Zeta)
            YZeta[m-1][j] = YZeta[0][j] = YZeta = (Y[1][j]-Y[m-2][j])/(2 * Zeta)
            XEta[m-1][j] = XEta[0][j] = XEta = (X[0][j+1]-X[0][j-1])/(2 * Eta) 
            YEta[m-1][j] = YEta[0][j] = YEta = (Y[0][j+1]-Y[0][j-1])/(2 * Eta)
            Alphf[m-1][j] = Alphf[0][j] = Alphf =  XEta**2 + YEta**2
            Beta[m-1][j] = Beta[0][j] = Beta = XZeta * XEta + YZeta * YEta
            gamma[m-1][j] = gamma[0][j] = gamma = XZeta ** 2 + YZeta**2
            X[m-1][j]=X[0][j] = (Alphf*(X[1][j]+X[m-2][j])+gamma*(X[0][j+1]+X[0][j-1])-(Beta*(X[1][j+1]+X[m-2][j-1]-X[1][j-1]-X[m-2][j+1]))/2)/(2*(Alphf+gamma))
            Y[m-1][j]=Y[0][j] = (Alphf*(Y[1][j]+Y[m-2][j])+gamma*(Y[0][j+1]+Y[0][j-1])-(Beta*(Y[1][j+1]+Y[m-2][j-1]-Y[1][j-1]-Y[m-2][j+1]))/2)/(2*(Alphf+gamma))
            for i in range(1,m-1):
                XZeta[i][j] = XZeta = (X[i+1][j]-X[i-1][j])/(2 * Zeta)
                YZeta[i][j] = YZeta = (Y[i+1][j]-Y[i-1][j])/(2 * Zeta)
                XEta[i][j] = XEta = (X[i][j+1]-X[i][j-1])/(2 * Eta) 
                YEta[i][j] = YEta = (Y[i][j+1]-Y[i][j-1])/(2 * Eta)
                Alphf[i][j] = Alphf =  XEta**2 + YEta**2
                Beta[i][j] = Beta = XZeta * XEta + YZeta * YEta
                gamma[i][j] = gamma = XZeta ** 2 + YZeta ** 2
                X[i][j] = (Alphf*(X[i+1][j]+X[i-1][j])+gamma*(X[i][j+1]+X[i][j-1])-(Beta*(X[i+1][j+1]+X[i-1][j-1]-X[i+1][j-1]-X[i-1][j+1]))/2)/(2*(Alphf+gamma))
                Y[i][j] = (Alphf*(Y[i+1][j]+Y[i-1][j])+gamma*(Y[i][j+1]+Y[i][j-1])-(Beta*(Y[i+1][j+1]+Y[i-1][j-1]-Y[i+1][j-1]-Y[i-1][j+1]))/2)/(2*(Alphf+gamma))        
        errorX = X0 - X
        errorY = Y0 - Y
        maxX = np.fabs(errorX).max()
        maxY = np.fabs(errorY).max()
        k += 1
        if ((maxX<=1e-7) and (maxY<=1e-7)):
            print(maxX,maxY,k)
            break
    XZeta[m-1][0] = XZeta[0][0] = (X[1][0] - X[m-2][0])/(2*Zeta)
    YZeta[m-1][0]=YZeta[0][0] = (Y[1][0] - Y[m-2][0])/(2*Zeta)
    XEta[m-1][0] = XEta[0][0] = (4*X[0][1]-X[0][2]-3*X[0][0])/(2*Eta)
    YEta[m-1][0] = YEta[0][0] = (4*Y[0][1]-Y[0][2]-3*Y[0][0])/(2*Eta)
    Alphf[0][0] = Alphf[m-1][0] = XEta[0][0] ** 2 + YEta[0][0] ** 2
    Beta[0][0] = Beta[m-1][0] = XZeta[0][0]*XEta[0][0] +YZeta[0][0]*YEta[0][0]
    gamma[0][0] = gamma[m-1][0] = XZeta[0][0] ** 2 + YZeta[0][0] ** 2
    for i in range(1,m-1):
        XZeta[i][0] = (X[i+1][0] - X[i-1][0])/(2*Zeta)
        YZeta[i][0] = (Y[i+1][0] - Y[i-1][0])/(2*Zeta)
        XEta[i][0] = (4*X[i][1] - X[i][2] - 3*X[i][0])/(2*Eta)
        YEta[i][0] = (4*Y[i][1] - Y[i][2] - 3*Y[i][0])/(2*Eta)
        Alphf[i][0] = XEta[i][0] ** 2 + YEta[i][0] ** 2
        Beta[i][0] = XZeta[i][0] * XEta[i][0] + YZeta[i][0] * YEta[i][0]
        gamma[i][0] = XZeta[i][0] ** 2 + YZeta[i][0] ** 2
    return X,Y,Alphf,Beta,gamma,XZeta,XEta,YZeta,YEta

def FlowField(X,Y,V,m,n,Alphf,Beta,gamma):
    k = 0
    flowfield = np.ones((m,n))
    flowfield = V * X
    while True:
        flowfield0 = flowfield.copy()
        flowfield[m-1][0] = flowfield[0][0] = (4*gamma[0][0]*flowfield[0][1]-gamma[0][0]*flowfield[0][2]-Beta[0][0]*flowfield[1][0]+Beta[0][0]*flowfield[m-2][0])/(3*gamma[0][0])
        for i in range(1,m-1):
            flowfield[i][0] = (4*gamma[i][0]*flowfield[i][1]-gamma[i][0]*flowfield[i][2]-Beta[i][0]*flowfield[i+1][0]+Beta[i][0]*flowfield[i-1][0])/(3*gamma[i][0])
        for j in range(1,n-1):
            flowfield[m-1][j]=flowfield[0][j] = (Alphf[0][j]*(flowfield[1][j]+flowfield[m-2][j])+gamma[0][j]*(flowfield[0][j+1]+flowfield[0][j-1])-(Beta[0][j]*(flowfield[1][j+1]+flowfield[m-2][j-1]-flowfield[1][j-1]-flowfield[m-2][j+1]))/2)/(2*(Alphf[0][j]+gamma[0][j]))
            for i in range(1,m-1):
                flowfield[i][j] = (Alphf[i][j]*(flowfield[i+1][j]+flowfield[i-1][j])+gamma[i][j]*(flowfield[i][j+1]+flowfield[i][j-1])-(Beta[i][j]*(flowfield[i+1][j+1]+flowfield[i-1][j-1]-flowfield[i+1][j-1]-flowfield[i-1][j+1]))/2)/(2*(Alphf[i][j]+gamma[i][j]))
        errorFlow = flowfield0 - flowfield
        maxFlowError = np.fabs(errorFlow).max()
        k += 1
        if (maxFlowError <= 1e-7):
            print(maxFlowError,k)
            break
    return flowfield
def ParameterSolve(Phi,XZeta,XEta,YZeta,YEta,V):
    Zeta = 0.025;Eta = 0.025
    Jobic = np.ones((m,n))
    Jobic = XZeta*YEta-XEta*YZeta
    Jobic[:,n-1] = 1
    PhiZeta = np.ones((m,n));PhiEta = np.ones((m,n))
    u = np.ones((m,n));v = np.ones((m,n));cp = np.ones((m,n));velocity = np.ones((m,n))
    for j in range(1,n-1):
        PhiZeta[m-1][j] = PhiZeta[0][j] = (Phi[1][j]-Phi[m-2][j])/(2 * Zeta)
        PhiEta[m-1][j] = PhiEta[0][j] = (Phi[0][j+1]-Phi[0][j-1])/(2 * Eta) 
        for i in range(1,m-1):
            PhiZeta[i][j] = (Phi[i+1][j]-Phi[i-1][j])/(2 * Zeta)
            PhiEta[i][j] = (Phi[i][j+1]-Phi[i][j-1])/(2 * Eta) 
    PhiZeta[m-1][0] = PhiZeta[0][0] = (Phi[1][0] - Phi[m-2][0])/(2*Zeta)
    PhiEta[m-1][0] = PhiEta[0][0] = (4*Phi[0][1]-Phi[0][2]-3*Phi[0][0])/(2*Eta)
    for i in range(1,m-1):
        PhiZeta[i][0] = (Phi[i+1][0] - Phi[i-1][0])/(2*Zeta)
        PhiEta[i][0] = (4*Phi[i][1] - Phi[i][2] - 3*Phi[i][0])/(2*Eta)
    u = (PhiZeta * YEta - PhiEta * YZeta)/Jobic
    v = (PhiEta * XZeta - PhiZeta * XEta)/Jobic
    u[:,n-1] = V;v[:,n-1] = 0
    velocity = u**2 + v**2
    cp = np.ones((m,n)) - (u**2 + v**2)/(V**2)
    return u,v,cp,velocity
#主程序
start = datetime.datetime.now()
X,Y,m,n = BoundrayCondition(15,4)
plt.plot(X, Y)
plt.show()  
X0,Y0,Alphf,Beta,gamma,XZeta,XEta,YZeta,YEta= iteration(X,Y,m,n)
end1 = datetime.datetime.now()
print(end1-start)
flowfield = FlowField(X0,Y0,20.0,m,n,Alphf,Beta,gamma)
U,V,Cp,Velocity = ParameterSolve(flowfield,XZeta,XEta,YZeta,YEta,20.0)
end = datetime.datetime.now()
print(end-end1)

X0.resize([m*n,1]);Y0.resize([m*n,1]);U.resize([m*n,1]);V.resize([m*n,1]);Cp.resize([m*n,1]);Velocity.resize([m*n,1])
dataX = pd.DataFrame(X0,columns=None,index=None);dataY = pd.DataFrame(Y0,columns=None,index=None)
dataU = pd.DataFrame(U,columns=None,index=None);dataV = pd.DataFrame(V,columns=None,index=None)
dataCp = pd.DataFrame(Cp,columns=None,index=None);dataVelocity = pd.DataFrame(Velocity,columns=None,index=None)
dataX.to_csv('dataX.csv');dataY.to_csv('dataY.csv')
dataU.to_csv('dataU.csv');dataV.to_csv('dataV.csv')
dataCp.to_csv('dataCp.csv');dataVelocity.to_csv('dataVelocity.csv')
