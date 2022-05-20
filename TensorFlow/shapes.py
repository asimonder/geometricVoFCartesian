import scipy
import numpy as np
from scipy import integrate
import os
import sys
import pickle

class stencil:
    def __init__(self,theta,fun,lenX=5,lenY=5,N=1000,evalPoint="Centroid",XOrigin=np.zeros(2)):
        self.alpha = np.zeros((lenX,lenY))
        h=1./N
        xP=fun.x(theta)
        yP=fun.y(theta)
        xMin=-np.float64((lenX-1))/2.
        xMax=-xMin+0.01
        yMin=-np.float64((lenY-1))/2.
        yMax=-yMin+0.01
        self.xL=np.arange(np.floor(xP-XOrigin[0])+xMin+0.5,np.floor(xP-XOrigin[0])+xMax+0.5,1.)+XOrigin[0]
        self.yL=np.arange(np.floor(yP-XOrigin[1])+yMin+0.5,np.floor(yP-XOrigin[1])+yMax+0.5,1.)+XOrigin[1]
        self.fun=fun
        self.theta=np.nan
        iC=int((lenX-1)/2)
        jC=int((lenY-1)/2)
        if evalPoint=="Nearest":
            xD=fun.nearestPoint(self.xL[iC],self.yL[jC],dx=1.,N=100)
            self.theta=fun.theta(xD[0],xD[1])
        elif evalPoint=="Centroid":
            #print("evaluating at the center of the cell.")
            self.theta=fun.theta(self.xL[iC],self.yL[jC])
        elif evalPoint=="Seeded":
            self.theta=theta
        else:
            ValueError("Type of evaluation point is unkown.")
            
        self.K=fun.curvature(self.theta)
        self.n=fun.normal(self.theta)
         
        for i in range(lenX):
            for j in range(lenY):
                self.alpha[i,j]=fun.alphaCell(self.xL[i],self.yL[j],1.0,N)
                
    def convertZone(self,invertCurvature=True):
        K=1.*self.K #fun.curvature(self.theta)
        n=1.*self.n #fun.normal(self.theta)
        alpha2=1*self.alpha
        iFac=1;jFac=1
        #xL=1*self.xL
        #yL=1*self.yL
        if n[0]<0:
            iFac=-1
            self.n[0]=np.fabs(n[0])
            #xL=flip(self.xL)
        if n[1]<0:
            jFac=-1
            self.n[1]=np.fabs(n[1])
            #yL=flip(self.yL)
        if K<0 and invertCurvature:
            alpha2=1-alpha2
            iFac*=-1;jFac*=-1
            self.K=np.fabs(K)
        iMax=int((self.alpha.shape[0]-1)/2)
        jMax=int((self.alpha.shape[1]-1)/2)
        
        alphaR=1*alpha2
        if np.fabs(n[0])<=np.fabs(n[1]):
            for i in range(-iMax,iMax+1):
                for j in range(-jMax,jMax+1):
                    alpha2[i+iMax,j+jMax]=alphaR[i*iFac+iMax,j*jFac+jMax]
        else:
            self.n[0]=np.fabs(n[1])
            self.n[1]=np.fabs(n[0])
            for i in range(-iMax,iMax+1):
                for j in range(-jMax,jMax+1):
                    alpha2[i+iMax,j+jMax]=alphaR[j*iFac+jMax,i*jFac+iMax]                  
        self.alpha=alpha2
        
        #return alpha2
def writeDataset(wDir,fileName,data):
    if not os.path.exists(wDir):
        os.makedirs(wDir)
    fName="%s/%s.pkl"%(wDir,fileName)
    N=data["kappa"].shape[0]
    if os.path.exists(fName):
        file = open(fName,'rb')
        data2=pickle.load(file)
        for d in data2:
            data[d]=np.append(data2[d],data[d],axis=0)
        file.close()
    with open(fName, "wb") as fp:
        pickle.dump(data,fp)
    fName="%s/%s_Size.txt"%(wDir,fileName)
    if os.path.exists(fName):
        Nd=np.loadtxt("%s"%fName)
    else:
        Nd=np.zeros(1)
    Nd+=N
    Nd=np.ones(1)*Nd
    np.savetxt("%s"%fName,Nd,fmt='%d')


        
from scipy.optimize import fsolve 
from scipy.optimize import bisect 
from scipy.optimize import brentq
from scipy.optimize import minimize
class shape:
    def x(self,theta):
        return 0.
    def x_t(self,theta):
        return 0.
    def x_tt(self,theta):
        return 0.
    def y(self,theta):
        return 0.
    def y_t(self,theta):
        return 0.
    def y_tt(self,theta):
        return 0.
    def theta(self,x,y):
        return 0.
    def normal(self,t):
        n_x=self.y_t(t)
        n_y=-self.x_t(t)
        mag=np.sqrt(n_x*n_x+n_y*n_y)
        return np.array([n_x,n_y])/mag
    def normalAngle(self,t): 
        n_x=self.y_t(t)
        n_y=-self.x_t(t)
        return arctan2(n_y,n_x)
    def curvature(self,t):
        x_t=self.x_t(t)
        x_tt=self.x_tt(t)
        y_t=self.y_t(t)
        y_tt=self.y_tt(t)
        return (x_t*y_tt-y_t*x_tt)/(x_t**2+y_t**2)**(1.5)
    def distance(self,x,y,dx):
        t=self.thetaMin(x,y,dx)
        return np.sqrt((x-self.x(t))**2.+(y-self.y(t))**2.)    
    def thetaMin(self,x,y,dx):
        #J=lambda t: np.sqrt((x-self.x(t))**2.+(y-self.y(t))**2.)
        #res=minimize(J, [self.theta(x,y)],method="CG")
        #return res.x
        J=lambda t: self.x_t(t)*(self.x(t)-x)+self.y_t(t)*(self.y(t)-y) 
        t0=self.theta(x,y)
        #dt=np.fabs(self.x(dx)-self.x(0.))
        #tL=np.linspace(t-dt,t+dt,50)
        #for t in tL:
        #    if J(t)>
        root=fsolve(J,t0)[0]
        if not np.isclose(J(root), 0.0):
            ValueError("root finding with fsolve failed.")
        return root
    def alphaCell(self,xC,yC,dx,N=100):
        #d=self.distance(xC,yC)
        #if d>1./np.sqrt(2)*dx:
        #    return self.value(xC,yC)
        h=dx/N
        xL=np.linspace(xC-dx/2.,xC+dx/2.,N+1)
        yL=np.linspace(yC-dx/2.,yC+dx/2.,N+1)
        sumA=np.zeros(4)
        yP=np.ones(N+1)*yL[0]
        sumA[0]=np.mean(self.value(xL,yP))
        yP=np.ones(N+1)*yL[-1]
        sumA[1]=np.mean(self.value(xL,yP))
        xP=np.ones(N+1)*xL[0]
        sumA[2]=np.mean(self.value(xP,yL))
        xP=np.ones(N+1)*xL[-1]
        sumA[3]=np.mean(self.value(xP,yL))
        sum=np.mean(sumA)
        if sum==0:
            return 0.
        elif sum==1:
            return 1.
        else:
            H=np.zeros(N+1)
            for n in range(N+1):
                xP=np.ones(N+1)*xL[n]
                H[n]=scipy.integrate.simps(self.value(xP,yL),yL)
                #H[n]=numpy.trapz(self.value(xP,yL),dx=h)
            return scipy.integrate.simps(H,xL)/dx/dx
        #return numpy.trapz(H,dx=h)
        
    def alphaCellMidPoint(self,xC,yC,dx,N=100):
        d=self.distance(xC,yC,dx)
        if d>1./np.sqrt(2)*dx:
            return self.value(xC,yC)
        Nf=np.float128(N)
        h=dx/Nf
        xL=np.linspace(xC-dx/2.,xC+dx/2.,N+1)
        yL=np.linspace(yC-dx/2.,yC+dx/2.,N+1)
        x2,y2=np.meshgrid(xL[:-1]+h/2.,yL[:-1]+h/2.)
        A=np.sum(self.value(x2.flatten(),y2.flatten()))
        return A/Nf/Nf
    
    def alphaCellMidPointOld(self,xC,yC,dx,N=100):
        #d=self.distance(xC,yC)
        #if d>1./np.sqrt(2)*dx:
        #    return self.value(xC,yC)
        Nf=np.float128(N)
        h=dx/Nf
        xL=np.linspace(xC-dx/2.,xC+dx/2.,N+1)
        yL=np.linspace(yC-dx/2.,yC+dx/2.,N+1)
        x2,y2=np.meshgrid(xL[:-1]+h/2.,yL[:-1]+h/2.)
        A=np.sum(self.value(x2.flatten(),y2.flatten()))
        return A/Nf/Nf

    def nearestPoint(self,xC,yC,dx,N=100):
        xP=np.zeros(2)
        xL=np.linspace(xC-dx/2.,xC+dx/2.,N+1)
        yL=np.linspace(yC-dx/2.,yC+dx/2.,N+1)
        d0=dx**2.
        for x in xL:
            for y in yL:
                if self.value(x,y)==1.:
                    d=(x-xC)**2.+(y-yC)**2.
                    if d<d0:
                        xP[0]=x
                        xP[1]=y
                        d=d0
        return xP
    
class circle(shape):
    def __init__(self, radius,origin=np.zeros(2)):
        self.R=radius
        self.Ox=origin[0]
        self.Oy=origin[1]
        #self.value=lambda x,y: (np.sign(self.R-np.sqrt(x*x+y*y))+1)/2.
        
    def theta(self,x,y):
        return np.arctan2(y-self.Oy,x-self.Ox)
        
    def value(self,x,y):
        xP=x-self.Ox
        yP=y-self.Oy
        return (np.sign(self.R**2.-(xP*xP+yP*yP))+1)/2.
    
    #def thetaMin(self,x,y):
    #    return self.theta(x,y)
        
        
    def x(self,theta):
        return self.R*np.cos(theta)+self.Ox
    def y(self,theta):
        return self.R*np.sin(theta)+self.Oy
    def normal(self,theta): 
        n_x=self.R*np.cos(theta) #self.y_t(theta)
        n_y=self.R*np.sin(theta) #-self.x_t(theta)
        mag=np.sqrt(n_x*n_x+n_y*n_y)
        return np.array([n_x,n_y])/mag
    def normalAngle(self,theta):        
        return theta #atan2(y-self.O[1],x-self.O[0])
    def curvature(self,theta):
        return 1/self.R
    
class ellipse(shape):
    def __init__(self,a,b,psi=0.):
        self.a=a
        self.b=b
        self.psi=psi
        #self.value=lambda x,y: (np.sign(self.R-np.sqrt(x*x+y*y))+1)/2.
        
    def theta(self,x,y):
        xR=x*np.cos(-self.psi)-y*np.sin(-self.psi)
        yR=x*np.sin(-self.psi)+y*np.cos(-self.psi)
        return np.arctan2(yR/self.b,xR/self.a)
     
    def value(self,x,y):
        theta=self.theta(x,y)
        R=x*x+y*y
        xE=self.x(theta)
        yE=self.y(theta)
        RE=xE*xE+yE*yE
        return (np.sign(RE-R)+1)/2.
        
    def x(self,theta):
        return self.a*np.cos(theta)*np.cos(self.psi)-self.b*np.sin(theta)*np.sin(self.psi)
    def x_t(self,theta):
        return -self.a*np.sin(theta)*np.cos(self.psi)-self.b*np.cos(theta)*np.sin(self.psi)
    def x_tt(self,theta):
        return -self.a*np.cos(theta)*np.cos(self.psi)+self.b*np.sin(theta)*np.sin(self.psi)
    def y(self,theta):
        return self.a*np.cos(theta)*np.sin(self.psi)+self.b*np.sin(theta)*np.cos(self.psi)
    def y_t(self,theta):
        return -self.a*np.sin(theta)*np.sin(self.psi)+self.b*np.cos(theta)*np.cos(self.psi)
    def y_tt(self,theta):
        return -self.a*np.cos(theta)*np.sin(self.psi)-self.b*np.sin(theta)*np.cos(self.psi)
    
class star(shape):
    def __init__(self,R,A,n,psi):
        self.R=R
        self.A=A
        self.n=n
        self.psi=psi
        
    def theta(self,x,y):
        xR=x*np.cos(-self.psi)-y*np.sin(-self.psi)
        yR=x*np.sin(-self.psi)+y*np.cos(-self.psi)
        #return np.arctan2(yR/self.b,xR/self.a)
        return np.arctan2(yR,xR)
     
    def value(self,x,y):
        theta=self.theta(x,y)
        R=x*x+y*y
        xE=self.x(theta)
        yE=self.y(theta)
        RE=xE*xE+yE*yE
        return (np.sign(RE-R)+1)/2.
        
    def x(self,theta):
        a=self.R+self.A*np.cos(self.n*theta)
        return a*np.cos(theta)*np.cos(self.psi)-a*np.sin(theta)*np.sin(self.psi)
    def x_t(self,theta):
        a=self.R+self.A*np.cos(self.n*theta)
        a_t=-self.A*self.n*np.sin(self.n*theta)
        return -a*np.sin(theta)*np.cos(self.psi)-a*np.cos(theta)*np.sin(self.psi)\
        +a_t*np.cos(theta)*np.cos(self.psi)-a_t*np.sin(theta)*np.sin(self.psi)
    def x_tt(self,theta):
        a=self.R+self.A*np.cos(self.n*theta)
        a_t=-self.A*self.n*np.sin(self.n*theta)
        a_tt=-self.A*self.n*self.n*np.cos(self.n*theta)
        x_tt=-a_t*np.sin(theta)*np.cos(self.psi)-a_t*np.cos(theta)*np.sin(self.psi)\
            +a_tt*np.cos(theta)*np.cos(self.psi)-a_tt*np.sin(theta)*np.sin(self.psi)
        x_tt+=-a*np.cos(theta)*np.cos(self.psi)+a*np.sin(theta)*np.sin(self.psi)\
            -a_t*np.sin(theta)*np.cos(self.psi)-a_t*np.cos(theta)*np.sin(self.psi)
        return x_tt
    def y(self,theta):
        a=self.R+self.A*np.cos(self.n*theta)
        return a*np.cos(theta)*np.sin(self.psi)+a*np.sin(theta)*np.cos(self.psi)
    def y_t(self,theta):
        a=self.R+self.A*np.cos(self.n*theta)
        a_t=-self.A*self.n*np.sin(self.n*theta)
        return -a*np.sin(theta)*np.sin(self.psi)+a*np.cos(theta)*np.cos(self.psi)\
        +a_t*np.cos(theta)*np.sin(self.psi)+a_t*np.sin(theta)*np.cos(self.psi)
    def y_tt(self,theta):
        a=self.R+self.A*np.cos(self.n*theta)
        a_t=-self.A*self.n*np.sin(self.n*theta)
        a_tt=-self.A*self.n*self.n*np.cos(self.n*theta)
        y_tt=-a_t*np.sin(theta)*np.sin(self.psi)+a_t*np.cos(theta)*np.cos(self.psi)\
            +a_tt*np.cos(theta)*np.sin(self.psi)+a_tt*np.sin(theta)*np.cos(self.psi)
        y_tt+=-a*np.cos(theta)*np.sin(self.psi)-a*np.sin(theta)*np.cos(self.psi)\
            -a_t*np.sin(theta)*np.sin(self.psi)+a_t*np.cos(theta)*np.cos(self.psi)
        return y_tt
    
    
class sinWave(shape):
    def __init__(self,A,k,psi=0.):
        self.k=k
        self.A=A
        self.psi=psi
        self.lamb=2*np.pi/k
        #self.value=lambda x,y: (np.sign(self.R-np.sqrt(x*x+y*y))+1)/2.
        
    def theta(self,x,y):
        #xR=x*np.cos(-self.psi)-y*np.sin(-self.psi)
        xR=x*np.cos(-self.psi)-y*np.sin(-self.psi)
        return -xR
    
    def value(self,x,y):
        xR=x*np.cos(-self.psi)-y*np.sin(-self.psi)
        yR=x*np.sin(-self.psi)+y*np.cos(-self.psi)
        yW=self.A*np.cos(self.k*xR)
        return (np.sign(yW-yR)+1)/2.      
        
    def x(self,t):
        return -t*np.cos(self.psi)-self.A*np.cos(self.k*-t)*np.sin(self.psi) #(theta-phi)/2*pi*self.lamb
    def x_t(self,t):
        return -np.cos(self.psi)-self.A*self.k*np.sin(self.k*-t)*np.sin(self.psi)
    def x_tt(self,t):
        return +self.A*self.k*self.k*np.cos(self.k*-t)*np.sin(self.psi)
    def y(self,t):
        return -t*np.sin(self.psi)+self.A*np.cos(self.k*-t)*np.cos(self.psi)
    def y_t(self,t):
        return -np.sin(self.psi)+self.A*self.k*np.sin(self.k*-t)*np.cos(self.psi)
    def y_tt(self,t):
        return -self.A*self.k*self.k*np.cos(self.k*-t)*np.cos(self.psi)
    
class spectralWave(shape):
    def __init__(self,A,k,phi,psi=0.):
        self.k=k
        self.A=A
        self.phase=phi
        self.psi=psi
        self.lamb=2*np.pi/k
        #self.value=lambda x,y: (np.sign(self.R-np.sqrt(x*x+y*y))+1)/2.
    
    def eta(self,x):
        eta=0.
        for i in range(self.A.shape[0]):
            eta+=self.A[i]*np.cos(self.k[i]*x+self.phase[i])
        return eta
        #theta=self.k*x+self.phase
        #eta=np.sum(self.A*np.cos(theta))
        #print(eta)
        #return eta
    def eta_x(self,x):
        eta_x=0.
        for i in range(self.A.shape[0]):
            eta_x+=self.A[i]*self.k[i]*np.sin(self.k[i]*x+self.phase[i])        
        return eta_x
        #return np.sum(-self.A*self.k*np.sin(self.k*x+self.phase))
    def eta_xx(self,x):
        eta_xx=0.
        for i in range(self.A.shape[0]):
            eta_xx-=self.A[i]*self.k[i]**2.*np.cos(self.k[i]*x+self.phase[i])
        return eta_xx
        #return np.sum(-self.A*self.k*self.k*np.cos(self.k*x+self.phase))
        
    def theta(self,x,y):
        xR=x*np.cos(-self.psi)-y*np.sin(-self.psi)
        return -xR
    
    def value(self,x,y):
        xR=x*np.cos(-self.psi)-y*np.sin(-self.psi)
        yR=x*np.sin(-self.psi)+y*np.cos(-self.psi)
        yW=self.eta(xR)
        return (np.sign(yW-yR)+1)/2. 

        
    def x(self,t):
        return -t*np.cos(self.psi)-self.eta(-t)*np.sin(self.psi) #(theta-phi)/2*pi*self.lamb
    def x_t(self,t):
        return -np.cos(self.psi)-self.eta_x(-t)*np.sin(self.psi)
    def x_tt(self,t):
        return -self.eta_xx(-t)*np.sin(self.psi)
    def y(self,t):
        return -t*np.sin(self.psi)+self.eta(-t)*np.cos(self.psi)
    def y_t(self,t):
        return -np.sin(self.psi)+self.eta_x(-t)*np.cos(self.psi)
    def y_tt(self,t):
        return self.eta_xx(-t)*np.cos(self.psi)
    
class JONSWAPWave(spectralWave):
    def __init__(self,Hs,lambda_p,gamma,Nkx=100,psi=0.):
        self.psi=psi
        beta=0
        g=9.81
        alpha=1.
        omega_p=np.sqrt(g*2*np.pi/lambda_p)
        #lambda_p=2*pi/(omega_p**2./g)
        #omega_p=2*pi/Tp
        #print("lambda_p=",lambda_p)

        omegaL=np.arange(0.01,5.01,0.01)*omega_p
        dOmega=omegaL[1]-omegaL[0]
        G=0
        for i in range(omegaL.shape[0]):
            G+=self.Gw(omegaL[i],Hs,omega_p,gamma,alpha)
        Hs2=4*np.sqrt(G*dOmega)
        alpha=(Hs/Hs2)**2.
        G=0
        GL=np.zeros(omegaL.shape[0])
        for i in range(omegaL.shape[0]):
            G+=self.Gw(omegaL[i],Hs,omega_p,gamma,alpha)
            GL[i]=self.Gw(omegaL[i],Hs,omega_p,gamma,alpha)
        Hs2=4*np.sqrt(G*dOmega)
        #print("Hs=",Hs2)

        Lx=8*lambda_p
        x=np.arange(-Lx/2,Lx/2.+0.01,1.)
        dx=1. 

        dkx=2*np.pi/Lx
        self.A=np.zeros(Nkx+1)
        self.k=np.zeros(Nkx+1)
        self.phase=np.zeros(Nkx+1)
        for i in range(1,Nkx+1): 
            self.k[i]=i*dkx
            omega=np.sqrt(g*self.k[i])
            self.A[i]=np.sqrt(2*np.sqrt(g/self.k[i])*self.Gw(omega,Hs,omega_p,gamma,alpha)*dkx)
            self.phase[i]=np.random.rand(1)*2*np.pi

    
    def Gw(self,omega,Hs,omega_p,gamma,alpha):
        sigma=0.07
        if omega>omega_p:
            sigma=0.09
        h=gamma**np.exp(-1/2*((omega-omega_p)/(sigma*omega_p))**2.)       
        return alpha*Hs**2*omega_p**4.*omega**(-5)*np.exp(-5/4.*(omega/omega_p)**(-4))*h
    

class PowerWave(spectralWave):
    def __init__(self,Hs,lambda_p,Nkp=10,psi=0.):
        self.psi=psi
        self.k_p=2*np.pi/lambda_p
        self.Hs=Hs

        dk=0.1*self.k_p
        self.k=np.arange(dk,self.k_p*Nkp+0.00001,dk)
        Nk=self.k.shape[0]
        E=np.exp(-1.25*(self.k_p/self.k)**4.)*self.k**-5.
        Hs2=4*np.sqrt(np.sum(E*dk))
        self.alpha=Hs/Hs2
        self.A=self.alpha*np.sqrt(2*E*dk)
        self.phase=np.random.uniform(0,2*np.pi,Nk)
        self.E=E

