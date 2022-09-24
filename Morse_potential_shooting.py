# Morse_potential_shooting of Hydrogen Molecule

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
from scipy.integrate import simps
from scipy.optimize import brentq
L=1.0
h=0.01
e=np.arange(-1000.,0.1,1000.)
all_zeros=[0]
D=(38292/60.8296) # Dissociation energy of Hydrogen molecule in scale of (h_bar**2/((2*M*(ro)**2))
alpha=1.440 # Hydrogen molecule parameter
# Function for determining eigen energy of harmonic oscillator in shooting method
def psiAtL(e):
    po,p1=0.0,1.0
    x=-L+h
    while x<=+L:
        po,p1=p1,(2-h**2*(e-D*(np.exp(-(2*alpha*x))-2*np.exp(-alpha*x))))*p1-po
        x+=h
        return p1

# plotting Given Function
plt.plot(e,psiAtL(e))
plt.grid(True)
plt.show()

# Find the zeros for the symmetrical case
s=np.sign(psiAtL(e))
for i in range(len(s)-1):
    if s[i]+s[i+1]==0:      # Find zeros of this crazy function
        zero=brentq(psiAtL,e[i],e[i+1])
        all_zeros.append(zero)
# Energy eigen values in scale (h_bar**2/((2*M*(ro)**2))
print (all_zeros)

# Energy in electron volt
all_zeros=np.asarray(all_zeros)
print (all_zeros*60.8296,'cm^-1')
print (all_zeros*60.8296*1.2398*10**(-4),'ev')
#RK4 for solving wave function
psi=0.;psid=1.;# initial psi and first derivative of psi
I=np.array([psi,psid])
xi=-0.5;xf=+0.5;h1=0.001;x1=xi
e=all_zeros[0]
xs=[];psi_s=[]

def RK4(x1,I):
    p=dI(x1,I)
    q=dI(x1+h1/2,I+h1/2*p)
    r=dI(x1+h1/2,I+h1/2*q)
    s=dI(x1+h1,I+h1*r)  
    I+=(p+2*q+2*r+s)/6*h1
    x1+=h1
    return x1,I
def dI(x1,I):
    I=psi,psid
    dphi_dx=psid
    dpsid_dx=-(e-D*(np.exp(-(2*alpha*x1))-2*np.exp(-alpha*x1)))*psi
    return np.array([dphi_dx,dpsid_dx])
while x1<xf:
    x1,I=RK4(x1,I)
    psi,psid=I
    xs.append(x1)
    psi_s.append(psi)
    arr_xs=np.asarray(xs)
    arr_psi=np.asarray(psi_s)
    sq_arr_psi=arr_psi**2
# Normalization constant using simpsons rule
print ("Normalization constant",simps(sq_arr_psi,dx=h1))
N=1./simps(sq_arr_psi,dx=h1)
arr_psi=N*arr_psi
plt.plot(arr_xs,arr_psi)
plt.show()

























    






























     

