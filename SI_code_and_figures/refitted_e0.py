# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 12:13:04 2023

@author: shreyas
"""

import math
import numpy as np
from numpy import genfromtxt
from numpy import linalg as LA
import matplotlib.pyplot as plt
import warnings
import time

t1 = time.time()
#suppress warnings
warnings.filterwarnings('ignore')

r16e0data = genfromtxt(r'C:\Users\shrey\Documents\Lab\NO-scattering\GHM_NAH model\Refitting\no_fitting_data\NOAu\r=1.6\e0.csv', delimiter=',')

#Define unit conversion constants.
htoev = 27.2114
btoang = 0.529177249
evtoh = 0.0367493
angtob = 1.8897259886
fstoau = 41.34 
autofs = 0.02419
kB = 1.380649e-23 
Jtoev = 6.242e18

#Define beta in au from T in K
def Ttobetaau(T):
    beta = (1.0/kB/T)/Jtoev/evtoh
    return(beta)
#Define parameters for the model.
#mN = 25536, mO = 29170, mNO -> (25536 + 29170), MuNO -> (1/25536 + 1/29170)^-1,
R0 = 1.151*angtob
a0 = 2.7968/angtob
D0 = 6.61*evtoh
b0 = 1.9535/angtob
Z0 = -0.26876*angtob
C0 = 6.5713*evtoh
a1 = 2.5194/angtob
R1 = 1.295*angtob
D1 = 4.1528*evtoh 
a2 = 1.0015/angtob
Z1 = 1.235*angtob
D2 = 2.4171*evtoh
C1 = 8.9587*evtoh
abar = 10*angtob
temp = 300
#Define morse potential
def Vmorse(R,d,a): 
    return(d*(math.exp(-2.0*a*R) - 2.0*math.exp(-a*R)))
#Define ground and excited states for NAH and coupling
def U0(R,Z): 
    return(Vmorse(R-R0,D0,a0) + math.exp(-b0*(Z - Z0)) + C0)

def U1(R,Z): 
    return(Vmorse(R-R1,D1,a1) + Vmorse(Z-Z1,D2,a2) + C1)

def Vk(Z,Gamma,BW,Nb): 
  return(math.sqrt(Gamma/2.0/math.pi)*math.sqrt(BW/Nb)*(1.0 - math.tanh(Z/abar)))

#Set up lists for plotting U0 and U1
zrange= np.arange(0.5,6.1,0.1)
rrange= np.arange(0.9,2.2,0.05)

#U0r117 = [U0(1.17*angtob,zk*angtob)*htoev for zk in zrange]
#U1r117 = [U1(1.17*angtob,zk*angtob)*htoev for zk in zrange]
#U0r16 = [U0(1.6*angtob,zk*angtob)*htoev for zk in zrange]
#U1r16 = [U1(1.6*angtob,zk*angtob)*htoev for zk in zrange]
#U0z16 = [U0(rk*angtob,1.6*angtob)*htoev for rk in rrange]
#U1z16 = [U1(rk*angtob,1.6*angtob)*htoev for rk in rrange]

#Set up diagonalization for NAH.
def fermi(ee):
    if np.exp(Ttobetaau(temp)*ee) > 1e-10 :
        return(1.0/(1.0 + np.exp(Ttobetaau(temp)*ee)))
    else:
        return(1.0)

#List of metal states:
def MetalEK(BW,Nb):
    Ek = [-BW/2.0 + k*BW/Nb for k in range(Nb)]
    return(Ek)

#Define potential matrix for NAH hamiltonian at R,Z 
def VMatrix(R,Z,Gamma,BW,Nb):
    Vmat = np.zeros((Nb+1,Nb+1))  #Set up
    Vmat[0,0] = U1(R*angtob,Z*angtob) - U0(R*angtob,Z*angtob) #molecule part of NAH potential
    Ek = MetalEK(BW*evtoh,Nb) 
    for k in range(1,Nb+1): 
        Vmat[k,k] = Ek[k-1]
        Vmat[0,k] = Vk(Z*angtob,Gamma*evtoh,BW*evtoh,Nb) #couplinh b/w metal and molecule
        Vmat[k,0] = Vk(Z*angtob,Gamma*evtoh,BW*evtoh,Nb) 
    return(Vmat)

#Function that calculates the ground state of the NAH hamiltonian.
def E0calc(R,Z,Gamma,BW,Nb):
    evals, evecs = LA.eig(VMatrix(R,Z,Gamma,BW,Nb))
    evals = np.sort(evals)
    E0 = U0(R*angtob,Z*angtob) + np.sum([evals[k]*fermi(evals[k]) for k in range(len(evals))])
    return(E0)

#Get E0 for Gamma = 1.5 eV, BW = 100eV, Nb = 200
E0BW100 = np.zeros((len(zrange),len(rrange)))
E0BW100 = [[E0calc(rk,zk,1.5,100,200)*htoev for rk in rrange] for zk in zrange]

#Set correct zero of energy for each calculation.
#Here r16e0data[-1,1] is DFT E0 for (z = 6.02 A, R = 1.6 A) and
# E0BW100[-1][14] is the E0 at the same point calculated here. 
E0BW100 = E0BW100 - E0BW100[-1][14] + r16e0data[-1,1] 

E0G3p5=np.zeros((len(zrange),len(rrange)))

Gamma  = 3.5
#Get E0 for Gamma = 3.5 eV, BW = 7eV, Nb = 200
E0G3p5 = [[E0calc(rk,zk,Gamma,7.0,200)*htoev for rk in rrange] for zk in zrange]
#Set correct zero of energy for each calculation.
#Here r16e0data[-1,1] is DFT E0 for (z = 6.02 A, R = 1.6 A) and
# E0BW100[-1][14] is the E0 at the same point calculated here. 
E0G3p5 = E0G3p5 - E0G3p5[-1][14] + r16e0data[-1,1]

t2 = time.time()

rg, zg = np.meshgrid(rrange, zrange)
np.shape
print("Time taken in secs", t2-t1)

fig,ax=plt.subplots(1,1)
cp = ax.contourf(rg, zg, E0BW100,np.arange(-1,5))
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title("Γ=1.5 eV, ΔE = 100 eV")
ax.set_xlabel('R (Å)')
ax.set_ylabel('Z (Å)')
ax.set_xlim(right=2.0)
plt.rc('axes', titlesize=20)     # fontsize of the axes title
plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
plt.rc('ytick', labelsize=16)    # fontsize of the tick labels
plt.rc('legend', fontsize=16)    # legend fontsize
plt.rc('figure', titlesize=22)  # fontsize of the figure title
plt.show()
fig,ax=plt.subplots(1,1)
cp = ax.contourf(rg, zg, E0G3p5,np.arange(-1,5))
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title("Γ=3.5 eV, ΔE = 7 eV")
ax.set_xlabel('R (Å)')
ax.set_ylabel('Z (Å)')
plt.rc('axes', titlesize=20)     # fontsize of the axes title
plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
plt.rc('ytick', labelsize=16)    # fontsize of the tick labels
plt.rc('legend', fontsize=16)    # legend fontsize
plt.rc('figure', titlesize=22)  # fontsize of the figure title
ax.set_xlim(right=2.0)
plt.show()






