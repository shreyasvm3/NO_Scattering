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

#Read in r=1.17 A data
r117e0data = genfromtxt(r'C:\Users\shrey\Documents\Lab\NO-scattering\GHM_NAH model\Refitting\no_fitting_data\NOAu\r=1.17\e0.csv', delimiter=',')
r117h00data = genfromtxt(r'C:\Users\shrey\Documents\Lab\NO-scattering\GHM_NAH model\Refitting\no_fitting_data\NOAu\r=1.17\h00.csv', delimiter=',')
r117h11data = genfromtxt(r'C:\Users\shrey\Documents\Lab\NO-scattering\GHM_NAH model\Refitting\no_fitting_data\NOAu\r=1.17\h11.csv', delimiter=',')

#Read in r=1.6 A data
r16e0data = genfromtxt(r'C:\Users\shrey\Documents\Lab\NO-scattering\GHM_NAH model\Refitting\no_fitting_data\NOAu\r=1.6\e0.csv', delimiter=',')
r16h00data = genfromtxt(r'C:\Users\shrey\Documents\Lab\NO-scattering\GHM_NAH model\Refitting\no_fitting_data\NOAu\r=1.6\h00.csv', delimiter=',')
r16h11data = genfromtxt(r'C:\Users\shrey\Documents\Lab\NO-scattering\GHM_NAH model\Refitting\no_fitting_data\NOAu\r=1.6\h11.csv', delimiter=',')

#Read in z=1.6 A data
z16e0data = genfromtxt(r'C:\Users\shrey\Documents\Lab\NO-scattering\GHM_NAH model\Refitting\no_fitting_data\NOAu\z=1.6\e0.csv', delimiter=',')
z16h00data = genfromtxt(r'C:\Users\shrey\Documents\Lab\NO-scattering\GHM_NAH model\Refitting\no_fitting_data\NOAu\z=1.6\h00.csv', delimiter=',')
z16h11data = genfromtxt(r'C:\Users\shrey\Documents\Lab\NO-scattering\GHM_NAH model\Refitting\no_fitting_data\NOAu\z=1.6\h11.csv', delimiter=',')

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

def Vk(Z,Gamma,BW,Nb,ab): 
  return(math.sqrt(Gamma/2.0/math.pi)*math.sqrt(BW/Nb)*(1.0 - math.tanh(Z/ab)))

#Set up lists for plotting U0 and U1
zrange= np.arange(0.5,6,0.01)
rrange= np.arange(0.9,2,0.01)
U0r117 = [U0(1.17*angtob,zk*angtob)*htoev for zk in zrange]
U1r117 = [U1(1.17*angtob,zk*angtob)*htoev for zk in zrange]
U0r16 = [U0(1.6*angtob,zk*angtob)*htoev for zk in zrange]
U1r16 = [U1(1.6*angtob,zk*angtob)*htoev for zk in zrange]
U0z16 = [U0(rk*angtob,1.6*angtob)*htoev for rk in rrange]
U1z16 = [U1(rk*angtob,1.6*angtob)*htoev for rk in rrange]

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
def VMatrix(R,Z,Gamma,BW,Nb,ab):
    Vmat = np.zeros((Nb+1,Nb+1))  #Set up
    Vmat[0,0] = U1(R*angtob,Z*angtob) - U0(R*angtob,Z*angtob) #molecule part of NAH potential
    Ek = MetalEK(BW*evtoh,Nb) 
    for k in range(1,Nb+1): 
        Vmat[k,k] = Ek[k-1]
        Vmat[0,k] = Vk(Z*angtob,Gamma*evtoh,BW*evtoh,Nb,ab*angtob) #couplinh b/w metal and molecule
        Vmat[k,0] = Vk(Z*angtob,Gamma*evtoh,BW*evtoh,Nb,ab*angtob) 
    return(Vmat)

#Function that calculates the ground state of the NAH hamiltonian.
def E0calc(R,Z,Gamma,BW,Nb,ab):
    evals, evecs = LA.eig(VMatrix(R,Z,Gamma,BW,Nb,ab))
    evals = np.sort(evals)
    E0 = U0(R*angtob,Z*angtob) + np.sum([evals[k]*fermi(evals[k]) for k in range(len(evals))])
    return(E0)

#Get E0 for Gamma = 1.5 eV, BW = 100eV, Nb = 200
E0r117BW100 = [E0calc(1.17,zk,1.5,100,200,10)*htoev for zk in zrange]
E0r16BW100 = [E0calc(1.6,zk,1.5,100,200,10)*htoev for zk in zrange]
E0z16BW100 = [E0calc(rk,1.6,1.5,100,200,10)*htoev for rk in rrange]
#Calculate reference energy at R = 1.6 A and Z = 6.02 A 
ErefBW100 = E0calc(1.6,r16e0data[-1,0],1.5,100,200,10)*htoev

#Set correct zero of energy for each calculation by setting Eref(calc) = ErefdFT
E0r117BW100 = [E0r117BW100[k] - ErefBW100 + r16e0data[-1,1] for k in range(len(E0r117BW100))]
E0r16BW100  = [E0r16BW100[k] - ErefBW100 + r16e0data[-1,1] for k in range(len(E0r16BW100))]
E0z16BW100  = [E0z16BW100[k] - ErefBW100 + r16e0data[-1,1] for k in range(len(E0z16BW100))]

#Get E0 for Gamma = 3.5 eV, BW = 7eV, Nb = 200, abar = 10 A^
E0r117 = [E0calc(1.17,zk,3.5,7.0,200,10)*htoev for zk in zrange]
E0r16  = [E0calc(1.6,zk,3.5,7.0,200,10)*htoev for zk in zrange]
E0z16  = [E0calc(rk,1.6,3.5,7.0,200,10)*htoev for rk in rrange]

#Calculate reference energy at R = 1.6 A and Z = 6.02 A 
Eref = E0calc(1.6,r16e0data[-1,0],3.5,7.0,200,10)*htoev

#Set correct zero of energy for each calculation.
E0r117 = [E0r117[k] - Eref + r16e0data[-1,1] for k in range(len(E0r117[:]))]
E0r16  = [E0r16[k] - Eref + r16e0data[-1,1] for k in range(len(E0r16[:]))]
E0z16  = [E0z16[k] - Eref + r16e0data[-1,1] for k in range(len(E0z16[:]))]

t2 = time.time()
print("Time taken in secs", t2-t1)
plt.plot(r117e0data[:,0],r117e0data[:,1],'.',markersize=10)
#plt.plot(r117h00data[:,0],r117h00data[:,1],'.')
#plt.plot(r117h11data[:,0],r117h11data[:,1],'.')
#plt.plot(zrange,U0r117,'-')
#plt.plot(zrange,U1r117,'-')
#for i in range(4):
plt.plot(zrange,E0r117,'-',linewidth='3.0') 
plt.plot(zrange,E0r117BW100,'-',linewidth='3.0')
plt.title("R = 1.17 Å")
plt.xlabel("Z (Å)")
plt.ylabel("Energy (eV)")
#plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=20)     # fontsize of the axes title
plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
plt.rc('ytick', labelsize=16)    # fontsize of the tick labels
plt.rc('legend', fontsize=16)    # legend fontsize
plt.rc('figure', titlesize=22)  # fontsize of the figure title
#plt.legend(["E0 DFT","Gamma=3.5 eV, BW = 7 eV","Gamma=1.5 eV, BW = 100 eV"],loc='upper right')
#plt.ylim([-0.2, 6.5])
plt.show()
plt.plot(r16e0data[:,0],r16e0data[:,1],'.',markersize=10)
#plt.plot(r16h00data[:,0],r16h00data[:,1],'.')
#plt.plot(r16h11data[:,0],r16h11data[:,1],'.')
#plt.plot(zrange,U0r16,'-')
#plt.plot(zrange,U1r16,'-')
#for i in range(4):
plt.plot(zrange,E0r16,'-',linewidth='3.0') 
plt.plot(zrange,E0r16BW100,'-',linewidth='3.0')
plt.legend(["cDFT Reference","Γ=3.5 eV, ΔE = 7 eV","Γ=1.5 eV, ΔE = 100 eV"],loc='upper right')
plt.title("R = 1.6 Å")
plt.xlabel("Z (Å)")
plt.ylabel("Energy (eV)")
plt.rc('axes', titlesize=20)     # fontsize of the axes title
plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
plt.rc('ytick', labelsize=16)    # fontsize of the tick labels
plt.rc('legend', fontsize=16)    # legend fontsize
plt.rc('figure', titlesize=22)  # fontsize of the figure title
#plt.ylim([2, 6.5])
plt.show()
plt.plot(z16e0data[:,0],z16e0data[:,1],'.',markersize=10)
#plt.plot(z16h00data[:,0],z16h00data[:,1],'.')
#plt.plot(z16h11data[:,0],z16h11data[:,1],'.')
#plt.plot(rrange,U0z16,'-')
#plt.plot(rrange,U1z16,'-')
#for i in range(4):
plt.plot(rrange,E0z16,'-',linewidth='3.0') 
plt.plot(rrange,E0z16BW100,'-',linewidth='3.0')
#plt.legend(["E0 DFT","Γ=3.5 eV, ΔE = 7 eV","Γ=1.5 eV, ΔE = 100 eV"],loc='upper right')
plt.xlabel("R (Å)")
plt.ylabel("Energy (eV)")
plt.title("Z = 1.6 Å")
plt.rc('axes', titlesize=20)     # fontsize of the axes title
plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
plt.rc('ytick', labelsize=16)    # fontsize of the tick labels
plt.rc('legend', fontsize=16)    # legend fontsize
plt.rc('figure', titlesize=22)  # fontsize of the figure title
#plt.ylim([-0.1, 6.5])
plt.show()



