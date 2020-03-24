# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 16:21:17 2020

Commentaire:
    FDTD
    Propagation d'onde dans le temps
"""

import numpy as np
import matplotlib.pyplot as plt

Lx = 10
Ly = 10

time = np.linspace(0,0.002,300)
dt = time[1] - time[0]

fact = 0.25
dx = 0.1*fact
Nx = int(np.round(Lx/dx)+1)
Ny = int(np.round(Ly/dx)+1)

x = np.linspace(0,Lx,Nx)
y = np.linspace(0,Ly,Ny)

alpha_b = 0 #3.21e-4
B_b = 10e9
rho_b = 640.72
alpha_e = 0 #1.18
B_e = 2.15e9
rho_e = 998.30

# Initialisation time
P = np.zeros((int(Ny),int(Nx),len(time)));

yy, xx = np.meshgrid(y, x, sparse=True)
P[:,:,1] = np.exp(-((xx-5)/dx)**2 - ((yy-5)/dx)**2)

P[:,:,2] = 0.99*P[:,:,1]

# Initialisation space
B=np.zeros((int(Ny),int(Nx))); alpha=np.zeros((int(Ny),int(Nx))); rho=np.zeros((int(Ny),int(Nx)));
for i in range(int(Ny)):
    yy=(i)*dx;
    for j in range(int(Nx)):
        xx=(j)*dx;
        
        # Délimitation des espaces
        if (xx>=6)&(xx<=7)&(yy>=3)&(yy<=4):
            # Bois
            rho[i,j] = rho_b;
            alpha[i,j] = alpha_b;
            B[i,j] = B_b;
        else:
            # eau
            rho[i,j] = rho_e;
            alpha[i,j] = alpha_e;
            B[i,j] = B_e;

P_now = np.zeros((int(Ny),int(Nx)));

for n in range(1,len(time)-1):
    
    # Core
#    for i in range(1,int(Ny)-1):
#        for j in range(1,int(Nx)-1):
#            P_now[i,j] = 1/(alpha[i,j]*B[i,j]*dt - 1) \
#                *( B[i,j-1]*dt**2./rho[i,j-1]/dx**2*P[i,j-1,n] \
#                + B[i,j+1]*dt**2/rho[i,j+1]/dx**2*P[i,j+1,n] \
#                + B[i-1,j]*dt**2/rho[i-1,j]/dx**2*P[i-1,j,n] \
#                + B[i+1,j]*dt**2/rho[i+1,j]/dx**2*P[i+1,j,n] \
#                - P[i,j,n]*(4.*B[i,j]*dt**2/rho[i,j] + 2) \
#                + P[i,j,n-1]*(alpha[i,j]*B[i,j]*dt + 1) )
                
    P_now[1:int(Ny)-1,1:int(Nx)-1] = 1/(alpha[1:int(Ny)-1,1:int(Nx)-1]*B[1:int(Ny)-1,1:int(Nx)-1]*dt - 1) \
        *( B[1:int(Ny)-1,0:int(Nx)-2]*dt**2./rho[1:int(Ny)-1,0:int(Nx)-2]/dx**2*P[1:int(Ny)-1,0:int(Nx)-2,n] \
        + B[1:int(Ny)-1,2:int(Nx)]*dt**2/rho[1:int(Ny)-1,2:int(Nx)]/dx**2*P[1:int(Ny)-1,2:int(Nx),n] \
        + B[0:int(Ny)-2,1:int(Nx)-1]*dt**2/rho[0:int(Ny)-2,1:int(Nx)-1]/dx**2*P[0:int(Ny)-2,1:int(Nx)-1,n] \
        + B[2:int(Ny),1:int(Nx)-1]*dt**2/rho[2:int(Ny),1:int(Nx)-1]/dx**2*P[2:int(Ny),1:int(Nx)-1,n] \
        - P[1:int(Ny)-1,1:int(Nx)-1,n]*(4.*B[1:int(Ny)-1,1:int(Nx)-1]*dt**2/rho[1:int(Ny)-1,1:int(Nx)-1] + 2) \
        + P[1:int(Ny)-1,1:int(Nx)-1,n-1]*(alpha[1:int(Ny)-1,1:int(Nx)-1]*B[1:int(Ny)-1,1:int(Nx)-1]*dt + 1) )    
    
    # Border
#    for i in range(1,int(Ny)-1):
#        # x=1
#        j = 0
#        P_now[i,j] = 1/(alpha[i,j]*B[i,j]*dt - 1) \
#            *( - 5.*B[i,j+1]*dt**2/rho[i,j+1]/dx**2*P[i,j+1,n] \
#            + 4.*B[i,j+2]*dt**2/rho[i,j+2]/dx**2*P[i,j+2,n] \
#            - B[i,j+3]*dt**2/rho[i,j+3]/dx**2*P[i,j+3,n] \
#            + B[i-1,j]*dt**2/rho[i-1,j]/dx**2*P[i-1,j,n] \
#            + B[i+1,j]*dt**2/rho[i+1,j]/dx**2*P[i+1,j,n] \
#            - P[i,j,n]*(2) \
#            + P[i,j,n-1]*(alpha[i,j]*B[i,j]*dt + 1) )
#        
#        # x=end
#        j = int(Nx)-1;
#        P_now[i,j] = 1/(alpha[i,j]*B[i,j]*dt - 1) \
#            *( - 5.*B[i,j-1]*dt**2/rho[i,j-1]/dx**2*P[i,j-1,n] \
#            + 4.*B[i,j-2]*dt**2/rho[i,j-2]/dx**2*P[i,j-2,n] \
#            - B[i,j-3]*dt**2/rho[i,j-3]/dx**2*P[i,j-3,n] \
#            + B[i-1,j]*dt**2/rho[i-1,j]/dx**2*P[i-1,j,n] \
#            + B[i+1,j]*dt**2/rho[i+1,j]/dx**2*P[i+1,j,n] \
#            - P[i,j,n]*(2) \
#            + P[i,j,n-1]*(alpha[i,j]*B[i,j]*dt + 1) )
    
    # x=1
    j=0  
    P_now[1:int(Ny)-1,j] = 1/(alpha[1:int(Ny)-1,j]*B[1:int(Ny)-1,j]*dt - 1) \
        *( - 5.*B[1:int(Ny)-1,j+1]*dt**2/rho[1:int(Ny)-1,j+1]/dx**2*P[1:int(Ny)-1,j+1,n] \
        + 4.*B[1:int(Ny)-1,j+2]*dt**2/rho[1:int(Ny)-1,j+2]/dx**2*P[1:int(Ny)-1,j+2,n] \
        - B[1:int(Ny)-1,j+3]*dt**2/rho[1:int(Ny)-1,j+3]/dx**2*P[1:int(Ny)-1,j+3,n] \
        + B[0:int(Ny)-2,j]*dt**2/rho[0:int(Ny)-2,j]/dx**2*P[0:int(Ny)-2,j,n] \
        + B[2:int(Ny),j]*dt**2/rho[2:int(Ny),j]/dx**2*P[2:int(Ny),j,n] \
        - P[1:int(Ny)-1,j,n]*(2) \
        + P[1:int(Ny)-1,j,n-1]*(alpha[1:int(Ny)-1,j]*B[1:int(Ny)-1,j]*dt + 1) )
    
    # x=end
    j = int(Nx)-1;
    P_now[1:int(Ny)-1,j] = 1/(alpha[1:int(Ny)-1,j]*B[1:int(Ny)-1,j]*dt - 1) \
        *( - 5.*B[1:int(Ny)-1,j-1]*dt**2/rho[1:int(Ny)-1,j-1]/dx**2*P[1:int(Ny)-1,j-1,n] \
        + 4.*B[1:int(Ny)-1,j-2]*dt**2/rho[1:int(Ny)-1,j-2]/dx**2*P[1:int(Ny)-1,j-2,n] \
        - B[1:int(Ny)-1,j-3]*dt**2/rho[1:int(Ny)-1,j-3]/dx**2*P[1:int(Ny)-1,j-3,n] \
        + B[0:int(Ny)-2,j]*dt**2/rho[0:int(Ny)-2,j]/dx**2*P[0:int(Ny)-2,j,n] \
        + B[2:int(Ny),j]*dt**2/rho[2:int(Ny),j]/dx**2*P[2:int(Ny),j,n] \
        - P[1:int(Ny)-1,j,n]*(2) \
        + P[1:int(Ny)-1,j,n-1]*(alpha[1:int(Ny)-1,j]*B[1:int(Ny)-1,j]*dt + 1) )    
    
#    for j in range(1,int(Nx)-1):
#        # y=1
#        i = 0;
#        P_now[i,j] = 1/(alpha[i,j]*B[i,j]*dt - 1) \
#            *( - 5.*B[i+1,j]*dt**2/rho[i+1,j]/dx**2*P[i+1,j,n] \
#            + 4.*B[i+2,j]*dt**2/rho[i+2,j]/dx**2*P[i+2,j,n] \
#            - B[i+3,j]*dt**2/rho[i+3,j]/dx**2*P[i+3,j,n] \
#            + B[i,j-1]*dt**2/rho[i,j-1]/dx**2*P[i,j-1,n] \
#            + B[i,j+1]*dt**2/rho[i,j+1]/dx**2*P[i,j+1,n] \
#            - P[i,j,n]*(2) \
#            + P[i,j,n-1]*(alpha[i,j]*B[i,j]*dt + 1) )
#        
#        # y=end
#        i = int(Ny)-1;
#        P_now[i,j] = 1/(alpha[i,j]*B[i,j]*dt - 1) \
#            *( - 5.*B[i-1,j]*dt**2/rho[i-1,j]/dx**2*P[i-1,j,n] \
#            + 4.*B[i-2,j]*dt**2/rho[i-2,j]/dx**2*P[i-2,j,n] \
#            - B[i-3,j]*dt**2/rho[i-3,j]/dx**2*P[i-3,j,n] \
#            + B[i,j-1]*dt**2/rho[i,j-1]/dx**2*P[i,j-1,n] \
#            + B[i,j+1]*dt**2/rho[i,j+1]/dx**2*P[i,j+1,n] \
#            - P[i,j,n]*(2) \
#            + P[i,j,n-1]*(alpha[i,j]*B[i,j]*dt + 1) )
            
    # y=1
    i = 0;
    P_now[i,1:int(Nx)-1] = 1/(alpha[i,1:int(Nx)-1]*B[i,1:int(Nx)-1]*dt - 1) \
        *( - 5.*B[i+1,1:int(Nx)-1]*dt**2/rho[i+1,1:int(Nx)-1]/dx**2*P[i+1,1:int(Nx)-1,n] \
        + 4.*B[i+2,1:int(Nx)-1]*dt**2/rho[i+2,1:int(Nx)-1]/dx**2*P[i+2,1:int(Nx)-1,n] \
        - B[i+3,1:int(Nx)-1]*dt**2/rho[i+3,1:int(Nx)-1]/dx**2*P[i+3,1:int(Nx)-1,n] \
        + B[i,0:int(Nx)-2]*dt**2/rho[i,0:int(Nx)-2]/dx**2*P[i,0:int(Nx)-2,n] \
        + B[i,2:int(Nx)]*dt**2/rho[i,2:int(Nx)]/dx**2*P[i,2:int(Nx),n] \
        - P[i,1:int(Nx)-1,n]*(2) \
        + P[i,1:int(Nx)-1,n-1]*(alpha[i,1:int(Nx)-1]*B[i,1:int(Nx)-1]*dt + 1) )
    
    # y=end
    i = int(Ny)-1;
    P_now[i,1:int(Nx)-1] = 1/(alpha[i,1:int(Nx)-1]*B[i,1:int(Nx)-1]*dt - 1) \
        *( - 5.*B[i-1,1:int(Nx)-1]*dt**2/rho[i-1,1:int(Nx)-1]/dx**2*P[i-1,1:int(Nx)-1,n] \
        + 4.*B[i-2,1:int(Nx)-1]*dt**2/rho[i-2,1:int(Nx)-1]/dx**2*P[i-2,1:int(Nx)-1,n] \
        - B[i-3,1:int(Nx)-1]*dt**2/rho[i-3,1:int(Nx)-1]/dx**2*P[i-3,1:int(Nx)-1,n] \
        + B[i,0:int(Nx)-2]*dt**2/rho[i,0:int(Nx)-2]/dx**2*P[i,0:int(Nx)-2,n] \
        + B[i,2:int(Nx)]*dt**2/rho[i,2:int(Nx)]/dx**2*P[i,2:int(Nx),n] \
        - P[i,1:int(Nx)-1,n]*(2) \
        + P[i,1:int(Nx)-1,n-1]*(alpha[i,1:int(Nx)-1]*B[i,1:int(Nx)-1]*dt + 1) )
    
    # Coins
    # x=1,y=1
    i=0;j=0;
    P_now[i,j] = 1/(alpha[i,j]*B[i,j]*dt - 1) \
        *( - 5.*B[i,j+1]*dt**2/rho[i,j+1]/dx**2*P[i,j+1,n] \
        + 4.*B[i,j+2]*dt**2/rho[i,j+2]/dx**2*P[i,j+2,n] \
        - B[i,j+3]*dt**2/rho[i,j+3]/dx**2*P[i,j+3,n] \
        - 5.*B[i+1,j]*dt**2/rho[i+1,j]/dx**2*P[i+1,j,n] \
        + 4.*B[i+2,j]*dt**2/rho[i+2,j]/dx**2*P[i+2,j,n] \
        - B[i+3,j]*dt**2/rho[i+3,j]/dx**2*P[i+3,j,n] \
        - P[i,j,n]*(-4.*B[i,j]*dt**2/rho[i,j] + 2) \
        + P[i,j,n-1]*(alpha[i,j]*B[i,j]*dt + 1) )
    
    # x=1,y=end
    i=int(Ny)-1;j=0;
    P_now[i,j] = 1/(alpha[i,j]*B[i,j]*dt - 1) \
        *( - 5.*B[i,j+1]*dt**2/rho[i,j+1]/dx**2*P[i,j+1,n] \
        + 4.*B[i,j+2]*dt**2/rho[i,j+2]/dx**2*P[i,j+2,n] \
        - B[i,j+3]*dt**2/rho[i,j+3]/dx**2*P[i,j+3,n] \
        - 5.*B[i-1,j]*dt**2/rho[i-1,j]/dx**2*P[i-1,j,n] \
        + 4.*B[i-2,j]*dt**2/rho[i-2,j]/dx**2*P[i-2,j,n] \
        - B[i-3,j]*dt**2/rho[i-3,j]/dx**2*P[i-3,j,n] \
        - P[i,j,n]*(-4.*B[i,j]*dt**2/rho[i,j] + 2) \
        + P[i,j,n-1]*(alpha[i,j]*B[i,j]*dt + 1) )
    
    # x=end,y=1
    i=0;j=int(Nx)-1;
    P_now[i,j] = 1/(alpha[i,j]*B[i,j]*dt - 1) \
        *( - 5.*B[i,j-1]*dt**2/rho[i,j-1]/dx**2*P[i,j-1,n] \
        + 4.*B[i,j-2]*dt**2/rho[i,j-2]/dx**2*P[i,j-2,n] \
        - B[i,j-3]*dt**2/rho[i,j-3]/dx**2*P[i,j-3,n] \
        - 5.*B[i+1,j]*dt**2/rho[i+1,j]/dx**2*P[i+1,j,n] \
        + 4.*B[i+2,j]*dt**2/rho[i+2,j]/dx**2*P[i+2,j,n] \
        - B[i+3,j]*dt**2/rho[i+3,j]/dx**2*P[i+3,j,n] \
        - P[i,j,n]*(-4.*B[i,j]*dt**2/rho[i,j] + 2) \
        + P[i,j,n-1]*(alpha[i,j]*B[i,j]*dt + 1) )
    
    # x=end,y=end
    i=int(Ny)-1;j=int(Nx)-1;
    P_now[i,j] = 1/(alpha[i,j]*B[i,j]*dt - 1) \
        *( - 5.*B[i,j-1]*dt**2/rho[i,j-1]/dx**2*P[i,j-1,n] \
        + 4.*B[i,j-2]*dt**2/rho[i,j-2]/dx**2*P[i,j-2,n] \
        - B[i,j-3]*dt**2/rho[i,j-3]/dx**2*P[i,j-3,n] \
        - 5.*B[i-1,j]*dt**2/rho[i-1,j]/dx**2*P[i-1,j,n] \
        + 4.*B[i-2,j]*dt**2/rho[i-2,j]/dx**2*P[i-2,j,n] \
        - B[i-3,j]*dt**2/rho[i-3,j]/dx**2*P[i-3,j,n] \
        - P[i,j,n]*(-4.*B[i,j]*dt**2/rho[i,j] + 2) \
        + P[i,j,n-1]*(alpha[i,j]*B[i,j]*dt + 1) )
    
    # Next step
    P[:,:,n+1] = P_now[:,:]

    #print(n);
    
    

plt.imshow(np.log(abs(P[:,:,-1])+1e-300), cmap="jet", alpha=1, interpolation='none',origin='lower')
    









