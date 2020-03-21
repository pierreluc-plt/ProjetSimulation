## Importationsdes fonctions
import numpy as np
from IPython.display import Image
import matplotlib.pyplot as plt
import scipy.sparse.csc, scipy.sparse.linalg
import time

# Importations des fonctions personalisées

from Functions import Bois,PML,Source,Coeff_Frontiere,Coeff_PML,p,Source_Cylindrique,Construction_Map,Construction_A,\
    Resolution, Plots_Results



# Fonction dépendant de Nx

# Nombre de points en x

Nx = 100


# Paramètres de la simulations
# Longueur en x (m)

Lx = 40
# Espace entre chaque noeud
dx = (Lx) / (Nx - 1)

# Nombre de points en y
Ny = Nx
# Longueur en x (m)
Ly = Lx
# Espace entre chaque noeud
dy = dx
h=dx
# Épaisseur (en points de la couche de PML)
N_PML = 5

# Emplacement du bois
centre_bois_x = 60
centre_bois_y = 60
# Longueur en x du bois (en points)
Nx_Bois = 10
Ny_Bois = 40

# Emplacement de la source
S_x = 30
S_y = 30

## Paramètres des milieux:

# Fréquence d'oscillation de la source
omega = 1e2
# Intensité de la source (arbitraire)
p_source = 1e2

# Eau
rho_eau = 998.3
alpha_eau = 1.18 * 0.0001
B_eau = 2.15e9

# Bois
rho_bois = 640.72
alpha_bois = 3.2e-4 * 0.1
B_bois = 10e9

# Paramètres calculés
k2_eau = rho_eau * (omega ** 2 / B_eau + 2j * omega * alpha_eau)
k2_bois = rho_bois * (omega ** 2 / B_bois + 2j * omega * alpha_bois)

gamma_eau = rho_eau * (alpha_eau * B_eau + 1j * omega)
gamma_bois = rho_bois * (alpha_bois * B_bois + 1j * omega)

## Paramètres modifiables pour l'exécution du code
forme = 'triangle'
coeff = 2*np.pi/4
#forme = 'cercle'
#coeff = 11

# Pour faire le code à 9 points ou pas
Neuf_points = True

# Décider si on utilise un source cylindrique
SourceCylindrique=False
Source_Map=np.ones([Nx,Ny])
# Main:



if __name__ == "__main__":


    if SourceCylindrique==True:
        Source_Map=Source_Cylindrique(Nx,Ny,S_x,S_y,dx,k2_eau,plot=False)

    Map,MapSB,Display_Map= Construction_Map(Nx,Ny,Nx_Bois,Ny_Bois,centre_bois_x, centre_bois_y,forme,coeff,S_x,S_y,dx,N_PML,plot=True)

    A,A_SB, b= Construction_A(Nx,Ny,dx,Neuf_points,k2_eau,k2_bois,gamma_eau,gamma_bois,rho_eau,p_source,SourceCylindrique,
                              Map,MapSB,Source_Map,coeff,centre_bois_x,centre_bois_y)




    MapSol,MapSolSB=Resolution(A,A_SB,b,Nx,Ny)

    Plots_Results(MapSol, MapSolSB, Display_Map, Interpolation="none")


