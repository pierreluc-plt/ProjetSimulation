## Importationsdes fonctions
import numpy as np
#from IPython.display import Image
#import matplotlib.pyplot as plt
#import scipy.sparse.csc, scipy.sparse.linalg
#import time

# Importations des fonctions personalisées
#Pas besoin de les importer*****************************
#from Functions import Bois,PML,Source,Coeff_Frontiere,Coeff_PML,p,Source_Cylindrique,Construction_Map,Construction_A,\
    #Resolution, Plots_Results

from Functions import Source_Cylindrique,Construction_Map,Construction_A,Resolution,Plots_Results,surface_directe


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
# Longueur en x,y du bois (en points)
Nx_Bois = 26
Ny_Bois = 40

# Emplacement de la source
S_x = 30
S_y = 30

# Emplacement du détecteur
D_x = 30
D_y = 60

## Paramètres des milieux:

# Fréquence d'oscillation de la source
    # Mandat demande entre 100 Hz et 10 kHz
omega = 1e2 

# Intensité de la source (arbitraire)
p_source = -1e10

# Eau
rho_eau = 998.3
alpha_eau = 1.18 * 0.0
B_eau = 2.15e9

# Bois
rho_bois = 640.72
alpha_bois = 3.2e-4 * 0.1
B_bois = 10e9

# Vitesse du son
v_eau = np.sqrt(B_eau/rho_eau)
v_bois = np.sqrt(B_bois/rho_bois)

# Paramètres calculés
k2_eau = rho_eau * (omega ** 2 / B_eau + 2j * omega * alpha_eau)
k2_bois = rho_bois * (omega ** 2 / B_bois + 2j * omega * alpha_bois)

gamma_eau = rho_eau * (alpha_eau * B_eau + 1j * omega)
gamma_bois = rho_bois * (alpha_bois * B_bois + 1j * omega)

## Paramètres modifiables pour l'exécution du code
forme = 'triangle'
# coeff doit être en radian et supérieur à 0 et inférieur à pi
coeff = np.pi/4
#forme = 'cercle'
# coeff dois être égal ou supérieur à Nx_Bois/2
#coeff = Nx_Bois/2

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
                              Map,MapSB,Source_Map,coeff,centre_bois_x,centre_bois_y,Nx_Bois,Ny_Bois)




    MapSol,MapSolSB,P_detecteur=Resolution(A,A_SB,b,Nx,Ny,D_x,D_y)

    Plots_Results(MapSol, MapSolSB, Display_Map, Interpolation="none")
    
    Surface = surface_directe(S_x, S_y, centre_bois_x, centre_bois_y, Nx_Bois, Ny_Bois, forme, coeff)


