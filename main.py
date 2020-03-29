## Importationsdes fonctions
import numpy as np
#from IPython.display import Image
import matplotlib.pyplot as plt
#import scipy.sparse.csc, scipy.sparse.linalg
#import time

# Importations des fonctions personalisées
#Pas besoin de les importer*****************************
#from Functions import Bois,PML,Source,Coeff_Frontiere,Coeff_PML,p,Source_Cylindrique,Construction_Map,Construction_A,\
    #Resolution, Plots_Results

from Functions import Source_Cylindrique,Construction_Map,Construction_A,Resolution,Plots_Results,surface_directe,Surface_equivalente,Construction_alpha_Map


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
N_PML = 20

# Emplacement du bois
centre_bois_x = 65
centre_bois_y = 60
# Longueur en x,y du bois (en points)
Nx_Bois = 20
Ny_Bois = 20

# Emplacement de la source
S_x = 50
S_y = 50

# Emplacement du détecteur
D_x = 23
D_y = 25

## Paramètres des milieux:

# Fréquence d'oscillation de la source
    # Mandat demande entre 100 Hz et 10 kHz

omega = 5e3


# Intensité de la source (arbitraire)
p_source = -1e12

# Eau
rho_eau = 998.3

alpha_eau = 1.18 *1e-7


B_eau = 2.15e9

# Bois
rho_bois = 640.72


alpha_bois = 3.2e-4 * 0.01

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
SourceCylindrique=True
Source_Map=np.ones([Nx,Ny])
# Main:

PML_mode=1  # Mode 2: PML avec le alpha map, Mode 1= PML classique
alpha_PML=5*alpha_eau


if __name__ == "__main__":


    if SourceCylindrique==True:
        Source_Map=Source_Cylindrique(Nx,Ny,S_x,S_y,dx,k2_eau,plot=False)


    Map,Display_Map= Construction_Map(Nx,Ny,Nx_Bois,Ny_Bois,centre_bois_x, centre_bois_y,forme,coeff,S_x,S_y,dx,N_PML,\
                                      plot=False,PML_mode=PML_mode, Bateau=True, Boisnez_bool=True)
    alpha_Map=Construction_alpha_Map(Nx,Ny,alpha_eau, alpha_PML,N_PML)
    #Temporaire
    SF_radius=10
    Q_map=np.ones([Nx,Ny])
    #for i in range(Nx):
    #    for j in range(Ny):
    #        if ((i-D_x)**2+(j-D_y)**2)<SF_radius**2:
    #            Q_map[i,j]=1
    Q_map[Display_Map==0]=0
    Q_map[Display_Map == 3] = 0




    A_sp,b_TFSF= Construction_A(Nx,Ny,dx,Neuf_points,k2_eau,k2_bois,gamma_eau,gamma_bois,rho_eau,p_source,SourceCylindrique,
                              Map,Source_Map ,Q_map,coeff,centre_bois_x,centre_bois_y,Nx_Bois,Ny_Bois, alpha_Map,omega,B_eau, PML_mode=PML_mode)


    MapSol_TFSF,P_detecteur=Resolution(A_sp, b_TFSF,Nx,Ny,D_x,D_y)

    ## Temporaire:

    fig,ax=plt.subplots(2,2,figsize=(16,8))
    ax[1][1].set_title("Scattered field seulement")
    SF_only=(MapSol_TFSF*Q_map)
    SF_only[SF_only==0]=np.nan
    ax[1][1].imshow(np.transpose((np.real(SF_only[N_PML:-N_PML,N_PML:-N_PML]))), alpha=1.0, cmap="jet")
    ax[1][0].set_title("Solution")
    ax[1][0].imshow(np.transpose((np.real(MapSol_TFSF))), alpha=1.0, cmap="jet")

    ax[0][0].set_title("Display Map")
    ax[0][0].imshow(np.transpose((Display_Map)), alpha=1.0, cmap="jet")

    ax[0][1].set_title("Région TF et Région SF en rouge")
    ax[0][1].imshow(np.transpose((Q_map)), alpha=1.0, cmap="jet")
    plt.show()


    # Plots_Results(MapSol, MapSolSB, MapSol_TFSF, Display_Map, Interpolation="none")
    

    #Surface = surface_directe(S_x, S_y, centre_bois_x, centre_bois_y, Nx_Bois, Ny_Bois, forme, coeff)
    
    #SER = Surface_equivalente(Nx_Bois,Ny_Bois,forme,coeff,P_incident,P_scattered,V_incident,Surface)



