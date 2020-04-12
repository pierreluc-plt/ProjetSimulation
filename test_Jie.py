## Importationsdes fonctions
import numpy as np
import math
#from IPython.display import Image
import matplotlib.pyplot as plt
#import scipy.sparse.csc, scipy.sparse.linalg
#import time

# Importations des fonctions personalisées
#Pas besoin de les importer*****************************
#from Functions import Bois,PML,Source,Coeff_Frontiere,Coeff_PML,p,Source_Cylindrique,Construction_Map,Construction_A,\
    #Resolution, Plots_Results

from Functions_Jie import Source_Cylindrique,Source_Lineaire,Source_Ponctuelle,Construction_Map,Construction_A,Resolution,\
    Plots_Results,surface_directe,Surface_equivalente,Construction_alpha_Map,Select_Source,SER2


# Fonction dépendant de Nx

# Nombre de points en x

Nx = 201


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
N_PML = 10


centre_bois_x = int(Nx / 2)  # 65
centre_bois_y = int(1.3 * Nx / 2)  # 600#60
# Longueur en x,y du bois (en points)
Nx_Bois = int(5 / dx)
Ny_Bois = int(10 / dx)

# Emplacement de la source
S_x = int(Nx / 5*2)
S_y = 50

# Emplacement du détecteur
D_x = 99
D_y = 99

## Paramètres des milieux:

# Fréquence d'oscillation de la source
    # Mandat demande entre 100 Hz et 10 kHz

omega = 1000


# Intensité de la source (arbitraire)
p_source = -1e12

# Eau
rho_eau = 998.3

alpha_eau = 0 # 1.18 *1e-7


B_eau = 2.15e9

# Bois
rho_bois = 640.72


alpha_bois = 3.2e-4 #* 0.0001

B_bois = 10e9

# Vitesse du son
v_eau = np.sqrt(B_eau/rho_eau)
v_bois = np.sqrt(B_bois/rho_bois)

# # Paramètres calculés
# k2_eau = rho_eau * (omega ** 2 / B_eau + 2j * omega * alpha_eau)
# k0_eau = rho_eau * (omega ** 2 / B_eau)
# k2_bois = rho_bois * (omega ** 2 / B_bois + 2j * omega * alpha_bois)
# n_eau = 1.33
#
# gamma_eau = rho_eau * (alpha_eau * B_eau + 1j * omega)
# gamma_bois = rho_bois * (alpha_bois * B_bois + 1j * omega)

## Paramètres modifiables pour l'exécution du code
forme = 'triangle'
# coeff doit être en radian et supérieur à 0 et inférieur à pi
coeff = np.pi/4
#forme = 'cercle'
# coeff dois être égal ou supérieur à Nx_Bois/2
#coeff = Nx_Bois/2

# Pour faire le code à 9 points ou pas
Neuf_points = True

# Décider type de source
Source = '' # Lineaire ou Cylindrique[Default]
theta = -5 # Angle en degrées

if Source=='Lineaire':
    SourceLineaire=True
    SourceCylindrique=False
    SourcePonctuelle=False
elif Source=='Ponctuelle':
    SourceLineaire=False
    SourceCylindrique=False
    SourcePonctuelle=True
else:
    SourceCylindrique=True
    SourceLineaire=False
    SourcePonctuelle=False
    
Source_Map=np.ones([Ny,Nx])
# Main:

PML_mode=1  # Mode 2: PML avec le alpha map, Mode 1= PML classique
alpha_PML=5*alpha_eau


SER_Array=[]
SER_Array_v2=[]
SF_only_Array=[]
if __name__ == "__main__":



    Map,Display_Map= Construction_Map(Nx,Ny,Nx_Bois,Ny_Bois,centre_bois_x, centre_bois_y,forme,coeff,S_x,S_y,dx,N_PML,\
                                      plot=False,PML_mode=PML_mode, Bateau=True, Boisnez_bool=True)
    alpha_Map=Construction_alpha_Map(Nx,Ny,alpha_eau, alpha_PML,N_PML)
    #Temporaire
    Q_map=np.ones([Ny,Nx])

    Q_map[Display_Map==0]=0
    Q_map[Display_Map == 3] = 0

    Surface = surface_directe(S_x, S_y, centre_bois_x, centre_bois_y, Nx_Bois, Ny_Bois, forme, coeff)


    omega_array=2*np.pi*np.array([omega])

    for omega in omega_array:
        # Paramètres calculés
        k2_eau = rho_eau * (omega ** 2 / B_eau + 2j * omega * alpha_eau)
        k0_eau = rho_eau * (omega ** 2 / B_eau)
        k2_bois = rho_bois * (omega ** 2 / B_bois + 2j * omega * alpha_bois)
        n_eau = 1.33
        gamma_eau = rho_eau * (alpha_eau * B_eau + 1j * omega)
        gamma_bois = rho_bois * (alpha_bois * B_bois + 1j * omega)

        Plot_Source = True
        Source_Map=Select_Source(SourceLineaire, SourceCylindrique, SourcePonctuelle, Nx, Ny, S_x, S_y, dx, k2_eau, k0_eau, n_eau,\
                      theta, Plot_Source)

        A_sp,b_TFSF= Construction_A(Nx,Ny,dx,Neuf_points,k2_eau,k2_bois,gamma_eau,gamma_bois,rho_eau,v_eau,p_source,SourceCylindrique,SourceLineaire,SourcePonctuelle,
                              Map,N_PML,Source_Map,Q_map,coeff,centre_bois_x,centre_bois_y,Nx_Bois,Ny_Bois, alpha_Map,omega,B_eau, PML_mode=PML_mode,TF_SF=True)


        MapSol_TFSF=Resolution(A_sp, b_TFSF,Nx,Ny,D_x,D_y)
        SF_only=(MapSol_TFSF)
        SF_only[SF_only==0]=np.nan
        SF_only_Array.append(SF_only)

    ## Temporaire:

    # fig,ax=plt.subplots(2,2,figsize=(16,8))
    # ax[1][1].set_title("Scattered field seulement")
    #
    # ax[1][1].imshow(np.transpose((np.real(SF_only[N_PML:-N_PML,N_PML:-N_PML]))), alpha=1.0, cmap="jet")
    # ax[1][0].set_title("Solution")
    # ax[1][0].imshow(np.transpose((np.real(MapSol_TFSF))), alpha=1.0, cmap="jet")
    #
    # ax[0][0].set_title("Source Map")
    # ax[0][0].imshow(np.transpose((abs(Source_Map))), alpha=1.0, cmap="jet")
    #
    # ax[0][1].set_title("Région TF et Région SF en rouge")
    # ax[0][1].imshow(np.transpose((Q_map)), alpha=1.0, cmap="jet")
    # plt.show()


    # Plots_Results(MapSol, MapSolSB, MapSol_TFSF, Display_Map, Interpolation="none")
        plt.figure()
        plt.imshow(np.transpose((np.real(SF_only[N_PML:-N_PML,N_PML:-N_PML]))), alpha=1.0, cmap="jet",interpolation="gaussian")
        plt.show()
        SERv2=SER2(S_x, S_y, MapSol_TFSF, Source_Map,p_source, centre_bois_x, centre_bois_y, 5, dx)
        SER_Array_v2.append(SERv2)

        SER = Surface_equivalente(S_x, S_y, p_source, Nx, Lx, Nx_Bois, Ny_Bois, forme, coeff, Source_Map, SF_only,
                                  Surface)
        SER_Array.append(SER)

        # plt.title("SF Only, omega={}".format(omega))
        # plt.show()


#    plt.figure()
#    plt.plot(omega_array,np.array(SER_Array)/SER_Array[0],label="V1")
#    plt.plot(omega_array, np.array(SER_Array_v2) / SER_Array_v2[0],label="V2")
#    plt.xlabel("omega")
#    plt.ylabel("SER")
#    plt.legend()
#    plt.show()
#    #np.save("Resultats/Carre1501_Source_f250Hz.npy", Source_Map)
#    np.save("Resultats/Carre1501_f5000Hz_SER" + str(SERv2) + ".npy", MapSol_TFSF)

### Test Gif



    # # fig, ax = plt.subplots(figsize=(12, 12))
    # Min = 0
    # Max = 0
    # for i in range(len(SF_only_Array)):
    #     S = np.real(SF_only_Array[i])
    #
    #     if Min > np.nanmin(S, ):
    #         Min = np.nanmin(S)
    #     if Max < np.nanmax(S):
    #         Max = np.nanmax(S)

    # def update(i):
    #     SF_only = SF_only_Array[i]
    #     ax.imshow(np.transpose((np.real(SF_only[N_PML:-N_PML,N_PML:-N_PML]))),interpolation="gaussian",vmin=Min,vmax=Max, cmap="jet")
    #     #ax.set_title("Omega: {}".format(omega_array[i]), fontsize=20)
    #     ax.set_axis_off()
    # anim = FuncAnimation(fig, update, frames=np.arange(0, len(SF_only_Array)), interval=100)
    # anim.save('colour_rotation.gif', dpi=40, writer='imagemagick')


