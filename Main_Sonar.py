import numpy as np
import math
import matplotlib.pyplot as plt
import time
import scipy.sparse.csc, scipy.sparse.linalg


## Importationsdes fonctions
import numpy as np
import math
# from IPython.display import Image
import matplotlib.pyplot as plt


from Functions_Sonar import Source_Cylindrique,Source_Lineaire,Source_Ponctuelle,Construction_Map,Construction_A,Resolution,\
    Plots_Results,surface_directe,Surface_equivalente,Construction_alpha_Map,Select_Source,SER3
# Paramètres de la simulations
Nx=401
forme="rectangle"
coeff=1.0
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
h = dx
# Épaisseur (en points de la couche de PML)
N_PML = 10

centre_bois_x = int(Nx / 2)  # 65
centre_bois_y = int(1.3 * Nx / 2)  # 600#60
# Longueur en x,y du bois (en points)
Nx_Bois = int(5 / dx)
Ny_Bois = int(10 / dx)

# Emplacement de la source
S_x = int(Nx / 2)
S_y = 75

# Emplacement du détecteur
D_x = 99
D_y = 99

## Paramètres des milieux:

# Fréquence d'oscillation de la source
# Mandat demande entre 100 Hz et 10 kHz

omega = 10000

# Intensité de la source (arbitraire)
p_source = 1e12

# Eau
rho_eau = 998.3

alpha_eau = 1.18 * 1e-7 # *1e-7 pour les tests

B_eau = 2.15e9

# Bois
rho_bois = 640.72

alpha_bois = 3.2e-4 * 0.01

B_bois = 10e9

# Vitesse du son
v_eau = np.sqrt(B_eau / rho_eau)
v_bois = np.sqrt(B_bois / rho_bois)


if forme == "cercle":
    coeff = coeff * Nx_Bois / 2

if forme == "rectangle":
    Pointe = False
    Nx_Bois = int(coeff * Nx_Bois)
else:
    Pointe = True

# Pour faire le code à 9 points ou pas
Neuf_points = True

# Décider type de source
Source_type = ''  # Lineaire ou Cylindrique[Default]
theta = -30  # Angle en degrées

if Source_type == 'Lineaire':
    SourceLineaire = True
    SourceCylindrique = False
    SourcePonctuelle = False
elif Source_type == 'Ponctuelle':
    SourceLineaire = False
    SourceCylindrique = False
    SourcePonctuelle = True
else:
    SourceCylindrique = True
    SourceLineaire = False
    SourcePonctuelle = False

# Main:

PML_mode = 1  # Mode 2: PML avec le alpha map, Mode 1= PML classique
alpha_PML = 5 * alpha_eau

SER_Array= []


if __name__ == "__main__":

    Map, Display_Map = Construction_Map(Nx, Ny, Nx_Bois, Ny_Bois, centre_bois_x, centre_bois_y, forme, coeff, S_x,
                                        S_y, dx, N_PML, \
                                        plot=False, PML_mode=PML_mode, Bateau=True, Boisnez_bool=Pointe)
    alpha_Map = Construction_alpha_Map(Nx, Ny, alpha_eau, alpha_PML, N_PML)
    # Temporaire
    Q_map = np.ones([Ny, Nx])

    Q_map[Display_Map == 0] = 0
    Q_map[Display_Map == 3] = 0



    omega_array = 2 * np.pi * np.linspace(100, 10000, 1)


    for omega in omega_array:
        # Paramètres calculés
        print("Fréquence={:.2f}".format(omega))
        k2_eau = rho_eau * (omega ** 2 / B_eau + 2j * omega * alpha_eau)
        k0_eau = rho_eau * (omega ** 2 / B_eau)
        k2_bois = rho_bois * (omega ** 2 / B_bois + 2j * omega * alpha_bois)
        n_eau = 1.33
        gamma_eau = rho_eau * (alpha_eau * B_eau + 1j * omega)
        gamma_bois = rho_bois * (alpha_bois * B_bois + 1j * omega)

        Plot_Source = False
        Source_Map = Select_Source(SourceLineaire, SourceCylindrique, SourcePonctuelle, Nx, Ny, S_x, S_y, dx,
                                   k2_eau, k0_eau, n_eau, \
                                   theta, Plot_Source)

        A_sp, b_TFSF = Construction_A(Nx, Ny, dx, Neuf_points, k2_eau, k2_bois, gamma_eau, gamma_bois, rho_eau,
                                      v_eau, p_source, SourceCylindrique, SourceLineaire, SourcePonctuelle,
                                      Map, N_PML, Source_Map, Q_map, coeff, centre_bois_x, centre_bois_y, Nx_Bois,
                                      Ny_Bois, alpha_Map, omega, B_eau, PML_mode=PML_mode, TF_SF=True)

        MapSol_TFSF = Resolution(A_sp, b_TFSF, Nx, Ny, D_x, D_y)
        SF_only = (MapSol_TFSF)
        SF_only[SF_only == 0] = np.nan


        fig, ax = plt.subplots(1, 1, figsize=(6, 6))
        ax.set_title("Scattered field seulement")

        ax.imshow(np.transpose((np.real(SF_only[N_PML:-N_PML, N_PML:-N_PML]))), alpha=1.0, cmap="jet")

        plt.show()


        SER_val = SER3(S_x, S_y, MapSol_TFSF, Source_Map, p_source, centre_bois_x, centre_bois_y, int(1 / dx), dx,
                     forme, coeff, Nx_Bois, Ny_Bois)
        SER_Array.append(SER_val)
