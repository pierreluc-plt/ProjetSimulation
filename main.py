## Importationsdes fonctions
import numpy as np
from IPython.display import Image
import matplotlib.pyplot as plt
import scipy.sparse.csc, scipy.sparse.linalg
import time

# Importations des fonctions personalisées

from Functions import Bois
from Functions import PML
from Functions import Source
from Functions import Coeff_Frontiere
from Functions import Coeff_PML

# Fonction dépendant de Nx

# Nombre de points en x
Nx = 100


# Position du point
def p(i, j):
    """Le noeud (i,j) correspond à la ligne L de la matrice A"""
    # J'ai modifié pour tenir compte de la convention de Python:
    L = i + (j) * Nx
    return L


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

## Paramètres modifiables pour l'exécution du code

# Pour afficher les schémas de la structure du maillage
Afficher_schema = True

# Pour faire le code à 9 points ou pas
Neuf_points = True

# Main:

if __name__ == "__main__":

    if Afficher_schema == True:
        Image("./Images/LegendGrille.jpg", width=200, height=200)
        Image("./Images/Grille.jpg", width=400, height=400)
        Image("./Images/situations.png", width=500, height=600)
        Image("./Images/PML2.png", width=300, height=300)

    Source_Map = np.zeros([Nx, Ny], dtype=np.complex)
    h = dx
    r_max = 8 * h
    ## Pourrait probablement être intégré dans la construction de A...
    for i in range(Nx):
        for j in range(Ny):
            r = np.sqrt((i - S_x) ** 2 + (j - S_y) ** 2) * h
            if r < r_max:
                if r == 0:
                    Source_Map[i, j] = 1 / np.sqrt(h / 2)
                else:
                    val = np.exp(1j * k2_eau * r) / np.sqrt(r)
                    # Vérifier si on devrait pas prendre la -1j*k2_eau*r
                    Source_Map[i, j] = np.exp(1j * np.real(k2_eau) * r) / np.sqrt(r)
    plt.imshow(np.transpose(np.abs(np.abs(Source_Map))))

    # Construction de la carte des noeuds
    Map = np.ones([Nx, Ny])
    # Création de la map avec les points
    Map = Bois(Map, centre_bois_x, centre_bois_y, Nx_Bois, Ny_Bois)
    Map = PML(Map, N_PML)
    Map = Source(Map, S_x, S_y)

    ## Map Sans bois
    MapSB = np.ones([Nx, Ny])
    MapSB = PML(MapSB, N_PML)
    MapSB = Source(MapSB, S_x, S_y)

    Display_Map = np.copy(Map)
    Display_Map[np.logical_and(Map > 10, Map < 19)] = 0
    Display_Map[np.logical_and(Map > 2, Map < 11)] = 3
    Display_Map[Map == 19] = 4

    cmap = plt.cm.get_cmap('jet', 19)
    plt.figure(figsize=(6, 6))
    plt.imshow(np.transpose(Map), cmap=cmap)
    cbar = plt.colorbar(ticks=np.arange(1, 20), norm=np.arange(1, 20))

    plt.title("La carte, chaque couleur=1 type de point", fontsize=14)
    plt.xlabel("x", fontsize=20)
    plt.ylabel("y", fontsize=20)
    plt.show()

    fig, ax = plt.subplots(figsize=(6, 6))
    cmap = plt.cm.get_cmap('jet', 5)
    plt.imshow(np.transpose(Display_Map), cmap=cmap)

    ["PML", "Eau", "Bois", "Frontières", "Source"]

    # **********************Construction de la matrice A************************

    # L'ordre des coefficients est toujours
    # [p(i-2,j),p(i-1,j) ,p(i,j-2),p(i,j-1),p(i,j),p(i+1,j),p(i+2,j),p(i,j+1),p(i,j+2)]

    # Cas 1:
    if Neuf_points == True:
        Coeff1 = [0, 1, 0, 1, -(4 - k2_eau * h ** 2), 1, 0, 1, 0]
    else:
        # Version à 9 points
        # [p(i-1,j-1),p(i-1,j) ,p(i-1,j+1),p(i,j-1),p(i,j),p(i,j+1),p(i+1,j-1),p(i+1,j),p(i+1,j+1)]
        Coeff1 = [1, 4, 1, 4, -20 + 6 * h ** 2 * k2_eau, 4, 1, 4, 1]

    # Cas 2:
    if Neuf_points == True:
        Coeff2 = [0, 1, 0, 1, -(4 - k2_bois * h ** 2), 1, 0, 1, 0]
    else:

        # Version à 9 points
        # [p(i-1,j-1),p(i-1,j) ,p(i-1,j+1),p(i,j-1),p(i,j),p(i,j+1),p(i+1,j-1),p(i+1,j),p(i+1,j+1)]
        Coeff2 = [1, 4, 1, 4, -20 + 6 * h ** 2 * k2_bois, 4, 1, 4, 1]

    # Cas 3 à 10:
    gamma_eau = rho_eau * (alpha_eau * B_eau + 1j * omega)
    gamma_bois = rho_bois * (alpha_bois * B_bois + 1j * omega)

    Coeff3 = Coeff_Frontiere(gamma_eau, gamma_bois, -1 / np.sqrt(2), -1 / np.sqrt(2))
    Coeff4 = Coeff_Frontiere(gamma_eau, gamma_bois, 0, -1)
    Coeff5 = Coeff_Frontiere(gamma_bois, gamma_eau, 1 / np.sqrt(2), -1 / np.sqrt(2))  # -ny
    Coeff6 = Coeff_Frontiere(gamma_bois, gamma_eau, 1, 0)
    Coeff7 = Coeff_Frontiere(gamma_bois, gamma_eau, 1 / np.sqrt(2), 1 / np.sqrt(2))
    Coeff8 = Coeff_Frontiere(gamma_bois, gamma_eau, 0, 1)
    Coeff9 = Coeff_Frontiere(gamma_eau, gamma_bois, -1 / np.sqrt(2), 1 / np.sqrt(2))
    Coeff10 = Coeff_Frontiere(gamma_eau, gamma_bois, -1, 0)

    # Cas 11 à 18 (PML): À compléter !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    # Cas 19 (source): Option 2
    Coeff19 = [0, 1, 0, 1, -(4 - k2_eau * h ** 2), 1, 0, 1, 0]

    Dict_Coeff = {1: Coeff1, 2: Coeff2, 3: Coeff3, 4: Coeff4, 5: Coeff5, 6: Coeff6, 7: Coeff7, 8: Coeff8, 9: Coeff9,
                  10: Coeff10, 19: Coeff19}

    A = np.zeros([Nx * Ny, Nx * Ny], dtype=complex)
    b = np.zeros([Nx * Ny], dtype=complex)

    for i in range(Nx):
        for j in range(Ny):
            L = p(i, j)

            Type = Map[i, j]
            if np.logical_and(Type >= 11, Type <= 18):
                Coefficient = Coeff_PML(Type, i, j, h, Nx, Ny, k2_eau)
            else:
                Coefficient = Dict_Coeff[Type]

            if np.logical_and(np.logical_or(Type == 1, Type == 2), Neuf_points == True):
                Position = [p(i - 1, j - 1), p(i - 1, j), p(i - 1, j + 1), p(i, j - 1), p(i, j), p(i, j + 1),
                            p(i + 1, j - 1), p(i + 1, j), p(i + 1, j + 1)]
            else:
                Position = [p(i - 2, j), p(i - 1, j), p(i, j - 2), p(i, j - 1), p(i, j), p(i + 1, j), p(i + 2, j),
                            p(i, j + 1), p(i, j + 2)]

            for k, pos in enumerate(Position):
                if np.logical_and(pos >= 0, pos < (Nx * Ny)):
                    A[L, int(pos)] = Coefficient[k]
            #         if Type==19:
            #             b[L]=h**2*rho_eau*p_source
            b[L] = Source_Map[i, j] * h ** 2 * rho_eau * p_source

    # Matrice sans bois
    A_SB = np.zeros([Nx * Ny, Nx * Ny], dtype=complex)

    for i in range(Nx):
        for j in range(Ny):
            L = p(i, j)

            Type = MapSB[i, j]
            if np.logical_and(Type >= 11, Type <= 18):
                Coefficient = Coeff_PML(Type, i, j, h, Nx, Ny, k2_eau)
            else:
                Coefficient = Dict_Coeff[Type]

            if np.logical_and(np.logical_or(Type == 1, Type == 2), Neuf_points == True):
                Position = [p(i - 1, j - 1), p(i - 1, j), p(i - 1, j + 1), p(i, j - 1), p(i, j), p(i, j + 1),
                            p(i + 1, j - 1), p(i + 1, j), p(i + 1, j + 1)]
            else:
                Position = [p(i - 2, j), p(i - 1, j), p(i, j - 2), p(i, j - 1), p(i, j), p(i + 1, j), p(i + 2, j),
                            p(i, j + 1), p(i, j + 2)]
            for k, pos in enumerate(Position):
                if np.logical_and(pos >= 0, pos < (Nx * Ny)):
                    A_SB[L, int(pos)] = Coefficient[k]
            #         if Type==19:
            #             b[L]=h**2*rho_eau*p_source
            b[L] = Source_Map[i, j] * h ** 2 * rho_eau * p_source

    ## Calcul du temps d'inversion
    t0 = time.perf_counter()
    # sol=np.linalg.solve(A,b)
    sol = scipy.sparse.linalg.spsolve(scipy.sparse.csc_matrix(A), b)
    solSB = scipy.sparse.linalg.spsolve(scipy.sparse.csc_matrix(A_SB), b)
    t = time.perf_counter() - t0
    print("Temps pour inverser les deux matrices: {:.3f} s.".format(t))

    # Création map de solution
    MapSol = np.zeros([Nx, Ny], dtype=complex)
    MapSolSB = np.zeros([Nx, Ny], dtype=complex)

    for i in range(Nx):
        for j in range(Ny):
            MapSol[i, j] = sol[int(p(i, j))]
            MapSolSB[i, j] = solSB[int(p(i, j))]

    ## Création des figures pour la distribution
    fig, ax = plt.subplots(1, 4, figsize=(15, 5))
    ax[1].imshow(np.transpose(np.log(abs((MapSol)) + 1e-300)), cmap="gist_heat", alpha=1, interpolation="none")
    ax[0].imshow(np.transpose(Display_Map), alpha=1, cmap="jet")

    ax[2].imshow(np.transpose((np.log(abs(MapSolSB)) + 1e-250)), alpha=1.0, cmap="gist_heat", interpolation="none")
    # ax[3].imshow(np.log(abs(MapSolSB)-abs(MapSol)+1e-300),alpha=1.0,cmap="coolwarm")
    ax[0].set_title("La map physique")
    ax[1].set_title("Distribution  avec bois")
    ax[2].set_title("Distribution  sans bois")

    # ax[3].imshow(np.transpose(np.log(abs(MapSol)+1e-300))-np.transpose(np.log(abs(MapSolSB)+1e-300)),
    # cmap="gist_heat",alpha=1,interpolation="gaussian")

    Diff = abs(np.real(MapSol)) - abs(np.real(MapSolSB))
    Diff = Diff[(S_x - 15):(S_x + 15), (S_y - 15):(S_y + 15)]
    Diff = Diff + 2 * abs(np.min(Diff))
    Diff2 = Display_Map[(S_x - 15):(S_x + 15), (S_y - 15):(S_y + 15)]
    Diff = np.log(Diff)

    Cs = ax[3].imshow(np.transpose(Diff), cmap="gist_heat", alpha=1, interpolation="gaussian")
    ax[3].imshow(np.transpose(Diff2), cmap="binary", alpha=0.2, interpolation="none")

    ax[3].set_title("Différence de pression\n autour de la source")
    # ax[3].imshow(np.transpose(np.log(abs(MapSol)-abs(MapSolSB)+1e-20)),cmap="gist_heat",alpha=1,
    # interpolation="gaussian")
    # fig.colorbar(Cs,ax=ax[3])
    plt.show()

    abs(MapSol[S_x + 4, S_y + 4]) - abs(MapSolSB[S_x + 4, S_y + 4])

    ## Sauvegarde des résultats
    np.save("Resultats/MapSol_cylindrique_9points.npy", MapSol)
    np.save("Resultats/MapSolSB_cylindrique_9points.npy", MapSolSB)

    ## load des résultats
    Cylindrique_5points = np.load("Resultats/MapSol_cylindrique_5points.npy")
    CylindriqueSB_5points = np.load("Resultats/MapSolSB_cylindrique_5points.npy")
    Cylindrique_9points = np.load("Resultats/MapSol_cylindrique_9points.npy")
    CylindriqueSB_9points = np.load("Resultats/MapSolSB_cylindrique_9points.npy")
    Ponctuelle_5points = np.load("Resultats/MapSol_ponctuelle_5points.npy")
    PonctuelleSB_5points = np.load("Resultats/MapSolSB_ponctuelle_5points.npy")
    Ponctuelle_9points = np.load("Resultats/MapSol_ponctuelle_9points.npy")
    PonctuelleSB_9points = np.load("Resultats/MapSolSB_ponctuelle_9points.npy")

    ## 4 différentes méthodes de représentation
    fig, ax = plt.subplots(2, 2, figsize=(10, 10))
    ax[0, 0].imshow(np.transpose(np.log(abs(Cylindrique_5points))), interpolation="none")
    ax[0, 0].set_title("Cylindrique 5 points")
    ax[1, 0].imshow(np.transpose(np.log(abs(Cylindrique_9points))), interpolation="none")
    ax[1, 0].set_title("Cylindrique 9 points")

    ax[0, 1].imshow(np.transpose(np.log(abs(Ponctuelle_5points))), interpolation="none")
    ax[0, 1].set_title("Ponctuelle 5 points")
    ax[1, 1].imshow(np.transpose(np.log(abs(Ponctuelle_9points))), interpolation="none")
    ax[1, 1].set_title("Ponctuelle 9 points")

    plt.savefig("Resultats/Comparaison")

    ##
    fig, ax = plt.subplots(1, 2, figsize=(15, 5))
    ax[0].imshow(np.log(abs(MapSol) + 1e-300), cmap="coolwarm", alpha=1, interpolation="none")
    ax[0].set_xlim(20, 60)
    ax[0].set_ylim(20, 60)

    ax[1].imshow(abs(Map) == 5, cmap="coolwarm", alpha=1)
    ax[1].set_xlim(20, 60)
    ax[1].set_ylim(20, 60)

    abs(MapSol[45, 37]) - abs(MapSol[42, 38])_Graphique():