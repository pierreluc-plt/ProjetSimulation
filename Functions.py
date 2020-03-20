import numpy as np
import matplotlib.pyplot as plt
import time
import scipy.sparse.csc, scipy.sparse.linalg

# Importer ce qui est nécessaire pour les fonctions

def p(i, j,Nx):
    """Le noeud (i,j) correspond à la ligne L de la matrice A"""
    # J'ai modifié pour tenir compte de la convention de Python:
    L = i + (j) * Nx
    return L


# Emplacement du bois:
def Bois(Map, centre_bois_x, centre_bois_y, Nx_Bois, Ny_Bois):
    """Cette fonction va placer les points pour le bois et les frontières du bois. Lorsqu'il s'agit de bois pur\
    (pas une frontière), le # du point est 2. Les autres points aux frontières sont numérotés de 3 à 10 selon\
    la convention exprimé ci-haut,
    Nx et Ny doivent être des multiples de 2\
    Il faudra modifier cette fonction lorsqu'on voudra ajouter des formes plus compliquées"""
    # On met tout en bois
    Map[int(centre_bois_x - Nx_Bois / 2):int(centre_bois_x + Nx_Bois / 2) + 1,
    int(centre_bois_y - Ny_Bois / 2):int(centre_bois_y + Ny_Bois / 2) + 1] = 2
    # On traite les frontières
    # Point 3
    Map[int(centre_bois_x - Nx_Bois / 2), int(centre_bois_y - Ny_Bois / 2)] = 3
    # Point 4
    Map[int(centre_bois_x - Nx_Bois / 2) + 1:int(centre_bois_x + Nx_Bois / 2), int(centre_bois_y - Ny_Bois / 2)] = 4
    # Point 5
    Map[int(centre_bois_x + Nx_Bois / 2), int(centre_bois_y - Ny_Bois / 2)] = 5
    # Point 6
    Map[int(centre_bois_x + Nx_Bois / 2), int(centre_bois_y - Ny_Bois / 2) + 1:int(centre_bois_y + Ny_Bois / 2)] = 6
    # Point 7
    Map[int(centre_bois_x + Nx_Bois / 2), int(centre_bois_y + Ny_Bois / 2)] = 7
    # Point 8
    Map[int(centre_bois_x - Nx_Bois / 2) + 1:int(centre_bois_x + Nx_Bois / 2), int(centre_bois_y + Ny_Bois / 2)] = 8
    # Point 9
    Map[int(centre_bois_x - Nx_Bois / 2), int(centre_bois_y + Ny_Bois / 2)] = 9
    # Point 10
    Map[int(centre_bois_x - Nx_Bois / 2), int(centre_bois_y - Ny_Bois / 2) + 1:int(centre_bois_y + Ny_Bois / 2)] = 10
    return Map


# Emplacement des PML
def PML(Map, N_PML):
    """Cette fonction vient placer les points de PML selon la convention décrite ci haut"""
    # Point 11
    Map[0:N_PML, -N_PML:] = 11
    # Point 12
    Map[N_PML:-N_PML, -N_PML:] = 12
    # Point 13
    Map[-N_PML:, -N_PML:] = 13
    # Point 14
    Map[-N_PML:, N_PML:-N_PML] = 14
    # Point 15
    Map[-N_PML:, 0:N_PML] = 15
    # Point 16
    Map[N_PML:-N_PML, 0:N_PML] = 16
    # Point 17
    Map[0:N_PML, 0:N_PML] = 17
    # Point 18
    Map[0:N_PML, N_PML:-N_PML] = 18
    return Map


# Emplacement de la source
def Source(Map, S_x, S_y):
    Map[S_x, S_y] = 19
    return Map





# Coefficients aux frontières
def Coeff_Frontiere(gamma1, gamma2, nx, ny):
    return [nx / gamma1, -4 * nx / gamma1, ny / gamma1, -4 * ny / gamma1,
            (nx + ny) * (gamma1 + gamma2) / (gamma1 + gamma2) \
        , -4 * nx / gamma2, nx / gamma2, -4 * ny / gamma2, ny / gamma2]  # Version avec gradient


# [nx/gamma1,-2*nx/gamma1,ny/gamma1,-2*ny/gamma1,(nx+ny)*gamma1*gamma2/(gamma1-gamma2),2*nx/gamma2,-nx/gamma2,
# 2*ny/gamma2,-2*ny/gamma2]


# Coefficients pour les PML
def Coeff_PML(Type, i, j, h, Nx, Ny, k2_eau):
    k = np.sqrt(k2_eau)
    x = i * h + h / 2
    y = j * h + h / 2  # Pour éviter les divisions par 0 ?

    if Type == 11:
        x0 = 0
        y0 = h * Ny  # Ny ou Ny-1 ???
        Beta_x = 1j / ((x - x0) * (k * abs(x - x0) + 1j))
        Beta_y = 1j / ((y - y0) * (k * abs(y - y0) + 1j))
        Gamma_x = 1 + 1j / k / (abs(x0 - x))
        Gamma_y = 1 + 1j / k / (abs(y0 - y))
        Coeff = [0, 0, 1 / Gamma_y ** 2, (-h * Beta_y - 2) / Gamma_y ** 2,
                 k ** 2 * h ** 2 + ((1 + h * Beta_y) / Gamma_y ** 2) + ((1 - h * Beta_x) / Gamma_x ** 2), \
                 (h * Beta_x - 2) / Gamma_x ** 2, 1 / Gamma_x ** 2, 0, 0]
    if Type == 12:
        y0 = h * Ny  # Ny ou Ny-1 ???
        Beta_y = 1j / ((y - y0) * (k * abs(y - y0) + 1j))
        Gamma_y = 1 + 1j / k / (abs(y0 - y))
        Coeff = [0, 1, 1 / Gamma_y ** 2, -(2 + h * Beta_y) / Gamma_y ** 2,
                 k ** 2 * h ** 2 + ((1 + h * Beta_y) / Gamma_y ** 2) - 2, \
                 1, 0, 0, 0]

    if Type == 13:
        x0 = h * Nx  # Nx ou Nx-1 ???
        y0 = h * Ny  # Ny ou Ny-1 ???
        Beta_x = 1j / ((x - x0) * (k * abs(x - x0) + 1j))
        Beta_y = 1j / ((y - y0) * (k * abs(y - y0) + 1j))
        Gamma_x = 1 + 1j / k / (abs(x0 - x))
        Gamma_y = 1 + 1j / k / (abs(y0 - y))
        Coeff = [1 / Gamma_x ** 2, -(h * Beta_x + 2) / Gamma_x ** 2, 1 / Gamma_y ** 2, -(h * Beta_y + 2) / Gamma_y ** 2, \
                 k ** 2 * h ** 2 + ((1 + h * Beta_y) / Gamma_y ** 2) + ((1 + h * Beta_x) / Gamma_x ** 2), 0, 0, 0, 0]
    if Type == 14:
        x0 = h * Nx  # Nx ou Nx-1 ???
        Beta_x = 1j / ((x - x0) * (k * abs(x - x0) + 1j))
        Gamma_x = 1 + 1j / k / (abs(x0 - x))
        Coeff = [1 / Gamma_x ** 2, -(h * Beta_x + 2) / Gamma_x ** 2, 0, 1, \
                 k ** 2 * h ** 2 - 2 + ((1 + h * Beta_x) / Gamma_x ** 2), 0, 0, 1, 0]

    if Type == 15:
        x0 = h * Nx  # Nx ou Nx-1 ???
        y0 = 0
        Beta_x = 1j / ((x - x0) * (k * abs(x - x0) + 1j))
        Beta_y = 1j / ((y - y0) * (k * abs(y - y0) + 1j))
        Gamma_x = 1 + 1j / k / (abs(x0 - x))
        Gamma_y = 1 + 1j / k / (abs(y0 - y))
        Coeff = [-1 / Gamma_x ** 2, -(h * Beta_x + 2) / Gamma_x ** 2, 0, 0, \
                 k ** 2 * h ** 2 + ((1 - h * Beta_y) / Gamma_y ** 2) + ((1 + h * Beta_x) / Gamma_x ** 2), 0, 0,
                 (h * Beta_y - 2) / Gamma_y ** 2, 1 / Gamma_y ** 2]
    if Type == 16:
        y0 = 0
        Beta_y = 1j / ((y - y0) * (k * abs(y - y0) + 1j))
        Gamma_y = 1 + 1j / k / (abs(y0 - y))
        Coeff = [0, 1, 0, 0, \
                 k ** 2 * h ** 2 + ((1 + h * Beta_y) / Gamma_y ** 2) - 2, 1, 0, (h * Beta_y - 2) / Gamma_y ** 2,
                 1 / Gamma_y ** 2]

    if Type == 17:
        x0 = 0
        y0 = 0
        Beta_x = 1j / ((x - x0) * (k * abs(x - x0) + 1j))
        Beta_y = 1j / ((y - y0) * (k * abs(y - y0) + 1j))
        Gamma_x = 1 + 1j / k / (abs(x0 - x))
        Gamma_y = 1 + 1j / k / (abs(y0 - y))
        Coeff = [0, 0, 0, 0, k ** 2 * h ** 2 + ((1 - h * Beta_y) / Gamma_y ** 2) + ((1 - h * Beta_x) / Gamma_x ** 2), \
                 (h * Beta_x - 2) / Gamma_x ** 2, 1 / Gamma_x ** 2, (h * Beta_y - 2) / Gamma_y ** 2, 1 / Gamma_y ** 2]

    if Type == 18:
        x0 = 0
        Beta_x = 1j / ((x - x0) * (k * abs(x - x0) + 1j))
        Gamma_x = 1 + 1j / k / (abs(x0 - x))
        Coeff = [0, 0, 0, 1, k ** 2 * h ** 2 - 2 + ((1 - h * Beta_x) / Gamma_x ** 2), (h * Beta_x - 2) / Gamma_x ** 2,
                 1 / Gamma_x ** 2, 1, 0]
    return Coeff


def Construction_Map(Nx,Ny,Nx_Bois,Ny_Bois,centre_bois_x, centre_bois_y,S_x,S_y,dx,N_PML,plot=True):
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

    if plot==True:

        fig,ax=plt.subplots(1,2,figsize=(20,8))
        cmap = plt.cm.get_cmap('jet', 19)

        ax[0].imshow(np.transpose(Map), cmap=cmap)
        #cbar = plt.colorbar(ticks=np.arange(1, 20), norm=np.arange(1, 20),ax=ax[0])
        ax[0].set_title("La carte, chaque couleur=1 type de point", fontsize=12)
        ax[0].set_xlabel("x", fontsize=20)
        ax[0].set_ylabel("y", fontsize=20)

        #cmap = plt.cm.get_cmap('jet', 5)
        ax[1].imshow(np.transpose(Display_Map), cmap=cmap)
        ax[1].set_title("Map Physique", fontsize=12)
        ax[1].set_xlabel("x", fontsize=20)
        ax[1].set_ylabel("y", fontsize=20)
        plt.show()

    return Map, MapSB,Display_Map



def Source_Cylindrique(Nx,Ny,S_x,S_y,dx,k2_eau,plot=False):
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

        if plot==True:
            fig, ax = plt.subplots(1, 1, figsize=(20, 8))
            cmap = plt.cm.get_cmap('jet', 19)

            ax.imshow(np.transpose(abs(Source_Map)), cmap=cmap)
            ax.set_title("La  source cylindrique", fontsize=12)
            ax.set_xlabel("x", fontsize=20)
            ax.set_ylabel("y", fontsize=20)
            plt.show()
        return Source_Map


def Construction_A(Nx,Ny,dx,Neuf_points,k2_eau,k2_bois,gamma_eau,gamma_bois,rho_eau,p_source,SourceCylindrique,Map,MapSB,Source_Map):
    h=dx
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
            L = p(i, j, Nx)

            Type = Map[i, j]
            if np.logical_and(Type >= 11, Type <= 18):
                Coefficient = Coeff_PML(Type, i, j, h, Nx, Ny, k2_eau)
            else:
                Coefficient = Dict_Coeff[Type]

            if np.logical_and(np.logical_or(Type == 1, Type == 2), Neuf_points == True):
                Position = [p(i - 1, j - 1, Nx), p(i - 1, j, Nx), p(i - 1, j + 1, Nx), p(i, j - 1, Nx), p(i, j, Nx),
                            p(i, j + 1, Nx),
                            p(i + 1, j - 1, Nx), p(i + 1, j, Nx), p(i + 1, j + 1, Nx)]
            else:
                Position = [p(i - 2, j, Nx), p(i - 1, j, Nx), p(i, j - 2, Nx), p(i, j - 1, Nx), p(i, j, Nx),
                            p(i + 1, j, Nx), p(i + 2, j, Nx),
                            p(i, j + 1, Nx), p(i, j + 2, Nx)]

            for k, pos in enumerate(Position):
                if np.logical_and(pos >= 0, pos < (Nx * Ny)):
                    A[L, int(pos)] = Coefficient[k]

                    if Type==19:
                        b[L]=h**2*rho_eau*p_source
            if SourceCylindrique==True:
                b[L] = Source_Map[i, j] * h ** 2 * rho_eau * p_source

    # Matrice sans bois
    A_SB = np.zeros([Nx * Ny, Nx * Ny], dtype=complex)

    for i in range(Nx):
        for j in range(Ny):
            L = p(i, j, Nx)

            Type = MapSB[i, j]
            if np.logical_and(Type >= 11, Type <= 18):
                Coefficient = Coeff_PML(Type, i, j, h, Nx, Ny, k2_eau)
            else:
                Coefficient = Dict_Coeff[Type]

            if np.logical_and(np.logical_or(Type == 1, Type == 2), Neuf_points == True):
                Position = [p(i - 1, j - 1, Nx), p(i - 1, j, Nx), p(i - 1, j + 1, Nx), p(i, j - 1, Nx), p(i, j, Nx),
                            p(i, j + 1, Nx), p(i + 1, j - 1, Nx), p(i + 1, j, Nx), p(i + 1, j + 1, Nx)]
            else:
                Position = [p(i - 2, j, Nx), p(i - 1, j, Nx), p(i, j - 2, Nx), p(i, j - 1, Nx), p(i, j, Nx),
                            p(i + 1, j, Nx), p(i + 2, j, Nx),
                            p(i, j + 1, Nx), p(i, j + 2, Nx)]
            for k, pos in enumerate(Position):
                if np.logical_and(pos >= 0, pos < (Nx * Ny)):
                    A_SB[L, int(pos)] = Coefficient[k]
                    if Type==19:
                        b[L]=h**2*rho_eau*p_source
            if SourceCylindrique==True:
                b[L] = Source_Map[i, j] * h ** 2 * rho_eau * p_source

    return A, A_SB,b

def Resolution(A,A_SB,b,Nx,Ny):

    ## Calcul du temps d'inversion
    MapSol = np.zeros([Nx, Ny], dtype=complex)
    MapSolSB = np.zeros([Nx, Ny], dtype=complex)

    t0 = time.perf_counter()

    solSB = scipy.sparse.linalg.spsolve(scipy.sparse.csc_matrix(A_SB), b)
    sol = scipy.sparse.linalg.spsolve(scipy.sparse.csc_matrix(A), b)
    t = time.perf_counter() - t0
    print("Temps pour inverser les deux matrices: {:.3f} s.".format(t))
    # Création map de solution

    for i in range(Nx):
        for j in range(Ny):
            MapSol[i, j] = sol[int(p(i, j,Nx))]
            MapSolSB[i, j] = solSB[int(p(i, j,Nx))]

    return MapSol, MapSolSB


def Plots_Results(MapSol,MapSolSB,Display_Map,Interpolation="none"):
    ## Création des figures pour la distribution
    fig, ax = plt.subplots(2, 2, figsize=(25, 25))

    ax[1][0].imshow(np.transpose(Display_Map), alpha=1, cmap="jet", interpolation="none")
    ax[1][0].set_title("La map physique")

    ax[0][0].imshow(np.transpose(np.log(abs((MapSol)) + 1e-300)), cmap="gist_heat", alpha=1, interpolation=Interpolation)
    ax[0][0].set_title("Distribution  avec bois")

    ax[0][1].imshow(np.transpose((np.log(abs(MapSolSB)) + 1e-250)), alpha=1.0, cmap="gist_heat", interpolation=Interpolation)
    ax[0][1].set_title("Distribution  sans bois")

    Diff = abs(np.real(MapSol)) - abs(np.real(MapSolSB))
    #Diff = Diff[(S_x - 15):(S_x + 15), (S_y - 15):(S_y + 15)]
    Diff = Diff + 2 * abs(np.min(Diff))
    Diff = np.log(Diff)

    ax[1][1].imshow(np.transpose(Diff), cmap="gist_heat", alpha=1, interpolation=Interpolation)
    ax[1][1].imshow(np.transpose(Display_Map), cmap="binary", alpha=0.1, interpolation="none")
    plt.show()

