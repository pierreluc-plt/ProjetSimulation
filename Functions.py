import numpy as np


# Importer ce qui est nécessaire pour les fonctions

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


# Position du point
def p(i, j):
    """Le noeud (i,j) correspond à la ligne L de la matrice A"""
    # J'ai modifié pour tenir compte de la convention de Python:
    L = i + (j) * Nx
    return L


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
