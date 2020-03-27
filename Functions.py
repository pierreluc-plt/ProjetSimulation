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

# Emplacement du nez du bois:
def Boisnez(Map, centre_bois_x, centre_bois_y, Nx_Bois, Ny_Bois, forme, coeff):
    # forme peut-etre soit un triangle ou un cercle
    # coeff est un angle en radian de la pointe pour le triangle ou le rayon du cercle
    if forme == 'triangle':
        hauteur = Nx_Bois/(2*np.tan(coeff/2))
        larg = Nx_Bois
        i=0
        up=1
        larg2 = larg
        while int(larg) > 1:
            if up == 1:
                # On met toute la pointe en bois
                Map[int(centre_bois_x-larg/2):int(centre_bois_x+larg/2)+1,int(centre_bois_y-Ny_Bois/2)-i]=2
                # On ajoute le contour
                Map[int(centre_bois_x-larg/2),int(centre_bois_y-Ny_Bois/2)-i]=20
                Map[int(centre_bois_x+larg/2),int(centre_bois_y-Ny_Bois/2)-i]=21
            else:
                Map[int(centre_bois_x-larg/2),int(centre_bois_y-Ny_Bois/2)-i]=20
                Map[int(centre_bois_x+larg/2),int(centre_bois_y-Ny_Bois/2)-i]=21
            larg1 = larg +1
            larg2 = larg
            i+=1
            larg = round(2*(hauteur-i)*np.tan(coeff/2))
            if (larg == larg2-2) or (larg == larg2-1) or (larg == larg2):
                up=1
            else:
                i-=1
                larg = larg2 - 1
                up=0
                
        # Pour ajouter la derniere pointe
        if up == 0:
            i+=1
        if larg1 == 3:
            Map[int(centre_bois_x-(larg1-1)/2+1),int(centre_bois_y-Ny_Bois/2)-i]=4
            #Map[int(centre_bois_x-larg/2),int(centre_bois_y-Ny_Bois/2)-i]=4
        else:
            Map[int(centre_bois_x-(larg1-2)/2):int(centre_bois_x+(larg1-2)/2),int(centre_bois_y-Ny_Bois/2)-i]=4
            
    elif forme == 'cercle':
        # hauteur est la distance entre le centre du cercle et le cote du rectangle
        hauteur = np.sqrt(coeff**2-(Nx_Bois/2)**2)
        centre_y = centre_bois_y - Ny_Bois/2 + hauteur
        x_val = np.arange(centre_bois_x-(Nx_Bois/2),centre_bois_x+Nx_Bois/2+1)
        for x in x_val:
            pos_y = round(-np.sqrt(int((coeff)**2 - (x-centre_bois_x)**2)) + centre_y)
            # Pour faire la frontiere du demi-cercle
            Map[int(x),int(pos_y)]=22
            if pos_y < (centre_bois_y-Ny_Bois/2+1):
                # Pour faire l'interieur du demi-cercle
                Map[int(x),(int(pos_y+1)):int(centre_bois_y-Ny_Bois/2+1)] = 2
        # Pour s'assurer qu'il n'y a pas de trou dans la frontiere du demi-cercle
        x_val1 = np.arange(centre_bois_x-Nx_Bois/2+1,centre_bois_x+Nx_Bois/2)
        for x in x_val1:
            pos_y = round(-np.sqrt(int((coeff)**2 - (x-centre_bois_x)**2)) + centre_y)
            for i in np.arange(1,5):
                if (Map[int(x),int(pos_y+i)] != 6 and Map[int(x-1),int(pos_y+i)] == 1 and Map[int(x-1),int(pos_y)] == 1 \
                    and x < centre_bois_x) or (Map[int(x),int(pos_y+i)] != 6 and Map[int(x+1),int(pos_y+i)] == 1 \
                                               and Map[int(x+1),int(pos_y)] == 1 and x > centre_bois_x):
                    Map[int(x),int(pos_y+i)]=22
    else:
        print('La variable de forme a ete mal defini.')
    return Map

# Calcul surface du demi-cercle
def surfdemicercle(S_x,S_y,centre_bois_x,centre_cercle_y,Nx_Bois,Ny_Bois,centre_bois_y,coeff):
    limite_pente = (-(-2*S_y*centre_bois_x+2*centre_bois_x*centre_cercle_y+2*S_x*(S_y-centre_cercle_y))+(np.sqrt((-2*S_y*centre_bois_x+2*centre_bois_x*centre_cercle_y+2*S_x*(S_y-centre_cercle_y))**2 - (4*(2*S_x*centre_bois_x+coeff**2-centre_bois_x**2-S_x**2)*(coeff**2-(S_y-centre_cercle_y)**2)))))/(2*(2*S_x*centre_bois_x+coeff**2-centre_bois_x**2-S_x**2))
    limite_b = S_y - limite_pente*S_x
    limite_x = (centre_bois_x-limite_pente*(limite_b - centre_cercle_y))/(limite_pente**2 + 1)
    limite_y = limite_x*limite_pente + limite_b
    surf = 2*coeff*np.arcsin((np.sqrt((limite_x - centre_bois_x-Nx_Bois/2)**2 + (centre_bois_y-Ny_Bois/2 - limite_y)**2))/(2*coeff))
    return surf

# Calcul coin_pente
def coinpente(S_y,Nx_Bois,Ny_Bois,centre_bois_y,coeff,cote):
    if coeff != Nx_Bois/2 and S_y > centre_bois_y - Ny_Bois/2:
        coin_pente = cote*(Nx_Bois/2)/(-np.sqrt(coeff**2 - (Nx_Bois/2)**2))
    elif coeff == Nx_Bois/2 and S_y > centre_bois_y - Ny_Bois/2:
        coin_pente = -999*cote
    elif coeff != Nx_Bois/2 and S_y <= centre_bois_y - Ny_Bois/2:
        coin_pente = cote*(Nx_Bois/2)/(np.sqrt(coeff**2 - (Nx_Bois/2)**2))
    else:
        coin_pente = 999*cote
    return coin_pente

# Surface en contact direct avec la source
def surface_directe(S_x, S_y, centre_bois_x, centre_bois_y, Nx_Bois, Ny_Bois, forme, coeff):
    if S_x < (centre_bois_x-Nx_Bois/2):
        if forme == 'triangle':
            y_min = centre_bois_y - Ny_Bois/2 + (centre_bois_x-Nx_Bois/2-S_x)/np.tan(coeff/2)
            y_max = centre_bois_y - Ny_Bois/2 - Nx_Bois/(2*np.tan(coeff/2)) - (centre_bois_x-S_x)/np.tan(coeff/2)
            if S_y <= y_max:
                surf = 2*(Nx_Bois/2)/(np.sin(coeff/2))
            elif S_y <= y_min and S_y > y_max:
                surf = (Nx_Bois/2)/(np.sin(coeff/2))
            elif S_y > y_min and S_y <= (centre_bois_y + Ny_Bois/2):
                surf = 0
            else:
                surf = Nx_Bois
        elif forme == 'cercle':
            centre_cercle_y = centre_bois_y - Ny_Bois/2 + np.sqrt(coeff**2-(Nx_Bois/2)**2)
            coin_pente = coinpente(S_y,Nx_Bois,Ny_Bois,centre_bois_y,coeff,1)
            
            if S_y > centre_bois_y - Ny_Bois/2:
                S_pente = (centre_bois_y-Ny_Bois/2 - S_y)/(centre_bois_x-Nx_Bois/2 - S_x)
                if S_pente > coin_pente:
                    #calcul arc
                    surf = surfdemicercle(S_x,S_y,centre_bois_x,centre_cercle_y,Nx_Bois,Ny_Bois,centre_bois_y,coeff)
                elif S_pente <= coin_pente and S_y <= (centre_bois_y + Ny_Bois/2):
                    surf = 0
                else:
                    surf = Nx_Bois
            elif S_y <= centre_bois_y - Ny_Bois/2 and S_y >= centre_cercle_y-coeff:
                surf = surfdemicercle(S_x,S_y,centre_bois_x,centre_cercle_y,Nx_Bois,Ny_Bois,centre_bois_y,coeff)
            elif S_y < centre_cercle_y-coeff:
                S_pente = (centre_bois_y-Ny_Bois/2 - S_y)/(centre_bois_x+Nx_Bois/2 - S_x)
                if S_pente < coin_pente:
                    surf = surfdemicercle(S_x,S_y,centre_bois_x,centre_cercle_y,Nx_Bois,Ny_Bois,centre_bois_y,coeff)
                else:
                    surf = 2*coeff*np.arcsin(Nx_Bois/(2*coeff))
        surface = Ny_Bois + surf
    elif S_x >= (centre_bois_x-Nx_Bois/2) and S_x <= (centre_bois_x+Nx_Bois/2):
        if forme == 'triangle':
            if S_x < centre_bois_x:
                y_max = centre_bois_y - Ny_Bois/2 - Nx_Bois/(2*np.tan(coeff/2)) - (centre_bois_x-S_x)/np.tan(coeff/2)
                if S_y < y_max:
                    surface = 2*(Nx_Bois/2)/(np.sin(coeff/2))
                elif S_y < (centre_bois_y-Ny_Bois/2) and S_y >= y_max:
                    surface = (Nx_Bois/2)/(np.sin(coeff/2))
                elif S_y > (centre_bois_y+Ny_Bois/2):
                    surface = Nx_Bois
            elif S_x == centre_bois_x:
                if S_y < (centre_bois_y-Ny_Bois/2):
                    surface = 2*(Nx_Bois/2)/(np.sin(coeff/2))
                elif S_y > (centre_bois_y+Ny_Bois/2):
                    surface = Nx_Bois
            elif S_x > centre_bois_x:
                y_max = centre_bois_y - Ny_Bois/2 - Nx_Bois/(2*np.tan(coeff/2)) - (S_x-centre_bois_x)/np.tan(coeff/2)
                if S_y < y_max:
                    surface = 2*(Nx_Bois/2)/(np.sin(coeff/2))
                elif S_y < (centre_bois_y-Ny_Bois/2) and S_y >= y_max:
                    surface = (Nx_Bois/2)/(np.sin(coeff/2))
                elif S_y > (centre_bois_y+Ny_Bois/2):
                    surface = Nx_Bois

        elif forme == 'cercle':
            centre_cercle_y = centre_bois_y - Ny_Bois/2 + np.sqrt(coeff**2-(Nx_Bois/2)**2)
            if S_y > (centre_bois_y+Ny_Bois/2):
                surface = Nx_Bois
            elif S_x < centre_bois_x:
                surface = surfdemicercle(S_x,S_y,centre_bois_x,centre_cercle_y,Nx_Bois,Ny_Bois,centre_bois_y,coeff)
            elif S_x == centre_bois_x:
                surface = 2*coeff*np.arcsin(Nx_Bois/(2*coeff))
            elif S_x > centre_bois_x:
                surface = 2*coeff*np.arcsin(Nx_Bois/(2*coeff)) - surfdemicercle(S_x,S_y,centre_bois_x,centre_cercle_y,Nx_Bois,Ny_Bois,centre_bois_y,coeff)
    elif S_x > (centre_bois_x+Nx_Bois/2):
        if forme == 'triangle':
            y_min = centre_bois_y - Ny_Bois/2 + (S_x - centre_bois_x-Nx_Bois/2)/np.tan(coeff/2)
            y_max = centre_bois_y - Ny_Bois/2 - Nx_Bois/(2*np.tan(coeff/2)) - (S_x-centre_bois_x)/np.tan(coeff/2)
            if S_y <= y_max:
                surf = 2*(Nx_Bois/2)/(np.sin(coeff/2))
            elif S_y <= y_min and S_y > y_max:
                surf = (Nx_Bois/2)/(np.sin(coeff/2))
            elif S_y > y_min and S_y <= (centre_bois_y + Ny_Bois/2):
                surf = 0
            else:
                surf = Nx_Bois
        elif forme == 'cercle':
            centre_cercle_y = centre_bois_y - Ny_Bois/2 + np.sqrt(coeff**2-(Nx_Bois/2)**2)
            coin_pente = coinpente(S_y,Nx_Bois,Ny_Bois,centre_bois_y,coeff,-1)
            
            if S_y > centre_bois_y - Ny_Bois/2:
                S_pente = (centre_bois_y-Ny_Bois/2 - S_y)/(S_x - centre_bois_x-Nx_Bois/2)
                if S_pente < coin_pente:
                    #calcul arc
                    surf = 2*coeff*np.arcsin(Nx_Bois/(2*coeff)) - surfdemicercle(S_x,S_y,centre_bois_x,centre_cercle_y,Nx_Bois,Ny_Bois,centre_bois_y,coeff)
                elif S_pente >= coin_pente and S_y <= (centre_bois_y + Ny_Bois/2):
                    surf = 0
                else:
                    surf = Nx_Bois
            elif S_y <= centre_bois_y - Ny_Bois/2 and S_y >= centre_cercle_y-coeff:
                surf = 2*coeff*np.arcsin(Nx_Bois/(2*coeff)) - surfdemicercle(S_x,S_y,centre_bois_x,centre_cercle_y,Nx_Bois,Ny_Bois,centre_bois_y,coeff)
            elif S_y < centre_cercle_y-coeff:
                S_pente = -(centre_bois_y-Ny_Bois/2 - S_y)/(S_x - centre_bois_x+Nx_Bois/2)
                if S_pente < coin_pente:
                    surf = 2*coeff*np.arcsin(Nx_Bois/(2*coeff)) - surfdemicercle(S_x,S_y,centre_bois_x,centre_cercle_y,Nx_Bois,Ny_Bois,centre_bois_y,coeff)
                else:
                    surf = 2*coeff*np.arcsin(Nx_Bois/(2*coeff))
        surface = Ny_Bois + surf
    return surface


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


def Construction_Map(Nx,Ny,Nx_Bois,Ny_Bois,centre_bois_x, centre_bois_y,forme,coeff,S_x,S_y,dx,N_PML,plot=True):
    Map = np.ones([Nx, Ny])
    # Création de la map avec les points
    Map = Bois(Map, centre_bois_x, centre_bois_y, Nx_Bois, Ny_Bois)
    Map=Boisnez(Map, centre_bois_x, centre_bois_y, Nx_Bois, Ny_Bois, forme, coeff)
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
        r_max = Nx * h
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
            ax.set_title("La source cylindrique", fontsize=12)
            ax.set_xlabel("x", fontsize=20)
            ax.set_ylabel("y", fontsize=20)
            plt.show()
        return Source_Map


def Construction_A(Nx,Ny,dx,Neuf_points,k2_eau,k2_bois,gamma_eau,gamma_bois,rho_eau,p_source,SourceCylindrique,Map,MapSB,Source_Map,coeff,centre_bois_x,centre_bois_y,Nx_Bois,Ny_Bois):
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

    # Cas 20 à 21 (triangle)
    # Cas 20
    Nx20 = -np.cos(coeff/2) #-
    Ny20 = -np.sin(coeff/2) # -
    Coeff20= Coeff_Frontiere(gamma_eau, gamma_bois, Nx20, Ny20)
    # Cas 21
    Nx21 = np.cos(coeff/2)
    Ny21 = -np.sin(coeff/2)
    Coeff21= Coeff_Frontiere(gamma_eau, gamma_bois, Nx21, Ny21)
    
    # Cas 22 (Cercle)    
    # Voir la boucle plus bas

    # Cas 11 à 18 (PML):Dans les fonctions suivantes

    # Cas 19 (source): Option 2
    Coeff19 = [0, 1, 0, 1, -(4 - k2_eau * h ** 2), 1, 0, 1, 0]

    Dict_Coeff = {1: Coeff1, 2: Coeff2, 3: Coeff3, 4: Coeff4, 5: Coeff5, 6: Coeff6, 7: Coeff7, 8: Coeff8, 9: Coeff9,
                  10: Coeff10, 19: Coeff19,20:Coeff20,21:Coeff21}

    A = np.zeros([Nx * Ny, Nx * Ny], dtype=complex)
    b = np.zeros([Nx * Ny], dtype=complex)
    b_TFSF = np.zeros([Nx * Ny], dtype=complex)

    for i in range(Nx):
        for j in range(Ny):
            L = p(i, j, Nx)

            Type = Map[i, j]
            if np.logical_and(Type >= 11, Type <= 18):
                Coefficient = Coeff_PML(Type, i, j, h, Nx, Ny, k2_eau)
            elif Type == 22:
                Nx22 = (i-centre_bois_x)/coeff
                # Coordonnées en y du centre du cercle
                centre_y = centre_bois_y - Ny_Bois/2 + np.sqrt(coeff**2-(Nx_Bois/2)**2)
                Ny22 = (j-centre_y)/coeff
                Coefficient = Coeff_Frontiere(gamma_eau, gamma_bois, Nx22, Ny22)
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
    Q=np.zeros([Nx * Ny, Nx * Ny], dtype=complex)


    # Mask function
    Q_map = Map - MapSB

    Q_map[Q_map > 0] = 1

    where_0 = np.where(Q_map == 0)
    where_1 = np.where(Q_map == 1)
    Q_map[where_0] = 1
    Q_map[where_1] = 0

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

            Q[L,L]=Q_map[i,j]

    # Mask function

    # Source vector for TF/SF
    A_sp = scipy.sparse.csc_matrix(A)
    Q_sp = scipy.sparse.csc_matrix(Q)
    b_TFSF = (Q_sp.dot(A_sp) - A_sp.dot(Q_sp)).dot(b)


    return A, A_SB, b, b_TFSF,Q

def Resolution(A,A_SB,b,b_TFSF,Nx,Ny,D_x,D_y):

    ## Calcul du temps d'inversion
    MapSol = np.zeros([Nx, Ny], dtype=complex)
    MapSolSB = np.zeros([Nx, Ny], dtype=complex)
    MapSol_TFSF = np.zeros([Nx, Ny], dtype=complex)

    t0 = time.perf_counter()

    solSB = scipy.sparse.linalg.spsolve(scipy.sparse.csc_matrix(A_SB), b)
    sol = scipy.sparse.linalg.spsolve(scipy.sparse.csc_matrix(A), b)
    sol_TFSF = scipy.sparse.linalg.spsolve(scipy.sparse.csc_matrix(A), b_TFSF)
    
    t = time.perf_counter() - t0
    print("Temps pour inverser les deux matrices: {:.3f} s.".format(t))
    # Création map de solution

    for i in range(Nx):
        for j in range(Ny):
            MapSol[i, j] = sol[int(p(i, j,Nx))]
            MapSolSB[i, j] = solSB[int(p(i, j,Nx))]
            MapSol_TFSF[i, j] = sol_TFSF[int(p(i, j,Nx))]
    
    #Intensite au detecteur
    P_detecteur = MapSol[D_x, D_y]

    return MapSol, MapSolSB, MapSol_TFSF, P_detecteur


def Plots_Results(MapSol,MapSolSB,MapSol_TFSF,Display_Map,Interpolation="none"):
    ## Création des figures pour la distribution
    fig, ax = plt.subplots(2, 2, figsize=(25, 25))

    ax[1][0].imshow(np.transpose(Display_Map), alpha=1, cmap="jet", interpolation="none")
    ax[1][0].set_title("La map physique")

    ax[0][0].imshow(np.transpose(np.log10(np.real(MapSol)-1.2*np.min(np.real(MapSol)))), cmap="jet", alpha=1, interpolation=Interpolation)
    ax[0][0].set_title("Distribution  avec bois")
    heatmap = ax[0][0].pcolor(np.transpose(np.log10(np.real(MapSol)-1.2*np.min(np.real(MapSol)))), cmap="jet")

    # legend
    cbar = plt.colorbar(heatmap)

    ax[0][1].imshow(np.transpose(np.log10(np.real(MapSolSB)-1.2*np.min(np.real(MapSolSB)))), alpha=1.0, cmap="jet", interpolation=Interpolation)
    ax[0][1].set_title("Distribution  sans bois")

    Diff = -((np.real(MapSol)) -(np.real(MapSolSB)))
    #Diff = Diff[(S_x - 15):(S_x + 15), (S_y - 15):(S_y + 15)]
    Diff = Diff + 2 * abs(np.min(Diff))
    Diff = np.log(Diff)

    ax[1][1].imshow(np.transpose(np.real(MapSol_TFSF)), cmap="jet", alpha=1, interpolation=Interpolation)
    ax[1][1].imshow(np.transpose(Display_Map), cmap="binary", alpha=0.1, interpolation="none")
    ax[1][1].set_title("Différence entre les deux distributions")
    plt.show()


def Surface_equivalente(Nx_Bois,Ny_Bois,forme,coeff,P_incident,P_scattered,V_incident,Surface):
    """
    P_incident: Pression initiale/ incident sur le bateau 
        (En réalité la puissance sonore donné par Aire * Pression * Vitesse.
        Puisque la vitesse est la même est on suppose que l'aire sur lequel ils parcourent
        sont à peu près la même chose, on dit que pression équivaut puissance sonore.)
    P_scattered: Calculé avec formulation TF/SF
    V_incident: Volume/Aire incident. Si source ponctuelle, juste prendre dx * dy.
    V_scattered: Volume/Aire du bateau
    Surface: Surface/Périmètre étant exposé au P_incident    
    """
    
    if forme == 'triangle':
        aire_nez = Nx_Bois**2/(2*np.tan(coeff/2))
    elif forme == 'cercle':
        angle = 2*np.arcsin(Nx_Bois/(2*coeff))
        aire_nez = 0.5*coeff**2*(angle-np.sin(angle))
        V_scattered = Nx_Bois*Ny_Bois + aire_nez
    
    SER = Surface*(P_scattered/V_scattered)/(P_incident/V_incident)
    
    return SER    
