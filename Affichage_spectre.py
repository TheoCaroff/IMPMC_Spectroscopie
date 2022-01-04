# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 18:51:51 2020

@author: Theo_C
Ce fichier contient toutes les chose utiles pour tracer et traiter les spectres
Les fonctions colorimétrique on été faites par DUNE ANDRE et l'affichage LAB3D par TANYA ASSERAF'
"""

import numpy as np
import os
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import ticker as ticker
from scipy import interpolate


from Lecture_input import mono2tab
from Lecture_input import Readspectre
from Lecture_input import normYminmax
from Lecture_input import nm2cm1
from Lecture_input import Gauss

from Nettoyage_spectre import baseline_Rubberband

try :
    import colour
    from colour import plotting as cplot
except ModuleNotFoundError :
    print('ATTENTION MODULE COLOR-SCIENCE PAS INSTALLER')


def set_graph(FIGSIZE=[12, 10], DPI = 120, GRAPHOUT=0, grid=True, mode='default', xylim=[-0.1, 0.8, 0, 1]):
    '''
    Paramètre la figure qui va servir à afficher les spectres

    Parameters
    ----------
    FIGSIZE : TYPE, optional
        DESCRIPTION. The default is [12, 6].
    DPI : TYPE, optional
        DESCRIPTION. The default is 120.
    grid : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    None.

    '''
    top=0.97
    bottom=0.12
    left=0.08
    right=0.985
    wspace=0
    hspace=0
    
    
    if GRAPHOUT == 1:
        FIGSIZE=[5, 3];
        bottom=0.15
        left=0.13
    
    if GRAPHOUT== 2:
        FIGSIZE=[10, 5]
        DPI = 100
    
    else:
        pass
    
    plt.style.use({'figure.figsize': FIGSIZE, 'figure.dpi': DPI})
#   plt.figure(figsize=(FIGSIZE), dpi=DPI)
    if mode == 'default':
        plt.figure()
        plt.grid(grid);
    elif mode == 'CIE1931':
        # Plotting the *CIE 1931 Chromaticity Diagram*.
        # The argument *standalone=False* is passed so that the plot doesn't get
        # displayed and can be used as a basis for other plots.
        cplot.plot_chromaticity_diagram_CIE1931(standalone=False,
                                                bounding_box=xylim, x_tighten=True, y_tighten=True)
    
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
    
    return(plt.gcf(), plt.gca())


def Sav_fig(Titre='Pouet', Repertoire='Graph', colorplot=False):
    '''
    Cette fonction permet de sauvegarer un graph matplotlib puis l'affiche'

    Parameters
    ----------
    Titre : TYPE, optional
        Titre du graph. The default is 'Pouet'.
    Repertoire : TYPE, optional
        Repertoire dans lequel les figures vont être sauvergarder. The default is 'Graph'.
    cplot : TYPE, optional
        si utilisé pour le module color science. The default is False.

    Returns
    -------
    None.

    '''
    try:
        os.mkdir(Repertoire)
    except OSError:
        pass
    plt.legend()
    plt.savefig(Repertoire+os.sep+Titre, bbox_inches='tight')

    if colorplot:
        cplot.render(standalone=True)
    else:
        plt.show()
        
def Affichage_abs(Liste, Legende, Autoaxe=True, Xlim=[4000, 35000], Ylim=[0, 1.5],
                  SecondAxe=True, TITRE='superposition', AdditionTr=0,
                  linewidth=1.5, valeurnorm=1, Modeaff='ABScm', modecouleurs='auto', optionplot='',
                  SHOW=True, GRAPHOUT=1,COUPURENORMminmax=[400, 2500], newgraph=True, SpectreVIS=True):
    '''
    Cette fonction affiche des spectes optique obtenu en %T en absorbance et
    sauvegarde le graph dans un dossier graph placé dans le repertoire courant
    Parameters
    ----------
    
    Liste : Tableau de string
        Contient le nom des fichier à afficher.
    Legende : Tableau de string
        Legende des fichiers.
    Autoaxe : bool
        Pour rélger les axes manuellement ou automatiquement
    Xlim : TYPE, optional
        DESCRIPTION. The default is [4000, 35000].
    Ylim : TYPE, optional
        DESCRIPTION. The default is [0, 1.5].
    SecondAxe : bool
        True si deuxième axe en nm.
    TITRE : TYPE, optional
        DESCRIPTION. The default is 'superposition'.
    AdditionTr: float
        Correction de transmitance
    linewidth : TYPE, optional
        DESCRIPTION. The default is 1.5.
    valeurnorm: float, optional
        Valeur qui sert à normaliser l'ABS dans les modes ABSnormep et Epsilon
    Modeaff : TYPE, optional
        Choix du mode, entre ABScm, Transmittance, reflectance, ABSnormep (normalisation épaisseur),
        Epsilon. The default is 'ABScm'.
    modecouleurs : string, optional
        Permet de choisir le mode d'affichage des couleurs (trois mode existe :
            'auto', bigdata', 'manuel'. The default is 'auto'.
    optionplot : string, optional
        Choix des options des matplolib à utiliser en mode manuel
    SHOW : Bool, optional
        Pour afficher et sauvergarder la figure. The default is True.
    SpectreVIS = Bool, optional
        Pour afficher le spectre visible. The default is True.
    Returns
    -------
        La figure courante

    '''
    
    #Liste=mono2tab(Liste, np.size(Liste))
    Legende=mono2tab(Legende, np.size(Liste))
    optionplot=mono2tab(optionplot, np.size(Liste))
    AdditionTr=mono2tab(AdditionTr, np.size(Liste))
    valeurnorm=mono2tab(valeurnorm, np.size(Liste))
    
    if GRAPHOUT == 1:
        FIGSIZE=[15, 10]
        DPI = 120
        tailletexte=14
    elif GRAPHOUT == 2 :
        FIGSIZE=[10, 5]
        DPI = 100
        tailletexte=10
    
    
    if newgraph : set_graph(FIGSIZE, DPI) # Cf fonction, on créer la figure à la bonne taille/résolution

    

    colorsbigdata = plt.cm.nipy_spectral(np.linspace(0,1,np.size(Liste))) # On fixe la gamme de couleur utilisé

    if Autoaxe:
        pass
    else:
        plt.xlim(Xlim);
        plt.ylim(Ylim);
    

    for i in np.arange(0, np.size(Liste)): # on lit tous les fichier  
        Fichier = Liste[i];
        
        Xnm, Ytr = Readspectre(Fichier)
        Ytr= Ytr + AdditionTr[i]
        X = 1/(Xnm*1E-7);
        Y = -np.log10(Ytr);
        plt.xlabel("Nombre d'onde ($cm^{-1}$)");
        plt.ylabel('Absorbance')
            
        if Modeaff == 'ABScm':
            RAJOUT= ''
        
        elif Modeaff == 'ABSnm':
            X = Xnm;
            plt.xlabel("Longueur d'onde (en nm)");
            RAJOUT= '_ABSnm'
        
        elif Modeaff == 'ABSnorm_min_max':
            Y = normYminmax(Xnm, Y, COUPURENORMminmax)
            plt.ylabel('Absorbance normalisé (u.a)')
            RAJOUT = '_normminmax'
            
        elif Modeaff == 'Reflectance':
            Y = 1-Ytr;
            plt.ylabel('Reflectance')
            RAJOUT = '_Reflexion'
       
        elif Modeaff == 'ABSnormep':
            Y = Y/valeurnorm[i];
            plt.ylabel('Absorbance linéaire ($cm^{-1}$)') 
            RAJOUT = '_norm_ep'
        
        elif Modeaff == 'Epsilon':
            Y = Y/valeurnorm[i];
            plt.ylabel('$\\varepsilon (L.mol^{-1}.cm^{-1})$') 
            RAJOUT = '_Epsilon'
            
        elif Modeaff == 'ABSnorm_min_ep':
            Y = Y/valeurnorm[i];
            Y=Y-np.min(Y)
            plt.ylabel('Absorbance normalisé (u.a)')
            RAJOUT = '_normmin_ep'
        
        elif Modeaff == 'SubBaseline':
            INDEXUV=Xnm>355
            X=X[INDEXUV]
            Y=Y[INDEXUV]
            
            Xbaseline, Ybaseline = baseline_Rubberband(X, Y)
            baseline=np.interp(X, Xbaseline, Ybaseline)
            Y=Y-baseline
            # #Y=baseline
            # Y = Y/valeurnorm[i];
            RAJOUT = '_subRubber_normep'

            Y = normYminmax(Xnm[INDEXUV], Y, COUPURENORMminmax)  
            #RAJOUT = '_subRubber'
            
            plt.ylabel('Absorbance normalisé (u.a)')# à l''épaisseur soustrait ligne de base Rubberband (u.a)')
        
        elif Modeaff == 'Gradient':
            INDEXUV=Xnm>350
            X=X[INDEXUV]
            Y=Y[INDEXUV]
            Y=np.gradient(Y, X)
            #Y=np.gradient(Y, X)
            plt.ylabel('dA/dx (u.a)')
            RAJOUT = '_gradient'
        
        else:
            X = Xnm;
            Y = Ytr;
            plt.xlabel("Longueur d'onde ($nm$)");
            plt.ylabel('Transmittance') 
            RAJOUT = '_Tr'
        
        if modecouleurs == 'auto':
            plt.plot(X, Y, linewidth=linewidth, label=Legende[i])
        
        elif modecouleurs == 'bigdata':
            plt.plot(X, Y, color=colorsbigdata[i], linewidth=linewidth, label=Legende[i])
            
        elif modecouleurs == 'manuel':
            if  optionplot[i] == '':
                plt.plot(X, Y, linewidth=linewidth, label=Legende[i])
            else:
                plt.plot(X, Y, optionplot[i], linewidth=linewidth, label=Legende[i])
            
    #plt.legend(Legende, loc="upper left");
    #plt.legend(Legende, bbox_to_anchor = [0.5, 0.2])
    plt.legend()
      
    TITRE=TITRE+RAJOUT
    
    
    ax1=plt.gca()
    fig=plt.gcf()   
    ax1.xaxis.set_minor_locator(ticker.AutoMinorLocator())

    
    if SecondAxe and not (Modeaff == 'ABSnm' or Modeaff == 'Transmittance') : # Parti double échelle mettre false pour avoir uniquement en cm^-1
        def cm12um(X):
            return(nm2cm1(X)*1/1000)
        majortickman=[2.5, 2, 1.5, 1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.33, 0.3, 0.2]
        
        secax = ax1.secondary_xaxis('top', functions=(cm12um, cm12um))
        # Créaction d'un axes secondaire pour afficher les nm.
        secax.set_xlabel('Longueur d\'onde (µm)')
        secax.xaxis.set_major_locator(ticker.FixedLocator(majortickman))
        secax.xaxis.set_minor_locator(ticker.FixedLocator(np.arange(0.02, 1, 0.02)))

        # secax.xaxis.set_major_locator(ticker.AutoLocator())
        # secax.xaxis.set_minor_locator(ticker.AutoMinorLocator())

    if SpectreVIS :
        #Représentation du spectre Visible en échelle linéaire en cm-1.
        xlim=ax1.get_xlim()
        ylim=ax1.get_ylim()
        taillespectre=0.05*(ylim[1]-ylim[0]); # Calcul de la hauteur de l'arc en ciel 
        
        if xlim[0]>nm2cm1(850) or xlim[1]<nm2cm1(350) :   
            axins=ax1.inset_axes([xlim[0],ylim[1]-(taillespectre+0.02*taillespectre),
                                  xlim[1]-xlim[0],taillespectre],transform=ax1.transData)
            #Création d'un graph suplémentaire dans la figure principale.
            #Le spectre est ajusté aux valeurs limites de x et de y.
            clim=(xlim[0],xlim[1])
           
        
        else :    
            axins=ax1.inset_axes([nm2cm1(850),ylim[1]-(taillespectre+0.02*taillespectre),
                                  nm2cm1(350)-nm2cm1(850),taillespectre],transform=ax1.transData)              
            clim=(nm2cm1(850),nm2cm1(350))
            #Le spectre tiendra dans une boîte de 350 à 780nm mais linéaire en cm-1
        
        norm = plt.Normalize(*clim)
        wn = np.linspace(clim[0],clim[1],2000) #sélection des nombres d'onde
        wl=nm2cm1(wn)
        
        colorlist = list(zip(norm(wn),[wavelength_to_rgb(w) for w in wl]))
        #on associe pour chaque nombre d'onde la valeur en RGB
        spectralmap = matplotlib.colors.LinearSegmentedColormap.from_list("spectrum", colorlist)
        # conversion en colormap
        
        y = np.linspace(0, 6, 100)
        X,Y = np.meshgrid(wn, y) # création du maillage pour remplir l'image
        
        extent=(np.min(wn), np.max(wn), np.min(y), np.max(y))
        axins.imshow(X, interpolation='none', clim=clim, extent=extent, cmap=spectralmap,
                     aspect='auto') #remplissage l'axe avec les couleurs
        axins.tick_params(axis='both', which='both', bottom=False, left=False,
                          labelbottom=False, labelleft=False) #aucun tick
        axins.axis('off')
        axins.text(nm2cm1(790),2, 'IR',fontsize=tailletexte) 
        axins.text(nm2cm1(370),2, 'UV',fontsize=tailletexte)

#Agrandir les légendes
    # params = {'legend.fontsize': tailletexte,'legend.handlelength': 1, 'font.size': tailletexte}
    # plt.rcParams.update(params)
    
    # ax1.tick_params(which='major',length=5)
    # ax1.tick_params(which='minor',length=5)
    # secax.tick_params(which='major',length=5)
    # secax.tick_params(which='minor',length=5)
    
    if SHOW : Sav_fig(TITRE)
    
    return(fig)

def Affichage_gauss(Liste, Legende, TITRE, SHOW=True, optplt='--') :
    '''
    Lit le fichier parametre_fit_gauss.csv et trace les gaussienne correspondante
    si le nom du fit match avec les strings dans la variable legende

    Parameters
    ----------
    Liste : TYPE
        DESCRIPTION.
    Legende : TYPE
        DESCRIPTION.
    TITRE : TYPE
        Nom de sauvegarde du fichier.
    SHOW : TYPE, optional
        True pour sauvegarder et afficher la legende. The default is True.
    optplt : TYPE, optional
        DESCRIPTION. The default is '--'.

    Returns
    -------
    None.

    '''
    Data=np.genfromtxt('./parametre_fit_gauss.csv', comments='#', skip_header=1, delimiter=';', dtype='str' )
    
    X = np.linspace(25000, 45000, 200)
    
    for i in np.arange(0, np.size(Data[:, 0])):
        if Data[i, 0] in Legende :
            Y = Gauss(X, float(Data[i,1]), float(Data[i,2]), float(Data[i,3]), float(Data[i,4]))
            plt.plot(X, Y, optplt, label='fit'+Data[i, 0])
        else:
            pass
    if SHOW : Sav_fig(TITRE)
    return()

def Calcul_XYZ(Liste, Legende):
    '''
    Parameters
    ----------
    Liste : chemin vers les fichier contenant les datas
        DESCRIPTION.
    Returns
    -------
    XYZ : tableau des XYZ correspondant à la liste d'entrée.
    '''
    XYZ=np.ones((np.size(Liste), 3))
    for i in np.arange(0, np.size(Liste)):
        Fichier=Liste[i]

        Xnm, Ytr = Readspectre(Fichier)
        INDEX=~np.isnan(Ytr) # supression des nan
        Ytr=Ytr[INDEX]
        Xnm=Xnm[INDEX]
        
        # Ytr=-np.log10(Ytr) #sibesoin de corriger le spectre
        # Ytr=np.power(10, -Ytr)
        
        nm=np.arange(int(np.min(Xnm))+1 , int(np.max(Xnm))-1) # Création des valeur en nm à interpoler (pas de 1, 5, 10 ou 20 nm)
        
        fsample = interpolate.interp1d(Xnm, Ytr) # On interpole le fichier pour avoir des valeur espacée de 1nm  (norme CIE)
        Tsample = np.nan_to_num(fsample(nm));
        
        Data_dict={}
        for compteur in range(np.size(nm)): Data_dict[nm[compteur]] = Tsample[compteur] # Mise des données sous forme de dictionnaire pour utilisation avec le module Colour

        sd = colour.SpectralDistribution(Data_dict, name=Legende[i]) # Créaction de l'objet densité spectral sample
        cmfs = colour.colorimetry.MSDS_CMFS_STANDARD_OBSERVER['CIE 1931 2 Degree Standard Observer'] # Réfence pour le calcul tu trimultus
        illuminant = colour.SDS_ILLUMINANTS['D65']
        
        
        #cplot.plot_single_sd(sd) # Affichage de la densité si besoin, possibilité d'en afficher plusieurs avec la bonne fonction 
    
        XYZ[i] = colour.sd_to_XYZ(sd, cmfs, illuminant) # Passage de la densité spectral au trismultus
    return(XYZ)
    
def AffichageCIE1931(Liste, Legende, TITRE='CIE1931', Marqueur='', Fleche=True,
                         xylim=[-0.1, 0.8, 0, 1], show=True):
    '''
    Cette fonction sert à déterminer les coordonnées colorimétrique des spectres et les affichers sur le diagrame xyY
    Il se base sur le module colour science

    Parameters
    ----------
    Liste : TYPE
        Chemin relatif ou absolu vers les fichiers.
    Legende : TYPE
        Nom à afficher sur le graph
    TITRE : TYPE, optional
        Titre du graph sauvegarder. The default is 'CIE1931'.
    Fleche : Bool, optional
        Pour afficher les noms directements sur le graphe ou avec une legende classique. The default is True.

    Returns
    -------
    None.

    '''
    FIGSIZE = (6, 3)
    DPI = 150
    
    xtable=[]
    ytable=[]

    
    set_graph(FIGSIZE=FIGSIZE, DPI=DPI, mode='CIE1931', xylim=xylim)

    if TITRE == 'CIE1931':
        pass
    else:
        TITRE = TITRE + '_CIE1931'
    
    Marqueur=mono2tab(Marqueur, np.size(Liste))    
    tristimulus=Calcul_XYZ(Liste, Legende)
    
    for i in np.arange(0, np.size(Liste)):
        XYZ=tristimulus[i]
        #print(XYZ)
        #RGB = colour.XYZ_to_sRGB(XYZ / 100) # Passage en RGB.
        x, y =  colour.XYZ_to_xy(XYZ) # passage en xyY pour le diagramme fer à Cheval CIE1931
        xtable.append(x)
        ytable.append(y)
        #cplot.plot_single_colour_swatch(cplot.ColourSwatch('Sample', RGB), text_parameters={'size': 'x-large'}) # Afficage de la couleur RGB
    
        if Marqueur[i]:
            #plt.plot(x, y, marker=Marqueur[i], color='white', label=Legende[i])
            plt.plot(x, y, marker=Marqueur[i], label=Legende[i])
            
        elif Fleche:
        #Plotting the *CIE xy* chromaticity coordinates.
            plt.plot(x, y, 'x', color='white')
        # Annotating the plot. 
            plt.annotate(Legende[i],
                            xy=(x, y),
                            xytext=(-50, 30),
                            textcoords='offset points',
                            arrowprops=dict(arrowstyle='->', connectionstyle='arc3, rad=-0.2'))
            
        else :
                plt.plot(x, y, 'x', label=Legende[i])
        
    plt.legend()

    if show: Sav_fig(TITRE, colorplot=True)
    
    return(xtable, ytable)


def Affichage_Lab2D(Liste, Legende, TITRE='Lab', Marqueur='', SHOW=True):
    '''
    Cette fonction sert à afficher des LAB sur deux graph 2D

    Parameters
    ----------
    Liste : TYPE
        DESCRIPTION.
    Legende : TYPE
        DESCRIPTION.
    TITRE : TYPE, optional
        DESCRIPTION. The default is 'Lab'.
    Marqueur : TYPE, optional
        DESCRIPTION. The default is ''.
    SHOW : TYPE, optional
        DESCRIPTION. The default is 'True'.

    Returns
    -------
    None.

    '''
    FIGSIZE=(5,5)
    DPI=120
    Marqueur=mono2tab(Marqueur, np.size(Liste))    

    tristimulus=Calcul_XYZ(Liste, Legende)
    
    if TITRE == 'Lab':
        pass
    else:
        TITRE = TITRE + '_Lab'


    plt.style.use({'figure.figsize': FIGSIZE, 'figure.dpi': DPI})
    fig, axs = plt.subplots(2, 1, sharex=True)
    fig.subplots_adjust(left=0.15, bottom=0.1, right=0.99, top=0.99, wspace=0, hspace=0)
    ax1 = axs[0]
    ax2 = axs[1]
    Lab_liste=[]
    
    for i in np.arange(0, np.size(Liste)):
        XYZ=tristimulus[i]
        
        Lab=colour.XYZ_to_Lab(XYZ);
        Lab_liste.append(Lab)
        
        if Marqueur[i]:
            ax1.plot(Lab[1], Lab[0],'', marker=Marqueur[i], label=Legende[i])
            ax2.plot(Lab[1], Lab[2], '', marker=Marqueur[i], label=Legende[i])
        else:
            ax1.plot(Lab[1], Lab[0], 'x', label=Legende[i])
            ax2.plot(Lab[1], Lab[2], 'x', label=Legende[i])
            
    ax1.set_ylabel('L')
    ax1.grid()
    ax1.legend()
    
    ax2.set_ylabel('b')
    ax2.set_xlabel('a')
    ax2.grid()
 
    #fig.tight_layout()
    if SHOW : Sav_fig(TITRE)
    return(Lab_liste)
        

def Affichage_Lab3D(Liste, Legende, TITRE='Lab3D', Marqueur='', SHOW=True):
    '''
    Cette fonction afficer les LAB sur un graph 3D des fichier dans LISTE.

    Parameters
    ----------
    Liste : TYPE
        DESCRIPTION.
    Legende : TYPE
        DESCRIPTION.
    TITRE : TYPE, optional
        DESCRIPTION. The default is 'Lab'.
    Marqueur : TYPE, optional
        DESCRIPTION. The default is ''.
    SHOW : TYPE, optional
        DESCRIPTION. The default is 'True'.

    Returns
    -------
    None.

    '''
    
    FIGSIZE=(5,5)
    DPI=120
    Marqueur=mono2tab(Marqueur, np.size(Liste))    

    tristimulus=Calcul_XYZ(Liste, Legende)
    Lab_liste=[]
    
    if TITRE == 'Lab3D':
        pass
    else:
        TITRE = TITRE + '_Lab3D'

    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
        
    for i in np.arange(0, np.size(Liste)):
        XYZ=tristimulus[i]
        
        Lab=colour.XYZ_to_Lab(XYZ);
        Lab_liste.append(Lab)
        
        L=Lab[0]
        a=Lab[1]
        b=Lab[2]
        
               
        if Marqueur[i]:
            ax.scatter(a, b, L, marker=Marqueur[i], label=Legende[i]) 
        else:
            ax.scatter(a, b, L, marker='o', label=Legende[i])

                
    ax.set_xlabel('a')
    ax.set_ylabel('b')
    ax.set_zlabel('L')
    # ax.set_xlim([-300, 300])
    # ax.set_ylim([-300, 300])
    # ax.set_zlim([-350,350])
 
    #fig.tight_layout()
    if SHOW : Sav_fig(TITRE)
    return(Lab_liste)

def Calcul_BeerLambert_XYZ(Fichier, coeflim=[0.1, 5]):
    '''
    Cette fonction calcule l'évolution des paramètres XYZ en faisant varier le paramètre
    lc dans la loi de Beer-Lambert pour un spectre donné.
    

    Parameters
    ----------
    Fichier : string
        Chemin vers le fichier contenant les données du spectre.
    coeflim : liste de float, optional
        Valeurs minimale et maximale des paramètres lc. The default is [0.1, 5].

    Returns
    -------
    Paramètres XYZ, Légende, XYZ du spectre donné.

    '''
    Xnm, T =Readspectre(Fichier)

    INDEX=~np.isnan(T) # supression des nan
    T=T[INDEX]
    Xnm=Xnm[INDEX]
   
    coef = np.arange(coeflim[0], coeflim[1], 0.1)
    XYZ=np.ones((np.size(coef), 3))

    i=0
    Legende=[]
    
    for coef in coef:
        # print(coef)
        Legende.append('x'+str(coef))
    
        Ytr=np.power(T, coef)
    
        nm=np.arange(int(np.min(Xnm))+1 , int(np.max(Xnm))-1) # Création des valeur en nm à interpoler (pas de 1, 5, 10 ou 20 nm)
        
        fsample = interpolate.interp1d(Xnm, Ytr) # On interpole le fichier pour avoir des valeur espacée de 1nm  (norme CIE)
        Tsample = np.nan_to_num(fsample(nm));
        
        Data_dict={}
        for compteur in range(np.size(nm)): Data_dict[nm[compteur]] = Tsample[compteur] # Mise des données sous forme de dictionnaire pour utilisation avec le module Colour

        sd = colour.SpectralDistribution(Data_dict, name=Legende[i]) # Créaction de l'objet densité spectral sample
        cmfs = colour.colorimetry.MSDS_CMFS_STANDARD_OBSERVER['CIE 1931 2 Degree Standard Observer'] # Réfence pour le calcul tu trimultus
        illuminant = colour.SDS_ILLUMINANTS['D65']
        
        XYZ[i] = colour.sd_to_XYZ(sd, cmfs, illuminant)
        
        if coef==1: ref_XYZ=XYZ[i]

        i+=1
        
    return(XYZ, Legende, ref_XYZ)
    
    
def Courbe_BeerLambert_XY(Fichier, Legende, TITRE='CIE1931', fondCIE=False,
                          xylim=[-0.1, 0.8, -0.1, 1], lc_lim=[0.1, 5], show=True):
    '''
    Cette fonction trace la courbe de Beer-Lambert en coordonnées colorimétriques xy 
    (i.e. évolution des coordonnées x et y avec le paramètre lc) à partir d'un spectre donné.

    Parameters
    ----------
    Fichier : string
        Chemin vers le fichier contenant les données du spectre.
    TITRE : string, optional
        Titre du graphique sauvegardé. The default is 'CIE1931'.
    fondCIE : Bool, optional
        Si vrai, affiche le fond CIE1931. The default is False.
    xylim : liste, optional
        Bornes du fond CIE1931. The default is [-0.1, 0.8, -0.1, 1].
    lc_lim : liste, optional
        Valeurs minimale et maximale du paramètre lc. The default is [0.1, 5].
    show : Bool, optional
        Si vrai, sauvegarde et affiche le graphique. The default is True.

    Returns
    -------
    L'ensemble des coordonnées x et y.

    '''
    xtable=[]
    ytable=[]
    FIGSIZE = (6, 3)
    DPI = 150
    
    (cheminfichier, nomfichier) = os.path.split(Fichier)
    Verre = Legende
    
    if TITRE == 'CIE1931':
        pass
    else:
        TITRE = TITRE + '_CIE1931'

    XYZ, Legende ,ref_XYZ = Calcul_BeerLambert_XYZ(Fichier, coeflim=lc_lim)


    if fondCIE:
        set_graph(FIGSIZE=FIGSIZE, DPI=DPI, mode='CIE1931', xylim=xylim)
        ref_x, ref_y = colour.XYZ_to_xy(ref_XYZ)    
        plt.plot(ref_x,ref_y,'x',label=Verre)
    
    else:
        pass
    
    for i in np.arange(0,np.size(XYZ, 0)):                
        x, y =  colour.XYZ_to_xy(XYZ[i]) # passage en xyY pour le diagramme fer à Cheval CIE1931
        xtable.append(x)
        ytable.append(y)
        

    plt.plot(xtable, ytable, label='Effet lc ' + Verre)

    #plt.legend(loc='upper right')
    plt.legend(bbox_to_anchor = [1, 1])

    if show: Sav_fig(TITRE, colorplot=True)
       
    return(xtable, ytable)
    
    
def Courbe_BeerLambert_Lab3D(Fichier, Legende='', TITRE='Lab', lc_lim=[0.1, 5], show=True, newfig=False):
    '''
    Cette fonction trace la courbe 3D de Beer-Lambert en coordonnées colorimétriques Lab
    (i.e. évolution des coordonnées Lab avec le paramètre lc) à partir d'un spectre donné.

    Parameters
    ----------
    Fichier : string
        Chemin vers le fichier contenant les données du spectre.
    TITRE : string, optional
        Titre du graphique sauvegardé. The default is 'Lab'.
    lc_lim : liste, optional
        Valeurs minimale et maximale du paramètre lc. The default is [0.1, 5].
    show : Bool, optional
        Si vrai, sauvegarde et affiche le graphique. The default is True.

    Returns
    -------
    L'ensemble des valeurs des coordonnées Lab.

    '''
    
    (cheminfichier, nomfichier) = os.path.split(Fichier)
    if Legende == '':
        Verre = nomfichier[:-19]
    else :
        Verre = Legende
    
    XYZ,Legende,ref_XYZ=Calcul_BeerLambert_XYZ(Fichier, coeflim=lc_lim)
    
    if newfig:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        Ref_L, ref_a, ref_b = colour.XYZ_to_Lab(ref_XYZ)

        ax.scatter(ref_a, ref_b, Ref_L, marker='x', label=Verre)
        ax.set_xlabel('a')
        ax.set_ylabel('b')
        ax.set_zlabel('L')
        # ax.set_xlim([-350, 350])
        # ax.set_ylim([-350, 350])
        # ax.set_zlim([-350, 350])
        
    L=[]
    a=[]
    b=[]
    
    if TITRE == 'Lab3D':
        pass
    else:
        TITRE = TITRE + '_Lab3D'

    for i in np.arange(0,np.size(XYZ, 0)):                
        L_temp, a_temp, b_temp =  colour.XYZ_to_Lab(XYZ[i]) # passage en xyY pour le diagramme fer à Cheval CIE1931
        L.append(L_temp)
        a.append(a_temp)
        b.append(b_temp) 
    
    plt.plot(a,b, L, label='Effet lc ' + Verre)

    plt.legend()
 

    if show : Sav_fig(TITRE)
    return(L,a,b)  
  

def wavelength_to_rgb(wavelength, gamma=0.8):
    ''' taken from http://www.noah.org/wiki/Wavelength_to_RGB_in_Python
    This converts a given wavelength of light to an 
    approximate RGB color value. The wavelength must be given
    in nanometers in the range from 380 nm through 750 nm
    (789 THz through 400 THz).

    Based on code by Dan Bruton
    http://www.physics.sfasu.edu/astro/color/spectra.html
    Additionally alpha value set to 0.5 outside range
    '''
    wavelength = float(wavelength)
    if wavelength >= 400 and wavelength <= 750:
        A = 1.
    elif wavelength > 750 and wavelength <= 850:
        A=0.01*(850-wavelength)
    elif wavelength < 400 and wavelength >= 350:
        A=0.02*(wavelength-350)
    else:
        A=0.0
    if wavelength < 380:
        wavelength = 380.
    if wavelength >750:
        wavelength = 750.
    if wavelength >= 380 and wavelength <= 440:
        attenuation = 0.3 + 0.7 * (wavelength - 380) / (440 - 380)
        R = ((-(wavelength - 440) / (440 - 380)) * attenuation) ** gamma
        G = 0.0
        B = (1.0 * attenuation) ** gamma
    elif wavelength >= 440 and wavelength <= 490:
        R = 0.0
        G = ((wavelength - 440) / (490 - 440)) ** gamma
        B = 1.0
    elif wavelength >= 490 and wavelength <= 510:
        R = 0.0
        G = 1.0
        B = (-(wavelength - 510) / (510 - 490)) ** gamma
    elif wavelength >= 510 and wavelength <= 580:
        R = ((wavelength - 510) / (580 - 510)) ** gamma
        G = 1.0
        B = 0.0
    elif wavelength >= 580 and wavelength <= 645:
        R = 1.0
        G = (-(wavelength - 645) / (645 - 580)) ** gamma
        B = 0.0
    elif wavelength >= 645 and wavelength <= 750:
        attenuation = 0.3 + 0.7 * (750 - wavelength) / (750 - 645)
        R = (1.0 * attenuation) ** gamma
        G = 0.0
        B = 0.0
    else:
        R = 0.0
        G = 0.0
        B = 0.0
    return (R,G,B,A)

def Affichage_ID20(Liste, Legende, Autoaxe=True, Xlim=[0, 800], Ylim=[0, 1],
                  TITRE='superposition', linewidth=1.5, Modeaff='default',
                  modecouleurs='auto', optionplot='', SHOW=True, GRAPHOUT=1,
                  newgraph=True):
 
    
    Legende=mono2tab(Legende, np.size(Liste))
    optionplot=mono2tab(optionplot, np.size(Liste))

    
    if newgraph : set_graph(GRAPHOUT=GRAPHOUT) # Cf fonction, on créer la figure à la bonne taille/résolution


    colorsbigdata = plt.cm.nipy_spectral(np.linspace(0,1,np.size(Liste))) # On fixe la gamme de couleur utilisé

    if Autoaxe:
        pass
    else:
        plt.xlim(Xlim);
        plt.ylim(Ylim);
    

    for i in np.arange(0, np.size(Liste)): # on lit tous les fichier  
        Fichier = Liste[i];
        
        X, Y = Readspectre(Fichier, skip_header=0, delimiter=' ')

        plt.xlabel("Perte d'énergie (en eV)");
        plt.ylabel('Intensité (U.A)')
            
        if Modeaff == 'default':
            RAJOUT= ''

        elif Modeaff == 'Gradient':
            Y=np.gradient(Y, X)
            #Y=np.gradient(Y, X)
            plt.ylabel('dA/dx (u.a)')
            RAJOUT = '_gradient'
        

        if modecouleurs == 'auto':
            plt.plot(X, Y, linewidth=linewidth, label=Legende[i])
        
        elif modecouleurs == 'bigdata':
            plt.plot(X, Y, color=colorsbigdata[i], linewidth=linewidth, label=Legende[i])
            
        elif modecouleurs == 'manuel':
            if  optionplot[i] == '':
                plt.plot(X, Y, linewidth=linewidth, label=Legende[i])
            else:
                plt.plot(X, Y, optionplot[i], linewidth=linewidth, label=Legende[i])
            
    #plt.legend(Legende, loc="upper left");
    #plt.legend(Legende, bbox_to_anchor = [0.5, 0.2])
    plt.legend()
      
    TITRE=TITRE+RAJOUT
    
    ax1=plt.gca()
    fig=plt.gcf()   

    if SHOW : Sav_fig(TITRE)
    
    return(fig)       
        
        
