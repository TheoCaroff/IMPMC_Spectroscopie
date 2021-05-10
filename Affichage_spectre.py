# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 18:51:51 2020

@author: Theo_C
Ce fichier contient toutes les chose utiles pour tracer et traiter les spectres
"""

import numpy as np
import os
from matplotlib import pyplot as plt
from scipy import interpolate


from Lecture_input import mono2tab
from Lecture_input import Readspectre

try :
    import colour
    from colour import plotting as cplot
except ModuleNotFoundError :
    print('ATTENTION MODULE COLOR-SCIENCE PAS INSTALLER')


def set_graph(FIGSIZE=[12, 6], DPI = 120, grid=True):
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
    plt.style.use({'figure.figsize': FIGSIZE, 'figure.dpi': DPI})
#   plt.figure(figsize=(FIGSIZE), dpi=DPI)
    plt.figure()
    plt.grid(grid);
    plt.subplots_adjust(left=0.05, bottom=0.1, right=0.95, top=0.9, wspace=0.1, hspace=0.1)


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
    
    plt.savefig(Repertoire+os.sep+Titre, bbox_inches='tight')

    if colorplot:
        cplot.render(standalone=True)
    else:
        plt.show()
        
def Affichage_abs(Liste, Legende, Autoaxe=True, Xlim=[4000, 35000], Ylim=[0, 1.5],
                  SecondAxe=True, TITRE='superposition', AdditionTr=0,
                  linewidth=1.5, valeurnorm=1, Modeaff='ABScm', modecouleurs='auto', optionplot='',
                  SHOW=True):
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
    Returns
    -------
    None.

    '''
    
    #Liste=mono2tab(Liste, np.size(Liste))
    Legende=mono2tab(Legende, np.size(Liste))
    optionplot=mono2tab(optionplot, np.size(Liste))
    AdditionTr=mono2tab(AdditionTr, np.size(Liste))
    valeurnorm=mono2tab(valeurnorm, np.size(Liste))
    
    set_graph() # Cf fonction, on créer la figure à la bonne taille/résolution

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
        
        if Modeaff == 'ABScm':
            X = 1/(Xnm*1E-7);
            Y = -np.log10(Ytr);
            plt.xlabel("Nombre d'onde ($cm^{-1}$)");
            plt.ylabel('Absorbance')
            RAJOUT= ''
        
        elif Modeaff == 'ABSnm':
            X = Xnm;
            Y = -np.log10(Ytr);
            plt.xlabel("Longueur d'onde (en nm)");
            plt.ylabel('Absorbance')
            RAJOUT= '_ABSnm'
        
        
        elif Modeaff == 'Reflectance':
            X = 1/(Xnm*1E-7);
            Y = 1-Ytr;
            plt.xlabel("Nombre d'onde ($cm^{-1}$)");
            plt.ylabel('Reflectance')
            RAJOUT = '_Reflexion'
       
        elif Modeaff == 'ABSnormep':
            X = 1/(Xnm*1E-7);
            Y = -np.log10(Ytr);
            Y = Y/valeurnorm[i];
            plt.xlabel("Nombre d'onde ($cm^{-1}$)");
            plt.ylabel('Absorbance normalisée ($cm^{-1}$)') 
            RAJOUT = 'norm_ep'
        
        
        elif Modeaff == 'Epsilon':
            X = 1/(Xnm*1E-7);
            Y = -np.log10(Ytr)/valeurnorm[i];
            plt.xlabel("Nombre d'onde ($cm^{-1}$)");
            plt.ylabel('$\\varepsilon (L.mol^{-1}.cm^{-1})$') 
            RAJOUT = '_Epsilon'
        
        else:
            X = Xnm;
            Y = Ytr;
            plt.xlabel("Longueur d'onde ($nm$)");
            plt.ylabel('Transmittance') 
            RAJOUT = '_Tr'
        
        if modecouleurs == 'auto':
            plt.plot(X, Y, linewidth=linewidth)
        
        elif modecouleurs == 'bigdata':
            plt.plot(X, Y, color=colorsbigdata[i], linewidth=linewidth)
            
        elif modecouleurs == 'manuel':
            if  optionplot[i] == '':
                plt.plot(X, Y, linewidth=linewidth)
            else:
                plt.plot(X, Y, optionplot[i], linewidth=linewidth)
            
    #plt.legend(Legende, loc="upper left");
    #plt.legend(Legende, bbox_to_anchor = [0.5, 0.2])
    plt.legend(Legende);
    
    TITRE=TITRE+RAJOUT
    
    
    ax1=plt.gca()
    fig=plt.gcf()   
        
    if SecondAxe and (Modeaff == 'ABScm') : # Parti double échelle mettre false pour avoir uniquement en cm^-1
        fig.canvas.draw()
        ax2 = ax1.twiny()
        axmin, axmax = ax1.get_xlim()
        ax2.set_xlim(axmin, axmax)
    
        # Calculate Major Ticks
        ax2_labels = []
        ax2_labels.append(float('inf')) #créaction de la première tique manuelle sinon ça merde, division par 0
        for item in ax1.get_xticklabels()[1:]:
        #for item in ax1.get_xticklabels():
            l = 1/float(item.get_text())*1E7
            l = "{:3.0f}".format(l)
            ax2_labels.append(l)
        ax2.set_xticklabels(ax2_labels)
        #ax2.set_xlabel('Wavelenght (in nm)');
        ax2.set_xlabel('Longueur d\'onde (nm)');
    
    if SHOW : Sav_fig(TITRE)

    
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
                         xylim=[-0.1, 0.8, 0, 1], show='True'):
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
    
    plt.style.use({'figure.figsize': FIGSIZE, 'figure.dpi': DPI})
    xtable=[]
    ytable=[]

    
    if TITRE == 'CIE1931':
        pass
    else:
        TITRE = TITRE + '_CIE1931'
    
    Marqueur=mono2tab(Marqueur, np.size(Liste))    

    # Plotting the *CIE 1931 Chromaticity Diagram*.
    # The argument *standalone=False* is passed so that the plot doesn't get
    # displayed and can be used as a basis for other plots.
    cplot.plot_chromaticity_diagram_CIE1931(standalone=False, bounding_box=xylim, x_tighten=True,
        y_tighten=True)
    
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
            plt.plot(x, y, marker=Marqueur[i], color='white', label=Legende[i])
        elif Fleche:
        #Plotting the *CIE xy* chromaticity coordinates.
            plt.plot(x, y, 'x-', color='white')
        # Annotating the plot. 
            plt.annotate(Legende[i],
                            xy=(x, y),
                            xytext=(-50, 30),
                            textcoords='offset points',
                            arrowprops=dict(arrowstyle='->', connectionstyle='arc3, rad=-0.2'))
            
        else :
                plt.plot(x, y, 'x-', label=Legende[i])
        
    plt.legend()

    if show: Sav_fig(TITRE, colorplot=True)
    
    return(xtable, ytable)



def plt_BeerLambert_xy(Fichier, Legende, optionplot='', show='True', TITRE='CIE1931'):
    '''
    

    Parameters
    ----------
    Fichier : string
        Chemin vers le data à afficher.
    Legende : TYPE
        DESCRIPTION.
    optionplot : TYPE, optional
        DESCRIPTION. The default is ''.
    show : TYPE, optional
        DESCRIPTION. The default is 'True'.
    TITRE : TYPE, optional
        DESCRIPTION. The default is 'CIE1931'.

    Returns
    -------
    None.

    '''
    
    if TITRE == 'CIE1931':
        pass
    else:
        TITRE = TITRE + '_CIE1931'

    x, y = Readspectre(Fichier)
    
    if optionplot == '':
        plt.plot(x,y, '--', label=Legende)
    else:
        plt.plot(x,y, optionplot, label=Legende)
    
    plt.legend()
    
    if show : Sav_fig(TITRE, colorplot=True);

def Affichage_Lab(Liste, Legende, TITRE='Lab', Marqueur='', SHOW='True'):
    
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
        
        
        
        
        
        