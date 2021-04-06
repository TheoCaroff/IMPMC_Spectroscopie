# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 18:51:51 2020

@author: Theo_C
Ce fichier contient toutes les chose utiles pour tracer et traiter les spectres
"""

import numpy as np
import os
from matplotlib import pyplot as plt
import colour
from scipy import interpolate
from colour import plotting as cplot
import math

def Affichage_abs(Liste, Legende, Autoaxe=True, nbonde_min=4000, nbonde_max=35000, Amin=0,
                  Amax=1.5, SecondAxe=True, TITRE='superposition', AdditionTr=0,
                  linewidth=1.5, ABScm=True, modecouleurs='auto', optionplot=''):
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
    nbonde_min : float
        borne basse de l'affichage des X
    nbonde_max : float
        borne haute de l'affichage des X
    Amin : float
        borne basse de l'affichage des Y.
    Amax : float
         borne haute de l'affichage des Y..
    SecondAxe : bool
        True si deuxième axe en nm.
    TITRE : TYPE, optional
        DESCRIPTION. The default is 'superposition'.
    AdditionTr: float
        Correction de transmitance
    ABScm : bool
        Choisir de plot en absorbance/nombre d'onde ou transmittance/nm
    linewidth : TYPE, optional
        DESCRIPTION. The default is 1.5.
    modecouleurs : string, optional
        Permet de choisir le mode d'affichage des couleurs (trois mode existe :
            'auto', bigdata', 'manuel'. The default is 'auto'.
    optionplot : string, optional
        Choix des options des matplolib à utiliser en mode manuel

    Returns
    -------
    None.

    '''

    plt.figure(figsize=([12,6]), dpi=120)
    #plt.figure(figsize=(10,6), dpi=140)
    plt.subplots_adjust(left=0.05, bottom=0.1, right=0.9, top=0.9, wspace=0.1, hspace=0.1)

    colorsbigdata = plt.cm.nipy_spectral(np.linspace(0,1,np.size(Liste))) # On fixe la gamme de couleur utilisé

    plt.grid(True);
    #plt.xlabel("Wavenumber ($cm^{-1}$)");
    #plt.ylabel('Absorbance [log(%T)]');

    #plt.ylabel('Absorbance [log(%T)]')

    if (np.size(optionplot) == 1): # Si pas d'argument un marqueur unique pour toute les données
        optionplot_ref = optionplot
        optionplot=[]
        for i in np.arange(0, np.size(Liste)):
            optionplot.append(optionplot_ref)

    if Autoaxe:
        pass
    else:
        plt.xlim(nbonde_min, nbonde_max);
        plt.ylim(Amin, Amax);
    
    if not ABScm : TITRE = TITRE + 'Tr'
    
    i=0;
    for Fichier in Liste: # on lit tous les fichier  
        try:
            Data = np.genfromtxt(Fichier, skip_header=2, delimiter=';'); # On récupère l'intensité de la référence
        except UnicodeDecodeError:
            Data = np.genfromtxt(Fichier, skip_header=2, delimiter=';', encoding='latin-1'); # On récupère l'intensité de la référence
         
        if ABScm:
            X = 1/(Data[:, 0]*1E-7);
            Y = -np.log10(Data[:, 1]+AdditionTr);
            plt.xlabel("Nombre d'onde ($cm^{-1}$)");
            plt.ylabel('Absorbance') 
        else:
            X = Data[:, 0];
            Y = Data[:, 1];
            plt.xlabel("Longueur d'onde ($nm$)");
            plt.ylabel('Transmitance (%)') 
        
        if modecouleurs == 'auto':
            plt.plot(X, Y, linewidth=linewidth)
        
        elif modecouleurs == 'bigdata':
            plt.plot(X, Y, color=colorsbigdata[i], linewidth=linewidth)
            
        elif modecouleurs == 'manuel':
            if  optionplot[i] == '':
                plt.plot(X, Y, linewidth=linewidth)
            else:
                plt.plot(X, Y, optionplot[i], linewidth=linewidth)
        
        i=i+1;
            
    #plt.legend(Legende, loc="upper left");
    #plt.legend(Legende, bbox_to_anchor = [0.5, 0.2])
    plt.legend(Legende);
    
    ax1=plt.gca()
    fig=plt.gcf()   
        
    if SecondAxe and ABScm: # Parti double échelle mettre false pour avoir uniquement en cm^-1
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
    
   
    try:
        os.mkdir('Graph')
    except OSError:
        pass
    plt.savefig('Graph'+os.sep+TITRE, bbox_inches='tight')
    plt.show()

def Affichage_epsi(nbonde_min, nbonde_max, epsimin, epsimax, SecondAxe, Liste,
                   Legende, I100, concentration, Epaisseur, Element):
    '''

      nbonde_min : float
        borne basse de l'affichage des X
    nbonde_max : float
        borne haute de l'affichage des X
    Amin : float
        borne basse de l'affichage des Y.
    Amax : float
         borne haute de l'affichage des Y.
    SecondAxe : bool
        True si deuxième axe en nm.
    Liste : Tableau de string
        Contient le nom des fichier à afficher.
    Legende : Tableau de string
        Legende des fichiers.
    I100 : float
        signal blanc de référence.
    concentration : tableau de float
        concentration de l'élément pour appliquer Beer-Lambert.
    Epaisseur : tableau de float
        Pour appliquer Beer-Lambert.
    Element : string
        Nom de l'élement étudié.

    Returns
    -------
    None.

    '''
    
    plt.figure(figsize=(5,3), dpi=120)
    #plt.figure(figsize=(10,6), dpi=140)
    plt.subplots_adjust(left=0.05, bottom=0.1, right=0.9, top=0.9, wspace=0.1, hspace=0.1)
    
    i=0;

    plt.grid(True);
    plt.xlabel("Nombre d'onde ($cm^{-1}$)");
    plt.ylabel('$\\varepsilon$ ($L.mol^{-1}.cm^{-1})$');
    plt.title(Element+'\n\n')
    plt.xlim(nbonde_min, nbonde_max);
    plt.ylim(epsimin, epsimax);
    
    for Fichier in Liste: # on lit tous les fichier       
        try:
            Data = np.genfromtxt(Fichier, skip_header=2, delimiter=';'); # On récupère l'intensité de la référence
        except UnicodeDecodeError:
            Data = np.genfromtxt(Fichier, skip_header=2, delimiter=';', encoding='latin-1'); # On récupère l'intensité de la référence
         
        X = 1/(Data[:, 0]*1E-7); # 
        #Y  = (Data[:, 1]); #Si data en ABS
        Y = -np.log10(Data[:, 1]); #Si data en %Tx
        
        Y=Y-I100;
        Y=Y/(concentration[i]*Epaisseur[i])
        
        #plt.plot(X, Y, color='tab:orange'); # On trace le nombre d'onde versus le coef d'extenction
        plt.plot(X, Y); # On trace le nombre d'onde versus le coef d'extenction
        i=i+1;
            
    plt.legend(Legende);
    ax1=plt.gca()
    fig=plt.gcf()
    
        
    if SecondAxe: # Parti double échelle mettre false pour avoir uniquement en cm^-1
    
        fig.canvas.draw()
        ax2 = ax1.twiny()
        axmin, axmax = ax1.get_xlim()
        ax2.set_xlim(axmin, axmax)
    
        # Calculate Major Ticks
        ax2_labels = []
        #♦ax2_labels.append(float('inf')) #créaction de la première tique manuelle sinon ça merde, division par 0
        ax2_labels.append(1/nbonde_min*1E7)
        for item in ax1.get_xticklabels()[1:]:
            l = 1/float(item.get_text())*1E7
            l = "{:3.0f}".format(l)
            ax2_labels.append(l)
        ax2.set_xticklabels(ax2_labels)
        ax2.set_xlabel('Longueur d\'onde (nm)');
     
    try:
        os.mkdir('Graph')
    except OSError:
        pass
    plt.savefig('Graph'+os.sep+'epsilon' + Element + '.png', bbox_inches='tight')
    
    plt.show()
    
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
        try:
             Data = np.genfromtxt(Fichier, skip_header=2, delimiter=';');
        except UnicodeDecodeError:
            Data = np.genfromtxt(Fichier, skip_header=2, delimiter=';', encoding='latin-1'); 
        
        nm=np.arange(int(Data[0,0])+1,int(Data[-1,0])-1) # Création des valeur en nm à interpoler (pas de 1, 5, 10 ou 20 nm)

        if np.size(nm) == 0:
            nm=np.arange(int(Data[-1,0])+1,int(Data[0,0])-1)

        fsample = interpolate.interp1d(Data[:, 0], Data[:, 1]) # On interpole le fichier pour avoir des valeur espacée de 1nm  (norme CIE)
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
    plt.style.use({'figure.figsize': (6, 3), 'figure.dpi': 150})
    xtable=[]
    ytable=[]
    try:
        os.mkdir('Graph')
    except OSError:
        pass
    
    if TITRE == 'CIE1931':
        pass
    else:
        TITRE = TITRE + '_CIE1931'
    

    if (np.size(Marqueur) == 1): # Si pas d'argument un marqueur unique pour toute les données
        refmark=Marqueur
        Marqueur=[]
        for i in np.arange(0, np.size(Liste)):
            Marqueur.append(refmark)

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
            
        # elif not RGBValueerror:
        #         plt.plot(x, y, 'x-', color=1-RGB, label=Legende[i])
        else :
                plt.plot(x, y, 'x-', label=Legende[i])
        
    plt.legend()
    
    plt.savefig('Graph'+os.sep+TITRE, bbox_inches='tight')

    # Displaying the plot.
    if show:
        cplot.render(standalone=True);

    
    return(xtable, ytable)



def plt_BeerLambert_xy(Fichier, Legende, optionplot='', show='True'):

    # if (np.size(optionplot) == 1): # Si pas d'argument un marqueur unique pour toute les données
    #     optionplot_ref = optionplot
    #     optionplot=[]
    #     for i in np.arange(0, np.size(Liste)):
    #         optionplot.append(optionplot_ref)


    try:
        Data = np.genfromtxt(Fichier, skip_header=2, delimiter=';');
    except UnicodeDecodeError:
        Data = np.genfromtxt(Fichier, skip_header=2, delimiter=';', encoding='latin-1');
    x=Data[:, 0]
    y=Data[:, 1]
    
    if optionplot == '':
        plt.plot(x,y, '--', label=Legende)
    else:
        plt.plot(x,y, optionplot, label=Legende)
    plt.legend()
    
    if show : cplot.render(standalone=True);


