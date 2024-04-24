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
import pandas as pd

import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
from matplotlib.widgets import RectangleSelector
from tkinter.simpledialog import askstring
from matplotlib.patches import FancyArrowPatch

from matplotlib.widgets import CheckButtons

from matplotlib import pyplot as plt
from matplotlib import ticker as ticker
from scipy import interpolate
from scipy.signal import savgol_filter

from Lecture_input import mono2tab
from Lecture_input import Readspectre
from Lecture_input import normYminmax
from Lecture_input import nm2cm1
from Lecture_input import eV2nm
from Lecture_input import Gauss
from Lecture_input import saveDATACIE_xy
from Lecture_input import saveDATACIE_LAB
from Lecture_input import saveDATACIE_XYZ
from Lecture_input import saveDATACIE_XYZ
from Lecture_input import Sav_fig
from Lecture_input import removeInfNan
from Lecture_input import  saveDataFrame
from Lecture_input import integrate2curve
from Nettoyage_spectre import baseline_Rubberband

try :
    import colour
    from colour import plotting as cplot
except ModuleNotFoundError :
    print('ATTENTION MODULE COLOR-SCIENCE PAS INSTALLER')


#%% General

def set_graph(GRAPHOUT=None, FIGSIZE=[12, 10], DPI = 120, grid=True, mode='default', xylim=[-0.1, 0.8, 0, 1],
              num = None):
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

    if GRAPHOUT:
        FIGSIZE, DPI = GRAPHOUT_mode(GRAPHOUT)
    
    plt.style.use({'figure.figsize': FIGSIZE, 'figure.dpi': DPI})
#   plt.figure(figsize=(FIGSIZE), dpi=DPI)
    if mode == 'default':
        fig = plt.figure(num=num)
        plt.grid(grid)
        ax = plt.gca()
        
    elif mode == 'CIE1931':
        # Plotting the *CIE 1931 Chromaticity Diagram*.
        # The argument *standalone=False* is passed so that the plot doesn't get
        # displayed and can be used as a basis for other plots.
        if xylim == 'bleu':
            xylim = [0.05, 0.25, 0, 0.20]
            
        elif xylim == 'pourpre':
            xylim = [0.30, 0.50, 0.26, 0.42]
            
        cplot.plot_chromaticity_diagram_CIE1931(standalone=False,
                                                bounding_box=xylim, x_tighten=True, y_tighten=True)
        fig = plt.gcf()
        ax = plt.gca()
    
    #plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
    
    return(fig, ax)


def GRAPHOUT_mode(number):
    if number == 1:
        FIGSIZE=[15, 10]
        DPI = 120

    elif number == 2 :
        FIGSIZE=[10, 5]
        DPI = 100

    elif number == 3 :
        FIGSIZE=[7, 3.5]
        DPI = 100

    elif number == 4 :
        FIGSIZE=[10, 3.5]
        DPI = 200
        
    elif number == 5 :
        # FIGSIZE=[7, 4.5]
        FIGSIZE=[5, 3.5]
        DPI = 200

    elif number == 6 :
        FIGSIZE=[7, 3.5]
        DPI = 200

    elif number == 7 :
        FIGSIZE=[6, 3.5]
        DPI = 200
    
    elif number == 8 :
        FIGSIZE=[9, 21/4]
        DPI = 200
    
    elif number == 9 :
        FIGSIZE=[7, 2.5]
        DPI = 200
        
    elif number == 11:
        FIGSIZE=[13, 10]
        DPI = 120

    return(FIGSIZE, DPI)



#%% Détaille affichage spectre

def zones_IR_UV(ax) :
    '''
    Dessine des rectangles dans les zones IR et UV du spectre tracé. 

    Parameters
    ----------
    ax : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    L = ax.get_xlim()
    H = ax.get_ylim()
    X_IR = nm2cm1(780)
    X_UV = nm2cm1(380)
    
    # Define custom colormap
    custom_cmap = mcolors.LinearSegmentedColormap.from_list('custom', ['#F8F8F8', '#D3D3D3']) #nuances de gris du plus clair au plus foncé

    # Create the gradient array
    gradient = np.linspace(0, 1, 200)

    # Plotting the gradient
    ax.imshow([gradient], aspect='auto', cmap=custom_cmap, extent=[X_IR,L[0], H[0], H[1]]) 
    ax.imshow([gradient], aspect='auto', cmap=custom_cmap, extent=[X_UV, L[1], H[0], H[1]])

    # Add text to the plot
    position_x_IR = nm2cm1(880)  # X-coordinates for the text
    position_x_UV = nm2cm1(360)
    position_y_txt = H[1]/2.5  # Y-coordinate for the text
    text_IR = 'Infrarouge'  # Texts to be displayed
    text_UV = 'Ultraviolet'
    ax.text(position_x_IR, position_y_txt, text_IR, fontsize=9, ha='center', va='center', fontstyle='italic', color = 'gray')
    ax.text(position_x_UV, position_y_txt, text_UV, fontsize=9, ha='center', va='center', fontstyle='italic', color = 'gray')
    
    ax.set_ylim(H[0], H[1])
    ax.set_xlim(L[0], L[1])
    
def Transm_window(X, Y, ax, lissage=True) :
    """
    Montre le centre de la fenêtre de transmission d'un spectre.
    
    Parameters
    ----------
    X : TYPE
        DESCRIPTION.
    Y : TYPE
        DESCRIPTION.
    ax : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    if lissage : Y = savgol_filter(Y, 31, 3)
        
    L = ax.get_xlim()
    H = ax.get_ylim()
    
    INDEX = np.logical_and(X > nm2cm1(780), X < nm2cm1(380))
    X=X[INDEX]
    Y=Y[INDEX]

    indice_Y_min = np.argmin(Y)
    
    X_min = X[indice_Y_min]
    
    ax.text(L[1]/2.05, H[1]*0.85, 'Centre de la fenêtre de transmission', backgroundcolor='white', fontsize=9, color = 'green')
    
    ax.axvline(x=X_min, linestyle='--', color='green')
    
    ax.set_ylim(H[0], H[1])
    ax.set_xlim(L[0], L[1])
    
    
def warning_light(ax):
    '''
    Dessine un rectangle dans la zone de faiblesse de la lampe
    (ou n'importe où ailleurs en modifiant les valeurs de x).
     Zone étendue définie entre 20490 et 23760cm-1.
     
    Parameters
    ----------
    ax : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    x_min = 22540
    x_max = 23300
  
    H = ax.get_ylim()

    ax.fill_betweenx(H, x_min, x_max, color='gray', alpha=0.2)
    ax.set_ylim(H[0], H[1])
    
def getylabel(ax) :
    return(ax.yaxis.get_label_text())

def SecondAxe_g(ax1, Langue='anglais', GRAPHOUT=None, unite = 'nm'):
    
    def mT2g(X):
        X=X*10
        X = (6.626e-34*9.5*1e9)/(9.274e-24*X*1e-4)

        # print(X)
        return(X)
    
        
 
    minortickman = np.arange(0.02, 1, 0.02)
    majortickman=[1, 2, 3, 4.3, 6, 11]

     # Créaction d'un axes secondaire pour afficher les nm.
    secax = ax1.secondary_xaxis('top', functions=(mT2g, mT2g))
 
    if Langue == 'anglais' : 
        secax.set_xlabel('Wavelenght (µm)')
    else : 
        secax.set_xlabel("Facteur de Landé effectif g")
     
    
    secax.xaxis.set_major_locator(ticker.FixedLocator(majortickman))
    secax.xaxis.set_minor_locator(ticker.FixedLocator(minortickman))
            
      
     # secax.xaxis.set_major_locator(ticker.AutoLocator())
     # secax.xaxis.set_minor_locator(ticker.AutoMinorLocator())

    return(secax)


def SecondAxe_nm(ax1, Langue='anglais', GRAPHOUT=None, unite = 'nm'):
    
    if unite == 'um' :
        
        def cm12um(X):
         return(nm2cm1(X)*1/1000)
     
        minortickman = np.arange(0.02, 1, 0.02)
        
        if GRAPHOUT in [3, 6] :
            majortickman=[2, 1.5, 1, 0.8, 0.7, 0.6, 0.5, 0.4, 0.33, 0.3, 0.2]
        elif GRAPHOUT in [5] :
            majortickman=[2, 1, 0.7, 0.6, 0.5, 0.4, 0.33, 0.3, 0.2]
            minortickman = np.arange(0.02, 1, 0.02)
        else:
            majortickman=[2.5, 2, 1.5, 1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.33, 0.3, 0.2]
    
         # Créaction d'un axes secondaire pour afficher les nm.
        secax = ax1.secondary_xaxis('top', functions=(cm12um, cm12um))
     
        if Langue == 'anglais' : 
            secax.set_xlabel('Wavelenght (µm)')
        else : 
            secax.set_xlabel("Longueur d'onde (µm)")
     
    elif unite == 'nm' :
       minortickman = np.arange(20, 1000, 20)
       
       if GRAPHOUT in [1, 6] :
           majortickman = [2000, 1500, 1000, 800, 700, 600, 500, 400, 333, 300, 200]
       else :
           majortickman = [2000, 1000, 800, 700, 600, 500, 400, 333, 300, 200]
    
       secax = ax1.secondary_xaxis('top', functions=(nm2cm1, nm2cm1))
       
       if Langue == 'anglais' : 
            secax.set_xlabel('Wavelenght (nm)')
       else : 
            secax.set_xlabel("Longueur d'onde (nm)")
    
    secax.xaxis.set_major_locator(ticker.FixedLocator(majortickman))
    secax.xaxis.set_minor_locator(ticker.FixedLocator(minortickman))
            
       
     # secax.xaxis.set_major_locator(ticker.AutoLocator())
     # secax.xaxis.set_minor_locator(ticker.AutoMinorLocator())

    return(secax)

def SpectreVIScm(ax1, GRAPHOUT=None, mode = 'cm'):
    #Représentation du spectre Visible en échelle linéaire en cm-1.
        xlim=ax1.get_xlim()
        ylim=ax1.get_ylim()
        taillespectre=0.05*(ylim[1]-ylim[0]); # Calcul de la hauteur de l'arc en ciel 
        
        limiteIR = 850
        limiteUV = 350
        
        if mode ==  'cm' :
            limiteIR = nm2cm1(limiteIR)
            limiteUV = nm2cm1(limiteUV)
            
        limite = np.array([limiteUV, limiteIR])
        # if xlim[0]>nm2cm1(850) or xlim[1]<nm2cm1(350) :  
        if xlim[0]>limite.min() or xlim[1]<limite.max():    
            axins=ax1.inset_axes([xlim[0],ylim[1]-(taillespectre+0.02*taillespectre),
                                  xlim[1]-xlim[0],taillespectre],transform=ax1.transData)
            #Création d'un graph suplémentaire dans la figure principale.
            #Le spectre est ajusté aux valeurs limites de x et de y.
            clim=(xlim[0],xlim[1])
           

        # else :    
        #     axins=ax1.inset_axes([nm2cm1(850),ylim[1]-(taillespectre+0.02*taillespectre),
        #                           nm2cm1(350)-nm2cm1(850),taillespectre],transform=ax1.transData)              
        #     clim=(nm2cm1(850),nm2cm1(350))
        else :    
            axins=ax1.inset_axes([limite.min(),ylim[1]-(taillespectre+0.02*taillespectre),
                                  limite.max()-limite.min(),taillespectre],transform=ax1.transData)              
            clim=(limite.min(),limite.max())
            #Le spectre tiendra dans une boîte de 350 à 780nm mais linéaire en cm-1
        
        norm = plt.Normalize(*clim)
        wn = np.linspace(clim[0],clim[1],2000) #sélection des nombres d'onde

        if mode == 'cm':
            wl=nm2cm1(wn)
        else :
            wl=wn
        
        colorlist = list(zip(norm(wn),[wavelength_to_rgb(w) for w in wl]))
        #on associe pour chaque nombre d'onde la valeur en RGB
        spectralmap = matplotlib.colors.LinearSegmentedColormap.from_list("spectrum", colorlist)
        # conversion en colormap
        
        ymaxinset = 6
        y = np.linspace(0, ymaxinset, 100)
        X,Y = np.meshgrid(wn, y) # création du maillage pour remplir l'image
        
        extent=(np.min(wn), np.max(wn), np.min(y), np.max(y))
        axins.imshow(X, interpolation='none', clim=clim, extent=extent, cmap=spectralmap,
                     aspect='auto', alpha=0.8) #remplissage l'axe avec les couleurs
        axins.tick_params(axis='both', which='both', bottom=False, left=False,
                          labelbottom=False, labelleft=False) #aucun tick
        axins.axis('off')
        


            
        if mode ==  'cm' :
            ytextinset = ymaxinset/3
            
            if GRAPHOUT in [2,3,4,5,6,7,8]:
                ytextinset = 1    
                
                if not GRAPHOUT in [5, 6, 7, 8] :
                   if axins.get_xlim()[0] < nm2cm1(800) : axins.text(nm2cm1(790), ytextinset, 'IR')#,fontsize=tailletexte) 
                   if axins.get_xlim()[1] > nm2cm1(370) : axins.text(nm2cm1(370), ytextinset, 'UV')#,fontsize=tailletexte)
        
    
#Agrandir les légendes
    # params = {'legend.fontsize': tailletexte,'legend.handlelength': 1, 'font.size': tailletexte}
    # plt.rcParams.update(params)
    
    # ax1.tick_params(which='major',length=5)
    # ax1.tick_params(which='minor',length=5)
    # secax.tick_params(which='major',length=5)
    # secax.tick_params(which='minor',length=5)

def Add_checkbutton(line_liste, fig=None, ax=None):
    figclic, rax = plt.subplots(num='Choix affichage')
    rax.clear()
    # fig = ax1.figure
    # rax = fig.add_axes([0.05, 0.4, 0.1, 0.15])
    
    lines_by_label = {l.get_label(): l for l in line_liste}
    line_colors = [l.get_color() for l in lines_by_label.values()]

    check = CheckButtons(
        ax=rax,
        labels=lines_by_label.keys(),
        actives=[l.get_visible() for l in lines_by_label.values()],
        # label_props={'color': line_colors},
        # frame_props={'edgecolor': line_colors},
        # check_props={'facecolor': line_colors},
    )
    def callback(label):
        # nonlocal figclic
        print(label)
        ln = lines_by_label[label]
        ln.set_visible(not ln.get_visible())
        ln.figure.canvas.draw_idle()
        # figclic.canvas.draw_idle()
    
    check.on_clicked(callback)
    return(check)

def ModeX_select(X, ModeX, ax1, valeurnormX, Langue = 'Anglais',  RAJOUT='') :
    if ModeX == 'cm' :
        if Langue == 'anglais' :
            ax1.set_xlabel("Wavenumber ($cm^{-1}$)")
        else :
            ax1.set_xlabel("Nombre d'onde ($cm^{-1}$)");
        RAJOUT = RAJOUT + '_cm'
        
    elif ModeX == 'nm':
        if Langue == 'anglais' :
            ax1.set_xlabel("Wavelenght ($nm$)");
        else :
            ax1.set_xlabel("Longueur d'onde ($nm$)");
        RAJOUT = RAJOUT + '_nm';
        
    elif ModeX == 'eV':
         if Langue == 'anglais' :
             ax1.set_xlabel("Energy ($eV$)");
         else :
             ax1.set_xlabel("Energie ($eV$)");
        
    elif ModeX == 'eVloose':
          if Langue == 'anglais' :
              ax1.set_xlabel("Energy loss ($eV$)");
          else :
              ax1.set_xlabel("Perte d'énergie ($eV$)");   
    
    elif  ModeX == 'g' :
        X = (6.626e-34*valeurnormX*1e9)/(9.274e-24*X*1e-4)
        #X = (6.626e-34*9.437463*1e9)/(9.274e-24*X*1e-4)
        print(X)
        print("pouet + modeXfonctioin")
        ax1.set_xlabel("g")
        RAJOUT = RAJOUT + '_g';
    
    elif  ModeX == 'mT' :
        if Langue == 'anglais' :
            ax1.set_xlabel("Magnetic field (mT)");
        else :
            ax1.set_xlabel("Champs magétique (mT)")
        RAJOUT = RAJOUT + '_champ';
    
    elif ModeX == 'mT_norm':
        # X = (6.626e-34*valeurnormX*1e9)/(9.274e-24*X*1e-4)
        X = (valeurnormX/9.5) * 1/10 * X
        
        if Langue == 'anglais' :
            ax1.set_xlabel("Magnetic field (mT)");
        else :
        
            ax1.set_xlabel("Champs magétique (mT)")
       
        RAJOUT = RAJOUT + '_mT_norm';
        
        
    else :
            ax1.set_xlabel(ModeX);
    return(X,RAJOUT)
    
def ModeY_select(X, Y, ModeY, ax1,  Langue = 'Anglais', valeurnorm='1', RAJOUT='',
                 COUPURENORMminmax=[0, 1], firstpassage=False) :
    
    if ModeY == 'ABS' :
        if Langue == 'anglais' :
            ax1.set_ylabel('Absorbance')
        else : 
              ax1.set_ylabel('Absorbance')
        RAJOUT=RAJOUT+'_ABS'
    
    elif ModeY == 'Int' :
        if Langue == 'anglais' :
              ax1.set_ylabel('Intensity (u.a)')
        else : 
              ax1.set_ylabel('Intensité (u.a)')
        RAJOUT=RAJOUT+'_Intens'
    
    elif ModeY == 'Coup' :
        if Langue == 'anglais' :
              ax1.set_ylabel('Intensity (count)')
        else : 
              ax1.set_ylabel('Intensité (coup)')
        RAJOUT=RAJOUT+'_Intens'
    
    elif ModeY == 'Int RPE' :
        if Langue == 'anglais' :
              ax1.set_ylabel('EPR signal (u.a)')
        else : 
              ax1.set_ylabel('signal RPE (u.a)')
            
    elif ModeY == 'Tr':
        Y = Y*100;
        ax1.set_ylabel('Transmittance (%)') 
        RAJOUT = RAJOUT + '_Tr'
    
    elif ModeY == 'Reflectance':
        Y = Y*100;
        ax1.set_ylabel('Reflectance (%)') 
        RAJOUT = RAJOUT + '_Reflectance'
        
    elif ModeY == 'Kubelka':
        Y = (1-Y)**2/(2*Y);
        ax1.set_ylabel('Approxiation Kubelka-Munk (u.A)')
        RAJOUT = RAJOUT + '_Reflexion'
   
    elif ModeY == 'ABSnormep':
        Y = Y/valeurnorm;
        if Langue == 'anglais' :
            ax1.set_ylabel('Linear absorbance ($cm^{-1}$)') 
        else :
            ax1.set_ylabel('Absorbance linéaire ($cm^{-1}$)') 
        RAJOUT = RAJOUT + '_norm_ep'
    
    elif ModeY == 'Epsilon':
        Y = Y/valeurnorm;
        ax1.set_ylabel('$\\varepsilon (L.mol^{-1}.cm^{-1})$') 
        RAJOUT = RAJOUT + '_Epsilon'
    
    elif ModeY == 'SubRubber':

        INDEXUV = X < (np.max(COUPURENORMminmax))#if COUPURENORMminmax[0] < 360 else 360)
        X=X[INDEXUV]
        Y=Y[INDEXUV]
         
        INDEXIR = X > (np.min(COUPURENORMminmax))#if COUPURENORMminmax[0] < 360 else 360)
        X=X[INDEXIR]
        Y=Y[INDEXIR]
       
        Xbaseline, Ybaseline = baseline_Rubberband(X, Y)
        baseline=np.interp(X, Xbaseline, Ybaseline)
        Y=Y-baseline 
        RAJOUT = RAJOUT + '_subRubber'
        
        if firstpassage :
            ax1.set_ylabel(getylabel(ax1)+
                   ' soustraction ligne de base Rubberband (u.a)')
    
    elif ModeY == 'FondRubber':
        INDEXUV = X < (np.max(COUPURENORMminmax))#if COUPURENORMminmax[0] < 360 else 360)
        X=X[INDEXUV]
        Y=Y[INDEXUV]
         
        INDEXIR = X > (np.min(COUPURENORMminmax))#if COUPURENORMminmax[0] < 360 else 360)
        X=X[INDEXIR]
        Y=Y[INDEXIR]
       
        Xbaseline, Ybaseline = baseline_Rubberband(X, Y)
        baseline=np.interp(X, Xbaseline, Ybaseline)
        Y=baseline 
        RAJOUT = RAJOUT + 'baseline_rubber'
        
        if firstpassage :
            ax1.set_ylabel(getylabel(ax1)+
                   'ligne de base Rubberband (u.a)')
    
    elif ModeY == 'Gradient':
        INDEXUV = X < np.min(COUPURENORMminmax)
        
        # X=X[INDEXUV]
        # Y=Y[INDEXUV]
        Y=np.gradient(Y, X)
        ax1.set_ylabel('dA/dx (u.a)')
        RAJOUT = RAJOUT + '_gradient'
        
    elif ModeY == 'DoubleGradient':
        Y=np.gradient(Y, X)
        Y=np.gradient(Y, X)
        plt.ylabel('ddA/ddx (u.a)')
        RAJOUT = '_doublegradient'
         
    elif ModeY == 'Integrate' :
        X, Y = integrate2curve(X,Y)
        ax1.set_ylabel('$\int f(x) $(u.a)')
        RAJOUT = RAJOUT + '_integrate'
        
    
    else :
        ax1.set_ylabel(ModeY)
    
    
    return(X, Y, RAJOUT)

def NormMode_select(X,Y, normmode, ax1, Langue='anglais', valeurnorm='1',
                    RAJOUT='', COUPURENORMminmax=[0,1], firstpassage=False) :
    
    if normmode == 'min_max' :
        Y = normYminmax(X, Y, COUPURENORMminmax)
        
        if firstpassage :
            ax1.set_ylabel(getylabel(ax1)+
                   ' normalisation min_max')
        RAJOUT = RAJOUT + '_normminmax'
     
    elif normmode == 'ep' :
         if firstpassage :
             if Langue == 'anglais': 
                 ax1.set_ylabel(getylabel(ax1)+
                            ' normalised to thickness')
             else :
                 ax1.set_ylabel(getylabel(ax1)+
                            ' normalisé épaisseur')
         Y = Y/valeurnorm;
         RAJOUT = RAJOUT + 'normep'
    
    elif normmode == 'area' :
        
        INDEX = np.logical_and(X>min(COUPURENORMminmax), X<max(COUPURENORMminmax))   
        Y = Y-Y[INDEX].min()
        norm = np.trapz(Y[INDEX], X[INDEX])
        Y = Y/np.absolute(norm)
        
    return(X, Y, RAJOUT)

def ReadData_spectrotype(Fichier, spectromode='', AdditionY=0, ModeX='cm', ModeY = 'Tr', Limite=None):
    
    if np.size(Limite) == 1 :
        if Limite == None : Limite = np.array([-np.inf, np.inf])

    if spectromode == 'OAS' :
        try : X, Y = Readspectre(Fichier)
        except IndexError:
            try : X, Y = Readspectre(Fichier, delimiter='\t')
            except IndexError : X, Y = Readspectre(Fichier, delimiter=',')
        
        
        Y = Y + AdditionY
 
        if ModeX == 'cm' : # Passage en cm-1
            X = 1/(X*1E-7);
            Limite = nm2cm1(Limite) # Conversion de la limite en cm-1
        elif ModeX == 'eV':
            X = eV2nm(X)
            Limite = eV2nm(Limite)
            
        if not (ModeY == 'Tr' or ModeY == 'Kubelka') : Y = -np.log10(Y);
        
    elif spectromode == 'RAMAN':
        try : X, Y = Readspectre(Fichier, delimiter='\t')
        except IndexError: X, Y = Readspectre(Fichier, delimiter=' ')
    
    elif spectromode == 'RPE':
        try : X, Y = Readspectre(Fichier,skip_header=3, delimiter='\t')
        except IndexError:
            try : X, Y = Readspectre(Fichier, skip_header=3, delimiter='\t')  
            except : X, Y = Readspectre(Fichier, skip_header=0, delimiter=' ')  

        if ModeX == 'mT':
            X=X/10
    
    elif spectromode == 'FLUOOPT' :
        try :
            df = pd.read_csv(Fichier, delimiter='\t', header=None, skiprows=2, usecols = [6, 7], index_col=0)
            X = df.index.to_numpy()
            Y = df.iloc[:, 0].to_numpy()
        except ValueError :
            try :
                df = pd.read_csv(Fichier, delimiter=',', header=None,
                        skiprows=23, usecols = [0, 1], index_col=0)
                X = df.index.to_numpy()
                Y = df.iloc[:, 0].to_numpy()
            except ValueError :
                try : X, Y = Readspectre(Fichier)
                except IndexError: X, Y = Readspectre(Fichier, delimiter='\t')
            

        if ModeX == 'cm' :
            X = 1/(X*1E-7);
            Limite = nm2cm1(Limite)
        elif ModeX == 'eV':
            X = eV2nm(X)
            Limite = eV2nm(Limite)
    else :
        print('\n\n LECTURE PAR DEFAUT \n\n')
        try : X, Y = Readspectre(Fichier, delimiter='\t')
        except IndexError:
            try : X, Y = Readspectre(Fichier, delimiter=';')
            except IndexError: X, Y = Readspectre(Fichier, delimiter=' ')
            
    INDEX =  np.logical_and(X > Limite.min(), X<Limite.max())
    if INDEX.any() : # On vérifie qu'il y est qqchose à afficher
        X=X[INDEX]
        Y=Y[INDEX]
    
    return(X,Y)

#%%  Main affichage
def Affichage_spectre(Liste, Legende, Autoaxe=True, Xlim=[4000, 35000], Ylim=[0, 1.5],
                  TITRE='superposition', ModeX='cm', ModeY='ABS', normmode='',
                  spectromode='OAS', modecouleurs='auto', optionplot='',
                  colormap = 'nipy', GRAPHOUT=1, legendin=True,
                  ax1=None, newgraph=True, SHOW=True, SAVE_HTML = True, markevery=None,
                  COUPURENORMminmax=[400, 2500], 
                  Sousmin=False,  etagage=False, 
                  Langue = 'anglais',valeurnorm=1, valeurnormX = 0, ANNOT = False, annotation='',
                  AdditionY=0, SecondAxe=False, SpectreVIS=True, Dosage=False,
                  LWARN = False, legncol=1, Lissage=False, LimiteAff=None,
                  Z_dim = None , third_dim = False, Z_Label='',
                  ZONES_IRUV = False, TRANSM_WIND = False, SUBPLOT_TRANSITION=False,
                  ListeElement=None, CHECKBUTTON = False) :
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
    AdditionY: float
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
    SAVE_HTML: Bool, optional
        Pour sauvegarder en HTML.     
    SpectreVIS : Bool, optional
        Pour afficher le spectre visible. The default is True.
    Sousmin : Bool, optional
        Pour soustraite par le mininum d'abs. The default is False
    legendin : Bool, optional
        Pour mettre la legende dans la figure où à coté. The default is True
    etagage : Bool, optional
        Pour Etager les coubres. The default is False
    markevery : float
        Option de pyplot pour les marqueurs. The default is 1.
    ax1 : ax
        Ax pyplot sur lequel tracé les coubres. The default is None
    Z_dim = valeur de la troisième dim
    third_dim = activé ou non la troisième dim
    Z_Label= label de l'axe Z,
    
    Returns
    -------
        La figure courante

    '''
    
    AdditionY=mono2tab(AdditionY, np.size(Liste))
    valeurnorm=mono2tab(valeurnorm, np.size(Liste))
    annotation=mono2tab(annotation, np.size(Liste))
    valeurnormX=mono2tab(valeurnormX, np.size(Liste))
    Z_dim = mono2tab(Z_dim, np.size(Liste))
    
    if np.size(LimiteAff) == 1 : LimiteAff=mono2tab(LimiteAff, np.size(Liste))
    
    
    margin=0.07
    etage=0;
    firstpassage=True;
    RAJOUT = ''
    
    FIGSIZE, DPI = GRAPHOUT_mode(GRAPHOUT)
    
    df=pd.DataFrame()
    
    otherpara={}
    line_liste = [] # variable qui servira a stocker les objects line2D
    
    if colormap == 'nipy': colormap = 'nipy_spectral' # Pour comptabilité avec les anciennes versions du code
    try :
        cmap  = matplotlib.colormaps[colormap]
        colorsbigdata = cmap(np.linspace(0,1,np.size(Liste)+2))
    except :
        print("METRE A JOUR MATPLOTLIB POUR UTILISER LE MODE bigdata")

    # if ax1 != None :
    #     if third_dim :
    #         plt.gca().remove() # pour suprimmer l'axes 2D créer lors de set_graph
    #         ax1 = fig.add_subplot(projection='3d')
    #         #ax1.view_init(elev=90, azim=-90, roll=0)
    #     else : ax1 = plt.gca()
    
    #### RAJOUT
    if newgraph and (not ax1): # A améliorier pour bien prendre en compte le subplot
        fig, ax1 = set_graph(FIGSIZE = FIGSIZE, DPI=DPI) # Cf fonction, on créer la figure à la bonne taille/résolution

        if SUBPLOT_TRANSITION : # Uniquement si on crée un graph, alors
            if not ListeElement : ListeElement = read_annot(annotation[0])
            
            if GRAPHOUT == 5 : nrowgrid_base=11
            if GRAPHOUT == 6 : nrowgrid_base=10
            else : nrowgrid_base=10
            
            nrowsub = len(ListeElement)
    
            dim_fig = fig.get_size_inches()
    
            nfig = plt.figure(figsize=(dim_fig[0], dim_fig[1]+0.05*nrowsub), layout="constrained")
            plt.close(fig)
            
            gs = gridspec.GridSpec(nrowgrid_base+nrowsub, 1)   # 2 rows: 1 for existing ax, 1 for the new subfigure
            ax1 = nfig.add_subplot(gs[0:nrowgrid_base-2, :])
            axsub = [nfig.add_subplot(gs[nrowgrid_base+i, :]) for i in range(nrowsub)]
            ax1.grid()
            fig = nfig
        else :
            fig=plt.gcf()
    
        if third_dim :
            margin=0.2
            plt.gca().remove() # pour suprimmer l'axes 2D créer lors de set_graph
            ax1 = fig.add_subplot(projection='3d')
            RAJOUT = RAJOUT + '_3D'
            #ax1.view_init(elev=90, azim=-90, roll=0)
        elif not SUBPLOT_TRANSITION : ax1 = plt.gca()
    
    if Autoaxe:
        pass
    else:
        ax1.set_xlim(Xlim);
        if Ylim : ax1.set_ylim(Ylim);
    
    
    for i in np.arange(0, np.size(Liste)): # on lit tous les fichier
    
        kwarg={}
        Fichier = Liste[i];
        
        X, Y = ReadData_spectrotype(Fichier, spectromode=spectromode,
                                    AdditionY=AdditionY[i], ModeX=ModeX, ModeY = ModeY, Limite=LimiteAff[i])
        
        if Lissage :
            Y = savgol_filter(Y, 51, 3)
            plt.title('ATTENTION ! Lissage!')
            RAJOUT = RAJOUT + '_lisser'
            
        X, RAJOUT = ModeX_select(X, ModeX, ax1, valeurnormX[i], Langue)
        
        X, Y, RAJOUT = ModeY_select(X, Y, ModeY, ax1, Langue, valeurnorm[i],
                                    RAJOUT, COUPURENORMminmax, firstpassage)
        
        X,Y, RAJOUT = NormMode_select(X, Y, normmode, ax1, Langue, valeurnorm[i],
                                    RAJOUT, COUPURENORMminmax, firstpassage)

        if Sousmin :
            if not Autoaxe:
                limsousmin = ax1.get_xlim()
                INDEXsousmin = np.logical_and(X < np.max(limsousmin),
                                              X > np.min(limsousmin));
                Y = Y-np.nanmin(Y[INDEXsousmin])
                
            else : Y = Y-np.nanmin(Y)
            plt.title('ATTENTION ! Soustraction minimum !')
            RAJOUT = RAJOUT + '_sousmin'
        
        Y = Y + etage;
        if etagage :
            INDEX = np.logical_and(X > np.min(COUPURENORMminmax),
                                   X < np.max(COUPURENORMminmax))
            etage = np.max(Y[INDEX]) * 1.2
        
        if optionplot[i] == 'nan':
            optionplot[i] = ''
       
        
        try : test =  float(Legende[i]) # castage de la legende en float, si jamais c'est un nan 
        except : test = 0
        
        if Legende[i] == '' or np.isnan(test) or Legende[i] == np.nan : # vérification que la légende est vide 
            kwarg['label'] = None
        
        else :
            kwarg['label'] = Legende[i]
       
        
        if markevery != None :
            kwarg['markevery'] = markevery
        
        if modecouleurs == 'manuel':
            fmt = optionplot[i]
        elif modecouleurs == 'bigdata':
            fmt = optionplot[i]
            kwarg['color'] = colorsbigdata[i]
        else:
            fmt = ''
            
        try :
            if np.isnan(optionplot[i]) :
                fmt = ''
        except :
            pass
        
            
            
            
        ################# OPTION CON 1 : deux echelle différente ##############################
        DoublePlot_pasmemeechele = True
        if DoublePlot_pasmemeechele == True :
            if Legende[i] == '$Fe^{3+}$' :
                ax1.set_ylabel('Absorbance linéaire ($cm^{-1}$)')
                # ax1.set_ylabel('$\\varepsilon (L.mol^{-1}.cm^{-1})$') 
                ax2 = ax1
                ax2.tick_params(axis='y', labelcolor='C0')
                leg = ax2.legend(loc='upper left')
                
                leg.set_draggable(True)
                ax1 = ax2.twinx()
                ax1.tick_params(axis='y', labelcolor='C8')
        ###############################
        
        ##################### AFFICHAGE
        if third_dim : 
             # print('OUIIIIIIIIIIIIIIIIIIIIIII')
             if Z_dim[i] == None : # SI on ne donne pas la valeur à mettre en 3ème dimnesion.
                     Z = mono2tab(i, np.size(X))
             else : Z = mono2tab(Z_dim[i], np.size(X))
             ax1.set_zlabel(getylabel(ax1)) #Changement nom axe z/Y
             ax1.set_ylabel(Z_Label)
             
             ax1.plot(X, Z, Y, fmt, **kwarg)
             
        else :
            try :
                ltemp, = ax1.plot(X, Y, fmt,  **kwarg)
                line_liste.append(ltemp)
            except :
                ax1.plot(X, Y, fmt,  **kwarg)
                
        ################## OPTION CON 2 : verre fer avec fill between             
        if Legende[i] == 'NCS05' and False:
            print('MODE VERRE VINCENT ON')
            X_OMCT = np.linspace(20000, 30000, 500)
            p = [102.827886466, 38852.5703782, 2930.82158635, 0.0514]
            Y_OMCT = Gauss(X_OMCT, *p)
            
            ax1.fill_between(X, Y, Y.min(),  where=X<20000, alpha = 0.5, facecolor='C0')
            ax1.fill_between(X, Y, Y.min(), where=X>20000, alpha = 0.5, color='y')
            ax1.fill_between(X_OMCT, Y_OMCT,Y_OMCT.min(), alpha = 0.8, color='k')
      ##########################################################""
        
        if ANNOT : 
            try : OAS_annotate(X, Y, annotation[i], ax1, lissage = False)
            except : print('annot fail') #pass
        
        if Dosage :
            try : dtemp = max2C(X,Y, annotation[i], Legende[i])
            except UnboundLocalError : dtemp = pd.DataFrame()
            
            df = pd.concat([df, dtemp])
        firstpassage=False;
        
    #ax1.legend()
    if third_dim :
        RAJOUT = RAJOUT + '_3D'
    TITRE=TITRE+RAJOUT

    ax1.xaxis.set_minor_locator(ticker.AutoMinorLocator())


    if spectromode == 'OAS' or spectromode == 'FLUOOPT' :
        if LWARN : warning_light(ax1)
        
        if ZONES_IRUV : zones_IR_UV(ax1)
        
        if TRANSM_WIND : Transm_window(X, Y, ax1)
        
        if SecondAxe and (ModeX == 'cm') : # Parti double échelle mettre false pour avoir uniquement en cm^-1
           SecondAxe_nm(ax1, Langue, GRAPHOUT)   
    
        if SpectreVIS and (ModeX == 'cm' or ModeX == 'nm'):
            if Autoaxe : ax1.margins(y=margin)
            SpectreVIScm(ax1, GRAPHOUT, mode = ModeX)
    
    if SecondAxe and (ModeX == 'mT_norm'):
        print('test')
        SecondAxe_g(ax1, Langue, GRAPHOUT)
        
    if SUBPLOT_TRANSITION :
        plot_transition_subplot(ax1, axsub, ListeElement)
        fig.subplots_adjust(top=0.875, bottom=0.015, left=0.09, right=0.99, hspace=0.1, wspace=0.0)
    # plt.title(TITRE)
    
    if Autoaxe : ax1.autoscale_view()    

    if CHECKBUTTON :
       check = Add_checkbutton(line_liste)
    else :
        check = None

    if SHOW : Sav_fig(TITRE, legendin=legendin, SAVEHTML=SAVE_HTML, legncol=legncol, fig=fig, ax=ax1)
    if Dosage : saveDataFrame(df, TITRE+'_dosage', 'Dosage')
    
    def on_resize(event):
        fig.tight_layout()
        fig.canvas.draw()
        
    if newgraph : fig.canvas.mpl_connect('resize_event', on_resize)
    else : fig = plt.gcf()

    otherpara['checkbox'] = check
    otherpara['line_liste'] = line_liste
    
    return(fig, ax1, TITRE, otherpara) #ON DOIT ABSOLUMENT GARDER UNE REFERENCE DE CHECK, sinon on pert l'interactivité


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
    Data=np.genfromtxt('./parametre_fit_gauss Fe-S.csv', comments='#', skip_header=1, delimiter=';', dtype='str' )
    
    X = np.linspace(3000, 45000, 200)
    plt.style.use({'figure.facecolor': 'white', 'savefig.transparent' : True, 'lines.markersize': 6.0})
    for i in np.arange(0, np.size(Data[:, 0])):
        if Data[i, 0] in Legende :
            Y1 = Gauss(X, float(Data[i,1]), float(Data[i,2]), float(Data[i,3]), float(Data[i,4]))
            Y2 = Gauss(X, float(Data[i,5]), float(Data[i,6]), float(Data[i,7]), float(Data[i,8]))
            #plt.plot(X, Y1, optplt, label='fit 1 '+Data[i, 0])
            #plt.plot(X, Y2, optplt, label='fit 2 '+Data[i, 0])
            plt.plot(X, Y1, optplt)
            plt.plot(X, Y2, optplt)
            Y_tot = Y1+Y2
            #plt.plot(X, Y_tot, optplt, label='fit total '+Data[i, 0])
            plt.plot(X, Y_tot, optplt)
        else:
            pass
    if SHOW : Sav_fig(TITRE, legendin=False)
    return()

#%% Colorimetrie

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
    
        XYZ[i] = colour.sd_to_XYZ(sd, cmfs, illuminant)/100 # Passage de la densité spectral au trismultus
    return(XYZ)
    
def AffichageCIE1931(Liste, Legende, TITRE='CIE1931', Marqueur='', Fleche=True,
                         xylim=[-0.1, 0.8, 0, 1], show=True, legendin=True,
                         GRAPHOUT = 2, full=False, inset=False):
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
    
    if GRAPHOUT == 1 :
        FIGSIZE = (9,5)
        DPI = 500
        
    elif GRAPHOUT == 2 : 
        FIGSIZE = (6, 3.2)
        DPI = 150
        
    elif GRAPHOUT == 3 : 
        FIGSIZE = (5, 5)
        DPI = 200
        
    elif GRAPHOUT == 4 :
        FIGSIZE = (4.5, 2.5)
        DPI = 200
        
    elif GRAPHOUT == 5 :
        FIGSIZE=[5, 3.5]
        DPI = 200
        
    xtable=[]
    ytable=[]

    
    set_graph(FIGSIZE=FIGSIZE, DPI=DPI, mode='CIE1931', xylim=xylim)

    if TITRE == 'CIE1931':
        pass
    else:
        TITRE = TITRE + '_CIE1931'
    
    Marqueur=mono2tab(Marqueur, np.size(Liste))    
    tristimulus=Calcul_XYZ(Liste, Legende)
    saveDATACIE_XYZ(TITRE[:-8], Legende, tristimulus[:, 0], tristimulus[:, 1],
                    tristimulus[:, 2])
    
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
            #plt.plot(x, y, marker=Marqueur[i], label=Legende[i])
            plt.plot(x, y, Marqueur[i], label=Legende[i])
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

    if full :
        plt.box(False)
        plt.xticks([])
        plt.yticks([])
        plt.xlabel('')
        plt.ylabel('')
        
    if inset :
        plt.title('')
        
    if show: Sav_fig(TITRE, colorplot=True, legendin=legendin)
    saveDATACIE_xy(TITRE[:-8], Legende, xtable, ytable)
    
    return(xtable, ytable)


def Affichage_Lab2D(Liste, Legende, TITRE='Lab', Marqueur='', SHOW=True,
                    legendin=True):
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
    fig.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.95, wspace=0, hspace=0)
    ax1 = axs[0]
    ax2 = axs[1]
    Lab_liste=[]
    
    for i in np.arange(0, np.size(Liste)):
        XYZ=tristimulus[i]
        
        Lab=colour.XYZ_to_Lab(XYZ);
        Lab_liste.append(Lab)
        
        if Marqueur[i]:
            ax1.plot(Lab[1], Lab[0], Marqueur[i], label=Legende[i])
            ax2.plot(Lab[1], Lab[2], Marqueur[i], label=Legende[i])
        else:
            ax1.plot(Lab[1], Lab[0], 'x', label=Legende[i])
            ax2.plot(Lab[1], Lab[2], 'x', label=Legende[i])
            
    ax1.set_ylabel('L')
    ax1.grid()
    if legendin : ax1.legend()
    else : ax1.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    
    ax2.set_ylabel('b')
    ax2.set_xlabel('a')
    ax2.grid()
 
    #fig.tight_layout()
    if SHOW : Sav_fig(TITRE, Tight=False, legendin=legendin, legend=False)
    return(Lab_liste)
        

def Affichage_Lab3D(Liste, Legende, TITRE='Lab3D', Marqueur='', SHOW=True,
                    legendin=True, GRAPHOUT=1):
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
    if GRAPHOUT == 1 :
        FIGSIZE=(8,8)
        DPI=120
    
    elif GRAPHOUT == 2 :
        FIGSIZE=(5,5)
        DPI=120
    
    
    Marqueur=mono2tab(Marqueur, np.size(Liste))    

    tristimulus=Calcul_XYZ(Liste, Legende)
    Lab_liste=[]
    
    if TITRE == 'Lab3D':
        pass
    else:
        TITRE = TITRE + '_Lab3D'

    
    fig = plt.figure(figsize=FIGSIZE, dpi=DPI)
    ax = fig.add_subplot(111, projection='3d')
    
        
    for i in np.arange(0, np.size(Liste)):
        XYZ=tristimulus[i]
        
        Lab=colour.XYZ_to_Lab(XYZ);
        Lab_liste.append(Lab)
        
        L=Lab[0]
        a=Lab[1]
        b=Lab[2]
        
               
        if Marqueur[i] == '' :            
            ax.scatter(a, b, L, marker='x', label=Legende[i])
        else:
            ax.plot(a, b, L, Marqueur[i], label=Legende[i])

                
    ax.set_xlabel('a')
    ax.set_ylabel('b')
    ax.set_zlabel('L')
    # ax.set_xlim([-300, 300])
    # ax.set_ylim([-300, 300])
    # ax.set_zlim([-350,350])
 
    #fig.tight_layout()
    if SHOW : Sav_fig(TITRE, legendin=legendin)
    
    L = [x[0] for x in Lab_liste] 
    a = [x[1] for x in Lab_liste] 
    b = [x[2] for x in Lab_liste] 
    saveDATACIE_LAB(TITRE[0:-6], Legende, L, a, b)
    
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
        
        XYZ[i] = colour.sd_to_XYZ(sd, cmfs, illuminant)/100
        
        if coef==1: ref_XYZ=XYZ[i]

        i+=1
        
    return(XYZ, Legende, ref_XYZ)
    
    
def Courbe_BeerLambert_XY(Fichier, Legende, TITRE='CIE1931', fondCIE=False,
                          xylim=[-0.1, 0.8, -0.1, 1], lc_lim=[0.1, 5],
                          legendin=True, show=True, style=''):
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
        

    if style == '' : plt.plot(xtable, ytable, label='Effet lc ' + Verre)
    else :  plt.plot(xtable, ytable, style, label='Effet lc ' + Verre)
    #plt.legend(loc='upper right')
    plt.legend(bbox_to_anchor = [1, 1])

    if show: Sav_fig(TITRE, colorplot=True, legendin=legendin)
    
    return(xtable, ytable)
    
    
def Courbe_BeerLambert_Lab3D(Fichier, Legende='', TITRE='Lab', lc_lim=[0.1, 5],
                             legendin=True, show=True, newfig=False):
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
 
    if show : Sav_fig(TITRE, legendin=legendin)
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

#%% ID 20
def Affichage_ID20(Liste, Legende, Autoaxe=True, Xlim=[0, 800], Ylim=[0, 1],
                  TITRE='superposition', linewidth=1.5, Modeaff='default',
                  modecouleurs='auto', optionplot='', SHOW=True, GRAPHOUT=1,
                  newgraph=True, REPsave='Graph'):
 
    
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
        
        elif Modeaff == 'DoubleGradient':
            Y=np.gradient(Y, X)
            Y=np.gradient(Y, X)
            plt.ylabel('ddA/ddx (u.a)')
            RAJOUT = '_doublegradient'
            
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

    if SHOW : Sav_fig(TITRE, REPsave)
    
    return(fig)       
        
#%% Annotate et récupération epsilon

def annot_max(x, y, text='pouet', ax=None, xlim=None, mode=1, xytext=None):
    '''
    Cette fonction sert à annoter une courbe au niveau de son max.

    Parameters
    ----------
    x : TYPE
        tableeau de donnée.
    y : TYPE
        tableau de donnée.
    text : TYPE, optional
        text a placer sur la figure. The default is 'pouet'.
    ax : TYPE, optional
        figure à annoter, si None récupère la figure courrante. The default is None.
    xlim : TYPE, optional
        pour sélectionner une limite en X pour trouve le y max . The default is None.
    mode : TYPE, optional
        DESCRIPTION. The default is 1.
    xytext : TYPE, optional
        Position du texte par rapport au x/y du maximum. The default is None.

    Returns
    -------
    None.

    '''
    if xytext == None:
        xytext=(10, 20)
    
    
    if not (xlim == None) :
        limite= np.logical_and(x > xlim[0], x < xlim[1])
        x=x[limite]
        y=y[limite]
        
    xmax = x[np.argmax(y)]
    ymax = y.max()
    
    if ax == None :
        ax=plt.gca()
    
    if mode == 1 :
        annotation = ax.annotate(text,
                          xy=(xmax, ymax),
                          xytext=xytext,
                          textcoords='offset points',
                          arrowprops=dict(arrowstyle='->', connectionstyle='arc3, rad=0.2'))
    if mode == 2 :
        bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
        arrowprops=dict(arrowstyle="->",connectionstyle="angle,angleA=0,angleB=60")
        kw = dict(xycoords='data',textcoords="axes fraction",
                  arrowprops=arrowprops, bbox=bbox_props, ha="right", va="top")
        annotation = ax.annotate(text, xy=(xmax, ymax), xytext=(0.94,0.96), **kw)
    
    annotation.draggable()
    
    return()


def plot_transition_subplot(mainax, axsub, Listelement) :

    barX = np.linspace(np.min(mainax.get_xlim()), np.max(mainax.get_xlim()), 1000)
    extent=(np.min(mainax.get_xlim()), np.max(mainax.get_xlim()), 0, 1)

    for Element, axgs in zip(Listelement, axsub) :
        custom_cmap = element2epsilon(Element)['cmap'].to_list()[0]
        
        #Limite = np.array([lim for lim in element2epsilon(Element)['Xlim']])
        barY = np.zeros(barX.shape)
        DfEl = element2epsilon(Element)
        name = DfEl['name'].to_list()[0]
  
        for transition in DfEl.index : 
            # Limite = np.array([lim for lim in element2epsilon(Element)['Xlim']])
            lim = np.array(DfEl.loc[transition, 'Xlim'])
            epsi = DfEl.loc[transition, 'epsilon']

            barY = barY + epsi*np.logical_and(barX>lim.min(), barX<lim.max()).astype(int)

            coordinence = DfEl.loc[transition, 'coordinence']

            if coordinence == 5:
                hatch_pattern = "/"
              
                rect = Rectangle((np.min(lim), 0), np.max(lim) - np.min(lim), 1, hatch=hatch_pattern, fill=False)
                axgs.add_patch(rect)
            
        bbarXX, bbarYY = np.meshgrid(barX, barY)
        bbarYY.astype(int)
        axgs.imshow(bbarYY.T, aspect='auto', cmap=custom_cmap.cmap, extent = extent) 
        axgs.text(-0.01, 0.5, name, va='center', ha='right', fontsize=10,
                 transform=axgs.transAxes)
        axgs.yaxis.set_visible(False)
        axgs.xaxis.set_visible(False)

def element2epsilon(Element):
    if Element == 'Co2+':
        #Données pour le Co2+ issu de la thèse de Myrtille (page 144, section 4.3 Relation entre composition et spéciation du Co2+)
        #Données range : page 94 thèse de Myrtille
        dfCo=pd.DataFrame(columns = ['Xlim', 'epsilon', 'coordinence'],
                          index=['4A2->4T1(4F)1', '4A2->4T1(4F)2', '4A2->4T1(4F)3',
                                 '4A2->2T1(2G)', '4A2->4T1(P)', '4A2->2T2(2G)']) # Association des transitions selon (1) Keppler, H. Crystal Field Spectra and Geochemistry of Transition Metal Ions in Silicate Melts and Glasses. American Mineralogist 1992, 77 (1–2), 62–75.

        dfCo.loc['4A2->4T1(4F)1', 'Xlim'] = [5809.968847, 6183.800623]
        dfCo.loc['4A2->4T1(4F)1', 'epsilon'] = 29.5 #Lecture graphique p139 NS2
        dfCo.loc['4A2->4T1(4F)1', 'coordinence'] = 4
        
        dfCo.loc['4A2->4T1(4F)2', 'Xlim'] = [6619.937695, 6993.76947]
        dfCo.loc['4A2->4T1(4F)2', 'epsilon'] = 35 #Lecture graphique p139 NS2
        dfCo.loc['4A2->4T1(4F)2', 'coordinence'] = 4
        
        dfCo.loc['4A2->4T1(4F)3', 'Xlim'] = [7492.211838, 7897.196262]
        dfCo.loc['4A2->4T1(4F)3', 'epsilon'] = 31 #Lecture graphique p139 NS2
        dfCo.loc['4A2->4T1(4F)3', 'coordinence'] = 4
        
        dfCo.loc['4A2->2T1(2G)', 'Xlim'] = [15529.59502, 15934.57944]
        dfCo.loc['4A2->2T1(2G)', 'epsilon'] = 100 #pour une matrice NaCa
        dfCo.loc['4A2->2T1(2G)', 'coordinence'] = 4
        
        dfCo.loc['4A2->4T1(P)', 'Xlim'] = [16526.47975, 16993.76947]
        dfCo.loc['4A2->4T1(P)', 'epsilon'] = 106 #pour une matrice NaCa
        dfCo.loc['4A2->4T1(P)', 'coordinence'] = 4
        
        dfCo.loc['4A2->2T2(2G)', 'Xlim'] = [18084.11215, 18676.01246]
        dfCo.loc['4A2->2T2(2G)', 'epsilon'] = 71 #pour une matrice NaCa
        dfCo.loc['4A2->2T2(2G)', 'coordinence'] = 4
        
        # dfCo.loc['', 'Xlim'] = [6619.937695, 8115.264798]
        # dfCo.loc['', 'epsilon'] = np.nan
        # dfCo.loc['', 'coordinence'] = 6
        
        # dfCo.loc['', 'Xlim'] = [11105.919, 12258.56698]
        # dfCo.loc['', 'epsilon'] = np.nan
        # dfCo.loc['', 'coordinence'] = 6
        
        # dfCo.loc['', 'Xlim'] = [17336.4486, 18800.62305]
        # dfCo.loc['', 'epsilon'] = np.nan
        # dfCo.loc['', 'coordinence'] = 6
        
        # dfCo.loc['', 'Xlim'] = [5062.305296, 6526.479751]
        # dfCo.loc['', 'epsilon'] = np.nan
        # dfCo.loc['', 'coordinence'] = 5
        
        dfCo.loc['T1(5)', 'Xlim'] = [8146.417445, 10358.25545]
        dfCo.loc['T1(5)', 'epsilon'] = 8 #Lecture graphique p139 NS2 
        dfCo.loc['T1(5)', 'coordinence'] = 5
        
        # dfCo.loc['T2(5)', 'Xlim'] = [15000, 17180.68536]
        # dfCo.loc['T2(5)', 'epsilon'] = np.nan
        # dfCo.loc['T2(5)', 'coordinence'] = 5
        
        dfCo.loc['T3(5)', 'Xlim'] = [18707.16511, 21012.46106]
        dfCo.loc['T3(5)', 'epsilon'] = 40 #Lecture graphique p139 NS2 
        dfCo.loc['T3(5)', 'coordinence'] = 5
        
        custom_cmap = mcolors.LinearSegmentedColormap.from_list('custom', ['#FFFFFF', '#173f87'])
        
        norm = matplotlib.colors.Normalize(vmin=0, vmax=np.nanmax(dfCo['epsilon']))
        cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=custom_cmap)

        dfCo['cmap'] = cmap
        dfCo['name'] = '$Co^{2+}$'

        df=dfCo
        
        df=dfCo
        
    elif Element == 'Cu2+'  :
        dfCu=pd.DataFrame(columns = ['Xlim', 'epsilon', 'coordinence'],
                           index=['2Eg(D)->2T2g(D)'])

        dfCu.loc['2Eg(D)->2T2g(D)', 'Xlim'] = [8000, 17000]
        dfCu.loc['2Eg(D)->2T2g(D)', 'epsilon'] = 30 #Matrice NaCa; Möncke et al. /Journal of Archaeological Science (2014)
        dfCu.loc['2Eg(D)->2T2g(D)', 'coordinence'] = 6 #Oh with tetragonal distorsion (Jahn-Teller d9)
        
        # df_temp = dfCu.filter(like='epsilon', level=1, axis=0)
        
        #Creation of the colormap
        custom_cmap = mcolors.LinearSegmentedColormap.from_list('custom', ['#FFFFFF', '#1C8FA3'])
        
        norm = matplotlib.colors.Normalize(vmin=0, vmax=dfCu['epsilon'].max())
        cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=custom_cmap)

        dfCu['cmap'] = cmap
        dfCu['name'] = '$Cu^{2+}$'
        df=dfCu
    
        
    elif Element == 'Cu+'  :
        dfCumono=pd.DataFrame(columns = ['Xlim', 'epsilon', 'coordinence'],
                           index=['intervalence_possible'])

        dfCumono.loc['intervalence_possible', 'Xlim'] = [18000, 30000]
        dfCumono.loc['intervalence_possible', 'epsilon'] = 200 #au pif
        dfCumono.loc['intervalence_possible', 'coordinence'] = 1 # pas etudier
        
        # df_temp = dfCu.filter(like='epsilon', level=1, axis=0)
        
        #Creation of the colormap
        custom_cmap = mcolors.LinearSegmentedColormap.from_list('custom', ['#FFFFFF', '#2ca02c'])
        
        norm = matplotlib.colors.Normalize(vmin=0, vmax=dfCumono['epsilon'].max())
        cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=custom_cmap)

        dfCumono['cmap'] = cmap
        dfCumono['name'] = '$Cu^{+}-Cu^{2+}$'
        df=dfCumono
    
    elif Element == 'Mn3+' : 
        dfMn=pd.DataFrame(columns = ['Xlim', 'epsilon', 'coordinence'],
                           index=['5Eg(D)->5T2g(D)'])
        
        dfMn.loc['5Eg(D)->5T2g(D)', 'Xlim'] = [14000, 25000] 
        dfMn.loc['5Eg(D)->5T2g(D)', 'epsilon'] = 135 #135 d'après le Handbook of Glass | 130 d'après la thèse de Natan p154 | 28,3 d'après Nelson & White, Geochemica et Cosmochimica Acta (1980)
        dfMn.loc['5Eg(D)->5T2g(D)', 'coordinence'] = 6 ##Oh with tetragonal distorsion (d4)
                
        #Creation of the colormap
        custom_cmap = mcolors.LinearSegmentedColormap.from_list('custom', ['#FFFFFF', '#AA1675'])
        
        norm = matplotlib.colors.Normalize(vmin=0, vmax=dfMn['epsilon'].max())
        cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=custom_cmap)

        dfMn['cmap'] = cmap
        dfMn['name'] = '$Mn^{3+}$'
        
        df=dfMn
        
    elif Element == 'Fe2+':
        dfFe=pd.DataFrame(columns = ['Xlim', 'epsilon', 'coordinence'],
                           index=['5T2g(D)->5Eg(D)'])
        
        dfFe.loc['5T2g(D)->5Eg(D)', 'Xlim'] = [4000, 15000] 
        dfFe.loc['5T2g(D)->5Eg(D)', 'epsilon'] = 30 #30 d'après lim inf Handbook of Glass
        dfFe.loc['5T2g(D)->5Eg(D)', 'coordinence'] = 6 ##Oh with tetragonal distorsion (d4)
        
        # dfFe.loc['5T2g(D)->5Eg(D)', 'Xlim'] = [20000, 22000] 
        # dfFe.loc['5T2g(D)->5Eg(D)', 'epsilon'] = np.nan() #d'après lim inf Handbook of Glass
        # dfFe.loc['5T2g(D)->5Eg(D)', 'coordinence'] = 6 ##Oh with tetragonal distorsion (d4)
                        
        #Creation of the colormap
        custom_cmap = mcolors.LinearSegmentedColormap.from_list('custom', ['#FFFFFF', '#1F77B4'])
        
        norm = matplotlib.colors.Normalize(vmin=0, vmax=dfFe['epsilon'].max())
        cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=custom_cmap)

        dfFe['cmap'] = cmap
        dfFe['name'] = '$Fe^{2+}$'
        
        df=dfFe

    elif Element == 'Fe3+':
        dfFe3=pd.DataFrame(columns = ['Xlim', 'epsilon', 'coordinence'],
                           index=['6A1(S)->4A1/4E', '6A1(S)->4T2(G)', '6A1(S)->4T1(G)'])
        
        dfFe3.loc['6A1(S)->4A1/4E', 'Xlim'] = [22750, 23200]
        dfFe3.loc['6A1(S)->4A1/4E', 'epsilon'] = 4 #30 d'après lim inf Handbook of Glass
        dfFe3.loc['6A1(S)->4A1/4E', 'coordinence'] = 6 ##Oh with tetragonal distorsion (d4)
        
        dfFe3.loc['6A1(S)->4T2(G)', 'Xlim'] = [23730, 24100]
        dfFe3.loc['6A1(S)->4T2(G)', 'epsilon'] = 4 #d'après lim inf Handbook of Glass
        dfFe3.loc['6A1(S)->4T2(G)', 'coordinence'] = 6 ##Oh with tetragonal distorsion (d4)
        
        dfFe3.loc['6A1(S)->4T1(G)', 'Xlim'] = [26170, 26500]
        dfFe3.loc['6A1(S)->4T1(G)', 'epsilon'] = 12 #d'après lim inf Handbook of Glass
        dfFe3.loc['6A1(S)->4T1(G)', 'coordinence'] = 6 ##Oh with tetragonal distorsion (d4)
                        
        #Creation of the colormap
        custom_cmap = mcolors.LinearSegmentedColormap.from_list('custom', ['#FFFFFF', '#BCBD22'])
        
        norm = matplotlib.colors.Normalize(vmin=0, vmax=dfFe3['epsilon'].max())
        cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=custom_cmap)

        dfFe3['cmap'] = cmap
        dfFe3['name'] = '$Fe^{3+}$'
        
        df=dfFe3
        
    elif Element == 'NiOh': # d'aprés Galoisy, L., Calas, G., 2022. Role of alkali field strength on the speciation of Ni2+ in alkali borate glasses: comparison with crystalline Ni-borates. Journal of Non-Crystalline Solids 577, 121320. https://doi.org/10.1016/j.jnoncrysol.2021.121320

        dfNiOh=pd.DataFrame(columns = ['Xlim', 'epsilon', 'coordinence'],
                           index=['T1', 'T2', 'T3', 'T4', 'T5'])
        
        dfNiOh.loc['T1', 'Xlim'] = [24000, 24500]
        dfNiOh.loc['T1', 'epsilon'] = 10             #ATTENTION epsilon au pif juste pour avoir les bonnes intensité dans le sousplot
        dfNiOh.loc['T1', 'coordinence'] = 6 
        
        dfNiOh.loc['T2', 'Xlim'] = [14500, 15000]
        dfNiOh.loc['T2', 'epsilon'] = 4 
        dfNiOh.loc['T2', 'coordinence'] = 6 
        
        dfNiOh.loc['T3', 'Xlim'] = [13000, 13500]
        dfNiOh.loc['T3', 'epsilon'] = 6 
        dfNiOh.loc['T3', 'coordinence'] = 6 

        dfNiOh.loc['T4', 'Xlim'] = [8200, 86000]
        dfNiOh.loc['T4', 'epsilon'] = 3 
        dfNiOh.loc['T4', 'coordinence'] = 6 
                        
        dfNiOh.loc['T5', 'Xlim'] = [6800, 7200]
        dfNiOh.loc['T5', 'epsilon'] = 3 
        dfNiOh.loc['T5', 'coordinence'] = 6 

                        
        #Creation of the colormap
        custom_cmap = mcolors.LinearSegmentedColormap.from_list('custom', ['#FFFFFF', 'g'])
        
        norm = matplotlib.colors.Normalize(vmin=0, vmax=dfNiOh['epsilon'].max())
        cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=custom_cmap)

        dfNiOh['cmap'] = cmap
        dfNiOh['name'] = '$^{[6]}Ni$ : Oh'
        
        df=dfNiOh
        
    elif Element == 'NiTd': # d'aprés Galoisy, L., Calas, G., 2022. Role of alkali field strength on the speciation of Ni2+ in alkali borate glasses: comparison with crystalline Ni-borates. Journal of Non-Crystalline Solids 577, 121320. https://doi.org/10.1016/j.jnoncrysol.2021.121320

        dfNiTd=pd.DataFrame(columns = ['Xlim', 'epsilon', 'coordinence'],
                           index=['T1', 'T2', 'T3', 'T4', 'T5', 'T6'])
        
        dfNiTd.loc['T1', 'Xlim'] = [19500, 20000] 
        dfNiTd.loc['T1', 'epsilon'] = 9            #ATTENTION epsilon au pif juste pour avoir les bonnes intensité dans le sous plot
        dfNiTd.loc['T1', 'coordinence'] = 4 
        
        dfNiTd.loc['T2', 'Xlim'] = [17500, 18000]
        dfNiTd.loc['T2', 'epsilon'] = 10
        dfNiTd.loc['T2', 'coordinence'] = 4
        
        dfNiTd.loc['T3', 'Xlim'] = [15200, 15700]
        dfNiTd.loc['T3', 'epsilon'] = 10 
        dfNiTd.loc['T3', 'coordinence'] = 4 

        dfNiTd.loc['T4', 'Xlim'] = [12500, 13000]
        dfNiTd.loc['T4', 'epsilon'] = 1 
        dfNiTd.loc['T4', 'coordinence'] = 4
                        
        dfNiTd.loc['T5', 'Xlim'] = [9000, 9500]
        dfNiTd.loc['T5', 'epsilon'] = 4 
        dfNiTd.loc['T5', 'coordinence'] = 4
  
        dfNiTd.loc['T6', 'Xlim'] = [4500, 5000]
        dfNiTd.loc['T6', 'epsilon'] = 3 
        dfNiTd.loc['T6', 'coordinence'] = 4
    
        #Creation of the colormap
        custom_cmap = mcolors.LinearSegmentedColormap.from_list('custom', ['#FFFFFF', 'b'])
        
        norm = matplotlib.colors.Normalize(vmin=0, vmax=dfNiTd['epsilon'].max())
        cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=custom_cmap)

        dfNiTd['cmap'] = cmap
        dfNiTd['name'] = '$^{[4]}Ni$ : Td'
        
        df=dfNiTd
        
    
    elif Element == 'Nisp': # d'aprés Galoisy, L., Calas, G., 2022. Role of alkali field strength on the speciation of Ni2+ in alkali borate glasses: comparison with crystalline Ni-borates. Journal of Non-Crystalline Solids 577, 121320. https://doi.org/10.1016/j.jnoncrysol.2021.121320

        dfNisp=pd.DataFrame(columns = ['Xlim', 'epsilon', 'coordinence'],
                           index=['T1', 'T2', 'T3', 'T4', 'T5'])
        
        dfNisp.loc['T1', 'Xlim'] = [23000, 23500]
        dfNisp.loc['T1', 'epsilon'] = 10             #ATTENTION epsilon au pif juste pour avoir les bonnes intensité dans le sousplot
        dfNisp.loc['T1', 'coordinence'] = 51
        
        dfNisp.loc['T2', 'Xlim'] = [19000, 19500]
        dfNisp.loc['T2', 'epsilon'] = 4
        dfNisp.loc['T2', 'coordinence'] = 51
        
        dfNisp.loc['T3', 'Xlim'] = [12000, 12500]
        dfNisp.loc['T3', 'epsilon'] = 6 
        dfNisp.loc['T3', 'coordinence'] = 51

        dfNisp.loc['T4', 'Xlim'] = [11000, 115000]
        dfNisp.loc['T4', 'epsilon'] = 3 
        dfNisp.loc['T4', 'coordinence'] = 51
                        
        dfNisp.loc['T5', 'Xlim'] = [5000, 5500]
        dfNisp.loc['T5', 'epsilon'] = 3 
        dfNisp.loc['T5', 'coordinence'] = 51

                        
        #Creation of the colormap
        custom_cmap = mcolors.LinearSegmentedColormap.from_list('custom', ['#FFFFFF', 'C1'])
        
        norm = matplotlib.colors.Normalize(vmin=0, vmax=dfNisp['epsilon'].max())
        cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=custom_cmap)

        dfNisp['cmap'] = cmap
        dfNisp['name'] = '$^{[5]}Ni$ : SP'
        
        df=dfNisp        
        
        

    elif Element == 'Nitbp': # d'aprés Galoisy, L., Calas, G., 2022. Role of alkali field strength on the speciation of Ni2+ in alkali borate glasses: comparison with crystalline Ni-borates. Journal of Non-Crystalline Solids 577, 121320. https://doi.org/10.1016/j.jnoncrysol.2021.121320

        dfNitpb=pd.DataFrame(columns = ['Xlim', 'epsilon', 'coordinence'],
                           index=['T1', 'T2', 'T3', 'T4', 'T5', 'T6'])
        
        dfNitpb.loc['T1', 'Xlim'] = [21500, 22000] 
        dfNitpb.loc['T1', 'epsilon'] = 10           #ATTENTION epsilon au pif juste pour avoir les bonnes intensité dans le sous plot
        dfNitpb.loc['T1', 'coordinence'] = 51
        
        dfNitpb.loc['T2', 'Xlim'] = [17600, 18100]
        dfNitpb.loc['T2', 'epsilon'] = 6
        dfNitpb.loc['T2', 'coordinence'] = 51
        
        dfNitpb.loc['T3', 'Xlim'] = [12200, 12800]
        dfNitpb.loc['T3', 'epsilon'] = 2
        dfNitpb.loc['T3', 'coordinence'] = 51

        dfNitpb.loc['T4', 'Xlim'] = [10100, 10500]
        dfNitpb.loc['T4', 'epsilon'] = 3 
        dfNitpb.loc['T4', 'coordinence'] =51
                        
        dfNitpb.loc['T5', 'Xlim'] = [9400, 9800]
        dfNitpb.loc['T5', 'epsilon'] = 3
        dfNitpb.loc['T5', 'coordinence'] = 51
  
        dfNitpb.loc['T6', 'Xlim'] = [5500, 6000]
        dfNitpb.loc['T6', 'epsilon'] = 3 
        dfNitpb.loc['T6', 'coordinence'] = 51
    
        #Creation of the colormap
        custom_cmap = mcolors.LinearSegmentedColormap.from_list('custom', ['#FFFFFF', 'r'])
        
        norm = matplotlib.colors.Normalize(vmin=0, vmax=dfNitpb['epsilon'].max())
        cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=custom_cmap)

        dfNitpb['cmap'] = cmap
        dfNitpb['name'] = '$^{[5]}Ni$ : TBP'
        
        df=dfNitpb


    return(df)

def read_annot(Listeannot):
    '''
    Convertie une liste d'annotation en Liste d'élément, sépare les éléments par virgule et supprime les espaces.

    Parameters
    ----------
    Listeannot : TYPE
        DESCRIPTION.

    Returns
    -------
    Listelement : TYPE
        DESCRIPTION.

    '''
    try : Listelement =  Listeannot.replace(' ', '').split(',')
    except AttributeError : 
        try : Listelement =  Listeannot[0].replace(' ', '').split(',')
        except AttributeError : Listelement = Listeannot[0][0].replace(' ', '').split(',')
    return Listelement   

def OAS_annotate (X,Y, Listeannot='', ax=None, mode='cm', lissage=True, Annot_type = 1):
        
    if lissage : Y = savgol_filter(Y, 31, 3)
    
    Listelement = read_annot(Listeannot)
    
    for Element in Listelement :
       if Annot_type == 1 : #Annotation du maximum
       
        Limite=element2epsilon(Element)
           
        if Element == 'Fe2+':
            annot_max(X, Y, '$Fe^{2+}$', ax, xlim=Limite.iloc[0,0])
                      
        elif Element == 'Fe3+':
            annot_max(X, Y, '$Fe^{3+}$', ax, xlim=Limite.iloc[0,0],
                      xytext=(-10, 30))
        
            annot_max(X, Y, '$Fe^{3+}$', ax, xlim=Limite.iloc[1,0],
                      xytext=(10, 20))
            
            annot_max(X, Y, '$Fe^{3+}$', ax, xlim=Limite.iloc[2,0],
                      xytext=(10, 30))
        
        elif Element == 'Cu2+' :
            annot_max(X, Y, '$^{[6]}Cu^{2+}$', ax, xlim=Limite.iloc[0,0])
            
        elif Element == 'Co2+':
            Limite=element2epsilon(Element)
            
            annot_max(X, Y, '$Co^{2+}$', ax, xlim=Limite.iloc[0, 0], #IR
                      xytext=(-10, 30))
        
            annot_max(X, Y, '$Co^{2+}$', ax, xlim=Limite.iloc[1, 0],
                      xytext=(10, 40))
            
            annot_max(X, Y, '$Co^{2+}$', ax, xlim=Limite.iloc[2, 0],
                      xytext=(20, 30))
            
            annot_max(X, Y, '$Co^{2+}$', ax, xlim=Limite.iloc[3, 0], #VIS
                      xytext=(-40, 20))
        
            annot_max(X, Y, '$Co^{2+}$', ax, xlim=Limite.iloc[4, 0],
                      xytext=(20, 20))
            
            annot_max(X, Y, '$Co^{2+}$', ax, xlim=Limite.iloc[5, 0],
                      xytext=(10, 20))
            
        elif Element == 'Mn3+':
            
            annot_max(X, Y, '$Mn^{3+}$', ax, xlim=Limite.iloc[0,0])
     
        else :
            pass


def max2C(X, Y, Listeannot, Legende):
    try : Listelement =  Listeannot.replace(' ', '').split(',')
    except AttributeError : Listelement =  Listeannot[0].replace(' ', '').split(',')
    
    print(Listelement)
    
    dfinal = pd.DataFrame(columns=['Coordinence', 'X (cm^-1)',
                                                 'Ymax (cm^-1)', 'Epsilon',
                                                 'Concentration (mol.L^-1)'])
    
    for Element in Listelement :
        
        df = element2epsilon(Element)
        
        indexori=df.index
        
        indexreturn = [[Legende for x in indexori], [Element for x in indexori], indexori]
        
        dfreturn = pd.DataFrame(columns=['Coordinence', 'X (cm^-1)',
                                                     'Ymax (cm^-1)', 'Epsilon',
                                                     'Concentration (mol.L^-1)'],
                                index = indexreturn)
        dfreturn.index.names=["Verre", "Element", "Transition"]
    
        dfreturn.loc[:, 'Coordinence'] = df.loc[:, 'coordinence'].to_numpy()
        dfreturn.loc[:, 'Epsilon'] = df.loc[:, 'epsilon'].to_numpy()
        
        for i, xlim in enumerate(df['Xlim']) :
            epsilontemp = df.loc[indexori[i], 'epsilon']

            if not np.isnan(epsilontemp) :
                x=X
                y=Y
                
                limite= np.logical_and(x > xlim[0], x < xlim[1])
                x=x[limite]
                y=y[limite]
                
                xmax = x[np.argmax(y)]
                ymax = y.max()
                indextemp = (Legende, Element, indexori[i])
                             
                dfreturn.loc[indextemp, 'X (cm^-1)']                = xmax
                dfreturn.loc[indextemp, 'Ymax (cm^-1)']             = ymax
                dfreturn.loc[indextemp, 'Concentration (mol.L^-1)'] = ymax/epsilontemp
                # print(dfreturn)
        dfinal = pd.concat([dfinal, dfreturn], axis=1)
        
    print('\n\n\n')
    print(dfinal)
    print('\n\n\n')
    # print(dfinal)          
    
    return(dfinal)


#%% Effet de zoom
from matplotlib.transforms import Bbox, TransformedBbox, blended_transform_factory

from mpl_toolkits.axes_grid1.inset_locator import BboxPatch, BboxConnector,\
    BboxConnectorPatch
    
def connect_bbox(bbox1, bbox2,
                 loc1a, loc2a, loc1b, loc2b,
                 prop_lines, prop_patches=None):
    if prop_patches is None:
        prop_patches = prop_lines.copy()
        prop_patches["alpha"] = prop_patches.get("alpha", 1)*0.2

    c1 = BboxConnector(bbox1, bbox2, loc1=loc1a, loc2=loc2a, **prop_lines)
    c1.set_clip_on(False)
    c2 = BboxConnector(bbox1, bbox2, loc1=loc1b, loc2=loc2b, **prop_lines)
    c2.set_clip_on(False)

    bbox_patch1 = BboxPatch(bbox1, **prop_patches)
    bbox_patch2 = BboxPatch(bbox2, **prop_patches)

    p = BboxConnectorPatch(bbox1, bbox2,
                           loc1a=loc1a, loc2a=loc2a, loc1b=loc1b, loc2b=loc2b,
                           **prop_patches)
    p.set_clip_on(False)

    return c1, c2, bbox_patch1, bbox_patch2, p


def zoom_effect02(ax1, ax2, **kwargs):
    """
    ax2 : the big main axes
    ax1 : the zoomed axes
    The xmin & xmax will be taken from the
    ax1.viewLim.
    """

    tt = ax1.transScale + (ax1.transLimits + ax2.transAxes)
    trans = blended_transform_factory(ax2.transData, tt)

    mybbox1 = ax1.bbox
    mybbox2 = TransformedBbox(ax1.viewLim, trans)

    prop_patches = kwargs.copy()
    prop_patches["ec"] = "none"
    prop_patches["alpha"] = 0.2

    c1, c2, bbox_patch1, bbox_patch2, p = \
        connect_bbox(mybbox1, mybbox2,
                     loc1a=2, loc2a=3, loc1b=1, loc2b=4, 
                     prop_lines=kwargs, prop_patches=prop_patches)

    ax1.add_patch(bbox_patch1)
    ax2.add_patch(bbox_patch2)
    ax2.add_patch(c1)
    ax2.add_patch(c2)
    ax2.add_patch(p)

    return c1, c2, bbox_patch1, bbox_patch2, p


#def Aff_ref(nom : str, )

#%% Big data


from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.decomposition import PCA


def Liste2PCA(Liste, Legende, Xmin=330, Xmax=2500, spectro='OAS'):
    '''
    Cette fonction réalise l'anayse en composante principal des spectres dans la Liste

    Parameters
    ----------
    Liste : TYPE
        DESCRIPTION.
    Legende : TYPE
        DESCRIPTION.
    Xmin : TYPE, optional
        DESCRIPTION. The default is 330.
    Xmax : TYPE, optional
        DESCRIPTION. The default is 2500.
    spectro : TYPE, optional
        DESCRIPTION. The default is 'OAS'.

    Returns
    -------
    None.

    '''

    if spectro == 'OAS' : 
        Xref = np.arange(Xmin, Xmax, 1)
        Bloc = np.ones([np.size(Liste), np.size(Xref)])
        
        for i, item in enumerate(Liste) :
            X, Y = Readspectre(item)
            
            if X[0] > X[1] : 
                X = np.flip(X)
                Y = np.flip(Y)
            
            Y = -np.log10(Y)
            X, Y = removeInfNan(X, Y)
            
            Bloc[i,:] = np.interp(Xref, X, Y)
            
        pca_kernel = PCA(n_components=3)
        pca_kernel.fit_transform(Bloc)

        plt.plot(np.arange(1,11,1),pca_kernel.explained_variance_ratio_,"s-")
        plt.xlabel("Number of components")
        plt.ylabel("Explained variance")

        print(pca_kernel.explained_variance_ratio_)
    
    return(pca_kernel)


def Liste2dendrogram(Liste, Legende, method='single', optimal_ordering=False,
                     orientation='left', spectro='OAS', Xmin=330, Xmax=2500):
    
    if spectro == 'OAS' : 
        Xref = np.arange(Xmin, Xmax, 1)
        Bloc = np.ones([np.size(Liste), np.size(Xref)])
        
        for i, item in enumerate(Liste) :
            X, Y = Readspectre(item)
            
            if X[0] > X[1] : 
                X = np.flip(X)
                Y = np.flip(Y)
            
            Y = -np.log10(Y)
            X, Y = removeInfNan(X, Y)
            
            Bloc[i,:] = np.interp(Xref, X, Y)

        # plt.figure()
        # plt.plot(X,Y, 'x', label = 'ori')
        # plt.plot(Xref,Bloc[i,:], '.', label = 'inter')
        

    Z = linkage(Bloc, method, optimal_ordering=optimal_ordering)
    dn = dendrogram(Z, labels=Legende, orientation=orientation)
    
    return(dn)

def IpVsMaxnu(Liste, Legende, optplt=None, Marker = 'x', GRAPHOUT=5, DEBUG=False, legncol=2):
    '''Cette fonction trace l'indice de polymérisation en fonction de la
    position du maxium d'intensité du massif d'élongation '''
    
    fig, ax1 = set_graph(GRAPHOUT)
    ax1.set_xlim([900, 1150])
    ax1.set_ylim([0, 4])
    
    if DEBUG :  fi2, ax2 = plt.figure()
    
    for i, Fichier in enumerate(Liste) :
        try : 
            if not optplt[i] == None : fmt = optplt[i]
            else : fmt = ''
        except :
            fmt = ''
        
        try : 
            if not Marker[i] == None : mark = Marker[i]
            else : mark = 'x'
        except :
            mark = 'x'
        
        X, Y  = ReadData_spectrotype(Fichier, spectromode = 'RAMAN')
        
        # if DEBUG : ax2.plot(X,Y, label=)
        
        INDEX_bending = X < 820
        
        A_bending = np.trapz(Y[INDEX_bending], X[INDEX_bending])
        A_stretching = np.trapz(Y[~INDEX_bending], X[~INDEX_bending])
        
        Numax = X[~INDEX_bending][Y[~INDEX_bending].argmax()]
        
        try : test =  float(Legende[i]) # castage de la legende en float, si jamais c'est un nan 
        except : test = 0
        if  not np.isnan(test) : legg = Legende[i]
        else : legg = None
        
        plt.plot(Numax, A_bending/A_stretching, fmt, marker = mark, label = legg)
    ax1.set_ylabel('Indice de polymérisation')
    ax1.set_xlabel('$ \\nu_{max} $ élongation Si-O')
    leg1 = ax1.legend(ncol=legncol)
    leg1.set_draggable(True)
    
    pass


def DecnuVsMaxnu(Liste, Legende, optplt=None, Marker = 'x', GRAPHOUT=5, DEBUG=False, legncol=2, molpercentage = None):
    fig, ax1 = set_graph(GRAPHOUT)
    fig2, ax2 = set_graph(GRAPHOUT)
    # ax1.set_xlim([900, 1150])
    # ax1.set_ylim([0, 4])
    
    if DEBUG :  fi2, ax2 = plt.figure()
    
    for i, Fichier in enumerate(Liste) :
        try : 
            if not optplt[i] == None : fmt = optplt[i]
            else : fmt = ''
        except  :
            fmt = ''

        try : 
            if not Marker[i] == None : mark = Marker[i]
            else : mark = 'x'
        except  :
            mark = 'x'
        
        X, Y  = ReadData_spectrotype(Fichier, spectromode = 'RAMAN')
        
        # if DEBUG : ax2.plot(X,Y, label=)
        
        INDEX_streaching = X > 820

        INDEX_bending = X < 820
        
        A_bending = np.trapz(Y[INDEX_bending], X[INDEX_bending])
        A_stretching = np.trapz(Y[~INDEX_bending], X[~INDEX_bending])
        
        Ip = A_bending/A_stretching
        
        Numax = X[~INDEX_bending][Y[~INDEX_bending].argmax()]


        Ystrech = Y[INDEX_streaching]
        Xstrech = X[INDEX_streaching]
        
        INDEX_part1bending = Xstrech < 1000
        
        numax2 = Xstrech[~INDEX_part1bending][Ystrech[~INDEX_part1bending].argmax()] # position du max de la première contribution 
        numax1 = Xstrech[INDEX_part1bending][Ystrech[INDEX_part1bending].argmax()] # Position du max de la seconde contribution 
        
        Deltanu = numax2-numax1
        
        ax1.plot(Deltanu, Numax , fmt, marker = mark, label = Legende[i])# Ecart des deux bosse vs IP 
        
        try : ax2.plot(Deltanu, molpercentage[i], fmt, marker = mark, label = Legende[i])
        except IndexError : pass
    
    ax1.set_xlabel('Delta nu')
    ax1.set_ylabel('$ \\nu_{max} $ élongation Si-O')      
        
    ax2.set_ylabel('Paramètre x de $0.6SiO2xK_2O(0.4-x)CaO$')
    ax2.set_xlabel('$ \\Delta\\nu $ massif élongation Si-O')
    
    leg1 = ax1.legend(ncol=legncol)
    leg1.set_draggable(True)
    
    fig.tight_layout()
    fig2.tight_layout()
    
    pass

#%% Fonction dynamique
from matplotlib.patches import Circle
from matplotlib.widgets import LassoSelector

def rectangle_plot(ax, color='y', alpha=1) :
    '''
    Trace un rectangle sur un axe

    Parameters
    ----------
    ax : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    def onselect(eclick, erelease):
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata
        width = x2 - x1
        height = y2 - y1
        ax.add_patch(Rectangle((x1, y1), width, height, fill=False, color=color, alpha=alpha))
        plt.draw()
    
    
    rs = RectangleSelector(ax, onselect)
    plt.show()
    return (rs)


def circle_plot(ax, color='red', alpha='1'):
    def onselect(verts):
        if len(verts) > 0:
            x, y = zip(*verts)
            x_center = sum(x) / len(verts)
            y_center = sum(y) / len(verts)
            radius = max([(x - x_center) ** 2 + (y - y_center) ** 2 for x, y in verts]) ** 0.5
            ax.add_patch(Circle((x_center, y_center), radius, fill=False, color=color, alpha=alpha))
            plt.draw()
    
    lasso = LassoSelector(ax, onselect)
    return(lasso)

def annotate_clic_text(fig, ax): 
    '''
    Créer un annotation à l'endroit du clic gauche. Clic droit pour supprimer'

    Parameters
    ----------
    fig : TYPE
        DESCRIPTION.
    ax : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    annotations = []

    def onclick(event):
        if event.inaxes:
            x, y = event.xdata, event.ydata
            if event.button == 1:  # Clic gauche pour ajouter une annotation
                annotation = askstring("Saisir une annotation", "Entrez votre annotation :")
                if annotation:
                    annotations.append((x, y, annotation))
                    ax.annotate(annotation, (x, y), textcoords="offset points", xytext=(0, 10), ha='center')
                    plt.draw()
                    print(f"ax.annotate('{annotation}', ({x:.2f}, {y:.2f}), textcoords='offset points', xytext=(0, 10), ha='center')")
            elif event.button == 3:  # Clic droit pour supprimer une annotation
                for annotation in annotations:
                    ann_x, ann_y, _ = annotation
                    x_margin = 0.02 * (ax.transData.transform((1, 1))[0] - ax.transData.transform((0, 0))[0])  # 2% de la largeur du graphique
                    y_margin = 0.02 * (ax.transData.transform((1, 1))[1] - ax.transData.transform((0, 0))[1])  # 2% de la hauteur du graphique
                    if abs(x - ann_x) < x_margin and abs(y - ann_y) < y_margin:
                        annotations.remove(annotation)
                        for txt in ax.texts:
                            if txt.get_text() == annotation[2]:
                                txt.remove()
                        plt.draw()
    
    fig.canvas.mpl_connect('button_press_event', onclick)
    
    plt.show()



def click_print_arrow(fig, ax, color='k'):
    drawing = False
    start_point = None
    arrow = None
    
    def onpick(event):
        nonlocal drawing, start_point, arrow
        if event.button == 1:  # Clic gauche pour ajouter une flèche
            if not drawing:
                start_point = (event.xdata, event.ydata)
                drawing = True
            else:
                end_point = (event.xdata, event.ydata)
                arrow = FancyArrowPatch(start_point, end_point, arrowstyle='->', mutation_scale=10, color='k')
                ax.add_patch(arrow)
                plt.draw()
                drawing = False
                print(f"ax.add_patch(FancyArrowPatch(({start_point[0]:.2f}, {start_point[1]:.2f}), ({end_point[0]:.2f}, {end_point[1]:.2f}), arrowstyle='->', mutation_scale=10, color='{color}'))")
        elif event.button == 3:  # Clic droit pour supprimer la dernière flèche
            if arrow:
                arrow.remove()
                plt.draw()
                arrow = None
    
    fig.canvas.mpl_connect('button_press_event', onpick)
    
    plt.show()


    