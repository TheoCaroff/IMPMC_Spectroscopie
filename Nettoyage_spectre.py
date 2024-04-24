# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 11:17:45 2020

@author: Théo Caroff
        Avec ajout de Dongxin sur la fonction filtrage_gauss
"""

import numpy as np
import pandas as pd
import os
try : import rampy
except : pass

from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d
from scipy.signal import savgol_filter
from scipy import interpolate
from IPython import display
from scipy.spatial import ConvexHull

from Lecture_input import Readspectre
from Lecture_input import Corr2Str_OAS
from Lecture_input import Corr2folder_OAS
from Lecture_input import Corr2Str_RAMAN
from Lecture_input import Corr2Str_RPE
from Lecture_input import Corr2folder_RAMAN
from Lecture_input import normYminmax
from Lecture_input import nm2cm1
from Lecture_input import Gauss
from Lecture_input import removeInfNan
from Lecture_input import Sav_fig
from Lecture_input import Releve_parametre

def show_interactif_loop(inline=False):
    '''
    équivalent de plt.show()
    Sensé permettre l'affichage à la matlab dans une boucle

    Returns
    -------
    None.

    '''
    
    if not inline :
        leg = plt.legend()
        leg.set_draggable(state=True)
        display.clear_output(wait=True)
        display.display(plt.gcf())
        plt.show(block=False)
        plt.pause(0.2)
        
    else :
        plt.show()


def baseline_linear(X,Y,Legend) :
    #on fait des tronçons de Y pour touver les minima locaux
    
    x_0 = min(X)
    x_1 = 5000
    x_2 = 14500
    x_3 = 20000
   
    X1 = np.absolute(x_1-X)
    X2 = np.absolute(x_2-X)
    X3 = np.absolute(x_3-X)
    
    p_0 = 0 #premier element de la liste X
    p_1 = X1.argmin()
    p_2 = X2.argmin()
    p_3 = X3.argmin()
   
    
    Y1 = Y[p_0:p_1]
    Y2 = Y[p_2:p_3]
    #X2 = X[p_2:p_3]
    
    min1 = np.min(Y1)
    #min2 n'est pas forcément le minimum car la diffusion rend la courbe toujours croissante : min2 = pt d'inflexion
    #deriv_sec = np.gradient(np.gradient(Y2,X2))
    #print(deriv_sec)
   # i_2 = np.where(deriv_sec==0)       
    min2 = np.min(Y2)
    i_1 = np.where(Y==min1)
    i_2 = np.where(Y==min2)
    
    Y_base = Y[i_1] + (X-X[i_1])*(( Y[i_2]-Y[i_1])/(X[i_2]-X[i_1]))
    
    plt.plot(X,Y)
    plt.plot(X, Y_base)

    Y_corr = Y-Y_base
    
    Y2=Y_corr[p_2:p_3]
    mini = np.min(Y2)
    pos_min = np.where(Y_corr==mini)
    #X_min = float(X[pos_min])
    #try :
       # text_file=open('Position_min.csv', "x")
        #ENTETE='Nom fichier;Position_min\n'
        #text_file.write(ENTETE)
    #except FileExistsError :
        #text_file=open('Position_min.csv', "a")
    
    #DATA= Legend + ';' + str(X_min) + '\n'
    
   # text_file.write(DATA)
    
    #text_file.close()
    
    return (X, Y_corr)
    

def baseline_Rubberband(X, Y): # Basé sur la méthode de rubberband (voir package hyperspectra pour R)
    # Find the convex hull$
    INDEXrub=~np.isnan(Y)
    X=X[INDEXrub]
    Y=Y[INDEXrub]
    
    array=np.concatenate((X[np.newaxis,:], Y[np.newaxis,:])).transpose()
    v = ConvexHull(array).vertices
    # Rotate convex hull vertices until they start from the lowest one
    v = np.roll(v, -v.argmin())
    # Leave only the ascending part
    v = v[:v.argmax()]

    # Create baseline using linear interpolation between vertices
    return (X[v], Y[v])

def fit_OMCT(X, Y, NOM, removenaninf=True): 
    '''
    Cette fonction sert à fiter un Gausienne
    Dans le cas d'un transfert de charge O-élément d.
    L'AFFICHAGE DE MATPLOLIB DOIT ETRE INLINE

    Parameters
    ----------
    X : float
        Nombre d'onde des données à ajuster.
    Y : float
        absorbance des donnés à ajuster.
    NOM : string
        nom de l'échantillon.

    Returns
    -------
    Y_coor, float qui est la soustration du Y original et du fit.

    '''
    
    if removenaninf : X,Y = removeInfNan(X, Y)
    
    nbonde_min_fit=27000; # limite basse de la partie des data qui vont être ajuster
    nbonde_max_fit=32000; # limite haute
    
    Intes_fit=10; # intensité de la gausienne pour l'ajustement
    Locali_fit=34000; #position central
    Sigma_fit = 1500; # écart type Gauss
    Shift_fit = 1; # shift de la gausienne
    
    X_ori=X;
    Y_ori=Y;
    
    #plt.ion() #Nécéssaire pour afficher les figures en %matplolib  
    while True:
        nbonde_min=nbonde_min_fit;
        nbonde_max=nbonde_max_fit;
        
        p=[Intes_fit, Locali_fit, Sigma_fit, Shift_fit] # on récupère les para de base du fit.
        
        limite_haute = (X_ori <= nbonde_max);
        X=X_ori[limite_haute]
        Y=Y_ori[limite_haute]
        
        #Y_min = np.min(Y); #On récupère la valeur minimum du spectre qui à un sens physique pour bien assoir la ligne de base.
        
        limite_basse = (X >= nbonde_min); ## On ne  sélection que les point dans la zone de fit.
        X=X[limite_basse]
        Y=Y[limite_basse]
        
        plt.figure(figsize=(5,3), dpi=120)
        
        plt.plot(X,Y, 'x',label='ref')
        
        plt.plot(X,Gauss(X, p[0], p[1], p[2], p[3]),'x',label='avant affinement')
        plt.legend();
        plt.grid();
        
        try:
            pop, cov = curve_fit(Gauss, X, Y, p0=p)
            
            Y_fit=Gauss(X, pop[0], pop[1], pop[2], pop[3]);
            
            plt.plot(X, Y_fit, '--', label='fit1G')
            plt.legend();
            
            plt.figure(figsize=(5,3), dpi=120)
            Y_fit=Gauss(X_ori, pop[0], pop[1], pop[2], 0); # on ne soustrait pas le shift
            #Y_fit=Gauss(X, pop[0], pop[1], pop[2], 0)
            
            Y_corr=Y_ori-Y_fit;
            plt.plot(X_ori,Y_ori,label='Origine')
            plt.plot(X_ori, Y_corr, '--', label='corrigé')
            plt.xlim(4000, 33000);
            plt.ylim(-0.1, 2)
            plt.legend();
            plt.grid();
            show_interactif_loop()
                
        except RuntimeError:
            show_interactif_loop()
                
            print('\nCurve fit n\'a vraisemblablement pas convergé\n'.format(RuntimeError));
            print('Essayer de changer les paramètres du fit\n')
            fit_OK=False;
        else:
            print(NOM + '  :   l\'ajustement est-il correct ?\n');
            try:
                fit_OK = int(input('Rentrez 0 si non, 1 si l\'ajustement est correct\n'))
            except ValueError:
                print('Vous n\'avez pas rentré un nombre, le fit est considéré comme mauvais')
                fit_OK=False;
        if fit_OK : break
        #print(pop)
        nbonde_min_fit = float(input('Valeur min nb onde en cm-1 (default = '+
                                 str(nbonde_min_fit) + ')\n') or nbonde_min_fit);
        nbonde_max_fit = float(input('Valeur max nb onde en cm-1 (default = '+
                                 str(nbonde_max_fit) + ')\n') or nbonde_max_fit);
        Intes_fit = float(input('Valeur intensité gauss (default = '+
                                 str(Intes_fit) + ')\n') or Intes_fit);
        Locali_fit = float(input('Valeur central gauss (default = '+
                                 str(Locali_fit) + 'cm-1)\n') or Locali_fit);
        Sigma_fit = float(input('Valeur sigma gauss (default = '+
                                 str(Sigma_fit) + 'cm-1 )\n') or Sigma_fit);
        Shift_fit = float(input('Valeur constance (default = '+
                                 str(Shift_fit) + 'cm-1 )\n') or Shift_fit);
        
    #return(Y_corr-Y_min) # On récupère la valeur minmum pour bien assoir la ligne de base 
    #! pas de sens physique, valeur de correction de reflexion à mesure
    print("Paramètre de fit  : ")
    print(pop)
    
    try :
        text_file=open('parametre_fit_gauss.csv', "x")
        ENTETE='Nom fichier;intensité gaussienne;centre gaussienne;sigma gaussienne;constante\n'
        text_file.write(ENTETE)
    except FileExistsError :
        text_file=open('parametre_fit_gauss.csv', "a")

    DATA=NOM + ';' + str(pop[0]) + ';' + str(pop[1]) + ';' + str(pop[2]) + ';' + str(pop[3])+'\n'
    
    text_file.write(DATA)
    
    text_file.close()
    return(X_ori, Y_corr) #On retourne le spectre corrigé

def lissage(X, Y, NOM='', sigma=2, selectmethode=2):
    '''
    Cette fonction sert à lisser des courbes, elle propose le choix entre un
    lissage gaussien (selectmethode=1) et Savgol (=2)

    Parameters
    ----------
    X : TYPE
        DESCRIPTION.
    Y : TYPE
        DESCRIPTION.
    NOM : TYPE, optional
        DESCRIPTION. The default is ''.
    sigma : TYPE, optional
        Valeur par defaut du filtre gaussien. The default is 2.
    selectmethode : TYPE, optional
        Pour choisir le type de lissage, gaussien =1, Savgol=2. The default is 2.

    Returns
    -------
    None.

    '''
    
    if np.nanmin(Y) < 0 : # Réglage des borne d'affichage en Y
        Amin = np.nanmin(Y)
        Amin = Amin + 0.2*Amin
    else :
        Amin = 0
    
    Amax = np.nanmax(X[X<33000])
    Amax = Amax+0.2*Amax
    if Amax == np.inf : Amax=4
    
    selectmethode = selectmethode;
    
    Xmin=4000
    Xmax=40000
    Zoom=0;
    zone_min=500;
    zone_max=10000;
    
    ordre_savgol = 3;
    nbpt_savgol = 31;
    
    X,Y = removeInfNan(X, Y)
    Y_fit = savgol_filter(Y, nbpt_savgol, ordre_savgol)
    Xfit = X
    
    #plt.ion() #Nécéssaire pour afficher les figures en %matplolib
    plt.figure(figsize=(5,3), dpi=180)
    plt.plot(X, Y, label='Original')
    plt.plot(X, Y_fit, ':', label='Lissé')
    plt.grid()
    plt.legend()
    plt.show()
    
    V1_OK = int(input('Est-ce que le lissage vous convient ? 1=oui, 0=non (default = 1) \n') or 1)
    
    if not V1_OK :
        selectzone = int(input('Voulez-vous selectionner une zone ou lisser toutes la courbe ?'
                               +'0 = une zone, 1 = toute la zone (default = 1) \n') or 1)
        while True: 
            selectmethode= int(input(' Selectionner le lissage 1=gaussien, 2= Savgol (default = '
                                     + str(selectmethode) +')\n') or selectmethode)
            
            if selectmethode == 1 :
                sigma = float(input('Rentrer la valeur sigma du filtre gaussien (default = '+
                                         str(sigma) + ')\n') or sigma);
            if selectmethode == 2 :
                ordre_savgol = int(input('Rentrer l\'ordre du polynome de lissage, 0 pour moyenne glissante (default = '+
                                         str(ordre_savgol) + ')\n') or ordre_savgol);
                nbpt_savgol = int(input('Rentrer le nombre de point de la zone de lissage glissante !IMPAIRE! (default = '+
                                         str(nbpt_savgol) + ')\n') or nbpt_savgol);
    
            if selectzone == 0 :
                zone_min=float(input('Nombre d\'onde: début de la zone (default='+ str(zone_min) + ') \n ') or zone_min);
                zone_max=float(input('Nombre d\'onde: fin de la zone (default='+ str(zone_max) + ') \n ') or zone_max);
                
                zone=X<=(zone_min) # on sépare les data inférieur à zone_min
                XR1=X[zone]
                YR1=Y[zone]  
                zone = np.logical_and(X > zone_min, X < zone_max) # on prend l'intervalle où on souhaite lisser
                XR2=X[zone]
                YR2=Y[zone]  
                zone=X>=(zone_max) # on sépare les data au dessus de zone_max
                XR3=X[zone]
                YR3=Y[zone]
                
                if selectmethode == 1 : Y_fit = gaussian_filter1d(YR2,sigma) # on filtre la partie qu'on souhaite
                if selectmethode == 2 : Y_fit = savgol_filter(YR2, nbpt_savgol, ordre_savgol)
                Xfit=XR2;
                
                
                X1 = np.concatenate([XR1,XR2]) # on rassemble les valeurs XR1 et XR2
                Y1 = np.concatenate([YR1,Y_fit])
                Y_fit = np.concatenate([Y1,YR3]) # on rassemble l'ensemble des valeurs Y
                Xfit = np.concatenate([X1,XR3]) # on rassemble l'ensemble des valeurs X
                
            if selectzone == 1 :
                if selectmethode==1 : Y_fit = gaussian_filter1d(Y,sigma) # on lisse sur l'ensemble 
                elif selectmethode == 2 : Y_fit = savgol_filter(Y, nbpt_savgol, ordre_savgol)
                
                Xfit=X;
                
            plt.figure(figsize=(5,3), dpi=180)
            plt.plot(X, Y, label='Original')
            plt.plot(Xfit, Y_fit, ':', label='Lissé')
            
            if Zoom == 1 :
                plt.ylim([Amin, Amax])
                plt.xlim([Xmin, Xmax])
            plt.grid()
            plt.legend()
            show_interactif_loop()
                    
            print(NOM + '  :   le filtrage est-il correct ?\n');
           
            try:
                fit_OK = int(input('Rentrez 0 si non, 1 si oui\n'))
            except ValueError:
                print('Vous n\'avez pas rentré un nombre, le fit est considéré comme mauvais')
                fit_OK=False;
            if fit_OK : break
           
            Zoom=int(input("Voulez vous zoomer sur une zone ? (0 non, 1 oui)\n") or 0)
           
            if Zoom :
                Xmin = float(input('Valeur min nb onde en cm-1 (default = '+
                                         str(Xmin) + ')\n') or Xmin);
                Xmax = float(input('Valeur max nb onde en cm-1 (default = '+
                                         str(Xmax) + ')\n') or Xmax);
                Amin = float(input('Valeur min abs (default = '+
                                         str(Amin) + ')\n') or Amin);
                Amax = float(input('Valeur max abs (default = '+
                                         str(Amax) + ')\n') or Amax);
                

    return(Xfit, Y_fit)


def Remontage_IR_VIS(X_BAS, Y_BAS, X_HAUT, Y_HAUT, mode='Perkin', X=0,Y=0,
                     NOM='pouet', VISsurNIR=False):
    '''
    Cette fonction sert à remonter deux partie de courbe au même niveau
    car issu de deux détecteur/spectro différent, par défaut on corrige l'IR sur le VIS
    Haut et Bas par rapport au nombre d'onde.
    Parameters
    ----------
    X_BAS : TYPE
        DESCRIPTION.
    Y_BAS : TYPE
        DESCRIPTION.
    X_HAUT : TYPE
        DESCRIPTION.
    Y_HAUT : TYPE
        DESCRIPTION.
    mode : TYPE, optional
        Perkin ou spectroportable (choix de limite suplémentaire). The default is 'Perkin'.
    X : TYPE, optional
        X originel si data coupée. The default is 0.
    Y : TYPE, optional
        Y originel si data coupée. The default is 0.
    NOM : TYPE, optional
        DESCRIPTION. The default is 'pouet'.

    Returns
    -------
    Y_BASDelta : TYPE
        Y remonté.
    interpo_limBAS : TYPE
        limite pour les IR.
    interpo_limHAUT : TYPE
        limite pour le VIS.

    '''
    #plt.ion() #Nécéssaire pour afficher les figures en %matplolib  
    FIGsize=(10,6)
    DPI=120
    

    Xmin=5000; # Affichage en x min en cm-1
    Xmax=20000;# Affichage en x max en cm-1

        
    if mode == 'PERKIN_micro':
        limIR= 9000;   # limite pour l'ajustement IR
        limVIS = 14000; # 
        point_jonction = 11000;
    
    elif mode == 'PERKIN_std':
        limIR= 11000;   # limite pour l'ajustement IR
        limVIS = 12000; # 
        point_jonction = 11627;
        Xmin = 9000; # Affichage en x min en cm-1
        Xmax = 14000;# Affichage en x max en cm-1
    
    elif mode == 'PERKIN_lampe': # A FINIRRRR
        limIR= 25000;   # limite pour l'ajustement IR
        limVIS = 28000 ; # 
        point_jonction = 26453;
        Xmin = 22000; # Affichage en x min en cm-1
        Xmax = 32000;# Affichage en x max en cm-1
        VISsurNIR = True;
        
    elif mode == 'PERKIN_filtre_319':
        limIR = 30000;   # limite pour l'ajustement IR
        limVIS = 33000; 
        Xmin = 28000; # Affichage en x min en cm-1
        Xmax = 34000;# Affichage en x max en cm-1 
        point_jonction = 31347;
        VISsurNIR = True;
        
    elif mode == 'PortableIR':
        limIR= 8000;   # limite pour l'ajustement IR
        limIR_haute = 10000;
        limVIS = 11000; # 
        limVIS_basse = 9800;
        X = np.concatenate([X_BAS, X_HAUT]);
        Y = np.concatenate([Y_BAS, Y_HAUT]);
        point_jonction = 9800;
   
    elif mode == 'PortableUV':
        limVIS_bas= 24000;   # limite pour l'ajustement UV
        limVIS_haut = 26000; # 
        limUV_bas = 25500;
        limUV_haut = 28000;
        point_jonction = 25500;
        X = np.concatenate([X_BAS, X_HAUT]);
        Y = np.concatenate([Y_BAS, Y_HAUT]);
        Xmin=20000; # Affichage en x min en cm-1
        Xmax=34000;# Affichage en x max en cm-1
    # print('\n\n\n X,Y = ')    
    # print(X,Y)    
    # print('\n\n\n')
    
    Borne=np.logical_and(X>Xmin, X<Xmax)
    
    if np.nanmin(Y[Borne]) < 0 : # Réglage des borne d'affichage en Y
        Ymin = np.nanmin(Y[Borne])
        Ymin = Ymin + 0.2*Ymin
    else : Ymin = 0
        
    Ymax = np.nanmax(Y[Borne])
    Ymax = Ymax+0.2*Ymax
    
    if Ymax == np.inf : Ymax=4
    
    plt.figure(figsize=FIGsize, dpi=DPI)
    plt.plot(X, Y, label=NOM)
    plt.xlim([Xmin, Xmax])
    plt.ylim([Ymin, Ymax])
    plt.grid()
    show_interactif_loop()
    
    while True: # Selection de la zone d'interpolation
        
        if (mode == 'PERKIN_std' or mode == 'PERKIN_micro' or mode == 'PERKIN_lampe'
            or mode =='PERKIN_filtre_319') :
            limIR=int(input('Rentrer la limite partie basse IR pour l ajustement (defaut ='
                               + str(limIR) +')\n') or limIR);
            limVIS=int(input('Rentrer la limite partie haute VIS pour l ajustement (defaut ='
                               + str(limVIS) +')\n') or limVIS);
            
            interpo_limBAS = X_BAS > limIR
            interpo_limHAUT = X_HAUT < limVIS 

        if mode == 'PortableIR' :
             limIR=int(input('Rentrer la limite partie basse IR pour l ajustement (defaut ='
                               + str(limIR) +')\n') or limIR);
             limIR_haute=int(input('Rentrer la limite partie HAUTE IR pour l ajustement (defaut ='
                           + str(limIR_haute) +')\n') or limIR_haute);
             limVIS_basse=int(input('Rentrer la limite partie BASSE VIS pour l ajustement (defaut ='
                           + str(limVIS_basse) +')\n') or limVIS_basse);
             limVIS=int(input('Rentrer la limite partie haute VIS pour l ajustement (defaut ='
                               + str(limVIS) +')\n') or limVIS);
             interpo_limBAS = np.logical_and(X_BAS > limIR,  X_BAS < limIR_haute)
             interpo_limHAUT = np.logical_and(X_HAUT < limVIS, X_HAUT > limVIS_basse)
             point_jonction = limVIS_basse
             
        if mode == 'PortableUV' :
             limUV_bas=int(input('Rentrer la limite partie basse UV pour l ajustement (defaut ='
                               + str(limUV_bas) +')\n') or limUV_bas);
             limUV_haute=int(input('Rentrer la limite partie haute UV pour l ajustement (defaut ='
                           + str(limUV_haut) +')\n') or limUV_haut);
             limVIS_basse=int(input('Rentrer la limite partie basse VIS pour l ajustement (defaut ='
                               + str(limVIS_bas) +')\n') or limVIS_bas);
             limVIS_haute=int(input('Rentrer la limite partie haute VIS pour l ajustement (defaut ='
                           + str(limVIS_haut) +')\n') or limVIS_haut);     
             
             interpo_limHAUT = np.logical_and(X_HAUT > limUV_bas, X_HAUT  < limUV_haute)
             interpo_limBAS = np.logical_and(X_BAS < limVIS_haute, X_BAS > limVIS_basse)
             point_jonction = limVIS_haute
             
        plt.figure(figsize=FIGsize, dpi=DPI)
        
        plt.plot(X, Y, label='Original')
        plt.plot(X_BAS[interpo_limBAS], Y_BAS[interpo_limBAS], label='Y_BAS_ajustement');
        plt.plot(X_HAUT[interpo_limHAUT], Y_HAUT[interpo_limHAUT], label='Y_HAUT_ajustement');
        plt.xlim([Xmin, Xmax])
        plt.ylim([Ymin, Ymax])
        plt.grid()
        plt.legend()
        show_interactif_loop()
                
        if(int(input('Est-ce correcte ?\n'
                     + '1 pour oui et continuer la correction, 0 pour non et recommencer \n') or 1)) : break
    
    Remonte_OK=0;
    Remonte = 3; # Par défaut fit des deux côtés
    
    # On mets la premièrer partie de la courbe au niveau de la seconde,
    # car le Photomultiplicateur du visible donne un valeur "vrai" et on corrige l'InGaAs dessus
    
    Delta=np.mean(Y_HAUT[:3])-np.mean(Y_BAS[-3:]) # Decallage entre les deux courbe pour le cas 2
    ordre_fit = 2

    while (Remonte_OK == 0):  #Partie ou on remonte la partie gauche (les IR) des données pour coller l'InGaAs sur le PM  
        
        Remonte = int(input('Selectionner la methode de remonté de l\InGaAs sur le PM :'+'(default'+str(Remonte)+')\n'
                      +'\n1 fit polynomiale à droite et moyenne 3 dernier point gauche\n'
                      +' 2 réhausser à la main\n'
                      +' 3 fit des deux cotés du même ordre (ordre 0 = valeur moyenne)\n')
                        or Remonte)
        
        if(Remonte == 1):
            try:                   
                ordre_fit=float(input('rentrer l\'ordre du polynome qui va servir à ajuster la remontée (default = '
                                      +str(ordre_fit)+')\n') or ordre_fit)
    
                fit_delta = np.polyfit(X_HAUT[interpo_limHAUT], Y_HAUT[interpo_limHAUT], ordre_fit) #fit pour remonter la courbe
                
                regY_HAUT= np.poly1d(fit_delta)
                
                #Delta=regY_HAUT(X_BAS[-1])-Y_BAS[-1]; #
                Delta=regY_HAUT(np.mean(X_BAS[-3:]))-np.mean(Y_BAS[-3:]); #
                
                #Y_BASDelta=Y_BAS+Delta
                plt.figure(figsize=FIGsize, dpi=DPI)
                plt.plot(X_HAUT, regY_HAUT(X_HAUT), '-', label='fit pour la remonté')
            
            except LinAlgError :
                print('Le fit ne converge pas')
        
        elif(Remonte == 2):

            Delta=float(input('Rentrer le Delta entre les deux partie (defaut =' + str(Delta) + ')\n') or Delta)
            #Y_BASDelta=Y_BAS+Delta
            plt.figure(figsize=FIGsize, dpi=DPI)
            
        elif(Remonte == 3) :
            try:
                ordre_fit=float(input('rentrer l\'ordre du polynome qui va servir à ajuster la remontée (default = '
                       +str(ordre_fit)+')\n') or ordre_fit)
 
                if ordre_fit !=0:
                    point_jonction=float(input('Rentrer la valeur en X du point de jonction (default = '
                                               + str(point_jonction)+ ')\n') or point_jonction)
                
                fit_VIS = np.polyfit(X_HAUT[interpo_limHAUT], Y_HAUT[interpo_limHAUT], ordre_fit) #fit pour remonter la courbe
                regY_HAUT= np.poly1d(fit_VIS)
                
                fit_IR = np.polyfit(X_BAS[interpo_limBAS], Y_BAS[interpo_limBAS], ordre_fit) #fit pour remonter la courbe
                regY_BAS= np.poly1d(fit_IR)
                
                if ordre_fit == 0:
                    Delta=regY_HAUT(X_BAS[-1])-regY_BAS(X_HAUT[-1]); #
                else:
                     Delta=regY_HAUT(point_jonction)-regY_BAS(point_jonction); #
              
                #Y_BASDelta=Y_BAS+Delta
                
                plt.figure(figsize=FIGsize, dpi=DPI)
                plt.plot(X_HAUT, regY_HAUT(X_HAUT), '--', label='fit VIS pour la remonté')
                plt.plot(X_BAS, regY_BAS(X_BAS), '--', label='fit IR pour la remonté')
            except:
                print("Ne marche pas, utilise une autre methodes\n")
                #Y_BASDelta=Y_BAS
        plt.plot(X, Y, label=NOM + 'Original')
        
        if VISsurNIR : # On remonte le VIS sur le NIR, pas commun
            Y_BASDelta=Y_BAS
            Y_HAUTDelta = Y_HAUT - Delta
            plt.plot(X_BAS, Y_BAS, '-', label='1ere partie')
            plt.plot(X_HAUT, Y_HAUTDelta, '-', label='2eme partie remontée')
        
        else:  # Classique, on corrige le NIR sur le VIS.
            Y_BASDelta=Y_BAS+Delta
            Y_HAUTDelta = Y_HAUT
            plt.plot(X_BAS, Y_BASDelta, '-', label='1ere partie remontée')
            plt.plot(X_HAUT, Y_HAUT, '-', label='2eme partie')


        plt.xlim([Xmin, Xmax])
        plt.ylim([Ymin, Ymax])
        plt.grid()
        plt.legend()
        show_interactif_loop()
        
        Remonte_OK = int(input('La correction est-elle correcte ? \n') or 1)
 
    return (Y_HAUTDelta, Y_BASDelta, interpo_limBAS, interpo_limHAUT)

def reconstruction_saut(X_IR, Y_IRDelta, X_VIS, Y_VIS, interpo_limIR,
                        interpo_limVIS, Xinterpo, X, Y, NOM, mode='Perkin_micro'):
    '''
    Permet de reconstruire des données trop bruitée et qui ont étés supprimées.

    Parameters
    ----------
    X_IR : TYPE
        DESCRIPTION.
    Y_IRDelta : TYPE
        DESCRIPTION.
    X_VIS : TYPE
        DESCRIPTION.
    Y_VIS : TYPE
        DESCRIPTION.
    interpo_limIR : TYPE
        DESCRIPTION.
    interpo_limVIS : TYPE
        DESCRIPTION.
    Xinterpo : TYPE
        DESCRIPTION.
    X : TYPE
        X originel.
    Y : TYPE
        Y originel.
    NOM : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    fmt = '-'
    
    if mode == 'Perkin_micro': 
        Xmin=5000; # Affichage en x min en cm-1
        Xmax=20000;# Affichage en x max en cm-1

    elif mode == 'Perkin_lampe':
        Xmin = 22000
        Xmax = 32000
        fmt = 'x'
        
    elif mode == 'PERKIN_filtre_319':
        Xmin = 28000; # Affichage en x min en cm-1
        Xmax = 34000;# Affichage en x max en cm-1 *
        fmt = '-x'
    Borne=np.logical_and(X>Xmin, X<Xmax)
     
    if np.nanmin(Y[Borne]) < 0 : # Réglage des borne d'affichage en Y
        Ymin = np.nanmin(Y[Borne])
        Ymin = Ymin + 0.2*Ymin
    else : Ymin = 0
        
    Ymax = np.nanmax(Y[Borne])
    Ymax = Ymax+0.2*Ymax
    
    if Ymax == np.inf : Ymax=4

    Forme_OK = 0;
    Forme = 1;
    #plt.ion() #Nécéssaire pour afficher les figures en %matplolib   

    Forme = int(input('Selectionner la methode d\'interpolation\n 1 pour un ajustement polynomiale\n'
                      +' 2 pour un ajustement poly sur le dernier point de chaque coubre\n')
                        or Forme)
    FIGsize=(10,6)
    DPI=120

    while (Forme_OK == 0): # Partie ou on selectionne l'interpolation
        if (Forme == 1): # si jamais besoin de mettre une boucle

            ordre_fit = 5
            
            # On mets la premièrer partie de la courbe au niveau de la seconde,
            # car le Photomultiplicateur du visible donne un valeur "vrai" et on corrige l'InGaAs dessus
            
            ordre_fit=float(input('rentrer l\'ordre du polynome qui va servir à faire l\'interpolation (defaut: '+str(ordre_fit)+'\n')
                                or ordre_fit)
   
            Xtemp=np.concatenate([X_IR[interpo_limIR], X_VIS[interpo_limVIS]])
            # On rassemble les deux jeux de donner pour interpoler
            Ytemp=np.concatenate([Y_IRDelta[interpo_limIR], Y_VIS[interpo_limVIS]]) 
            
            fit_interpo = np.polyfit(Xtemp, Ytemp, ordre_fit) # fit pour interpoler les data
            fonctionYinterpo = np.poly1d(fit_interpo)
            Yinterpo = fonctionYinterpo(Xinterpo);
            
        elif (Forme == 2): # meme que 1 mais on prend les 2 derniers points de chaque courbe

            ordre_fit = 0
            ordre_fit=float(input('rentrer l\'ordre du polynome qui va servir à faire l\'interpolation\n')
                                or ordre_fit)
   
            Xtemp=np.concatenate([X_IR[interpo_limIR][:-1], X_VIS[interpo_limVIS][:1]])
            # On rassemble les deux jeux de donner pour interpoler
            Ytemp=np.concatenate([Y_IRDelta[interpo_limIR][:-1], Y_VIS[interpo_limVIS][:1]]) 
            
            fit_interpo = np.polyfit(Xtemp, Ytemp, ordre_fit) # fit pour interpoler les data
            fonctionYinterpo = np.poly1d(fit_interpo)
            Yinterpo = fonctionYinterpo(Xinterpo);

        elif (Forme == 3):
           Xtemp=np.concatenate([X_IR[interpo_limIR], X_VIS[interpo_limVIS]])
           Ytemp=np.concatenate([Y_IRDelta[interpo_limIR], Y_VIS[interpo_limVIS]]) 
           
           def sinfit(x, a, b, c):
               return a * np.sin(b * x) + c

           p, params_covariance = curve_fit(sinfit, Xtemp, Ytemp)
                                               #,p0=[2, 2])

           Yinterpo = sinfit(Xinterpo, p[0], p[1], p[2])
           print(p)
           print(Yinterpo)
           Yfit = sinfit(Xtemp, p[0], p[1], p[2])
    
        else:
            Yinterpo=np.ones(np.size(Xinterpo))*-1;
            print("Correction non implémentée\n")
            break

        
        plt.figure(figsize=FIGsize, dpi=DPI)
        plt.plot(X, Y, label= NOM + 'Original')
        plt.plot(X_IR, Y_IRDelta, '-', label='1ere parti remonté')
        plt.plot(X_VIS, Y_VIS, '-', label='2eme parti')
        try : plt.plot(Xtemp, Yfit, '--', label='fit')
        except : pass
        plt.plot(Xinterpo, Yinterpo, fmt, label='interpolation')
        plt.xlim([Xmin, Xmax])
        plt.ylim([Ymin, Ymax])
        plt.grid()
        plt.legend()
        show_interactif_loop()
        
        if(int(input('L\'interpolation est-elle correcte ?\nO pour non, 1 pour oui\n') or 1)): break
    
    Ycorr=np.concatenate([Y_IRDelta, Yinterpo, Y_VIS])

    return(Ycorr)
    
def correction_saut_detect(X, Y, NOM, mode='PERKIN_std'):
    '''
    Cette fonction sert à corriger le saut de detecteurdes PERKIN,
    plusieurs type d'ajustement sont proposé'

    Parameters
    ----------
    X : tableau de float
        en nombre d'onde.
    Y : tableau de float
        En absorbance.
    NOM : str
        Nom de ce qui est corrigé.
    mode : str, PERKIN_std ou PERKIN_micro
        choix du mode de travail pour la correction du saut
        Pour la correction de la lampe l'interpolation se fait automatiquement en ordre 3'

    Returns
    -------
    Ycorr : TYPE
        Correction.
    '''
    kwarg = {}
    fmt=''
    
    if mode == 'PERKIN_std':
        Xmin = 9000; # Affichage en x min en cm-1
        Xmax = 14000;# Affichage en x max en cm-1
        Xdebcoup = 11615; # Nombre d'onde en cm^-1 à partir duquel on suprime les valeurs du saut
        Xfincoup = 11615; # Nombre d'onde en cm^-1 jusqu'auquel on suprime les valeurs du saut

        
    elif mode == 'PERKIN_micro':
        Xmin = 5000; # Affichage en x min en cm-1
        Xmax = 20000;# Affichage en x max en cm-1
        Xdebcoup = 11500; # Nombre d'onde en cm^-1 à partir duquel on suprime les valeurs du saut
        Xfincoup = 12500; # Nombre d'onde en cm^-1 jusqu'auquel on suprime les valeurs du saut
    
    elif mode == 'PERKIN_lampe':
        Xmin = 22000; # Affichage en x min en cm-1
        Xmax = 32000;# Affichage en x max en cm-1  
        Xdebcoup = 26453; # Nombre d'onde en cm^-1 à partir duquel on suprime les valeurs du saut
        Xfincoup = 26470; # Nombre d'onde en cm^-1 jusqu'auquel on suprime les valeurs du saut
        fmt = '-x'
        
    elif mode == 'PERKIN_filtre_319':
        Xmin = 28000; # Affichage en x min en cm-1
        Xmax = 34000;# Affichage en x max en cm-1  
        Xdebcoup = 31280; # Nombre d'onde en cm^-1 à partir duquel on suprime les valeurs du saut
        Xfincoup = 31620; # Nombre d'onde en cm^-1 jusqu'auquel on suprime les valeurs du saut
        fmt = 'x'
    Borne=np.logical_and(X>Xmin, X<Xmax)
    Ymin = 0
    Ymax = Y[Borne]
    Ymax = np.nanmax(Ymax[Ymax != np.inf])
    
    Ymax = Ymax+0.2*Ymax

    
    
    sep_OK = False;
    #plt.ion() #Nécéssaire pour afficher les figures en %matplolib  
    # FIGsize=(10,6)
    # DPI=120

    FIGsize=[10,3.5]
    DPI=200


    while True : # DansYR2[interpo_limVIS]] cette partie on suprime les DATA qui déconne et on sépare en deux le spectre
        while True : # Partie ou on suprimer le saut de decteur
            RES=X<(Xdebcoup) # On sépare les data  au dessus de Xdebcoup
            XR1=X[RES]
            YR1=Y[RES]
            
            RES_1=RES
            
            RES=X>(Xfincoup) # Jusqu'à Xfincoup
            XR2=X[RES]
            YR2=Y[RES]
            
            RESinter = np.invert(RES+RES_1);
            Xinterpo = X[RESinter] #valeur en X des points à recréer
            
            sep_OK = True;
            
            plt.figure(figsize=FIGsize, dpi=DPI)
            plt.plot(X, Y, 'C3', label='Original')
            plt.plot(XR1, YR1, '-C0', label='1ere parti')
            plt.plot(XR2, YR2, '-C2', label='2eme parti')
            plt.xlim([Xmin, Xmax])
            plt.ylim([Ymin, Ymax])
            plt.grid()
            plt.legend()
            plt.show()
            
            print(NOM + '  :   la séparation est-elle correct ?\n');
            try:
                sep_OK = int(input('Rentrez 0 si non, 1 si oui (default = 1)\n') or 1)
            except ValueError:
                print('Vous n\'avez pas rentré un nombre, la séparation est considéré comme mauvaise')
                sep_OK=False;
    
            if sep_OK : break
        
            Xdebcoup = float(input('Nombre d\'onde à partir duquel on suprime (default = '+
                                          str(Xdebcoup) + ')\n') or Xdebcoup);
        
            Xfincoup = float(input('Nombre d\'onde jusqu\'auquel on suprime (default = '+
                                          str(Xfincoup) + ')\n') or Xfincoup);

        if (sep_OK) : 
            Y_VISDelta, Y_IRDelta, interpo_limIR, interpo_limVIS = Remontage_IR_VIS(XR1, YR1, XR2, YR2,
                                                                        mode, X,Y, NOM)
            if mode == 'PERKIN_micro':
                Ycorr = reconstruction_saut(XR1, Y_IRDelta, XR2, Y_VISDelta, interpo_limIR,
                                     interpo_limVIS, Xinterpo, X, Y, NOM);
    
            elif mode == 'PERKIN_lampe' :
                Xtemp=np.concatenate([XR1[interpo_limIR], XR2[interpo_limVIS]]) # On rassemble les deux jeux de donner pour interpoler
                Ytemp=np.concatenate([Y_IRDelta[interpo_limIR], Y_VISDelta[interpo_limVIS]]) 

                fit_interpo = np.polyfit(Xtemp, Ytemp, 3) # fit pour interpoler les data
                fonctionYinterpo = np.poly1d(fit_interpo)
                Yinterpo = fonctionYinterpo(Xinterpo);
                
                print(Xinterpo)
                print('Y interpo =')
                print(Yinterpo)
                Ycorr=np.concatenate([Y_IRDelta, Yinterpo, Y_VISDelta])
                
                #pour reconstruction tordu
                # Ycorr = reconstruction_saut(XR1, Y_IRDelta, XR2, Y_VISDelta, interpo_limIR,
                                     # interpo_limVIS, Xinterpo, X, Y, NOM, mode='Perkin_lampe'); 
                
            elif mode == 'PERKIN_filtre_319' :
                 Xtemp=np.concatenate([XR1[interpo_limIR], XR2[interpo_limVIS]]) # On rassemble les deux jeux de donner pour interpoler
                 Ytemp=np.concatenate([Y_IRDelta[interpo_limIR], Y_VISDelta[interpo_limVIS]]) 

                 # fit_interpo = np.polyfit(Xtemp, Ytemp, 3) # fit pour interpoler les data
                 # fonctionYinterpo = np.poly1d(fit_interpo)
                 # Yinterpo = fonctionYinterpo(Xinterpo);
                 
                 # print(Xinterpo)
                 # print('Y interpo =')
                 # print(Yinterpo)
                 # Ycorr=np.concatenate([Y_IRDelta, Yinterpo, Y_VISDelta])
                 
                 # pour reconstruction tordu
                 Ycorr = reconstruction_saut(XR1, Y_IRDelta, XR2, Y_VISDelta, interpo_limIR,
                                       interpo_limVIS, Xinterpo, X, Y, NOM, mode='PERKIN_filtre_319'); 
                
            else:
                Ycorr = np.concatenate([Y_IRDelta, YR2])
                
        
        plt.figure(figsize=FIGsize, dpi=DPI)

        plt.plot(X, Y, label= NOM + 'Original')
        #plt.plot(Xtemp, Ytemp, label= NOM + 'temporaire')        
        plt.plot(X, Ycorr, label = NOM +'reconstruite')
        plt.xlim([Xmin, Xmax])
        plt.ylim([Ymin, Ymax])
        plt.grid()
        plt.legend()
        show_interactif_loop()
        if(int(input('La correction est-elle correcte ?\nO pour non et recommencer depuis le debut, 1 pour oui\n') or 1)): break
        

    return (Ycorr)

def Recollage_portable(X_IR, Y_IRDelta, X_VIS, Y_VIS, interpo_limIR, interpo_limVIS, NOM):
    '''
    Cette fonction premet de coller deux jeux de données sans trou.
    Utilisé pour le spectro portable
    A complété

    Parameters
    ----------
    X_IR : TYPE
        DESCRIPTION.
    Y_IRDelta : TYPE
        DESCRIPTION.
    X_VIS : TYPE
        DESCRIPTION.
    Y_VIS : TYPE
        DESCRIPTION.
    interpo_limIR : TYPE
        DESCRIPTION.
    interpo_limVIS : TYPE
        DESCRIPTION.
    NOM : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    Recol_OK=False;
    
    FIGsize=(10,6)
    DPI=120


    #plt.ion() #Nécéssaire pour afficher les figures en %matplolib  
    while not Recol_OK:

        
        Xtemp=np.concatenate([X_IR[interpo_limIR], X_VIS[interpo_limVIS]]) # Fusion des deux jeux de donnée
        Ytemp=np.concatenate([Y_IRDelta[interpo_limIR], Y_VIS[interpo_limVIS]])
        
        INDEX=np.argsort(Xtemp) # Trie pour avoir les nombres d'onde en ordre croissant.
        
        Xtemp=Xtemp[INDEX] # application du trie
        Ytemp=Ytemp[INDEX]
        
        INDEX_IR = X_IR < np.min(Xtemp)
        INDEX_VIS = X_VIS > np.max(Xtemp)
        
        # FiltrageIR=int(input('Voulez-vous appliquer un filtre gaussien sur la partie IR ? ') or 0)
        
        # if FiltrageIR:
        #     XIRtemp = np.concatenate([X_IR[INDEX_IR], Xtemp])
        #     YIRtemp = np.concatenate([Y_IRDelta[INDEX_IR], Ytemp])
            
        #     YIRtempcorr=filtrage_gauss(XIRtemp, YIRtemp, NOM, 'auto')
            
        #     plt.plot(XIRtemp, YIRtemp, label='avant = filtrage')
        #     plt.plot(XIRtemp, YIRtempcorr, label='aprés filtrage')
            
        #     Xcorr=np.concatenate([XIRtemp, X_VIS[INDEX_VIS]])
        #     Ycorr=np.concatenate([YIRtempcorr ,Y_VIS[INDEX_VIS]])
        
        
        # else :    
        #     Xcorr=np.concatenate([X_IR[INDEX_IR], Xtemp, X_VIS[INDEX_VIS]])
        #     Ycorr=np.concatenate([Y_IRDelta[INDEX_IR], Ytemp ,Y_VIS[INDEX_VIS]])
        
        
        Xcorr=np.concatenate([X_IR[INDEX_IR], Xtemp, X_VIS[INDEX_VIS]])
        Ycorr=np.concatenate([Y_IRDelta[INDEX_IR], Ytemp ,Y_VIS[INDEX_VIS]])
    
        INDEX=np.argsort(Xcorr) # Pour avoir un jeu de donnée continue
        Xcorr=Xcorr[INDEX]
        Ycorr=Ycorr[INDEX]
        
        # plt.figure(figsize=FIGsize, dpi=DPI)
        # plt.plot(Xcorr,Ycorr, label='Spectre joint')
        # plt.grid()
        # plt.legend()
        # show_interactif_loop()
            
        # Recol_OK = int(input('La jonction est-elle correcte ? \n') or 0)
        Recol_OK=True;
    return(Xcorr, Ycorr)

def Traitement_spectro_portable(CHEMIN_IR, CHEMIN_VIS, NOM='pouet',
                                Addition_Tr='0', VISsurNIR=False, mode='IR'):
    '''
    Cette fonction déclanche le traitement des spectres issu du spectro portable

    Parameters
    ----------
    CHEMIN_IR : TYPE
        DESCRIPTION.
    CHEMIN_VIS : TYPE
        DESCRIPTION.
    NOM : TYPE, optional
        DESCRIPTION. The default is 'pouet'.
    Addition_Tr : TYPE, optional
        DESCRIPTION. The default is '0'.
    mode : STRING, optional
        Choix de travailler en UV ou en IR, defaut : IR

    Returns
    -------
    None.

    '''
    #print(Addition_Tr)
    #plt.ion() #Nécéssaire pour afficher les figures en %matplolib  

    if mode == 'UV' : COUPURE_UV=280;
    else            : COUPURE_UV=300;
     
    FIGsize=(10,6)
    DPI=120

    SHIFT_NIRQUEST=0 # décallage en nm entre le NIRQUEST et le PERKIN
    # utiliser si spectro décalibré
    
    Trait_OK = False;
    
    
    Xnm_IR, Ytr_IR=Readspectre(CHEMIN_IR) # Lecture fichier IR
    
    # Xnm_IR=Xnm_IR - SHIFT_NIRQUEST # On corrige le NIRQUEST d'un mauvaise calibration.
    
    # if SHIFT_NIRQUEST:
    #     print('ATTENTION PARTIE NIR DECALLER DE ' + str(SHIFT_NIRQUEST)+ ' NM')
    
    X_IR = (1/(Xnm_IR*1E-7))[::-1]; #
    Y_IR = (- np.log10(Ytr_IR+Addition_Tr))[::-1];# Si data en %T
    X_IR, Y_IR = removeInfNan(X_IR,Y_IR)
    
    Xnm, YTr=Readspectre(CHEMIN_VIS) # Lecture fichier VIS

    INDEX_UV = Xnm > COUPURE_UV;
    
    Xnm=Xnm[INDEX_UV]
    YTr=YTr[INDEX_UV]
    
    X_VIS = (1/(Xnm*1E-7)); # [::-1] = inversion des valeurs dans le tableau
    Y_VIS = (-np.log10(YTr+Addition_Tr));# Si data en %T
    X_VIS, Y_VIS= removeInfNan(X_VIS,Y_VIS)
    
    
    while not Trait_OK:
        print(NOM)
        
        if mode == 'UV':
            Y_VISDelta, Y_IRDelta, interpo_limIR, interpo_limVIS = Remontage_IR_VIS(X_IR, Y_IR, X_VIS, Y_VIS,
                                                                        'PortableUV', X=0,Y=0,NOM=NOM,
                                                                        VISsurNIR=VISsurNIR)
        
        else :
            Y_VISDelta, Y_IRDelta, interpo_limIR, interpo_limVIS = Remontage_IR_VIS(X_IR, Y_IR, X_VIS, Y_VIS,
                                                                    'PortableIR', X=0,Y=0,NOM=NOM,
                                                                    VISsurNIR=VISsurNIR)

        Xcorr, Ycorr = Recollage_portable(X_IR, Y_IRDelta, X_VIS, Y_VISDelta,
                                          interpo_limIR, interpo_limVIS, NOM)
    
        plt.figure(figsize=FIGsize, dpi=DPI)
        
        plt.plot(X_IR, Y_IR, label='IR_ori');
        plt.plot(X_VIS, Y_VIS, label='VIS_ori');
        plt.plot(Xcorr,Ycorr, label='Spectre traité')
        plt.grid()
        plt.legend()
        show_interactif_loop()
            
        Trait_OK = int(input('Le traitement final est-il correcte ? \n') or 1)
        
        Xcorr = 1/(Xcorr*1E-7); # On repasse en nm
    return(Xcorr, Ycorr)

def Soustraction_reflex(Liste, I100):
    '''
    Le but de cette fonction est d'extraire la partie de l'absorbance qui est liée à la refléxion dans les verres.
    Pour utiliser cette fonction mettre uniquement deux verre de même épaisseur, l'un corrigé de la reflectivités, l'autre non
    Le non corrigé doit être placé en second 
    
    '''

    if (np.size(Liste) != 2):
        raise ValueError('Il faut juste 2 spectre')
    
    else:
        Xnm, Ytr= Readspectre(Liste[0])
        Xcorr = 1/(Xnm*1E-7); # 
        Ycorr  = - np.log10(Ytr)-I100; #
        
        
        Xnm, Ytr= Readspectre(Liste[0])
        
        #Xpascorr = 1/(Data[:, 0]*1E-7); # 
        Ypascorr  = - np.log10(Ytr)-I100; #
        
        Ydiff=Ypascorr-Ycorr;

        plt.figure(figsize=(5,3), dpi=120)
        plt.plot(Xcorr, Ydiff)
        plt.grid(True);
        plt.xlabel("Nombre d'onde ($cm^{-1}$)");
        plt.ylabel('Absorbance')  
        plt.xlim(5000, 25000);
        plt.ylim(-0.1, 0.2);
        plt.show
    return(Ydiff)

def Soustraction_interpolation(X1, Y1, X2, Y2, signe='-'):
    '''
    Cette fonction retourne Y2-Y1 sur les valeurs commune de X1 et X2 en interpollant linéairement les points
    de Y1 ou Y2.

    Parameters
    ----------
    X1 : TYPE
        DESCRIPTION.
    Y1 : TYPE
        DESCRIPTION.
    X2 : TYPE
        DESCRIPTION.
    Y2 : TYPE
        DESCRIPTION.

    Returns
    -------
    Ydiff :
        Soustraction ou addition entre Y2 et Y1.
    Xref:
        Gamme des X communes au deux courbe

    '''
    # plt.figure(figsize=([6,3]), dpi=120)
    # plt.plot(X1, Y1, '.', label='DATA 1', markersize=1)
    # plt.plot(X2, Y2, '.', label='DATA 2', markersize=3)
    # plt.legend()
    
    if (np.nanmin(X1)<np.nanmin(X2)):
        if (np.nanmax(X1)>np.nanmax(X2)): # Les valeur de X2 sont strictement comprise dans X1
            f = interpolate.interp1d(X1, Y1)
            Xref=X2;
            Y1=f(Xref);
            
        else: #X1 va plus bas en valeur de que X2 mais X2 va plus haut, on coupe donc X2
            coupure_sup = (X2 <= np.nanmax(X1))    
            X2=X2[coupure_sup]
            Xref=X2;
            
            Y2=Y2[coupure_sup]
            
            f = interpolate.interp1d(X1, Y1)
            Y1=f(Xref);

    else: # min X2 > min X1
        if (np.nanmax(X1)>np.nanmax(X2)): # 
            coupure_sup = (X1 <= np.nanmax(X2))    
            X1=X1[coupure_sup]
            Xref=X1;
            
            Y1=Y1[coupure_sup]
            f = interpolate.interp1d(X2, Y2)
            Y2=f(Xref)
            #print('pouet1')
            
        else: # Les valeurs de X1 sont comrpise dans X2.
            f = interpolate.interp1d(X2, Y2)
            Xref=X1;
            Y2=f(Xref);
            #print('pouet2')
    # print(Xref)
    # print('\n\n')
    # Xref=np.flipud(Xref)    
    # print(Xref)
    # print('\n\n')
    
    if signe == '-' :
        Y_diff=Y2-Y1;
        
    elif signe == '+' :
        print('ADDITION\n\n\n')
        Y_diff=Y2+Y1;   
        
    # #Y_diff=-Y_diff
    # plt.plot(Xref, Y_diff, '^', label="différence", markersize=1)
    # plt.legend()
    # plt.show
    
    return(Xref, Y_diff)

def Soustraction_norm(X1, Y1, X2, Y2, COUPUREminmax=[200, 3300]):
    '''
    Cette fonction soustrait en interpollant deux jeux de données en
    les normalisant entre 0 et 1 au préalable.
    Y2-Y1

    Parameters
    ----------
    X1 : TYPE
        DESCRIPTION.
    Y1 : TYPE
        DESCRIPTION.
    X2 : TYPE
        DESCRIPTION.
    Y2 : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    Y1 = normYminmax(X1, Y1, COUPUREminmax)
    
    Y2 = normYminmax(X2, Y2, COUPUREminmax)

    return(Soustraction_interpolation(X1, Y1, X2, Y2))


def SavCSV(X, Y, NOM, Legende='', Dossier='./Data_corr', ENTETE='\n X en nm; Y en %T',
           delimiter=';'):
    try:
        os.mkdir(Dossier)
    except OSError:
        pass
    
    HEADER =NOM+ ';' +Legende + ENTETE
    
    Save=np.array([X, Y]).transpose() # On met sous la forme XY en colonne. DEPRECIER ?

    #Save=np.concatenate((X[np.newaxis,:], Y[np.newaxis,:])).transpose() à vérifier la fiabilité
    np.savetxt(Dossier+os.sep+NOM, Save, delimiter=delimiter, header=HEADER, comments='');
    
def Reflexion(X,Y,Legend, SHOW=True)    :
    '''
    Cette fonction permet de supprimer les pertes de transmittance dues à la reflexion sur les 2 faces du verre 
    X tableau de float en longueur d'onde 
    

    Returns
    -------
    None.

    '''
    Xnm=nm2cm1(X)
    n = 1.517+6200E-14/(Xnm**2)+1000000E-28/(Xnm**4) #indice de réfraction
    R = np.power(((1-n)/(1+n)),2)#coeff de reflexion
    A_ref =-2*np.log10(1-R) #absorbance due à la reflexion et à soustraire à l'absorbance totale
    print('A_ref : ')
    print(A_ref)
    Ycorr  = Y-A_ref; 
    
    if SHOW == True :
        plt.plot(X, Ycorr, label = Legend + ' corrigé reflexion')
        plt.plot(X, Y, label= Legend + ' non corrigé')
        plt.grid(True);
        plt.xlabel("Nombre d'onde ($cm^{-1}$)");
        plt.ylabel('Absorbance')  
        plt.xlim(3000, 33000)
        plt.legend()
        plt.show()
        
    return(Ycorr)
        
def Sub_baseline (X,Y,Legend, lim=np.array([nm2cm1(2500), nm2cm1(360)]), retBaseline=False) :
    '''
    Cette fonction permet de calculer le spectre sans la ligne de base
    
    Returns
    -------
    None.
    '''
    INDEX = np.logical_and(X < lim.max(), X > lim.min())
    X=X[INDEX]
    Y=Y[INDEX]
    
    Xbaseline, Ybaseline = baseline_Rubberband(X, Y)
    baseline=np.interp(X, Xbaseline, Ybaseline)
    Ycorr=Y-baseline
    
    if retBaseline :
        Ycorr = baseline
        
    return(X, Ycorr);



def Nettoyage_spectre_OAS(Liste, Legende, Liste_ref, correction, Liste_refN='',
                      Dossier='',  Addition_Tr=0,
                      Borne=[nm2cm1(2500), nm2cm1(330)], valeurnorm=0, DataExcel = False):
    '''
    Cette fonction pilote les différents corrections et modifcations de spectre.

    Parameters
    ----------
    Liste : TYPE
        DESCRIPTION.
    Legende : TYPE
        DESCRIPTION.
    Liste_ref : TYPE
        DESCRIPTION.
    correction : INT
            Entier qui permet de choisir la correction,
            1 pour soustraire la référence,
            10 pour corriger avec un noir et un blanc,
            11 correction saut perkin,
            12 correction saut perkin microscope,
            13 correction saut de lampe PERKIN
            14 correction filtre PERKIN 319nm
            113 correction I0 + saut detecteurs std + saut lampe 
            123 correction I0 + saut detecteurs microscope + saut lampe 
            2 soustration,
            21 soustraction normalisé,
            22 soustraction de la bande UV,
            23 soustraction ligne de base,
            24 soustraction reflexion par le calcul,
25
26
27
            3 lissage (au choix gaussien au Savgol),
            4 enlevève inf et nan,
            41 sauvegarde normaliser par l'épaisseur
            411 sauvergarde normaliser par l'épaisseur et minimum soustrait
            4123 norm ep et baseline rubberband soustrait
            5 pour passer une fichier d'ABS en Tr,
            6 pour joindre les spectres issu du spectro portable,
            -6 pour VIS sur NIR,
            66 pour rajouter UV APRES correction 6 ou -6,
            7 addition d'un spectre en interpolant (pour DO spectro portable),
            9 extraction ligne de base.
    Liste_refN : TYPE, optional
        DESCRIPTION. The default is ''.
    Dossier : STR, optional
        Repertoir de sauvegarde des data_traité. The default is './Data_corr'.
    Addition_Tr : TYPE, optional
        DESCRIPTION. The default is 0.
    Borne : TYPE, optional
        Borne en cm-1. The default is [nm2cm1(2500), nm2cm1(330)].

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    La liste des fichier corriger.

    '''
    Liste_corr=[]
    ENTETE = '\n X en nm; Y en %T'
    RAJOUT = ''

    if os.path.exists("Position_min.csv"):
        os.remove("Position_min.csv")
    Borne = np.array(Borne)
    
    if np.size(Addition_Tr) == 1: Addition_Tr=(np.zeros(np.size(Liste)) + Addition_Tr)
    
    for i in np.arange(0, np.size(Liste), 1):
        #print(correction)
        
        if(correction[i]==0): # Si pas de correction on ne fait rien
            Liste_corr.append(Liste[i])
            
        else:
            Fichier=Liste[i];
            (cheminfichier, nomfichier) = os.path.split(Fichier)
           
            if DataExcel :  #Lecture de données de fichier excel
                pouet = pd.read_excel(Fichier)
                Xnm = pouet['nm'].to_numpy()
                Ytr = pouet['T'].to_numpy()
            else :
                try : Xnm, Ytr = Readspectre(Fichier)
                except IndexError : Xnm, Ytr = Readspectre(Fichier, delimiter=',')
                
            X = nm2cm1(Xnm); #
            
            Y = - np.log10(Ytr+Addition_Tr[i])# Si data en %T
            
            Xsave=Xnm;
            if Dossier == '' : Dossier = Corr2folder_OAS(correction[i])
            
            if(correction[i] == 1): # uniquement soustraction de la référence
        
               
                try :  Xnm_ref, Ytr_ref=Readspectre(Liste_ref[i])
                except IndexError :  Xnm_ref, Ytr_ref=Readspectre(Liste_ref[i], delimiter=',')
                
                I100 = -np.log10(Ytr_ref)

                Y_corr=Y-I100
                Fichier_corr=nomfichier[0:-4] + Corr2Str_OAS(correction[i])
                #Fichier_corr= Legende[i]+'_diff.csv'
              
            elif(correction[i] == 10): # Division blanc (tr) et soustraction noir (tr)
                Xnm_refN, Ytr_noir=Readspectre(Liste_refN[i])
                Xnm_ref, Ytr_blanc=Readspectre(Liste_ref[i])
                
                Ytr_blanc   = Ytr_blanc - Ytr_noir
                Ytr         = Ytr       - Ytr_noir
                
                Ytr_corr    = Ytr / Ytr_blanc
                
                Y_corr = - np.log10(Ytr_corr+Addition_Tr[i])
                
                Fichier_corr=nomfichier[0:-4] + Corr2Str_OAS(correction[i])
            
            elif(correction[i] == 11) : # correction saut std PERKIN
                Xnm_ref, Ytr_ref=Readspectre(Liste_ref[i])
                I100 = -np.log10(Ytr_ref)
                Y=Y-I100
            
                Y_corr= correction_saut_detect(X, Y, Legende[i], mode='PERKIN_std')
                Fichier_corr=nomfichier[0:-4] + Corr2Str_OAS(correction[i])
            
            elif(correction[i] == 12): # Saut de détecteur Perkin microscope
                Y_corr= correction_saut_detect(X, Y, Legende[i], mode='PERKIN_micro')
                Fichier_corr=nomfichier[0:-4] + Corr2Str_OAS(correction[i])
            
            elif(correction[i] == 13) : # correction saut std PERKIN
                 Xnm_ref, Ytr_ref=Readspectre(Liste_ref[i])
                 I100 = -np.log10(Ytr_ref)
                 Y=Y-I100
                
                 Y_corr= correction_saut_detect(X, Y, Legende[i], mode='PERKIN_lampe')
                 Fichier_corr=nomfichier[0:-4] + Corr2Str_OAS(correction[i])
                 
            elif(correction[i] == 14) : # correction saut filtre 319 PERKIN # Enfaite LAMPES a corriger
                 Y_corr= correction_saut_detect(X, Y, Legende[i], mode='PERKIN_filtre_319')
                 Fichier_corr=nomfichier[0:-4] + Corr2Str_OAS(correction[i])
                 
            elif(correction[i] == 113) : 
                Xnm_ref, Ytr_ref=Readspectre(Liste_ref[i])
                I100 = -np.log10(Ytr_ref)
                Y_corr=Y-I100
                
                Y_corr= correction_saut_detect(X, Y_corr, Legende[i], mode='PERKIN_std')
                Y_corr= correction_saut_detect(X, Y_corr, Legende[i], mode='PERKIN_lampe')
                Fichier_corr=nomfichier[0:-4] + Corr2Str_OAS(correction[i])
                
            elif(correction[i] == 123) : 
                Xnm_ref, Ytr_ref=Readspectre(Liste_ref[i])
                I100 = -np.log10(Ytr_ref)
                Y_corr=Y-I100
                
                Y_corr= correction_saut_detect(X, Y_corr, Legende[i], mode='PERKIN_micro')
                Y_corr= correction_saut_detect(X, Y_corr, Legende[i], mode='PERKIN_lampe')
                Fichier_corr=nomfichier[0:-4] + Corr2Str_OAS(correction[i])
            
            elif(correction[i] == 1234) :
                Xnm_ref, Ytr_ref=Readspectre(Liste_ref[i])
                I100 = -np.log10(Ytr_ref)
                Y=Y-I100
                
                Y_corr= correction_saut_detect(X, Y, Legende[i], mode='PERKIN_micro')
                Y_corr= correction_saut_detect(X, Y_corr, Legende[i], mode='PERKIN_lampe')
                Y_corr= correction_saut_detect(X, Y_corr, Legende[i], mode='PERKIN_filtre_319')
                Fichier_corr=nomfichier[0:-4] + Corr2Str_OAS(correction[i])
        
        
        
            elif(correction[i] == 2): # Soustration
                Xnm_ref, Ytr_ref=Readspectre(Liste_ref[i])
                YREF = -np.log10(Ytr_ref)
                XREF = 1/(Xnm_ref*1E-7);

                X, Y_corr = Soustraction_interpolation(XREF, YREF, X, Y)
                X, Y_corr = removeInfNan(X, Y_corr)
                
                Xsave=1/(X*1E-7);                
                Fichier_corr=nomfichier[0:-4] + Corr2Str_OAS(correction[i])
                Fichier_corr= Legende[i]+'_diff.csv'
            
            elif(correction[i] == 21): # Soustration normalisé
                Xnm_ref, Ytr_ref=Readspectre(Liste_ref[i])
                YREF = -np.log10(Ytr_ref)
                XREF = 1/(Xnm_ref*1E-7);

                Xsave, Y_corr = Soustraction_norm(XREF, YREF, X, Y,
                                                  COUPUREminmax=Borne)
                        
                Fichier_corr=nomfichier[0:-4] + Corr2Str_OAS(correction[i])
            
            elif(correction[i] == 22): # soustraction bande UV                
                Xsave, Y_corr=fit_OMCT(X,Y,Legende[i])
                
                Fichier_corr=nomfichier[0:-4]+ Corr2Str_OAS(correction[i])
           
            elif(correction[i] == 23) : # soustraction ligne de base
                X, Y_corr = Sub_baseline(X, Y, Legende[i], Borne)
                Xsave = nm2cm1(X)
                
                Fichier_corr=nomfichier[0:-4] + Corr2Str_OAS(correction[i])    

            elif(correction[i] == 231) : # extraction ligne de base
                X, Y_corr = Sub_baseline(X, Y, Legende[i], Borne, retBaseline=True)
                Xsave = nm2cm1(X)
                
                Fichier_corr = nomfichier[0:-4] + Corr2Str_OAS(correction[i])
                
            elif(correction[i]==2341) : # sub baseline puis normalise par l'épaisseur et

                X, Y_corr = Sub_baseline(X, Y, Legende[i], Borne)
                Y_corr = Y_corr/valeurnorm[i]
                
                Xsave = nm2cm1(X)
                
                Fichier_corr=nomfichier[0:-4] + Corr2Str_OAS(correction[i])    

            elif(correction[i] == 24) : # Soustration reflexion calcul
                 Y_corr = Reflexion(X, Y, Legende[i])
                 
                 Fichier_corr=nomfichier[0:-4] + Corr2Str_OAS(correction[i])
            
            elif(correction[i]==25) : # soustrait le minimum et normalise par l'épaisseur
                 Y_corr = Y-Y.min()    
                 Y_corr = Y_corr/valeurnorm[i]
 
                 Fichier_corr=nomfichier[0:-4] + Corr2Str_OAS(correction[i]) 
            
            elif(correction[i] == 26) : #soustrait les gaussiennes du fichier parametre_fit_gauss
                Data=np.genfromtxt('./parametre_fit_gauss.csv', comments='#', skip_header=1, delimiter=';', dtype='str' )
                if Data[i, 0] in Legende :
                    Y1 = Gauss(X, float(Data[i,1]), float(Data[i,2]), float(Data[i,3]), float(Data[i,4]))
                    Y2 = Gauss(X, float(Data[i,5]), float(Data[i,6]), float(Data[i,7]), float(Data[i,8]))
                Y_corr = Y-Y1-Y2
                Fichier_corr=nomfichier[0:-4] + Corr2Str_OAS(correction[i])      
  
            elif(correction[i] == 27) : # soustraction ligne de base lineaire
                X, Y_corr = baseline_linear(X, Y, Legende[i])
                Xsave = nm2cm1(X)
                
                Fichier_corr=nomfichier[0:-4] + Corr2Str_OAS(correction[i])  
                       
                
                
            elif(correction[i] == 3): #Lissage
  
                 Xsave, Y_corr= lissage (X, Y, Legende[i])
                 Xsave = 1/(Xsave*1E-7);
                        
                 Fichier_corr=nomfichier[0:-4]+ Corr2Str_OAS(correction[i])
            
            elif(correction[i]==4) : #enlève inf et nan
                Xsave, Y_corr = removeInfNan(Xnm,Y)
                      
                Fichier_corr=nomfichier[0:-4] + Corr2Str_OAS(correction[i])    
            
            elif(correction[i]==41) : #normalise par l'épaisseur
                Y_corr = Y/valeurnorm[i]
                Xsave, Y_corr = removeInfNan(Xnm,Y_corr)
                
                Fichier_corr=nomfichier[0:-4] + Corr2Str_OAS(correction[i])    
            


            elif(correction[i]==5): # passe en Tr, si Ttr en abosrbance il sera remis en %T lors de la sauv
                Y_corr  = Ytr
                Fichier_corr=nomfichier[0:-4] + Corr2Str_OAS(correction[i])
                Dossier = Corr2folder_OAS(correction[i])  
                
            elif(correction[i] == 6): # Spectro portable, pour le raccord bien penser à mettre le spectre IR dans la Liste et le VIS dans la ref   
                Dossier = Corr2folder_OAS(correction[i])      
                Xnm, Y_corr=Traitement_spectro_portable(Liste_ref[i], Liste[i],
                                                        Legende[i],Addition_Tr[i],
                                                        VISsurNIR=False, mode='IR')
                Xsave = Xnm;
                
                Fichier_corr=nomfichier[0:-4]+ Corr2Str_OAS(correction[i])
                
            elif(correction[i] == -6): # Idem 6 mais on remonte VIS sur NIR   
                Dossier = Corr2folder_OAS(correction[i])    
                Xnm, Y_corr=Traitement_spectro_portable(Liste_ref[i], Liste[i],
                                                        Legende[i], Addition_Tr[i],
                                                        VISsurNIR=True, mode='IR')
                Xsave = Xnm;
                
                Fichier_corr=nomfichier[0:-4]+ Corr2Str_OAS(correction[i])
                RAJOUT = 'VIS remonté sur NIR'
            
            elif(correction[i] == 66): # On colle UV sur VIS-NIR
                (cheminfichier, nomfichier) = os.path.split(Liste[i])
                Liste_temp = Corr2folder_OAS(correction[i]) + os.sep + nomfichier[:-4]+Corr2Str_OAS(6)
                # Recherche fichier joint NIR (correction 6 ou -6)
                try:
                    Xnm, Y_corr=Traitement_spectro_portable(Liste_temp, Liste_refN[i],
                                                        Legende[i], Addition_Tr[i],
                                                        VISsurNIR=True, mode='UV')
                except OSError :
                    raise ValueError('\n\n fais la correction 6 (jonction NIR) avant, bisous tocard \n\n')
                
                Xsave = Xnm;
                
                Fichier_corr=nomfichier[:-4]+ Corr2Str_OAS(correction[i])
                RAJOUT = 'Rajout ' + Liste_refN[i]
            
            elif(correction[i] == 7): # Ajout d'un spectre, utilise si DO mesure spectro portable
                Xnm_ref, Ytr_ref=Readspectre(Liste_ref[i])
                REF = -np.log10(Ytr_ref)

                #Y_corr=Y+REF
                Xsave, Y_corr = Soustraction_interpolation(Xnm_ref, REF, Xnm, Y,'+')
                Fichier_corr=nomfichier[0:-4] + Corr2Str_OAS(correction[i])
            
            
            elif(correction[i]==9): #extraction ligne de base
                INDEXUV=Xnm>Borne[1]
                
                Xsave, Y_corr = baseline_Rubberband(X[INDEXUV], Y[INDEXUV])
                Xsave=nm2cm1(Xsave)
                Fichier_corr=nomfichier[0:-4] + Corr2Str_OAS(correction[i])
                Dossier = Corr2folder_OAS(correction[i])  
                ENTETE = '\n X en nm; Y en %T'
           
            
            else:
                raise ValueError('La correction demandée n\'est pas implémentée, vérifier le fichier d\'input')
            
            Liste_corr.append(Dossier+os.sep+Fichier_corr);
            
            if not (correction[i] == 100000):
                Y_save=np.power(10, -Y_corr); # On repasse en %T pour assurer la comptabilité avec le reste des scripts

            # if DataExcel:
            #     save = pd.DataFrame(np.array([X, Y, -np.log(Y)]), columns = ['nm', 'T', 'A'])
            #     save.to_excel(Fichier_corr+'.xlsx')
                
            # else :
            SavCSV(Xsave, Y_save, Fichier_corr, Legende[i]+RAJOUT, Dossier, ENTETE);

    return(Liste_corr)

#%% RAMAN :

def Indice_polymerisation(Xsave,Ysave,Legende) :           
    x_1 = 304  #valeur à modifier en fonction de la ligne de base 
    x_2 = 735
    x_3 = 856
    x_4 = 1252
    
    X1 = np.absolute(x_1-Xsave)
    X2 = np.absolute(x_2-Xsave)
    X3 = np.absolute(x_3-Xsave)
    X4 = np.absolute(x_4-Xsave)
    
    p_1 = X1.argmin()
    p_2 = X2.argmin()
    p_3 = X3.argmin()
    p_4 = X4.argmin()
       
    
    Y_bend = Ysave[p_1:p_2]
    Y_stretch = Ysave[p_3:p_4]
    X_bend = Xsave[p_1:p_2]
    X_stretch = Xsave[p_3:p_4]
    
    i_max = Y_stretch.argmax()
    print(i_max)
    max_stretch = X_stretch[i_max]
    
    area_bend = np.trapz(Y_bend, X_bend)
    area_stretch = np.trapz(Y_stretch, X_stretch)
    Ip = area_bend/area_stretch
    
    try :
        text_file=open('Indice_polymerisation.csv', "x")
        ENTETE='Nom fichier;Aire bande bending; Aire bande stretching;Ip;Max stretching\n'
        text_file.write(ENTETE)
    except FileExistsError :
        text_file=open('Indice_polymerisation.csv', "a")
    
    DATA= Legende + ';' + str(area_bend) + ';' + str(area_stretch) + ';' + str(Ip) + ';'+ str(max_stretch) + '\n'
    
    text_file.write(DATA)
    
    text_file.close()
    
    return(Ip)      
            

def traitement_Raman_portable_spline(X, Y, Legende='pouet', borne=[80,1700],
                              roi = np.array([[80,500],[720,740],
                                              [810, 830],[1300,1700]])):
    #plt.style.use({'lines.linewidth': 1})  
    plt.figure(dpi=200)
    ax = plt.gca()
    ax.grid()
    
    Xcorr = np.arange(np.min(borne), np.max(borne), 0.5) # we generate the new X values with numpy.arange()
    Y_temp = rampy.resample(X, Y, Xcorr)
    
    Ycorr, base_gcvspl = rampy.baseline(Xcorr,Y_temp,roi,'gcvspline',s=0.1 ) # activate if you have installed gcvspline
    
    Ycorr=Ycorr[:, 0]
    # doing the figure

    ax.plot(X, Y, "-", label=Legende)
    ax.plot(Xcorr, Ycorr, "g-", label="extract_gcvspl")
    ax.plot(Xcorr,base_gcvspl,"-",color="red",label="gcvspline baseline") # activate if you have installed gcvspline
    
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    
    ax.legend()
    Sav_fig(Legende, './Graph/Ligne_de_base_spline')

    return(Xcorr, Ycorr)    


def traitement_Raman_portable_polyline(X, Y, Legende='pouet', borne=[80,1700],
                                       pt = np.array([80, 100, 200, 300, 736, 843, 1302, 1387, 1616, 1701])):
    #plt.style.use({'lines.linewidth': 1})  
    plt.figure(dpi=200)       
    ax = plt.gca()
    ax.grid()
    
    Xcorr = np.arange(np.min(borne), np.max(borne), 0.5) # we generate the new X values with numpy.arange()
    #Y_temp = rampy.resample(X, Y, Xcorr)
    Y_temp = np.interp(Xcorr, X, Y)
    
    background = interp1d(X, Y, kind='linear')
    background = background(pt)
    
    #baseline=np.interp(Xcorr, pt, background)
    baseline = interp1d(pt, background, kind='linear')
    baseline = baseline(Xcorr)
    
    Ycorr=Y_temp-baseline
    
    #Ycorr=Ycorr[:, 0]
    # doing the figure

    ax.plot(X, Y, "-", label=Legende)
    ax.plot(Xcorr, Ycorr, "g-", label="substract_polyline")
    ax.plot(Xcorr, baseline,"-",color="red",label="baseline") # activate if you have installed gcvspline
    
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    
    ax.legend()
    Sav_fig(Legende, './Graph/Ligne_de_base_polyline')

    return(Xcorr, Ycorr)

def traitement_Raman_portable_unispline(X, Y, Legende='pouet', borne=[80,1700],
                              roi = np.array([[80,500],[720,740],
                                              [810, 830],[1300,1700]])) : 

    plt.figure(dpi=200)
    ax = plt.gca()
    ax.grid()
    
    Xcorr = np.arange(np.min(borne), np.max(borne), 0.5) # we generate the new X values with numpy.arange()
    Y_temp = rampy.resample(X, Y, Xcorr)
    
    Ycorr, base_unispline = rampy.baseline(Xcorr, Y_temp,roi,'unispline', s=0.5    )

    Ycorr=Ycorr[:,0]

    ax.plot(X, Y, "-", label=Legende)
    ax.plot(Xcorr, Ycorr, "g-", label="extract_unispline")
    ax.plot(Xcorr,base_unispline,"-",color="red",label="unispline baseline")

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
        
    ax.legend()
    Sav_fig(Legende, './Graph/Ligne_de_base_unispline')

    return(Xcorr, Ycorr)

def baseline_rubberband_smooth_raman(X, Y, Legende, 
                                     ordre_savgol = 3, nbpt_savgol = 31,
                                     Dossier_baseline='baseline', Limite=[500, 1300]):
    '''
    Cette fonction extrait une ligne de base style rubber-band (fond convexe),
    la ligne de base est dÃ©terminÃ©e sur les donnÃ©es lissÃ©es via la mÃ©thode savgol.

    Parameters
    ----------
    X : TYPE
        DESCRIPTION.
    Y : TYPE
        DESCRIPTION.
    Legende : TYPE
        DESCRIPTION.
    ordre_savgol : TYPE, optional
        DESCRIPTION. The default is 3.
    nbpt_savgol : TYPE, optional
        DESCRIPTION. The default is 31.

    Returns
    -------
    Ysave : les donnÃ©es traitÃ©es.

    '''
    plt.figure(dpi=200)  
    ax = plt.gca()
    ax.grid()
    
    INDEX = np.logical_and(X > min(Limite), X< max(Limite))
    X=X[INDEX]
    Y=Y[INDEX]
    
    

    Y_liss = savgol_filter(Y, nbpt_savgol, ordre_savgol)

    X_baseline, Y_baseline = baseline_Rubberband(X, Y_liss)
    #Make sure the baseline goes through the last point of the interval
    X_last_pt=X[-1]
    Y_last_pt=Y[-1]
    
    # and the first point :
    X_first_pt=X[0]
    Y_first_pt=Y[0]
    
    X_baseline = np.concatenate([[X_first_pt], X_baseline, [X_last_pt]])
    Y_baseline = np.concatenate([[Y_first_pt], Y_baseline, [Y_last_pt]])
    
    baseline=np.interp(X, X_baseline, Y_baseline)
    Ysave=Y-baseline
    
    SavCSV(X_baseline, Y_baseline, 'baseline_rubberband_extraite_de_'+Legende, Legende, Dossier_baseline,
           ENTETE='', delimiter='\t')
    
    ax.plot(X, Y_liss, "-", label=Legende+'_  lissage')
    ax.plot(X, Ysave, "g-", label="extract_rubberband_smooth")
    ax.plot(X, baseline,"-",color="red",label="rubberband_baseline")
    ax.minorticks_on()
    
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
        
    ax.legend()
    Sav_fig(Legende, './Graph/Ligne_de_base_rubberband')
    
    return(X, Ysave)

def modified_z_score(Y): # Pour extraction spike
    # First we calculate âˆ‡x(i):
    dist = 0
    delta_intensity = [] 
    for i in np.arange(len(Y)-1):
        dist = Y[i+1] - Y[i]
        delta_intensity.append(dist)
       
    # Alternatively to the for loop one can use: 
    # delta_int = np.diff(intensity)
        
    delta_int = np.array(delta_intensity)
    median_int = np.median(Y)
    mad_int = np.median([np.abs(Y - median_int)])
    modified_z_scores = 0.6745*(delta_int - median_int) / mad_int
    
    return modified_z_scores

def fixer(X, Y, m=5):  # Pour extraction spike
    threshold = float(input("Choisir la valeur seuil du Z score. \n Celle-ci doit correspondre à une valeur SUPERIEURE au Z-score moyen. \n"))
    
    while not threshold == 0.0 :
        Y_threshold = threshold*np.ones(len(X))
        Y_modified_z_score = abs(np.array(modified_z_score(Y)))   
     
        plt.plot(X, Y_threshold, label = 'threshold', color = 'orange') #Superpose le threshold au Z-score pour check
        plt.plot(X[1:], Y_modified_z_score)
        plt.xlabel('Raman shift (cm$^{-1}$)')
        plt.ylabel('Modified Z-Score of âˆ‡x(i)')
        show_interactif_loop()
    
        spikes = abs(np.array(modified_z_score(Y))) > threshold
        Y_out = Y.copy()

        for i in np.arange(len(spikes)):
            if spikes[i] != 0: # If we have a spike at position i
                w = np.arange(i-m, i+1+m) # we select 2 m + 1 points around our spike
                w2 = w[spikes[w] == 0] # From such interval, we choose the ones which are not spikes
                Y_out[i] = np.mean(Y[w2]) # and we average their values
        
        return Y_out

def suppression_spikes(X, Y, Legende):
    '''Adapté de l'algorithme de Whitaker-Hayes sous R et de la version python de Coca.
    On utilise le score z modifiÃ© pour detecter les spikes.
    Le score z indique de combien d'Ã©cart-types la valeur regardÃ©e diffÃ¨re de la moyenne
    (z-score = (x-Âµ)/sigma
     modified z-score = 0.6745(x â€“ median of the dataset M) / median absolute deviation of the dataset = median(|x-M|))
    Whitaker & Hayes proposent de profiter de la forte intensité et de la faible largeur des spikes 
    et donc de remplacer x par âˆ‡x = x(i)-x(i-1) : z = 0.6745 (âˆ‡x(i)-M) / MAD. 
    Ce calcul est effectué par la fonction modified_z_score.
    Le spectre est ensuite corrigÃ© par la fonction fixer.'''
    
    Trait_OK = False #Pour permettre de recommencer le traitement si non satisfaisant
    i = 0
    Despiking = ''
    
    while not Trait_OK :
        plt.plot(X, Y) #Plot les data brutes
        plt.grid()
        plt.minorticks_on()
        plt.xlabel('Raman shift (cm$^{-1}$)')
        plt.ylabel('Intensity (a.u.)')
        show_interactif_loop()
    
        despicable = int(input('Voulez-vous procéder à la suppression des spikes ? \n (1 pour oui [default], 0 pour non) \n') or 1)
        i = i + 1
        
        if despicable == 1 :
            Despiking = '_despiked'
            Y_modified_z_score = abs(np.array(modified_z_score(Y)))
            
            plt.plot(X[1:], Y_modified_z_score) #Plot le Z-score pour permettre de trouver le threshold
            plt.grid()
            plt.minorticks_on()
            plt.xlabel('Raman shift (cm$^{-1}$)')
            plt.ylabel('Modified Z-Score of âˆ‡x(i)')
            show_interactif_loop()
            
            try : 
                Y_out = fixer(X, Y, m=5)
                
            except : 
                print('La valeur choisie pour le seuil ne convient pas car elle inclut un trop grand nombre de points. \n Try again !' )

            iteration = str(i)
            plt.figure(dpi=200)     
            plt.plot(X, Y, label = 'Original data')
            plt.plot(X, Y_out, label = 'Despiked spectrum')
            plt.xlabel('Raman shift (cm$^{-1}$)')
            plt.ylabel('Intensity (a.u.)')
            Sav_fig(Legende+'_iteration_'+iteration, './Graph/Despiked')
    
        else :
            return (Y, Despiking)
          
        Trait_OK = int(input('La suppression des spikes est-elle correcte ? \n 1 pour oui [default] \n 3 pour poursuivre le traitement \n') or 1)  #DOUTE #fin de la boucle de traitement ?
        
        if Trait_OK == 3:
            Y=Y_out
            Trait_OK = 0
        
    return (Y_out, Despiking)

def normalise_spectra(X, Y, Legende, Normalisation_Mode=None) :
    """
    Utilise les options de normalisation de rampy (modes 1, 3 & 4) et permet de normaliser
    à l'aire d'une partie du spectre (mode 2).
    """
    
    if Normalisation_Mode == None:
        Normalisation_Mode = int(input("Veuillez choisir un mode de normalisation. \n 1 pour une normalisation à l'aire totale [default] \n 2 pour une normalisation à l'aire d'une bande (à  définir) \n 3 pour une normalisation à  l'intensité \n 4 pour une normalisation minmax \n") or 1) 
        
    if Normalisation_Mode == 1 or Normalisation_Mode == 'area':
        Normalisation_type = 'area'
        try :
            Ysave = rampy.normalise(Y,X,method="area")
        except : 
            Ysave = Y/np.trapz(Y, X)
    
    if Normalisation_Mode == 2 : 
        Normalisation_type = 'zone_area'
        
        Borne_inf = float(input("Borne inférieure de la zone de normalisation : "))
        Borne_sup = float(input("Borne supérieure de la zone de normalisation : "))
        INDEX = np.logical_and(X > Borne_inf, X < Borne_sup)
        
        plt.figure(dpi=200)  
        ax = plt.gca()
        plt.style.use({'legend.loc': 'upper left'})
        ax.plot(X, Y, label = 'data_before_normalization')
        ax.fill_betweenx(plt.ylim(), Borne_inf, Borne_sup, color='green', alpha=0.2, label='zone de normalisation')
        ax.text(Borne_inf, 250, [Borne_inf, Borne_sup])
        ax.set_xlabel('Raman shift (cm$^{-1}$)')
        ax.set_ylabel('Intensity (a.u.)')
        Sav_fig(Legende+'_parameters_for_normalization_to_zone_area', './Graph/Normalized')
        
        Zone_Area = np.trapz(Y[INDEX], X[INDEX], axis=0)
        
        Normalisation_Factor = 1/Zone_Area
        Ysave = Y * Normalisation_Factor
        
        print('Aire de la zone : ')
        print(np.trapz(Ysave[INDEX], X[INDEX], axis=0))
        print('Aire totale : ')
        print(np.trapz(Ysave, X, axis = 0))
        
    if Normalisation_Mode == 3 :
        Normalisation_type = 'intensity'
        Ysave = rampy.normalise(Y,method="intensity")
        
    if Normalisation_Mode == 4 :
        Normalisation_type = 'minmax'
        Ysave = rampy.normalise(Y,method="minmax")

    return(X, Ysave, Normalisation_type)

def mean_spectra(X, Y, Legende, Liste_ref, Liste_refN):
    plt.figure(dpi=200)  
    ax = plt.gca()
    ax.grid()
    
    Xsave = X
    
    Xnm_ref, Ytr_ref=Readspectre(Liste_ref, delimiter='\t')
    
    try : 
        Xnm_refN, Ytr_refN=Readspectre(Liste_refN, delimiter = '\t')
        Ysave = (Y + Ytr_ref + Ytr_refN) / 3
        
        ax.plot(X, Y, "-", label='spot 1')
        ax.plot(Xnm_ref, Ytr_ref, "-", label='spot 2')
        ax.plot(Xnm_ref, Ytr_ref, "-", label='spot 3')
        
    except :
        Ysave = (Y + Ytr_ref) / 2
    
        ax.plot(X, Y, label='spot 1')
        ax.plot(Xnm_ref, Ytr_ref, label='spot 2')
    
    ax.plot(Xsave, Ysave,linestyle = "--",color="red",label="moyenne")
    
    return(Xsave, Ysave)

def traitement_Raman_IMPMC (X, Y, Legende, Despik = True):
    '''Cette fonction est destinée à réaliser le traitement de spectres Raman 
    acquis sur le dispositif situé en 23-13 3ème à l'IMPMC.
    1) Elle permet de supprimer les spikes en utilisant la fonction suppression_spikes
    (voir documentation).
    2) Elle permet de soustraire une ligne de base et de normer le spectre.
    Plusieurs types de baseline sont proposés : 
        1 - Polyline 
        2 - Spline (ne fonctionne que si gcvspline est actif)
        3 - Unispline (peut remplacer approximativement spline)
        4 - Rubberband [default]. Cette option enregistre aussi la ligne de base 
        déterminée dans un fichier .csv séparé.
    
    Pensez à optimiser les valeurs des points d'accroche pour polyline et les roi 
    pour unispline si besoin.
    '''
    
    if Despik :
        Y_despiked, Despiking = suppression_spikes(X, Y, Legende)
    
    else :
        Y_despiked = Y ; Despiking=''
    
    Baseline_mode = int(input('Nous allons soustraire la ligne de fond.\n Pour choisir un traitement "polyline", tapez 1. \n Pour un traitement "spline", tapez 2. \n Pour "unispline", tapez 3. \n Pour "rubberband", tapez 4. [default]\n') or 4)    
    
    if Baseline_mode == 1 : 
        Baseline_type = 'polyline'
        Xcorr, Ycorr = traitement_Raman_portable_polyline(X, Y_despiked, Legende, borne=[400,1240], 
                                       pt = np.array([[400, 500, 750, 800, 1210, 1240]]))
    elif Baseline_mode == 2 :
        Baseline_type = 'gcvspline'
        Xcorr, Ycorr = traitement_Raman_portable_spline(X, Y_despiked, Legende, borne=[400,1240], roi = np.array([[400,540],[750, 780], [810, 820], [1210,1240]]))
    elif Baseline_mode == 3 :
        Baseline_type = 'unispline'
        Xcorr, Ycorr = traitement_Raman_portable_unispline(X, Y_despiked, Legende, borne=[400,1275],
                                      roi = np.array([[400, 450], [490, 540], [740, 750], [820, 840], [1210,1275]]))
    elif Baseline_mode == 4 : 
        Baseline_type = 'rubberband'
        Dossier_baseline = 'Extraction_baseline'
        Xcorr, Ycorr = baseline_rubberband_smooth_raman(X, Y_despiked, Legende,
                                                Dossier_baseline=Dossier_baseline)
        
    Xsave, Ysave, Normalisation_type = normalise_spectra(Xcorr, Ycorr, Legende)
    
    return(Xsave, Ysave, Despiking, Normalisation_type, Baseline_type)

def Nettoyage_spectre_Raman(Liste, Legende, Liste_ref, correction,
                      Dossier='./Data_corr_baseline',  AdditionY=0,
                      Borne=[304, 1400], valeurnorm=0, Indice_polym = False,
                      Liste_refN=''):
    """
    Parameters
    ----------
    Liste : TYPE
        DESCRIPTION.
    Legende : TYPE
        DESCRIPTION.
    Liste_ref : TYPE
        DESCRIPTION.
    correction : int
        0  Pas de correction
        1  Spectro portable baseline polyline + normalisation aire
        11 Spectro portable baseline polyline verre altÃ©rÃ© + normalisation aire
        2  Spectro portable baseline spline (si gcvspline activÃ©) + normalisation aire
        3  Spectro portable rubberband + normalisation aire
        31 Spectro portable rubberband pour altération + normalisation aire
        6  Permet de supprimer les spikes (optimisÃ© pour le spectro 23-13 3Ã¨me IMPMC)
        66 Permet de traiter les spectres acquis sur le spectro 23-13 3Ã¨me IMPMC
        (suppression des spikes puis suppression d'une baseline au choix 
         puis normalisation')
        99 Permet de faire la moyenne de 2 ou 3 spectres.
    Dossier : TYPE, optional
        DESCRIPTION. The default is './Data_corr_baseline'.
    AdditionY : TYPE, optional
        DESCRIPTION. The default is 0.
    Borne : TYPE, optional
        DESCRIPTION. The default is [304, 1400].
    valeurnorm : TYPE, optional
        DESCRIPTION. The default is 0.
    Indice_polym : TYPE, optional
        DESCRIPTION. The default is False.
    Liste_refN : TYPE, optional
        DESCRIPTION. The default is ''.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    None.

    """
    import os
    if os.path.exists("Indice_polymerisation.csv"):
        os.remove("Indice_polymerisation.csv")
   
    Liste_corr=[]
    ENTETE = '\n X en cm^-1; Y intensité (u.a)'
    RAJOUT = ''
    
    if np.size(AdditionY) == 1: AdditionY=(np.zeros(np.size(Liste)) + AdditionY)
    
    for i in np.arange(0, np.size(Liste), 1):
        #print(correction)
        
        if(correction[i]==0): # Si pas de correction on ne fait rien
            Liste_corr.append(Liste[i])
            
        else:
            if Dossier == '' : Dossier = Corr2folder_RAMAN(correction[i])
            
            Fichier=Liste[i];
            (cheminfichier, nomfichier) = os.path.split(Fichier)
            Fichier_corr=nomfichier[0:-4] + Corr2Str_RAMAN(correction[i])
            
            if(correction[i] == 6)  :
                X, Y = Readspectre(Fichier, delimiter=',', skip_header=1)
                Ysave, Despiking = suppression_spikes(X, Y, Legende[i])
                Xsave = X;
                
            elif(correction[i] == 66):
                try :
                    X, Y = Readspectre(Fichier, delimiter=',', skip_header=2)
                except :
                    X, Y = Readspectre(Fichier, delimiter='\t', skip_header=2)
                #delimiter = ',' -> Fonctionne sur les fichier .dat en sortie de spectro.
                #delimiter ='\t' -> Pour des fichiers déjà  despikées (6) par ex.
                Xsave, Ysave, Despiking, Normalisation_type, Baseline_type = traitement_Raman_IMPMC(X, Y, Legende[i], Despik=False)
                Fichier_corr=nomfichier[0:-4] + Corr2Str_RAMAN(correction[i])
                RAJOUT = Despiking+'_baseline_'+Baseline_type+'_normalise_'+Normalisation_type
            
            else :
                X, Y = Readspectre(Fichier, delimiter='\t')
                if(correction[i] == 1) :
                    Xcorr, Ycorr = traitement_Raman_portable_polyline(X, Y, Legende[i], Borne)
                    Xsave, Ysave = normalise_spectra(Xcorr, Ycorr, Legende[i])
                    
                elif(correction[i] == 11) :
                    Xcorr, Ycorr = traitement_Raman_portable_polyline(X, Y, Legende[i],
                                borne = [304, 1600], pt = np.array(
                                    [80, 100, 200, 300, 736, 843, 961, 1250, 1302, 1387, 1616, 1701]))
                    Xsave, Ysave = normalise_spectra(Xcorr, Ycorr, Legende[i])

                elif(correction[i] == 2) : 
                    Xcorr, Ycorr = traitement_Raman_portable_spline(X, Y, Legende[i], borne=[80,1700],
                                                  roi = np.array([[80,500],[720,740],[810, 830],[1300,1700]]))
                    Xsave, Ysave = normalise_spectra(Xcorr, Ycorr, Legende[i])
                    
                    
                elif(correction[i] == 3) :                     
                    Xcorr, Ycorr = baseline_rubberband_smooth_raman(X, Y, Legende[i], Limite=[500, 1200])
                    Xsave, Ysave, Normalisation_type = normalise_spectra(Xcorr, Ycorr, Legende[i], Normalisation_Mode='area')
                    Fichier_corr=nomfichier[0:-4] + Corr2Str_RAMAN(correction[i])
                
                elif(correction[i] == 31) :                     
                    Xcorr, Ycorr = baseline_rubberband_smooth_raman(X, Y, Legende[i], Limite=[500, 1600])
                    Xsave, Ysave, Normalisation_type = normalise_spectra(Xcorr, Ycorr, Legende[i], Normalisation_Mode='area')
                    Fichier_corr=nomfichier[0:-4] + Corr2Str_RAMAN(correction[i])
                    
                    
                elif(correction[i] == 99) : #moyenne
                    try : Xsave, Ysave = mean_spectra(X, Y, Legende[i], Liste_ref[i], Liste_refN[i])
                    except IndexError :  Xsave, Ysave = mean_spectra(X, Y, Legende[i], Liste_ref[i], '')
                  
                    
                else:
                    raise ValueError('La correction demandÃ©e n\'est pas implÃ©mentÃ©e, vÃ©rifier le fichier d\'input')
                        
            Liste_corr.append('./' + Dossier +'/' + Fichier_corr)
            
            SavCSV(Xsave, Ysave, Fichier_corr, Legende[i]+RAJOUT, Dossier,
                   ENTETE, delimiter='\t');
            
            if Indice_polym == True:
                Ip = Indice_polymerisation(Xsave, Ysave, Legende[i])
                
    return(Liste_corr)


#%% RPE

from Lecture_input import integrate2curve
from matplotlib.widgets import SpanSelector

def polynomial_fit(x, *coeffs):
    return sum(c * x**i for i, c in enumerate(coeffs))



def Nettoyage_spectre_RPE(Liste, Legende, Liste_ref,Correction, Frequence = 0,
                        valeurnorm = 0, Dossier='Data_corr') :
    
   '''
   Cette fonction pilote les différents corrections et modifcations de spectre.

    Parameters
    ----------
    Liste : TYPE
        DESCRIPTION.
    Legende : TYPE
        DESCRIPTION.
    Liste_ref : TYPE
        DESCRIPTION.
    correction : INT
            Entier qui permet de choisir la correction,
            1 pour soustraire le signal de la cavité
            2 pour  normaliser par les différents parametres d'aquisition
            21 pour normaliser par les paramètres d'aquisation sauf g'
            3 pour soustraire la ligne de base
            31 pour soustraire la ligne de base sur le signal intégrer avec une fit polynomiale d'ordre 6
            4 pour lisser
    '''
   #dataset_filenames = Liste
   #importer_factory = aspecd.io.ImporterFactory()
   #baseline_subtraction = aspecd.processing.BaselineCorrection()
   #baseline_subtraction.parameters = {"kind": "polynomial", "order": 0}
   #plotter = aspecd.plotting.SinglePlotter1D()
   
    
   import os
          
   Liste_g = Facteur_g(Liste, Legende, Frequence)
   Liste_corr=[]
   ENTETE = '\n X en G; Y intensité (u.a)'
   RAJOUT = '' 
   
   #if (Correction ==3) :
       #for idx, dataset_source in enumerate(dataset_filenames):
           #dataset = aspecd.dataset.Dataset()
          # importer = importer_factory.get_importer(dataset_source)
           #dataset.import_from(importer)
           #dataset.process(baseline_subtraction)
           #plot = dataset.plot(plotter)
           #saver = aspecd.plotting.Saver()
           #saver.filename = figure_filenames[idx]
           #plot.save(saver)
           #Fichier_corr=nomfichier[0:-4] + '_sub_baseline.asc'
           #Xsave, Ysave = X, Y_base_corr
           #RAJOUT = 'sub_baseline'
           
   for i in np.arange(0, np.size(Liste), 1):
       
       if(Correction[i]==0): # Si pas de correction on ne fait rien
            Liste_corr.append(Liste[i])
            
       else:
             Fichier=Liste[i];
             (cheminfichier, nomfichier) = os.path.split(Fichier)
             X, Y = Readspectre(Fichier,skip_header=3, delimiter='\t')
            
             if (Correction[i] == 1) :
                X_cav,Y_cav = Readspectre('./Data_trait/cavite_vide_4.asc',skip_header=3, delimiter='\t')
                #Gain_cavite = Releve_parametre(CLEFDETRIE='*cavite_vide_gain5x10^3*',CHEMIN='.', DOSSIER = 'Data_parametre')
                #gain_cav = Gain_cavite[0]
                #Y_cav_norm = (Y_cav/1e5)*Gain[i]
                Y_cav_norm = (Y_cav/3.17e4)*1.002374e+004 # Diviser par gain cavité et multiplier par le gain
                Fichier_corr=nomfichier[0:-4] + Corr2Str_RPE(Correction[i])
                Xsave, Ysave = Soustraction_interpolation(X_cav, Y_cav_norm, X, Y)
                RAJOUT = 'sous_cavite'
        
             elif (Correction[i]== 2) :
                 g=Liste_g[i]
                 Y_norm = Y/(valeurnorm*g)
                 #Y_norm = Y/(Facteur_cavite*gain*m*am*ns*np.sqrt(p))
                 Fichier_corr=nomfichier[0:-4] + Corr2Str_RPE(Correction[i])
                 Xsave, Ysave = X, Y_norm
                 RAJOUT = '_normalise'
                 
             elif (Correction[i]== 3) :
                 X_corr, Y_corr = traitement_RPE_spline(X,Y,Legende[i])
                 Fichier_corr=nomfichier[0:-4] + Corr2Str_RPE(Correction[i])
                 Xsave, Ysave = X, Y_corr
                 RAJOUT = '_sub_baseline'
                 
                 
             elif (Correction[i]== 31) :
                 X_corr, Y_corr = remove_polynomiale_baseline(X,Y)
                 Fichier_corr=nomfichier[0:-4] + Corr2Str_RPE(Correction[i])
                 Xsave, Ysave = X_corr, Y_corr
                 RAJOUT = '_sub_baseline'
                 
             elif (Correction[i]== 32) : # Retrait spline ligne de base
                 norm1temp = np.max(Y)
                 Y = Y / norm1temp # On normalise à 1 pour que le choixd e la spline ne merde pas.
                 X, Y = integrate2curve(X, Y)
                 X_corr, Y_corr = remove_spline(X, Y, legend=Legende[i], gradient=True)
                 Fichier_corr=nomfichier[0:-4] + Corr2Str_RPE(Correction[i])
                 
                 Y_corr = Y_corr * norm1temp
                 # Y_corr = np.gradient(Y_corr)
                 
                 Xsave, Ysave = X_corr, Y_corr
                 RAJOUT = '_sub_baseline_spline'
                 
                 
             elif(Correction[i] == 4): #Lissage
                 Xsave, Ysave= lissage (X, Y, Legende[i])
                 Fichier_corr=nomfichier[0:-4]+ Corr2Str_RPE(Correction[i])
                 
             elif(Correction[i]== 5) : #calibration en champ avec DPPH
                 X_DPPH, Y_DPPH = Readspectre('./Data_trait/DPPH.asc',skip_header=3, delimiter='\t')
                 i_max,i_min = np.argmax(Y_DPPH), np.argmin(Y_DPPH)
                 x_max, x_min = X_DPPH[i_max],X_DPPH[i_min]
                 res = (x_max+x_min)/2
                 #Param_DPPH = Releve_parametre(CLEFDETRIE='*DPPH*',CHEMIN='.', DOSSIER = 'Data_parametre')
                 #Frequence_DPPH= Param_DPPH[1]
                 #Frequence_DPPH = re.sub("\n", '', Frequence_DPPH)
                 #res_DPPH_théo = (1e4*6.626e-34*Frequence_DPPH)/(2.0036*9.274e-24)
                 res_DPPH_théo = 3518.10925
                 decalage = res-res_DPPH_théo
                 print(decalage)
                 Fichier_corr=nomfichier[0:-4] + Corr2Str_RPE(Correction[i])
                 X_corr, Y_corr = X-decalage, Y
                 Xsave, Ysave = removeInfNan(X_corr, Y_corr)
                 RAJOUT = '_calibre'
             else:
                 raise ValueError('La correction demandée n\'est pas implémentée, vérifier le fichier d\'input')
          
             Liste_corr.append('./' + Dossier +'/' + Fichier_corr)
          
             SavCSV(Xsave, Ysave, Fichier_corr, Legende[i]+RAJOUT, Dossier,
                 ENTETE, delimiter='\t');
         
   return(Liste_corr)
    

def Facteur_g(Liste, Legende, Frequence =0) :
    
    Liste_g = []
    for i in np.arange(0, np.size(Liste), 1):
        
        Fichier=Liste[i];
        (cheminfichier, nomfichier) = os.path.split(Fichier)       
        X, Y = Readspectre(Fichier,skip_header=3, delimiter='\t')
        #removeInfNan(X,Y)
        pos_res = np.where(Y == np.nanmax(Y))
        pos_res = np.argmax(Y)
        res = X[pos_res]
        g = (6.626e-34*Frequence)/(9.274e-24*res*1e-4)
        Liste_g.append(g)
    
    return(Liste_g)
        
    

def Integrale(Liste, Legende, x_min, x_max, Gain=0,Frequence =0, 
                     Puissance =0, Masse=0, Hauteur=0, Dossier='Data_Integ') :
    
    import os
    ENTETE = '\n X en G; Y intensité (u.a)'
    Liste_Integrale=[]
    
    for i in np.arange(0, np.size(Liste), 1):
        Fichier = Liste[i];
        (cheminfichier, nomfichier) = os.path.split(Fichier)
        
        X, Y = Readspectre(Fichier,skip_header=3, delimiter='\t')
        
        p_1 = (np.abs(X - x_min)).argmin()
        p_2 =(np.abs(X - x_max)).argmin()
        Y_int = [0]*(p_2-p_1)
        Y_int = [0]*(p_2-p_1)
        X_int = X[p_1:p_2]

        Y_int[p_1-1]=0
               
        for i in range (np.size(X_int)):
            Y_int[i] = Y_int[i-1] + (X[i] - X[i-1]) * (Y[i-1] + Y[i]) / 2
            #Integrale.append(Y_int[i])
        #int_f = interp1d(X_integ, Integrale) 
        
        plt.plot(X_int,Y_int)
        
        Fichier_int=nomfichier[0:-4] + '_int.asc'
          
        Liste_Integrale.append('./' + Dossier +'/' + Fichier_int)
       
        SavCSV(X_int, Y_int, Fichier_int, Legende[i]+'_Int', Dossier,ENTETE, delimiter='\t');
        
    return(Liste_Integrale)


def traitement_RPE_spline(X, Y, Legende, borne=[0,3200],
                              roi = np.array([[0,730],[1320,1360],[2350,3200]])):
    plt.style.use({'lines.linewidth': 1})  
    plt.figure(dpi=200)
    ax = plt.gca()
    ax.grid()
    
    
    Xcorr = np.arange(np.min(borne), np.max(borne), 0.5) # we generate the new X values with numpy.arange()
    Y_temp = rampy.resample(X, Y, Xcorr, fill_value="extrapolate")
    print(np.shape(Xcorr))
    print(np.shape(Y_temp))
    
    Ycorr, base_gcvspl = rampy.baseline(Xcorr,Y_temp,roi,'unispline',s=0.5 ) # activate if you have installed gcvspline
    
    Ycorr=Ycorr[:, 0]
    # doing the figure

    ax.plot(X, Y, "-", label=Legende)
    ax.plot(Xcorr, Ycorr, "g-", label="extract_unispline")
    ax.plot(Xcorr,base_gcvspl,"-",color="red",label="unispline baseline") # activate if you have installed gcvspline
    
    
    
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    
    ax.legend()
    #Sav_fig(Legende,'./Graph/Ligne_de_base')

    #Ycorr = rampy.normalise(Ycorr,x=X,method="area")

    return(Xcorr, Ycorr)    


#%% FLUO

def Nettoyage_spectre_Fluo(Liste, Legende, Liste_ref, correction, Liste_refN='',
                      Dossier='',  Addition_Tr=0,
                      Borne=[nm2cm1(2500), nm2cm1(330)], valeurnorm=0, DataExcel = False):
    '''
    Cette fonction pilote les différents corrections et modifcations de spectre.

    Parameters
    ----------
    Liste : TYPE
        DESCRIPTION.
    Legende : TYPE
        DESCRIPTION.
    Liste_ref : TYPE
        DESCRIPTION.
    correction : INT
    Liste_refN : TYPE, optional
        DESCRIPTION. The default is ''.
    Dossier : STR, optional
        Repertoir de sauvegarde des data_traité. The default is './Data_corr'.
    Addition_Tr : TYPE, optional
        DESCRIPTION. The default is 0.
    Borne : TYPE, optional
        Borne en cm-1. The default is [nm2cm1(2500), nm2cm1(330)].

    Raises
    ------
    ValueError
        DESCRIPTION.
        cor 77

    Returns
    -------
    La liste des fichier corriger.

    '''
    Liste_corr=[]
    ENTETE = '\n X en nm; Y en %T'
    RAJOUT = ''


    Borne = np.array(Borne)
    
    if np.size(Addition_Tr) == 1: Addition_Tr=(np.zeros(np.size(Liste)) + Addition_Tr)
    
    for i in np.arange(0, np.size(Liste), 1):
        #print(correction)
        
        if(correction[i]==0): # Si pas de correction on ne fait rien
            Liste_corr.append(Liste[i])
            
        else:
            Fichier=Liste[i];
            (cheminfichier, nomfichier) = os.path.split(Fichier)
           
            try :
                df = pd.read_csv(Fichier, delimiter='\t', header=None, skiprows=2, usecols = [6, 7], index_col=0)
            except ValueError :
                df = pd.read_csv(Fichier, delimiter=',', header=None,
                            skiprows=23, usecols = [0, 1], index_col=0)
       
            Xnm = df.index.to_numpy()
            Ytr = df.iloc[:, 0].to_numpy()
            
            Y=Ytr
            X=nm2cm1(Xnm)
            Xsave=Xnm;
            
            if Dossier == '' : Dossier = Corr2folder_OAS(correction[i])
            
            if(correction[i] == 77): # Soustraction spline
                X, Y_corr = remove_spline(X, Y, legend=Legende[i], gradient=False)
 
                Fichier_corr=nomfichier[0:-4] + Corr2Str_OAS(correction[i])
                #Fichier_corr= Legende[i]+'_diff.csv'
           
            
            else:
                raise ValueError('La correction demandée n\'est pas implémentée, vérifier le fichier d\'input')
            
            Liste_corr.append(Dossier+os.sep+Fichier_corr);
            
            if not (correction[i] == 100000):
                Y_save = Y_corr
                # Y_save=np.power(10, -Y_corr); # On repasse en %T pour assurer la comptabilité avec le reste des scripts

            # if DataExcel:
            #     save = pd.DataFrame(np.array([X, Y, -np.log(Y)]), columns = ['nm', 'T', 'A'])
            #     save.to_excel(Fichier_corr+'.xlsx')
                
            # else :
            SavCSV(Xsave, Y_save, Fichier_corr, Legende[i]+RAJOUT, Dossier, ENTETE);

    return(Liste_corr)

#%% Fonction dynamique 
from scipy.interpolate import CubicSpline

def remove_spline(x, y, legend='', gradient=False):
    '''
    Cette fonction retire un spline cubic, le choix des points se fait à la souris
    Code en partie générer à l'aide de chatgpt

    Parameters
    ----------
    x : TYPE
        DESCRIPTION.
    y : TYPE
        DESCRIPTION.
    legend : TYPE, optional
        DESCRIPTION. The default is ''.
    
    gradient: true pour activer l'affichage de la dérivié du signal soustrait dela ligne de base
    Returns
    -------
    Xcorr, Ycorr : courbe avec la spline soustraite
    '''
    # Créer une figure
    if gradient : 
        fig, axs = plt.subplots(1, 2, sharey = False, figsize=(10, 6))
    else :
        fig, axs = plt.subplots(1, 2, sharey = True, figsize=(10, 6))
    
    fig.suptitle('Cliquez avec le bouton gauche pour ajouter un point, avec le bouton droit pour le retirer'
                +'\n q ou entrer pour valider la selection')
    
    ax1 = axs[0]
    ax2 = axs[1]
        
    
    # Créaction des variables partagée
    Xcorr = None
    Ycorr = None
    selection_done = False
    
    # Tracer la courbe originale sans marqueurs
    line, = ax1.plot(x, y, 'b-', label=legend)
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')

    # Initialiser une spline cubique à 0
    x_cs = np.linspace(np.min(x), np.max(x), 1000)
    y_cs = np.zeros_like(x_cs)  # Spline initialisée à 0
    spline, = ax1.plot(x_cs, y_cs, 'r-', label='Spline cubique (ligne de base)')
    linecorr, = ax2.plot(x_cs, y_cs, label = legend + 'corrigée')


    # Liste pour stocker les indices des points à conserver
    selected_indices = []
    
    # Marqueurs pour les points de la spline sur la courbe originale
    spline_points, = ax1.plot([], [], 'ro', label='Points de la spline')
    
    # Fonction pour ajouter ou retirer un point et mettre à jour la spline
    def update_point(event):
        nonlocal Xcorr, Ycorr
        
        if event.inaxes == ax1 : 
            x_data, y_data = line.get_data()
            x_click, y_click = event.xdata, event.ydata
        
            if event.button == 1:  # Clic gauche pour ajouter un point
                # Vérifier si le point a déjà été ajouté
                if np.argmin(np.abs(x_data - x_click)) not in selected_indices:
                    # Trouver l'indice le plus proche du point cliqué
                    distances = np.sqrt((x_data - x_click) ** 2 + (y_data - y_click) ** 2)
                    index_to_add = np.argmin(distances)
        
                    # Ajouter le point sélectionné à la liste
                    selected_indices.append(index_to_add)
        
            elif event.button == 3:  # Clic droit pour retirer un point
                if len(selected_indices) > 0:
                    # Trouver l'indice le plus proche du point cliqué
                    distances = np.sqrt((x_data[selected_indices] - x_click) ** 2 + (y_data[selected_indices] - y_click) ** 2)
                    index_to_remove = selected_indices[np.argmin(distances)]
        
                    # Retirer le point sélectionné de la liste
                    selected_indices.remove(index_to_remove)
            # Affichage point cliqué
            spline_points.set_data(x_data[selected_indices], y_data[selected_indices])
    
            # Recalculer la spline cubique mise à jour
            if len(selected_indices) >= 2:
                selected_x = x_data[selected_indices]
                selected_y = y_data[selected_indices]
                # Trier les points par ordre croissant de x
                sorted_indices = np.argsort(selected_x)
                selected_x = selected_x[sorted_indices]
                selected_y = selected_y[sorted_indices]
                cs = CubicSpline(selected_x, selected_y)
                y_cs = cs(x_cs)
        
                # Mettre à jour les données affichées
                # line.set_data(x_data, y_data)
            
                # Mettre à jour la ligne de base (spline cubique)
                spline.set_data(x_cs, y_cs)
                Xcorr = x_data
                Ycorr = y_data-cs(x_data)
                if gradient :
                    Ycorr = np.gradient(Ycorr)
                    ax2.set_ylim([np.min(Ycorr)*1.2, np.max(Ycorr)*1.2,])
                
                linecorr.set_data(Xcorr, Ycorr)
                # ax2.autoscale()
                
                # Mettre à jour les points de la spline sur la courbe originale
    
            plt.draw()
        
    def on_key(event): # Ce qui se passe quand on click un truc sur le clavier
        nonlocal selection_done
        if event.key in ['enter', 'q']:
            # Arrêter la boucle d'événements
            # plt.close()
            selection_done = True
            plt.close()
    
            
    def on_close(event):# NE FONCTIONNE PAS, à investiguer
        nonlocal selection_done
        selection_done = True  
        
    fig.canvas.mpl_connect('key_press_event', on_key)
    fig.canvas.mpl_connect('close_event', on_close) # NE FONCTIONNE PAS, à investiguer

    
    # Connecter la fonction update_point à l'événement de clic
    fig.canvas.mpl_connect('button_press_event', update_point)
    # Afficher la légende
    ax1.legend()
    ax2.legend()
    plt.show()
     
    while not selection_done : # Blocage car dans spyder block = True ça fait n'impt
        # print('a')
        fig.waitforbuttonpress()  # attente pour que matplolib soit le process principal
        
    return(Xcorr, Ycorr)

def remove_polynomiale_baseline(x_data, y_data, integrate=True) :
    '''
    Cette fonction retire un polynome en selectionnant la zone de fit à la souris.

    Parameters
    ----------
    x_data : TYPE
        DESCRIPTION.
    y_data : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    
    ycorr = None
    
    if integrate :
        x_data, y_data = integrate2curve(x_data, y_data)
    
    # Création de la figure initiale
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 6))
    axmain = ax1.twinx() # POur le span selector, l'axe de selection doit être dessus, sinon pas d'accés
    # axmain sert juste pour le span selector, c'est un axe fantome.
    
    axmain.plot(x_data, y_data, 'w', label='Données', visible=False) # Créaction d'un ligne invisible qui servira a voir le pan selector au premier plan sans que le clear ax le suprimer
    axmain.set_xlabel('X')
    axmain.set_ylabel('Y')
    
    fig.suptitle('Sélection de deux zones de fit à la souris (appuyer sur enter pour finir, r pour clean subplot 2)')

    
    # Initialisation des listes pour stocker les zones sélectionnées
    selected_zones = []
    selection_done = False
    
    # Fonction pour ajuster les deux zones sélectionnées
    def increment_counter(reset=False): # Fonction qui count de manière persistente
        if not hasattr(increment_counter, "counter") or reset :
            increment_counter.counter = 0
        increment_counter.counter += 1
        return str(increment_counter.counter)
        
    
    def fit_selected_zones(xmin, xmax):
        nonlocal selected_zones
        nonlocal selection_done
        nonlocal ycorr
        
        selected_zones.append((xmin, xmax)) 
     
        if len(selected_zones) == 2 :
            zone1 = selected_zones[0]
            zone2 = selected_zones[1]
    
            x_selected1 = x_data[(x_data >= zone1[0]) & (x_data <= zone1[1])]
            y_selected1 = y_data[(x_data >= zone1[0]) & (x_data <= zone1[1])]
    
            x_selected2 = x_data[(x_data >= zone2[0]) & (x_data <= zone2[1])]
            y_selected2 = y_data[(x_data >= zone2[0]) & (x_data <= zone2[1])]
    
            x_all = np.concatenate([x_selected1, x_selected2])
            y_all = np.concatenate([y_selected1, y_selected2])
            # Ajustement polynomial d'ordre 6 pour les deux zones
            coeffs1, _ = curve_fit(polynomial_fit, x_all, y_all, p0=np.ones(7))
    
    
            # Génération des courbes ajustées pour les deux zones
            y_fit1 = polynomial_fit(x_data, *coeffs1)
    
            # Affichage des zones sélectionnées et des ajustements
            count = increment_counter() # Count le nombre de fois ou un ajustement à été fais
            ax1.plot(x_selected1, y_selected1, '--r', label='Zone 1 sélectionnée')
            ax1.plot(x_selected2, y_selected2, '--b', label='Zone 2 sélectionnée')
            ax1.plot(x_data, y_fit1, '--', label='Ajustement'+count)
            ax1.legend()
            
            ycorr = y_data - y_fit1
            if integrate : ycorr = np.gradient(ycorr) # Si on retire la ligne de base sur la courbe dérivé
            ax2.plot(x_data, ycorr, label=count)
            ax2.legend()
            plt.draw()
    
            selected_zones =[]

    def start():
        '''
        Fonction de base qui affiche les données sur ax1 après avoir clear (axmain n'est pas touché')

        Returns
        -------
        None.

        '''
        ax1.clear()
        ax2.clear()
        ax1.plot(x_data, y_data, label='Données')
        increment_counter(True)
        plt.plot()
        
    def on_key(event): # Ce qui se passe quand on click un truc sur le clavier
        nonlocal selection_done
        if event.key in ['enter', 'q']:
            # Arrêter la boucle d'événements
            # plt.close()
            selection_done = True
            plt.close()
        elif event.key == 'r':
            start()

            
    def on_close(event):# NE FONCTIONNE PAS, à investiguer
        nonlocal selection_done
        selection_done = True

    
    fig.canvas.mpl_connect('key_press_event', on_key)
    fig.canvas.mpl_connect('close_event', on_close) # NE FONCTIONNE PAS, à investiguer

    start() # Initialisation du graph
    # Création d'un objet SpanSelector pour la sélection de zones
    span_selector = SpanSelector(axmain, fit_selected_zones,
                                 'horizontal', useblit=False, props=dict(alpha=0.5, facecolor='blue'),
                                 interactive=False)
    
    while not selection_done : # Blocage car dans spyder block = True ça fait n'impt
        # print('a')
        fig.waitforbuttonpress()  # attente pour que matplolib soit le process principal

    return(x_data, ycorr)


