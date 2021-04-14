# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 11:17:45 2020

@author: Theo_C
"""

import numpy as np
import os
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d
from scipy import interpolate
from IPython import display

from Affichage_spectre import Readspectre

def show_interactif_loop():
    '''
    équivalent de plt.show()
    Sensé permettre l'affichage à la matlab dans une boucle

    Returns
    -------
    None.

    '''
    # display.clear_output(wait=True)
    # display.display(plt.gcf())
    plt.show(block=True)
      
def fit_OMCT(X, Y, NOM): 
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

    nbonde_min_fit=27000; # limite basse de la partie des data qui vont être ajuster
    nbonde_max_fit=31000; # limite haute
    
    Intes_fit=500; # intensité de la gausienne pour l'ajustement
    Locali_fit=42000; #position central
    Sigma_fit = 3500; # écart type Gauss
    Shift_fit = 0.1; # shift de la gausienne
    
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
        
        def Gauss(X, A1, pos1, sigma1, b): # Note H largeur à mis hauteurs = 2.3548*sigma
            RES=A1*np.exp(-(X-pos1)**2/(2*sigma1**2))+b
            return(RES)
        
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
            
            Y_corr=Y_ori-Y_fit;
            plt.plot(X_ori,Y_ori,label='Origine')
            plt.plot(X_ori, Y_corr, '--', label='corrigé')
            plt.xlim(4000, 31000);
            plt.ylim(-0.1, 1)
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
    #return(Y_corr-Y_min) # On récupère la valeur minmum pour bien assoir la ligne de base 
    #! pas de sens physique, valeur de correction de reflexion à mesure
    return(Y_corr) #On retourne le spectre corrigé

def filtrage_gauss(X, Y, NOM, zoomfig ='auto', sigma=2):
    '''
    Cette fonction sert à débruité un signal avec un filtrage gaussien (convolution par une petite gausienne,
            ce qui équivaux dans l'espace de fourrier à une multiplication par une large gausienne')
    ==> on élimine le bruit haute fréquence.

    Parameters
    ----------
    X : float
        en cm^-1.
    Y : float
        Signal en absorbance.
    NOM : string
        nom à afficher.
    Sigma : float
        Valeur par défaut du filtre gaussien

    Returns
    -------
    Ygauss : le signal filtré

    '''
    Amin=0.2
    Amax=1
    Xmin=4000
    Xmax=32000
    Zoom=0;
    
    #plt.ion() #Nécéssaire pour afficher les figures en %matplolib  
    while True:
        plt.figure(figsize=(5,3), dpi=180)
        Ygauss=gaussian_filter1d(Y, sigma)
    
        plt.plot(X, Y, label='Original')
        plt.plot(X, Ygauss, label='filtrée Gauss')
        if zoomfig != 'auto' :
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
       
        if Zoom and (zoomfig != 'auto') :
            Xmin = float(input('Valeur min nb onde en cm-1 (default = '+
                                     str(Xmin) + ')\n') or Xmin);
            Xmax = float(input('Valeur max nb onde en cm-1 (default = '+
                                     str(Xmax) + ')\n') or Xmax);
            Amin = float(input('Valeur min abs (default = '+
                                     str(Amin) + ')\n') or Amin);
            Amax = float(input('Valeur max abs (default = '+
                                     str(Amax) + ')\n') or Amax);
            
        sigma = float(input('Rentrer la valeur sigma du filtre gaussien (default = '+
                                 str(sigma) + ')\n') or sigma);
        
    
    return(Ygauss)

def Remontage_IR_VIS(X_IR, Y_IR, X_VIS, Y_VIS, mode='Perkin', X=0,Y=0,NOM='pouet'):
    '''
    Cette fonction sert à remonter deux partie de courbe au même niveau
    car issu de deux détecteur/spectro différent 

    Parameters
    ----------
    X_IR : TYPE
        DESCRIPTION.
    Y_IR : TYPE
        DESCRIPTION.
    X_VIS : TYPE
        DESCRIPTION.
    Y_VIS : TYPE
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
    Y_IRDelta : TYPE
        Y remonté.
    interpo_limIR : TYPE
        limite pour les IR.
    interpo_limVIS : TYPE
        limite pour le VIS.

    '''
    #plt.ion() #Nécéssaire pour afficher les figures en %matplolib  
    
    Xmin=5000; # Affichage en x min en cm-1
    Xmax=20000;# Affichage en x max en cm-1
    point_jonction = 11000;
    
    if mode == 'Perkin':
        limIR= 9000;   # limite pour l'ajustement IR
        limVIS = 14000; # 
    elif mode == 'Portable':
        limIR= 9000;   # limite pour l'ajustement IR
        limVIS = 12000; # 
        limIR_haute = 12000;
        limVIS_basse = 10000;



    if mode == 'Portable':
        X = np.concatenate([X_IR, X_VIS]);
        Y = np.concatenate([Y_IR, Y_VIS]);
    
    FIGsize=(10,6)
    DPI=120
    
    plt.figure(figsize=FIGsize, dpi=DPI)
    plt.plot(X, Y, label=NOM)
    plt.xlim([Xmin, Xmax])
    plt.grid()
    show_interactif_loop()
    
    while True: # Selection de la zone d'interpolation
        limIR=int(input('Rentrer la limite partie basse IR pour l ajustement (defaut ='
                           + str(limIR) +')\n') or limIR);
        limVIS=int(input('Rentrer la limite partie haute VIS pour l ajustement (defaut ='
                           + str(limVIS) +')\n') or limVIS);
        
        interpo_limIR = X_IR > limIR
        interpo_limVIS = X_VIS < limVIS 

        
        if mode == 'Portable':
             limIR_haute=int(input('Rentrer la limite partie HAUTE IR pour l ajustement (defaut ='
                           + str(limIR_haute) +')\n') or limIR_haute);
             limVIS_basse=int(input('Rentrer la limite partie BASSE VIS pour l ajustement (defaut ='
                           + str(limVIS_basse) +')\n') or limVIS_basse);       
             interpo_limIR = np.logical_and(X_IR > limIR,  X_IR < limIR_haute)
             interpo_limVIS = np.logical_and(X_VIS < limVIS, X_VIS > limVIS_basse)

        plt.figure(figsize=FIGsize, dpi=DPI)
        
        plt.plot(X, Y, label='Original')
        plt.plot(X_IR[interpo_limIR], Y_IR[interpo_limIR], label='Y_IR_ajustement');
        plt.plot(X_VIS[interpo_limVIS], Y_VIS[interpo_limVIS], label='Y_VIS_ajustement');
        plt.xlim([Xmin, Xmax])
        plt.grid()
        plt.legend()
        show_interactif_loop()
                
        if(int(input('Est-ce correcte ?\n'
                     + '1 pour oui et continuer la correction, 0 pour non et recommencer \n') or 1)) : break
    
    Remonte_OK=0;
    Remonte = 1; # Par défaut fit linéaire
    
    # On mets la premièrer partie de la courbe au niveau de la seconde,
    # car le Photomultiplicateur du visible donne un valeur "vrai" et on corrige l'InGaAs dessus
    
    Delta=np.mean(Y_VIS[:3])-np.mean(Y_IR[-3:]) # Decallage entre les deux courbe pour le cas 2
    ordre_fit = 1

    while (Remonte_OK == 0):  #Partie ou on remonte la partie gauche (les IR) des données pour coller l'InGaAs sur le PM  
        
        Remonte = int(input('Selectionner la methode de remonté de l\InGaAs sur le PM :\n 1 fit polynomiale\n'
                      +' 2 réhausser à la main\n'
                      +' 3 fit des deux cotés du même ordre (ordre 0 = valeur moyenne) \n')
                        or Remonte)
        
        if(Remonte == 1):
            try:                   
                ordre_fit=float(input('rentrer l\'ordre du polynome qui va servir à ajuster la remontée\n') or ordre_fit)
    
                fit_delta = np.polyfit(X_VIS[interpo_limVIS], Y_VIS[interpo_limVIS], ordre_fit) #fit pour remonter la courbe
                
                regY_VIS= np.poly1d(fit_delta)
                
                #Delta=regY_VIS(X_IR[-1])-Y_IR[-1]; #
                Delta=regY_VIS(np.mean(X_IR[-3:]))-np.mean(Y_IR[-3:]); #
                
                Y_IRDelta=Y_IR+Delta
                plt.figure(figsize=FIGsize, dpi=DPI)
                plt.plot(X_VIS, regY_VIS(X_VIS), '-', label='fit pour la remonté')
            
            except LinAlgError :
                print('Le fit ne converge pas')
        
        elif(Remonte == 2):

            Delta=float(input('Rentrer le Delta entre les deux partie (defaut =' + str(Delta) + ')\n') or Delta)
            Y_IRDelta=Y_IR+Delta
            plt.figure(figsize=FIGsize, dpi=DPI)
            
        elif(Remonte == 3) :
            try:
                ordre_fit=float(input('rentrer l\'ordre du polynome qui va servir à ajuster la remontée\n') or ordre_fit)
                
                if ordre_fit !=0:
                    point_jonction=float(input('Rentrer la valeur en X du point de jonction (default = '
                                               + str(point_jonction)+ ')\n') or point_jonction)
                
                fit_VIS = np.polyfit(X_VIS[interpo_limVIS], Y_VIS[interpo_limVIS], ordre_fit) #fit pour remonter la courbe
                regY_VIS= np.poly1d(fit_VIS)
                
                fit_IR = np.polyfit(X_IR[interpo_limIR], Y_IR[interpo_limIR], ordre_fit) #fit pour remonter la courbe
                regY_IR= np.poly1d(fit_IR)
                
                if ordre_fit == 0:
                    Delta=regY_VIS(X_IR[-1])-regY_IR(X_VIS[-1]); #
                else:
                     Delta=regY_VIS(point_jonction)-regY_IR(point_jonction); #
              
                Y_IRDelta=Y_IR+Delta
                
                plt.figure(figsize=FIGsize, dpi=DPI)
                plt.plot(X_VIS, regY_VIS(X_VIS), '--', label='fit VIS pour la remonté')
                plt.plot(X_IR, regY_IR(X_IR), '--', label='fit IR pour la remonté')
            except:
                print("Ne marche pas, utilise une autre methodes\n")
                Y_IRDelta=Y_IR

        plt.plot(X, Y, label=NOM + 'Original')
        plt.plot(X_IR, Y_IRDelta, '-', label='1ere parti remonté')
        plt.plot(X_VIS, Y_VIS, '-', label='2eme parti')
        plt.xlim([Xmin, Xmax])
        #plt.ylim([Ymin, Ymax])
        plt.grid()
        plt.legend()
        show_interactif_loop()
        
        Remonte_OK = int(input('La correction est-elle correcte ? \n') or 0)
 
    return (Y_IRDelta, interpo_limIR, interpo_limVIS)

def reconstruction_saut(X_IR, Y_IRDelta, X_VIS, Y_VIS, interpo_limIR, interpo_limVIS, Xinterpo, X, Y, NOM):
    '''
    Permet de reconstruire des données trops bruité et donc supprimés.

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
    Xmin=5000; # Affichage en x min en cm-1
    Xmax=20000;# Affichage en x max en cm-1

    Forme_OK = 0;
    Forme = 1;
    #plt.ion() #Nécéssaire pour afficher les figures en %matplolib   

    Forme = int(input('Selectionner la methode d\'interpolation\n 1 pour un ajustement polynomiale\n')
                        or Forme)
    FIGsize=(10,6)
    DPI=120

    while (Forme_OK == 0): # Partie ou on selectionne l'interpolation
        if (Forme == 1): # si jamais besoin de mettre une boucle

            ordre_fit = 0
            
            # On mets la premièrer partie de la courbe au niveau de la seconde,
            # car le Photomultiplicateur du visible donne un valeur "vrai" et on corrige l'InGaAs dessus
            
            ordre_fit=float(input('rentrer l\'ordre du polynome qui va servir à faire l\'interpolation\n')
                                or ordre_fit)
   
            Xtemp=np.concatenate([X_IR[interpo_limIR], X_VIS[interpo_limVIS]])
            # On rassemble les deux jeux de donner pour interpoler
            Ytemp=np.concatenate([Y_IRDelta[interpo_limIR], Y_VIS[interpo_limVIS]]) 
            
            fit_interpo = np.polyfit(Xtemp, Ytemp, ordre_fit) # fit pour interpoler les data
            fonctionYinterpo = np.poly1d(fit_interpo)
            Yinterpo = fonctionYinterpo(Xinterpo);

        else:
            Yinterpo=np.ones(np.size(Xinterpo))*-1;
            print("Correction non implémentée\n")
            break
        
        plt.figure(figsize=FIGsize, dpi=DPI)
        plt.plot(X, Y, label= NOM + 'Original')
        plt.plot(X_IR, Y_IRDelta, '-', label='1ere parti remonté')
        plt.plot(X_VIS, Y_VIS, '-', label='2eme parti')
        plt.plot(Xinterpo, Yinterpo, label='interpolation')
        plt.xlim([Xmin, Xmax])
        #plt.ylim([Ymin, Ymax])
        plt.grid()
        plt.legend()
        show_interactif_loop()
        
        if(int(input('L\'interpolation est-elle correcte ?\nO pour non, 1 pour oui\n') or 0)): break
    
    Ycorr=np.concatenate([Y_IRDelta, Yinterpo, Y_VIS])
    return(Ycorr)
    
def correction_saut_detect(X, Y, NOM):
    '''
    Cette fonction sert à corriger le saut de detecteur des PERKIN, plusieurs type d'ajustement sont proposé'

    Parameters
    ----------
    X : tableau de float
        en nombre d'onde.
    Y : tableau de float
        En absorbance.
    NOM : str
        Nom de ce qui est corrigé.

    Returns
    -------
    Ycorr : TYPE
        Correction.

    '''
    
    Xmin=5000; # Affichage en x min en cm-1
    Xmax=20000;# Affichage en x max en cm-1
    
    Xdebcoup = 11500; # Nombre d'onde en cm^-1 à partir duquel on suprime les valeurs du saut
    Xfincoup = 12500; # Nombre d'onde en cm^-1 jusqu'auquel on suprime les valeurs du saut
    
    sep_OK = False;
    #plt.ion() #Nécéssaire pour afficher les figures en %matplolib  
    FIGsize=(10,6)
    DPI=120

    while True :
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
            plt.plot(X, Y, label='Original')
            plt.plot(XR1, YR1, '-', label='1ere parti')
            plt.plot(XR2, YR2, '-', label='2eme parti')
            plt.xlim([Xmin, Xmax])
            plt.grid()
            plt.legend()
            plt.show()
            
            print(NOM + '  :   la séparation est-elle correct ?\n');
            try:
                sep_OK = int(input('Rentrez 0 si non, 1 si oui\n'))
            except ValueError:
                print('Vous n\'avez pas rentré un nombre, la séparation est considéré comme mauvaise')
                sep_OK=False;
    
            if sep_OK : break
        
            Xdebcoup = float(input('Nombre d\'onde à partir duquel on suprime (default = '+
                                          str(Xdebcoup) + ')\n') or Xdebcoup);
        
            Xfincoup = float(input('Nombre d\'onde jusqu\'auquel on suprime (default = '+
                                          str(Xfincoup) + ')\n') or Xfincoup);

        if (sep_OK) : 
            Y_IRDelta, interpo_limIR, interpo_limVIS = Remontage_IR_VIS(XR1, YR1, XR2, YR2,
                                                                        'Perkin', X,Y, NOM)
            Ycorr=reconstruction_saut(XR1, Y_IRDelta, XR2, YR2, interpo_limIR,
                                     interpo_limVIS, Xinterpo, X, Y, NOM);
        
        
        
        plt.figure(figsize=FIGsize, dpi=DPI)

        plt.plot(X, Y, label= NOM + 'Original')
        plt.plot(X, Ycorr, label = NOM +'reconstruite')
        plt.xlim([Xmin, Xmax])
        #plt.ylim([Ymin, Ymax])
        plt.grid()
        plt.legend()
        display.clear_output(wait=True)
        display.display(plt.gcf())
        plt.show(block=True)
        if(int(input('La correction est-elle correcte ?\nO pour non et recommencer depuis le debut, 1 pour oui\n') or 1)): break
        
    #Ycorr=np.concatenate([YR1Delta, Yinterpo, YR2])
    return (Ycorr)

def Recollage_IR_VIS(X_IR, Y_IRDelta, X_VIS, Y_VIS, interpo_limIR, interpo_limVIS, NOM):
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
        
        plt.figure(figsize=FIGsize, dpi=DPI)
        plt.plot(Xcorr,Ycorr, label='Spectre joint')
        plt.grid()
        plt.legend()
        show_interactif_loop()
            
        Recol_OK = int(input('La jonction est-elle correcte ? \n') or 0)
    return(Xcorr, Ycorr)

def Traitement_spectro_portable(CHEMIN_IR, CHEMIN_VIS, NOM='pouet', Addition_Tr='0'):
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

    Returns
    -------
    None.

    '''
    #print(Addition_Tr)
    #plt.ion() #Nécéssaire pour afficher les figures en %matplolib  

    COUPURE_UV=350;
     
    FIGsize=(10,6)
    DPI=120

    Trait_OK = False;
    
    
    Xnm_IR, Ytr_IR=Readspectre(CHEMIN_IR) # Lecture fichier IR

    X_IR = (1/(Xnm_IR*1E-7))[::-1]; #
    Y_IR = (- np.log10(Ytr_IR+Addition_Tr))[::-1];# Si data en %T
    
    
    Xnm, YTr=Readspectre(CHEMIN_VIS) # Lecture fichier VIS

    INDEX_UV = Xnm > COUPURE_UV;
    
    Xnm=Xnm[INDEX_UV]
    YTr=YTr[INDEX_UV]
    
    X_VIS = (1/(Xnm*1E-7)); # [::-1] = inversion des valeurs dans le tableau
    Y_VIS = (-np.log10(YTr+Addition_Tr));# Si data en %T
    
    
    while not Trait_OK:
        Y_IRDelta, interpo_limIR, interpo_limVIS = Remontage_IR_VIS(X_IR, Y_IR, X_VIS, Y_VIS,
                                                                    'Portable', X=0,Y=0,NOM=NOM)
       
        Xcorr, Ycorr = Recollage_IR_VIS(X_IR, Y_IRDelta, X_VIS, Y_VIS, interpo_limIR, interpo_limVIS, NOM)
    
        plt.figure(figsize=FIGsize, dpi=DPI)
        
        plt.plot(X_IR, Y_IR, label='IR_ori');
        plt.plot(X_VIS, Y_VIS, label='VIS_ori');
        plt.plot(Xcorr,Ycorr, label='Spectre traité')
        plt.grid()
        plt.legend()
        show_interactif_loop()
            
        Trait_OK = int(input('Le traitement final est-il correcte ? \n') or 0)
        
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

def Soustraction_interpolation(X1, Y1, X2, Y2):
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
        Soustraction entre Y2 et Y1.
    Xref:
        Gamme des X communes au deux courbe

    '''
    plt.figure(figsize=([6,3]), dpi=120)
    plt.plot(X1,Y1, '.', label='DATA 1', markersize=1)
    plt.plot(X2, Y2, '.', label='DATA 2', markersize=1)
     
    if (min(X1)<min(X2)):
        if (max(X1)>max(X2)): # Les valeur de X2 sont strictement comprise dans X1
            f = interpolate.interp1d(X1, Y1)
            Xref=X2;
            Y1=f(Xref);
            
        else: #X1 va plus bas en valeur de que X2 mais X2 va plus haut, on coupe donc X2
            coupure_sup = (X2 <= max(X1))    
            X2=X2[coupure_sup]
            Xref=X2;
            
            Y2=Y2[coupure_sup]
            
            f = interpolate.interp1d(X1, Y1)
            Y1=f(Xref);
    
    else: # min X2 > min X1
        if (max(X1)>max(X2)): # 
            coupure_sup = (X1 <= max(X2))    
            X1=X1[coupure_sup]
            Xref=X1;
            
            Y1=Y1[coupure_sup]
            f = interpolate.interp1d(X2, Y2)
            Y2=f(Xref)
            
        else: # Les valeurs de X1 sont comrpise dans X2.
            f = interpolate.interp1d(X1, Y1)
            Xref=X1;
            Y1=f(Xref);
    
    Y_diff=Y2-Y1;
    #Y_diff=-Y_diff
    plt.plot(Xref, Y_diff, '^', label="différence", markersize=1)
    plt.legend()
    plt.show
    
    return(Y_diff, Xref)

def Corr2Str(Correction_number):
    '''
    Cette fonction sert à donner l'extension d'un fichier en fonction de la correction demandé

    Parameters
    ----------
    Correction_number : TYPE
        DESCRIPTION.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
        Nom de la correction

    '''
    # Correction_number = int(Correction_number)
    # print(Correction_number)
    if Correction_number == 1:
       Correction_NAME = '_cor_I100.csv'
       
    elif Correction_number == 2:
       Correction_NAME = '_cor_UV.csv'
       
    elif Correction_number == 3:
       Correction_NAME = '_cor_filtre_gauss.csv'
       
    elif Correction_number == 4:
       Correction_NAME = '_cor_I100_filtre_saut_filtre.csv'
       Correction_NAME = '_cor_I100_saut.csv'
   
    elif Correction_number == 5:
       Correction_NAME = '_Tr.csv'
       #Correction_NAME = '.csv'
    
    elif Correction_number == 6:
        Correction_NAME = '_jointVIS.csv'
    
    elif Correction_number == 7:
        Correction_NAME = '_ABScm.csv'
    
    elif Correction_number == 8:
        Correction_NAME = '_cor_DO.csv'
       
    elif Correction_number == 0:
        Correction_NAME=''
        
    else :
        raise ValueError("Pb la correction n'existe pas")
    
    return(Correction_NAME)

def SavCSV(X, Y, NOM, Legende, Dossier='./Data_corr', ENTETE='\n X en nm; Y en %T'):
    try:
        os.mkdir(Dossier)
    except OSError:
        pass
    
    HEADER =NOM+ ';' +Legende + ENTETE
    
    Save=np.array([X, Y]).transpose() # On met sous la forme XY en colonne.
    
    np.savetxt(Dossier+os.sep+NOM, Save, delimiter=';', header=HEADER, comments='');
      

def Nettoyage_spectre(Liste, Legende, Liste_ref, correction, Addition_Tr=0):
    '''
    Cette fonction permet de réaliser de traitement sur des spectres en nm/%T. 

    Parameters
    ----------
    Liste : TYPE
        Liste des fichier à traiter
    Legende : TYPE
    I100 : TYPE
        Référence EN ABS.
    correction : TYPE
        Entier qui permet de choisir la correction, 1 pour soustraire la référence,
        2 pour soustraitre UNIQUEMENT la bande UV (paramètre de base optimisé pour le fer),
        3 pour filtrer le bruit haute fréquence (filtre gaussien),
        4 pour effectué une correction du saut de détecteur du Perkin, avec filtrage et soustraction du blanc.
        5 pour passer une fichier d'ABS en Tr
        6 pour joindre les spectres issu du spectro portable
        7 pour sauvegarder les spectres en ABS et cm^-1
    Addition_Tr : float
        Valeur ajouté à la transmittance en cas de valeurs trop basse qui donne des nan suite au log
    Raises
    ------
    ValueError
        La correction demandée n'est pas implémenté.
        
    Returns
    -------
    Liste_corr : la liste des fichier corrigé

    '''
    Liste_corr=[]
    Dossier='./Data_corr'
    ENTETE = '\n X en nm; Y en %T'
    
    if np.size(Addition_Tr) == 1: Addition_Tr=(np.zeros(np.size(Liste)) + Addition_Tr)
    
    for i in np.arange(0, np.size(Liste), 1):
        #print(correction)
        
        if(correction[i]==0): # Si pas de correction on ne fait rien
            Liste_corr.append(Liste[i])
            
        else:
            Fichier=Liste[i];
            (cheminfichier, nomfichier) = os.path.split(Fichier)
              
            Xnm, Ytr = Readspectre(Fichier)
            X = 1/(Xnm*1E-7); #
            
            Y = - np.log10(Ytr+Addition_Tr[i])# Si data en %T
            
            Xsave=Xnm;
            
            if(correction[i] == 1): # uniquement soustraction de la référence
        
                Xnm_ref, Ytr_ref=Readspectre(Liste_ref[i])
                I100 = -np.log10(Ytr_ref)

                Y_corr=Y-I100
                Fichier_corr=nomfichier[0:-4] + Corr2Str(correction[i])
                
            elif(correction[i]==2): # fit bande UV                
                Y_corr=fit_OMCT(X,Y,Legende[i])    
                Fichier_corr=nomfichier[0:-4]+ Corr2Str(correction[i])
                
            elif(correction[i]==3): #Lissage
                Y_corr= filtrage_gauss(X, Y, Legende[i])
                Fichier_corr=nomfichier[0:-4]+ Corr2Str(correction[i])
                
            elif(correction[i]==4): # Saut de détecteur
                Xnm_ref, Ytr_ref=Readspectre(Liste_ref[i])
                I100 = -np.log10(Ytr_ref)

                Y_corr=Y-I100;
                #Y_corr= filtrage_gauss(X, Y_corr, Legende[i])
                Y_corr= correction_saut_detect(X, Y_corr, Legende[i])
                #Y_corr= filtrage_gauss(X, Y_corr, Legende[i] + 'corr', 2)
                Fichier_corr=nomfichier[0:-4] + Corr2Str(correction[i])
            
            elif(correction[i]==5): # passe en Tr, si Ttr en abosrbance il sera remis en %T lors de la sauv
                Y_corr  = Ytr
                Fichier_corr=nomfichier[0:-4] + Corr2Str(correction[i])
                Dossier = './Data_trait'
             
            elif(correction[i]==6): # Pour le raccord bien penser à mettre le spectre IR dans la Liste et le VIS dans la ref   
                Xnm, Y_corr=Traitement_spectro_portable(Liste_ref[i], Liste[i], Legende[i], Addition_Tr[i])
                Xsave = Xnm;
                
                Fichier_corr=nomfichier[0:-4]+ Corr2Str(correction[i])
            
            elif(correction[i]==7): #On sauverager en ABS / cm-1
                Y_save  = Y
                Xsave = X;
                Fichier_corr=nomfichier[0:-4] + Corr2Str(correction[i])
                Dossier = './Data_ABScm'
                ENTETE = '\n X en cm^-1; Y en Absorbance'
             
            elif(correction[i] == 8): # uniquement soustraction de la référence
                Xnm_ref, Ytr_ref=Readspectre(Liste_ref[i])
                REF = -np.log10(Ytr_ref)

                Y_corr=Y+REF
                Fichier_corr=nomfichier[0:-4] + Corr2Str(correction[i])
            
            else:
                raise ValueError('La correction demandée n\'est pas implémentée, vérifier le fichier d\'input')
            
            Liste_corr.append(Dossier+os.sep+Fichier_corr);
            
            if not (correction[i] == 7):
                Y_save=np.power(10, -Y_corr); # On repasse en %T pour assurer la comptabilité avec le reste des scripts

            SavCSV(Xsave, Y_save, Fichier_corr, Legende[i], Dossier, ENTETE);

    return(Liste_corr)
