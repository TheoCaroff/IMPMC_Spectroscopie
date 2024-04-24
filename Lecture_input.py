#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 14:16:08 2021

@author: theo
"""

from copy import copy
import numpy as np
import os
import warnings
import pandas as pd
import re
from tkinter.filedialog import askdirectory, askopenfilename
import tkinter as tk
import tempfile
import fnmatch
import ast
from matplotlib import pyplot as plt

try : import mpld3
except ModuleNotFoundError : pass

try : from colour import plotting as cplot
except ModuleNotFoundError : pass

from Calcul_masse_mol_concentration import Pmass2concentrationMol

#%% Utilitaire

def find_figsize_latexwidth(width_=320.52, xy_ratio=None):  
    """
        From a StackExchange answer  
        Set the size of the matplotlib figure in inches  
        according to the witdh of the LaTeX page  
        argument:  
        ---------  
        width_: the \textwidth of the LaTeX document, in pts  
                by default, 320.25 pts, probably the width of the standard  
                article document.find with  the  \the\textwidth latex command
        return:  
        (x_size, y_size): size of the matplotlib figure, in inches.  
            horizontally, it fills the width of the LaTeX page  
            vertically, the ratio y/x is set by xy_ratio, by default  
                    the golden ratio = 0.62  
    """

    inches_per_pt = 1 / 72.27
    
    if xy_ratio == None :
        golden_ratio = (5**.5-1)/2
        
        xy_ratio = golden_ratio
        
    return (width_*inches_per_pt, width_*inches_per_pt*xy_ratio)


def Sav_HTML(fig,Titre='Pouet', Repertoire='Graph_HTML') :
    '''
    Cette fonction permet de transformer un graph matplotlib en fichier html'

    Parameters
    ----------
    fig :  renvoyé par Affichage_spectre
    Titre : TYPE, optional
        Titre du graph. The default is 'Pouet'.
    Repertoire : TYPE, optional
        Repertoire dans lequel les figures vont être sauvergarder. The default is 'Graph_HTML'.

    Returns
    -------
    None.

    '''
    try:
        os.mkdir(Repertoire)
    except OSError:
        pass
    
    html_str = mpld3.fig_to_html(fig)
    titre_html = Repertoire+os.sep+str(Titre)+'.html'
    if os.path.exists (titre_html) : os. remove (titre_html)
    Html_file= open(titre_html,"w")
    Html_file.write(html_str)
    Html_file.close()


def Sav_fig(Titre='Pouet', Repertoire='Graph', colorplot=False, Tight=True,
            legendin=True, legend=True, SAVEHTML=False, legncol=1, save_format='.png',
            PDF=False, fig = None, ax=None, pad=0.001) :
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
    if not fig :
        fig = plt.gcf()
    if not ax :
        ax = plt.gca()
    
    if legend == True:
        if legendin == True : leg = ax.legend(ncol = legncol)
        else : leg = ax.legend(bbox_to_anchor=(1, 1), loc='upper left',
                                ncol=legncol)
        
        
        leg.set_draggable(state=True)
        
    
    if Tight : fig.tight_layout(pad=pad)     
    
    fig.savefig(Repertoire+os.sep+Titre+save_format) # bbox_inches='tight'
    if PDF and not (save_format == '.pdf'):
        fig.savefig(Repertoire+os.sep+Titre+'.pdf', bbox_inches='tight') 
    
    print(Repertoire+os.sep+Titre)
    

    
    if SAVEHTML : Sav_HTML(fig, Titre, Repertoire + '_HTML')
        
    if colorplot:
        cplot.render(standalone=True)
    else:
        plt.show()

def mono2tab(Value, size):
    '''
    Cette fonction sert à générer un tableau de la taille size, à partir d'un valeur

    Parameters
    ----------
    Value : TYPE
        DESCRIPTION.
    size : TYPE
        DESCRIPTION.

    Returns
    -------
        Liste de taille size avec Value dans chaque case.

    '''
    if (np.size(Value) == 1): # Si pas d'argument un marqueur unique pour toute les données
        Value_ref = Value
        Value = []
        for i in np.arange(0, size):
            Value.append(Value_ref)
    return(Value)


def list_folder_filter(CLEFDETRIE='', CHEMIN='.', DOSSIER='/Data_trait'):
    '''
    Renvoi une liste triée à l'humaine des fichiers contenant la clef
        de trie dans le repertoir DOSSIER

    Parameters
    ----------
    CLEFDETRIE : TYPE, optional
        DESCRIPTION. The default is ''.
    CHEMIN : TYPE, optional
        DESCRIPTION. The default is '.'.
    DOSSIER : TYPE, optional
        DESCRIPTION. The default is '/Data_trait'.

    Returns
    -------
    None.

    '''
    os.chdir(CHEMIN)
    CHEMIN = CHEMIN + os.sep + DOSSIER
    Liste = os.listdir(CHEMIN); #Récupère la liste des fichiers
    Liste = natural_sort(Liste)
    Liste = fnmatch.filter(Liste, CLEFDETRIE)
    return(Liste)

    
def Readspectre(Fichier, skip_header=2, delimiter=';', colonne=1):
    '''
    Cette fonction lis un CSV dont le chemin est renseigné dans la variable Fichier,
    il renvoit un valeur X correpsondant à la première colonne (0) et Y correspondant à la deuxième (1)

    Parameters
    ----------
    Fichier : string
        Chemin vers le fichier à lire.

    Returns
    -------
    None.

    '''
    
    try:
        Data = np.genfromtxt(Fichier, skip_header=skip_header, delimiter=delimiter); 
    except UnicodeDecodeError:
        Data = np.genfromtxt(Fichier, skip_header=skip_header, delimiter=delimiter, encoding='latin-1'); # si jamais l'encodage n'est pas UTF8
    X = Data[:, 0]
    Y = Data[:, colonne]/100

    return([X, Y])

def natural_sort(l):#Fonction qui sert à trier dans l'ordre naturel (ou humain) les listes.
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

def nm2cm1(X):
    return(1/(X*1E-7))

def eV2nm(X):
    return(1239.8/X)

def Htog(X, frequence):
    return((6.626e-34*frequence*1e9)/(9.274e-24*X*1e-4))
     

def normYminmax(X, Y, COUPURENORMminmax=[400, 2500], onlymin=False):
    '''
    Cette fonction normalise entre 0 et 1 le jeu de donnés Y dans la gamme des X compris dans la liste COUPURENORMminmax.
    

    Parameters
    ----------
    X : TYPE
        DESCRIPTION.
    Y : TYPE
        DESCRIPTION.
    COUPURENORMminmax : TYPE, optional
        DESCRIPTION. The default is [400, 2500].

    Returns
    -------
    None.

    '''
    # INDEX = np.logical_and(X > np.min(COUPURENORMminmax), X < np.max(COUPURENORMminmax))             
    # Y = Y - np.nanmin(Y[INDEX])
    # if not onlymin : Y = Y/np.nanmax(Y[INDEX])
    INDEX = np.logical_and(X>COUPURENORMminmax[0], X<COUPURENORMminmax[1])             
    Y = Y - np.min(Y[INDEX])
    if not onlymin : Y = Y/np.max(Y[INDEX])
    # print(Y)
    return(Y)

def integrate2curve(X, Y):
    integral_y = np.cumsum(Y[:-1] * np.diff(X))
    return (X[:-1], integral_y)
    

def Gauss(X, A1, pos1, sigma1, b): # Note H largeur à mis hauteurs = 2.3548*sigma
    '''
    calcul pour les valeur de X un gausienne + une constante
    
    Parameters
    ----------
    A1 : float
        Intensité de la gaussienne
    pos1 : float
        Centre de la gaussienne
    sigma1 : float
        sigma, FHW = 2.3548*sigma
    b : float
        constante pour shift en Y la gaussienne
        
    Returns
    -------
    None.
        
    '''
    RES=A1*np.exp(-(X-pos1)**2/(2*sigma1**2))+b;
    return(RES)

def removeInfNan(X,Y):
    '''
    Cette fonction retire les valeurs nan et inf en Y d'un jeux de données X/Y

    Parameters
    ----------
    X : TYPE
        DESCRIPTION.
    Y : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    Filtrenan = np.invert(np.isnan(Y))
    Filtreinf = np.invert(np.isinf(Y))
    
    
    FILTRE = np.logical_and(Filtrenan, Filtreinf)
        
    Xcorr=X[FILTRE]
    Ycorr=Y[FILTRE]
    
    return(Xcorr,Ycorr)

def Chemin2Liste(TITRE, Recuperer_nom_dossier_temporaire=False, CLEFDETRIE='', CHEMIN='.', LISTEDOSSIER=['Data_trait'], mode='OAS'):
    '''
    Cette fonction permet d'utiliser les codes écrit pour fonctionner avec un input en listant les fichier
    dans DOSSIER situé dans CHEMIN

    Parameters
    ----------
    TITRE : TYPE
        DESCRIPTION.
    Recuperer_nom_dossier_temporaire : TYPE, optional
        DESCRIPTION. The default is False.
    CLEFDETRIE : TYPE, optional
        DESCRIPTION. The default is ''.
    CHEMIN : TYPE, optional
        DESCRIPTION. The default is '.'.
    DOSSIER : TYPE, optional
        DESCRIPTION. The default is 'Data_trait'.

    Returns
    -------
    None.

    '''
    if CLEFDETRIE== '' : 
        CLEFDETRIE = '*'+TITRE+'*'
 
    if Recuperer_nom_dossier_temporaire:
        repspectro = open(tempfile.gettempdir()+os.sep+'repspectro', 'r');
        RACINE = repspectro.read();
        repspectro.close();

    else:
        RACINE = CHEMIN
        #RACINE = askdirectory()
    
    Liste=[]
    Legende=[]
    
    for DOSSIER in LISTEDOSSIER:
        Liste_temp = list_folder_filter(CLEFDETRIE, RACINE, DOSSIER);
        Legende=Legende+[x[:-4] for x in Liste_temp]
        Liste=Liste+[DOSSIER + os.sep + x for x in Liste_temp]
    
    if mode == 'OAS' :
        Liste_ref=''
        Liste_corr=''
        MarqueurCIE=''
        Addition_Tr=0
        valeurnorm=1
        Correction=0
        optplt=''
        Liste_refN=''
        Annot = ''
    
        return(Liste, Legende, Liste_ref, Liste_refN, Correction, optplt, MarqueurCIE, Addition_Tr, valeurnorm, Liste_corr, TITRE, Annot)

    if mode == 'RAMAN' :
        return(Liste, Legende)

#%% OAS

def Corr2Str_OAS(Correction_number):
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
    
    elif Correction_number == 10:
        Correction_NAME = '_cor_I100_I0.csv'
    
    elif Correction_number == 11:
        Correction_NAME = '_cor_saut_std.csv'
        
    elif Correction_number == 12:
        Correction_NAME = '_cor_saut_micro.csv'

    elif Correction_number == 13:
        Correction_NAME = '_cor_saut_lampe.csv'
        
    elif Correction_number == 14:
        Correction_NAME = '_cor_saut_filtre_319.csv'
            
    elif Correction_number == 113:
        Correction_NAME = '_cor_I00_detectstd_lampe.csv'
    
    elif Correction_number == 123:
        Correction_NAME = '_cor_I00_detect_lampe.csv'
            
    elif Correction_number == 1234:
        Correction_NAME = '_cor_detect_lampe_filtre319.csv'
    
    elif Correction_number == 2:
       Correction_NAME = '_soustrait.csv'
    
    elif Correction_number == 21:
       Correction_NAME = '_soustrait_norm.csv'
    
    elif Correction_number == 22:
          Correction_NAME = '_sous_UV.csv'
          
    elif Correction_number == 23:
          Correction_NAME = '_sous_baseline.csv'
          
    elif Correction_number == 231:
          Correction_NAME = '_baseline.csv'
          
    elif Correction_number == 2341:
        Correction_NAME = '_sous_baseline_norm_ep.csv'
        
    elif Correction_number == 24:
          Correction_NAME = '_sous_reflexion.csv'

    elif Correction_number == 25:
        Correction_NAME = '_sousmin_norm_ep.csv'

    elif Correction_number == 26 :
        Correction_NAME = '_sous_fit_gauss.csv'  

    elif Correction_number == 27:
        Correction_NAME = '_lin_baseline.csv'   
    elif Correction_number == 3:
       Correction_NAME = '_lisse.csv'
    
    elif Correction_number == 4:
        Correction_NAME = 'corr_NAN_INF.csv'      

    elif Correction_number == 41:
        Correction_NAME = '_norm_ep.csv'
    
        
    elif Correction_number == 5:
       #Correction_NAME = '_Tr.csv'
       Correction_NAME = '.csv'
    
    elif Correction_number == 6 or Correction_number == -6 :
        Correction_NAME = '_jointNIR.csv'
    
    elif Correction_number == 66 :
        Correction_NAME = '_jointNIR_jointUV.csv'
    
    elif Correction_number == 7:
        Correction_NAME = '_corr_DO.csv'
    
    elif Correction_number == 9:
        Correction_NAME = '_baseline.csv'
    
    elif Correction_number == 0:
        Correction_NAME=''
    

    elif Correction_number == 77:
        Correction_NAME='_remove_spline.csv'
 

    else :
        raise ValueError("Pb la correction n'existe pas")
    
    return(Correction_NAME)

def Corr2folder_OAS(Correction_number):
    '''
    Cette fonction sert à donner le repertoire de sauvegarde
    en fonction de la correction

    Parameters
    ----------
    Correction_number : int
        numéro de la correction.
    Returns
    -------
        Repertoire de la correction

    '''

    if Correction_number in [6, -6, 66]:
       Dossier = './Data_joint'
    
    elif Correction_number == 113 :
        Dossier = './Data_corr_I0_detect_lampe'

    elif Correction_number == 4 :
        Dossier = './Data_sans_inf_nan'

    elif Correction_number == 41 :
        Dossier = './Data_norm_ep'
    
    elif Correction_number == 25 :
        Dossier = './Data_sousmin_norm_ep'
        
    elif Correction_number == 231 :
        Dossier = './Baseline'
        
    elif Correction_number == 2341 :
        Dossier = './Data_convexhull_normep'
    
    elif Correction_number == 23 :
        Dossier = './Data_corr_convexhull'
    
    elif Correction_number == 5 :
        Dossier = './Data_traite'
    
    elif Correction_number == 9 :
        Dossier = './baseline'
    
    elif Correction_number == 77 :
        Dossier = './removespline'
    
    else :
        Dossier = './Data_corr'
    
    return(Dossier)

def Chemin2input_OAS(TITRE, CLEFDETRIE='', CHEMIN='.', LISTEDOSSIER=['Data_trait'], mode='portable', correction='0'):
    '''
    Cette fonction sert à créer un ficher d'input contenant le chemin relatif
    vers les fichiers matchant avec la clef de trie et étant dans les dossier 
    dans la liste LISTEDOSSIER
    
    Parameters
    ----------
    TITRE : string
        Titre.
    CLEFDETRIE : string, optional
        Chaine de caractère pour le trie, si vide le trie s'effectuer via le titre. The default is ''.
    CHEMIN : TYPE, optional
        DESCRIPTION. The default is '.'.
    LISTEDOSSIER : tab de string, optional
        Tableau des différents dossier à visiter. The default is ['Data_trait'].
    mode : string, optional
        choix du mode, PERKIN, portable ou ID20. The default is 'portable'.
    correction : TYPE, optional
        valeur de remplissage de la case correction. The default is '0'.
    
    Returns
    -------
    None.
    
    '''

    if CLEFDETRIE== '' : 
        CLEFDETRIE = '*'+TITRE+'*'
        
    Liste       =[]
    Liste_ref   =[]
    Legende     =[]
    
    for DOSSIER in LISTEDOSSIER:
       
        if mode == 'portable':
            Liste_temp       = list_folder_filter(CLEFDETRIE+'*VIS*', CHEMIN, DOSSIER);
            Liste_ref_temp   = list_folder_filter(CLEFDETRIE+'*NIR*', CHEMIN, DOSSIER);
            Legende_temp     = [re.sub("_Transmission.*csv", '', x) for x in Liste_temp]
            Legende_temp     = [re.sub("_T.*csv", '', x) for x in Legende_temp]
            Correction       = correction
        
        if mode == 'portable_trait':
            Liste_temp       = list_folder_filter(CLEFDETRIE+'*VIS*', CHEMIN, DOSSIER);
            Liste_ref_temp   = mono2tab('', np.size(Liste_temp))
            Legende_temp     = [re.sub("_T.*csv", '', x) for x in Liste_temp]
            Legende_temp     = [DOSSIER + ' ' + x for x in Legende_temp]
            Correction       = correction
        
        if mode == 'PERKIN':
            Liste_temp       = list_folder_filter(CLEFDETRIE, CHEMIN, DOSSIER);
            Liste_ref_temp   = mono2tab('', np.size(Liste_temp))
            Legende_temp     = [re.sub(".Sample.*csv", '', x) for x in Liste_temp]
            Correction       = correction

        if mode == 'Fluo':
            Liste_temp       = list_folder_filter(CLEFDETRIE, CHEMIN, DOSSIER);
            Liste_ref_temp   = mono2tab('', np.size(Liste_temp))
            Legende_temp     = [re.sub("_Graph.dat", '', x) for x in Liste_temp]
            Correction       = correction

        if mode == 'ID20' :
            Liste_temp       = list_folder_filter(CLEFDETRIE, CHEMIN, DOSSIER);
            Liste_ref_temp   = mono2tab('', np.size(Liste_temp))
            Legende_temp     = [re.sub(".dat", '', x) for x in Liste_temp]
            Correction       = correction
            
        
        Liste     = Liste     + ['.'+os.sep + DOSSIER + os.sep + x for x in Liste_temp]
        Liste_ref = Liste_ref + ['.'+os.sep + DOSSIER + os.sep + x for x in Liste_ref_temp]
        Legende   = Legende   + Legende_temp

    taille=np.size(Liste)
    df = pd.DataFrame(dict(
                      Nom_fichier = Liste,
                      Legende       = Legende,
                      Nom_ref       = Liste_ref,
                      Correction    = np.ones(taille)*Correction,
                      Matplot_opt   = mono2tab('', taille),
                      Marqueur_CIE  = mono2tab('', taille),
                      Epaisseur     = np.zeros(taille),
                      Concentration = np.ones(taille)*-1,
                      Tx_redox      = np.ones(taille),
                      Oxyde         = mono2tab('', taille),
                      Densite       = np.zeros(taille),
                      Addition_Tr   = np.zeros(taille),
                      Nom_refN      = mono2tab('', taille),
                      Annot         = mono2tab('', taille),
                      LimiteInf     = mono2tab('', taille),
                      LimiteSup     = mono2tab('', taille)
                      ))
    
    SAVE=df.to_csv(sep=';', index=False);
    
    print(SAVE)
    
    text_file = open('input_'+TITRE+'.csv', "w")

    text_file.write(SAVE)
    
    text_file.close()
    
    

def readinput_OAS(INPUTNAME='input.csv', mode='numpy', concentrationinput='none', comments = '#'):
    '''
    Cette fonction sert à lire les inputs de traitements et d'affichage des spectres'

    Parameters
    ----------
    INPUTNAME : TYPE, optional
        Nom du fichier d'input. The default is 'input.csv'.
    mode : TYPE, optional
        numpy ou pandas. The default is 'numpy'.
    concentrationinput : TYPE, optional
        none, epaisseur, molaire, massique. The default is 'none'. A adatper en fonction de l'input et de la normalisation souhaitée

    Returns
    -------
    None.

    '''
    
    TITRE = INPUTNAME[6:-4]

    if mode == 'numpy':
        try: # on ouvre en string pour récuperer le nom et le chemin des fichiers
            Data = np.genfromtxt(INPUTNAME, comments=comments, skip_header=1, delimiter=';', dtype='str'); 
        except UnicodeDecodeError:
            print("INPUT PAS UNICODE")
            Data = np.genfromtxt(INPUTNAME, comments=comments, skip_header=1, delimiter=';', dtype='str', encoding='latin-1');
    # Récupération des données de base du fichier d'input
        #print(Data)
        if np.size(Data.shape) < 2 : Data=np.array([Data]) # Si jamais le fichier ne contient qu'une seule ligne.
            
        Liste = (Data[:, 0].tolist())
        Legende = (Data[:, 1])
        Liste_ref = Data[:, 2] # Chemin de la référence, si vide le fichier à déja été traiter
        
        try:
            Liste_refN = Data[:,12]; # Chemin de la référence noir, si vide le noir n'est pas à soustraire
        except IndexError:
            Liste_refN = mono2tab('', np.size(Liste))
        
        try:
            Correction =  Data[:, 3].astype(float); 
        except ValueError:
            Correction = np.zeros(np.size(Liste))
        
        try:    optplt = Data[:, 4]
        except: optplt=''
        
        try :   MarqueurCIE = Data[:, 5]
        except: MarqueurCIE =''    
        
        
        try :       Addition_Tr = Data[:, 11].astype(float)
        except :    Addition_Tr = 0;
        
        try :       Epaisseur = Data[:, 6].astype(float)*1E-1 
        except :    Epaisseur = 0;
        
        try :       Concentration = Data[:, 7].astype(float)
        except :    Concentration=-1
        
        try :       Tx_redox = Data[:, 8].astype(float) # Le taux de rédox = [Fe2+]/[Fetot]
        except :    Tx_redox=1;
        
        try :       Oxyde = Data[:, 9] # Le taux de rédox = [Fe2+]/[Fetot]
        except :    Oxyde='';

        try :       Densite = Data[:, 10].astype(float) # Le taux de rédox = [Fe2+]/[Fetot]
        except :    Densite=0;
        
        try :       Annot = Data[:, 13] # Le taux de rédox = [Fe2+]/[Fetot]
        except :    Annot = '';
        
        try :       LimiteInf = Data[:, 14]
        except :    LimiteInf = -np.inf
        
        try :       LimiteSup = Data[:, 15]
        except :    LimiteSup = np.inf
        #Gestion des limites
        LimiteInf = mono2tab(LimiteInf, np.size(Liste)) # Création de la liste si valeur unique
        LimiteSup = mono2tab(LimiteSup, np.size(Liste))
        
    elif mode == 'pandas' :
        Data=pd.read_csv(INPUTNAME, delimiter=';', comment=comments, header=0)
            
        Liste =     Data.Nom_fichier.tolist()
        Legende =   Data.Legende.to_numpy()
        Liste_ref = Data.Nom_ref.tolist() # Chemin de la référence, si vide le fichier à déja été traiter
        
        try: Liste_refN = Data.Nom_refN.tolist(); 
        except : Liste_refN = mono2tab('', np.size(Liste))
        
        try:
            Correction = Data.Correction.to_numpy()  
        except AttributeError:
            Correction = np.zeros(np.size(Liste))
        
        try:        optplt = Data.Matplot_opt.tolist()
        except:     optplt=''
        
        try :       MarqueurCIE = Data.Marqueur_CIE.tolist()
        except:     MarqueurCIE ='';  
        
        
        try :       Addition_Tr = Data.Addition_Tr.tolist()
        except :    Addition_Tr = 0;

        try :       Epaisseur = Data.Epaisseur.to_numpy()*1E-1
        except :    Epaisseur = 0;
        
        try :       Concentration = Data.Concentration.to_list()
        except :    Concentration=-1
        
        try :       Tx_redox = Data.Tx_redox.to_list()# Le taux de rédox = [Fe2+]/[Fetot]
        except :    Tx_redox=1;
        
        try :       Oxyde = Data.Oxyde.to_list()
        except :    Oxyde = '';
        
        try :       Densite = Data.Densite.to_list()
        except :    Densite = 0
        
        try :       Annot = Data.Annot.to_list()
        except :    Annot = 0
        
        try :       LimiteInf = Data.LimiteSup.to_list()
        except :    LimiteInf = -np.inf
        
        try :       LimiteSup = Data.LimiteInf.to_list()
        except :    LimiteSup = np.inf
        
        
    print(concentrationinput)
    
    if concentrationinput == 'epaisseur':
        valeurnorm = Epaisseur

    elif concentrationinput == 'molaire' : # Récupération des information suplémentaire du fichier d'input
        valeurnorm=Concentration*Epaisseur*Tx_redox
    
    elif concentrationinput == 'massique' : # Récupération des information suplémentaire du fichier d'input  
        Concentration_molaire = Pmass2concentrationMol(Oxyde, Concentration, Densite)
        valeurnorm=Concentration_molaire*Epaisseur * Tx_redox
        
        print(valeurnorm)
        print('cm.g/mol')
        
        valeurnorm = valeurnorm.to_list()

    else:
        valeurnorm = 1

    
    Limite = np.zeros((np.size(Liste), 2)) # Créaction de la matrice limite
    LimiteInf = mono2tab(LimiteInf, np.size(Liste))
    LimiteSup = mono2tab(LimiteSup, np.size(Liste))
    
    
    for i, lim in enumerate(LimiteInf): # on parcous la liste pour mettre les trucs au bon endroit

        if LimiteInf[i] == '':
            Limite[i, 0] = -np.inf
        else :
            try : 
                Limite[i, 0]=LimiteInf[i]
            except :
                Limite[i, 0]=LimiteInf[i][0]
        
        if LimiteSup[i] == '' :
            Limite[i, 1]= np.inf
        else :
            try :
                Limite[i, 1]=LimiteSup[i]
            except :
                Limite[i, 1]=LimiteSup[i][0]
        
    
    try : # Récupération de la liste des fichier traité si applicable
        Liste_temp=copy(Liste);
        for i in np.arange(0, np.size(Liste), 1):
            if Correction[i] != 0:
                try:
                    (cheminfichier, nomfichier) = os.path.split(Liste[i])
                    DOSSIER = Corr2folder_OAS(Correction[i])
                    Liste_temp[i]=DOSSIER+os.sep+nomfichier[:-4]+Corr2Str_OAS(Correction[i])
                except :
                    warnings.warn('\nLe dossier '+DOSSIER+' est vide ou il y a trop de fichier par rapport au fichier d\'input\n')
            else:
                    Liste_temp[i]=Liste[i]            
    
        Liste_corr=np.array(Liste_temp);
    except FileNotFoundError:
        pass

    return(Liste, Legende, Liste_ref, Liste_refN, Correction, optplt,
           MarqueurCIE, Addition_Tr, valeurnorm, Liste_corr, TITRE,
           Annot, Limite)
#%% SAVE CIE
def saveDATACIE_xy(TITRE, Legende, x, y):
    DOSSIER='Colorimetrie'
    try :
        os.mkdir(DOSSIER)
    except :
        pass

    df = pd.DataFrame(dict(Nom_fichier  = Legende,
                           x            =  x,
                           y            = y))
    
    SAVE=df.to_csv(sep=';', index=False);
    
    print(SAVE)
    text_file = open('.' + os.sep + DOSSIER + os.sep + TITRE + '_CIE_xy.csv', "w")
    text_file.write(SAVE)
    text_file.close()
    
def saveDATACIE_LAB(TITRE, Legende, L, a, b):
    DOSSIER='Colorimetrie'
    try :
        os.mkdir(DOSSIER)
    except :
        pass
    
    df = pd.DataFrame(dict(Nom_fichier  = Legende,
                           L            =  L,
                           a            = a,
                           b            = b))
    
    SAVE=df.to_csv(sep=';', index=False);
    
    print(SAVE)
    text_file = open('.' + os.sep + DOSSIER + os.sep + TITRE + '_CIE_LAB.csv', "w")
    text_file.write(SAVE)
    text_file.close()

def saveDATACIE_XYZ(TITRE, Legende, X, Y, Z):
    DOSSIER='Colorimetrie'
    try :
        os.mkdir(DOSSIER)
    except :
        pass
    
    df = pd.DataFrame(dict(Nom_fichier  = Legende,
                           X            =  X,
                           Y            = Y,
                           Z            = Z))
    
    SAVE=df.to_csv(sep=';', index=False);
    
    print(SAVE)
    text_file = open('.' + os.sep + DOSSIER + os.sep + TITRE + '_CIE_XYZ.csv', "w")
    text_file.write(SAVE)
    text_file.close()
    
def saveDataFrame(df, TITRE='pouet', DOSSIER='.'):
    try :
        os.mkdir(DOSSIER)
    except :
        pass
    
    SAVE=df.to_csv(sep=';', index=True);
    
    text_file = open('.' + os.sep + DOSSIER + os.sep + TITRE + '.csv', "w")
    text_file.write(SAVE)
    text_file.close()

#%% Raman


def Corr2Str_RAMAN(Correction_number):
    '''
    Cette fonction sert Ã  donner l'extension d'un fichier en fonction de la correction demandÃ©

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
    Correction_number = int(Correction_number)
    # print(Correction_number)
    if Correction_number == 0:
        Correction_NAME=''
    
    elif Correction_number == 1:
       Correction_NAME = '_baseline_polyline_normalise_area.csv'
    
    elif Correction_number == 11:
          Correction_NAME = '_baseline_alteration_polyline_normalise_area.csv'

    elif Correction_number == 2 :
          Correction_NAME = '_baseline_gcvspline_normalise_area.csv'
          
    elif Correction_number == 3 :
          Correction_NAME = '_baseline_rubber_normalise_area.csv'
          
    elif Correction_number == 31 :
          Correction_NAME = '_baseline_rubber_alt_normalise_area.csv'
    #elif Correction_number == 10:
        #Correction_NAME = '_cor_I100_I0.csv'
    
    elif Correction_number == 6:
        Correction_NAME='_despiked.csv'
        
    elif Correction_number == 66:
        Correction_NAME = '_traitement_IMPMC.csv'
        
    elif Correction_number == 99 :
        Correction_NAME = "_averaged.csv"
        
    else :
        raise ValueError("Pb la correction n'existe pas")
    
    return Correction_NAME

def Corr2folder_RAMAN(Correction_number):
    '''
    Cette fonction sert à donner le repertoire de sauvegarde
    en fonction de la correction

    Parameters
    ----------
    Correction_number : int
        numéro de la correction.
    Returns
    -------
        Repertoire de la correction

    '''

    
    if Correction_number in [1, 2, 3]:
        Dossier = './Data_corr_baseline'
    else :
        Dossier = './Data_corr'
       
        
    return(Dossier)

def create_raman_txt(REP='.') :
    '''
    Créer un ficher txt du même nom que le fichier tsf (raman portable) dans le repertoir REP

    Parameters
    ----------
    REP : TYPE, optional
        DESCRIPTION. The default is '.'.

    Returns
    -------
    None.

    '''
    for file in list_folder_filter('*.tsf', '.', REP):
        os.system('touch ' + file[:-4] + '.txt')
            
    

def Chemin2input_RAMAN(TITRE, CLEFDETRIE='', CHEMIN='.', LISTEDOSSIER=['Data_trait'], mode='portable', correction='0'):
    '''
    Cette fonction sert à créer un ficher d'input contenant le chemin relatif
    vers les fichiers matchant avec la clef de trie et étant dans les dossier 
    dans la liste LISTEDOSSIER
    
    Parameters
    ----------
    TITRE : string
        Titre.
    CLEFDETRIE : string, optional
        Chaine de caractère pour le trie, si vide le trie s'effectuer via le titre. The default is ''.
    CHEMIN : TYPE, optional
        DESCRIPTION. The default is '.'.
    LISTEDOSSIER : tab de string, optional
        Tableau des différents dossier à visiter. The default is ['Data_trait'].
    mode : string, optional
        choix du mode, PERKIN, portable ou ID20. The default is 'portable'.
    correction : TYPE, optional
        valeur de remplissage de la case correction. The default is '0'.
    
    Returns
    -------
    None.
    
    '''

    if CLEFDETRIE== '' : 
        CLEFDETRIE = '*'+TITRE+'*'
        
    Liste       =[]
    Liste_ref   =[]
    Legende     =[]
    
    for DOSSIER in LISTEDOSSIER:
       
        if mode == 'RAMAN_portable':
            Liste_temp       = list_folder_filter(CLEFDETRIE+'*.txt', CHEMIN, DOSSIER);
            Legende_temp     = [re.sub("_300mW*.txt", '', x) for x in Liste_temp]
            #Legende_temp     = [re.sub("_T.*csv", '', x) for x in Legende_temp]
            Correction       = correction
        elif mode == 'RAMAN_IPG':
            Liste_temp       = list_folder_filter(CLEFDETRIE+'*.txt', CHEMIN, DOSSIER);
            Legende_temp     = [re.sub(".txt", '', x) for x in Liste_temp]
            #Legende_temp     = [re.sub("_T.*csv", '', x) for x in Legende_temp]
            Correction       = correction
        elif mode == 'RAMAN_CSV':
            Liste_temp       = list_folder_filter(CLEFDETRIE+'*.csv', CHEMIN, DOSSIER);
            #Legende_temp     = [re.sub("_despiked_baseline_normalise_area.csv", '', x) for x in Liste_temp]
            Legende_temp     = [re.sub(".csv", '', x) for x in Liste_temp]
            #Legende_temp     = [re.sub("_T.*csv", '', x) for x in Legende_temp]
            Correction       = correction

        elif mode == 'RAMAN_IMPMC':
            Liste_temp       = list_folder_filter(CLEFDETRIE, CHEMIN, DOSSIER);
            Legende_temp     = [re.sub(".dat", '', x) for x in Liste_temp]
            Legende_temp     = [re.sub("_r554_80p_vertical", '', x) for x in Legende_temp]
            Correction       = correction
        
        elif mode == 'RIXS' :
            Liste_temp       = list_folder_filter(CLEFDETRIE, CHEMIN, DOSSIER);
            Legende_temp     = [re.sub("RIXS_*.pkl", '', x) for x in Liste_temp]
            Correction       = correction
        
        Liste     = Liste     + ['.'+os.sep + DOSSIER + os.sep + x for x in Liste_temp]
        Liste_ref = Liste_ref + ['' for x in Liste_temp]
        Legende   = Legende   + Legende_temp

    taille=np.size(Liste)
    df = pd.DataFrame(dict(
                      Nom_fichier = Liste,
                      Legende       = Legende,
                      Nom_ref       = Liste_ref,
                      Correction    = np.ones(taille)*Correction,
                      Matplot_opt   = mono2tab('', taille),
                      Nom_refN      = mono2tab('', taille)
                      ))
    
    SAVE=df.to_csv(sep=';', index=False);
    
    print(SAVE)
    
    text_file = open('input_'+TITRE+'.csv', "w")

    text_file.write(SAVE)
    
    text_file.close()
    


def readinput_RAMAN(INPUTNAME='input.csv', comment = '#'):
    '''
    Cette fonction sert à lire les inputs de traitements et d'affichage des spectres'

    Parameters
    ----------
    INPUTNAME : TYPE, optional
        Nom du fichier d'input. The default is 'input.csv'.
    mode : TYPE, optional
        numpy ou pandas. The default is 'numpy'.
    concentrationinput : TYPE, optional
       

    Returns
    -------
    None.

    '''
    
    TITRE = INPUTNAME[6:-4]


    Data=pd.read_csv(INPUTNAME, delimiter=';', comment=comment, header=0)
        
    Liste =     Data.Nom_fichier.tolist()
    Legende =   Data.Legende.tolist()
    Liste_ref = Data.Nom_ref.tolist() # Chemin de la référence, si vide le fichier à déja été traiter
    
    
    try:
        Correction = Data.Correction.tolist()  
    except AttributeError:
        Correction = np.zeros(np.size(Liste))
    
    try:        optplt = Data.Matplot_opt.to_numpy(dtype='str')
    except:     optplt=''
    
    try: Liste_refN = Data.Nom_refN.to_list()
    except : Liste_refN = ''
    
    try:        Marker = Data.Marker.to_numpy(dtype='str')
    except:     Marker=''
    

    try : # Récupération de la liste des fichier traité si applicable
        Liste_temp=copy(Liste);
        for i in np.arange(0, np.size(Liste), 1):
            if Correction[i] != 0:
                
                try:
                    (cheminfichier, nomfichier) = os.path.split(Liste[i])
                    DOSSIER = Corr2folder_RAMAN(Correction[i])
                    Liste_temp[i]=DOSSIER+os.sep+nomfichier[:-4]+Corr2Str_RAMAN(Correction[i])
                except IndexError:
                    warnings.warn('\nLe dossier '+DOSSIER+' est vide ou il y a trop de fichier par rapport au fichier d\'input\n')
                except TypeError:
                    warnings.warn('\nL\'INPUT contient des lignes vide mal enterpréter, nettoyer tous les lignes qui semble vide\n')
            else:
                    Liste_temp[i]=Liste[i]            
    
        Liste_corr=np.array(Liste_temp);
    except FileNotFoundError:
        pass

    return(Liste, Legende, Liste_ref, Correction, optplt, Liste_corr, TITRE, Liste_refN, Marker)

#%% RPE

def Corr2Str_RPE(Correction_number):
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
   
    if Correction_number == 1:
       Correction_NAME = '_sous_cavite.asc'
    
    elif Correction_number == 2:
        Correction_NAME = '_normalise.asc'
    
    elif Correction_number == 3:
        Correction_NAME = '_sub_baseline.asc'
        
    elif Correction_number == 31:
        Correction_NAME = '_sub_baseline_IntPoly6.asc'
        
    elif Correction_number == 32:
        Correction_NAME = '_sub_baseline_cubicspline.asc'
    
    elif Correction_number == 4:
        Correction_NAME='_lisse.asc'
    
    elif Correction_number == 5:
        Correction_NAME='_calibre.asc'
        
    else :
        raise ValueError("Pb la correction n'existe pas")
    
    return(Correction_NAME)

def Releve_parametre(CLEFDETRIE='',CHEMIN='.', DOSSIER = 'Data_parametre'):
    '''
    Cette fonction sert à aller chercher les paramètre nécessaires dans le fichier de parametres
    '''
    Liste_parametre=[]
    Liste_Gain = []
    Liste_Frequence=[]
    Liste_Puissance = []
    Liste_Amplitude_Mod=[]
    Liste_Nb_scan=[]
    
    Liste_parametre_temp = list_folder_filter(CLEFDETRIE+'.par', CHEMIN, DOSSIER);
    Liste_parametre = Liste_parametre  + ['.'+os.sep + DOSSIER + os.sep + x for x in Liste_parametre_temp]
  
    for i in Liste_parametre :
        file_parametre = open(i)
        content = file_parametre.readlines()
        for line in content :
            if line.startswith('RRG') :
                Gain = re.sub('\n', '', re.sub("RRG ", '', line))
            if line.startswith('MF') :
                Frequence =  re.sub('\n', '', re.sub("MF  ", '', line))
            if line.startswith('MP  ') :
                Puissance =  re.sub('\n', '', re.sub("MP  ", '', line))
            if line.startswith('RMA ') :
                Amplitude_Mod = re.sub('\n', '', re.sub("RMA ", '', line))
            if line.startswith('JSD ') :
                Nb_scan =  re.sub('\n', '', re.sub("JSD ", '', line))
        Liste_Gain.append(Gain)
        Liste_Frequence.append(Frequence)
        Liste_Puissance.append(Puissance)
        Liste_Amplitude_Mod.append(Amplitude_Mod)
        Liste_Nb_scan.append(Nb_scan)
        
    return(Liste_Gain,Liste_Frequence,Liste_Puissance,Liste_Amplitude_Mod,Liste_Nb_scan )
          


    
def Chemin2input_RPE(TITRE, CLEFDETRIE='', CHEMIN='.',LISTEDOSSIER=['Data_trait'],mode='RPE',correction=1):
    '''
    Cette fonction sert à créer un ficher d'input contenant le chemin relatif
    vers les fichiers matchant avec la clef de trie et étant dans les dossier 
    dans la liste LISTEDOSSIER
    
    Parameters
    ----------
    TITRE : string
        Titre.
    CLEFDETRIE : string, optional
        Chaine de caractère pour le trie, si vide le trie s'effectuer via le titre. The default is ''.
    CHEMIN : TYPE, optional
        DESCRIPTION. The default is '.'.
    LISTEDOSSIER : tab de string, optional
        Tableau des différents dossier à visiter. The default is ['Data_trait'].
    mode : string, optional
        choix du mode, PERKIN, portable, RPE, RAMAN ou ID20. The default is 'portable'.
    correction : TYPE, optional
        valeur de remplissage de la case correction. The default is '0'.
    
    Returns
    -------
    None.
    
    '''

    if CLEFDETRIE== '' : 
        CLEFDETRIE = '*'+TITRE+'*'
        
    Liste       =[]
    Liste_ref   =[]
    Legende     =[]
    
    Liste_Gain  =[]
    Liste_Frequence = []
    Liste_Puissance = []
    Liste_Amplitude_Mod = []
    Liste_Nb_scan = []
    for DOSSIER in LISTEDOSSIER:
       
        if mode == 'RPE':
            Liste_temp       = list_folder_filter(CLEFDETRIE+'.asc', CHEMIN, DOSSIER);
            Legende_temp     = [re.sub(".asc", '', x) for x in Liste_temp]
            Correction       = correction
            Liste_Gain_temp, Liste_Frequence_temp, Liste_Puissance_temp, Liste_Amplitude_Mod_temp, Liste_Nb_scan_temp = Releve_parametre(CLEFDETRIE ,CHEMIN, DOSSIER)
        
            Liste     = Liste     + ['.'+os.sep + DOSSIER + os.sep + x for x in Liste_temp]
            Liste_ref = Liste_ref + ['' for x in Liste_temp]
            Legende   = Legende   + Legende_temp
            
            Liste_Gain  = Liste_Gain + Liste_Gain_temp
            Liste_Frequence = Liste_Frequence + Liste_Frequence_temp
            Liste_Puissance = Liste_Puissance + Liste_Puissance_temp
            Liste_Amplitude_Mod = Liste_Amplitude_Mod + Liste_Amplitude_Mod_temp
            Liste_Nb_scan = Liste_Nb_scan + Liste_Nb_scan_temp

        elif mode == 'fityk' :
            Liste_temp       = list_folder_filter(CLEFDETRIE+'.dat', CHEMIN, DOSSIER);
            Legende_temp     = [re.sub(".dat", '', x) for x in Liste_temp]
            Correction       = correction
        
            Liste     = Liste     + ['.'+os.sep + DOSSIER + os.sep + x for x in Liste_temp]
            Liste_ref = Liste_ref + ['' for x in Liste_temp]
            Legende   = Legende   + Legende_temp
            
            Liste_Gain  = 0
            Liste_Frequence = 0
            Liste_Puissance = 0
            Liste_Amplitude_Mod = 0
            Liste_Nb_scan = 0


    taille=np.size(Liste)
    df = pd.DataFrame(dict(
                      Nom_fichier   = Liste,
                      Legende       = Legende,
                      Nom_ref       = Liste_ref,
                      Correction    = np.ones(taille)*Correction,
                      Matplot_opt   = mono2tab('', taille),
                      Gain          = Liste_Gain,
                      Frequence     = Liste_Frequence,
                      Puissance     = Liste_Puissance,
                      Amplitude_Mod = Liste_Amplitude_Mod,
                      Nb_scan       = Liste_Nb_scan,
                      Hauteur       = np.zeros(taille),
                      Masse         = np.zeros(taille) ))

    SAVE=df.to_csv(sep=';', index=False);
    
    print(SAVE)
    
    text_file = open('input_'+TITRE+'.csv', "w")

    text_file.write(SAVE)
    
    text_file.close()
    
def readinput_RPE(INPUTNAME='input.csv'):
    '''
    Cette fonction sert à lire les inputs de traitements et d'affichage des spectres'

    Parameters
    ----------
    INPUTNAME : TYPE, optional
        Nom du fichier d'input. The default is 'input.csv'.
    mode : TYPE, optional
        numpy ou pandas. The default is 'numpy'.
    concentrationinput : TYPE, optional
       

    Returns
    -------
    None.

    '''
    
    TITRE = INPUTNAME[6:-4]
    
    Data=pd.read_csv(INPUTNAME, delimiter=';', comment='#', header=0)
    Liste =     Data.Nom_fichier.tolist()
    Legende =   Data.Legende.tolist()
    Liste_ref = Data.Nom_ref.tolist() # Chemin de la référence, si vide le fichier à déja été traiter
    
    
    try:
        Correction = Data.Correction.tolist()  
    except AttributeError:
        Correction = np.zeros(np.size(Liste))
    
    try:        optplt = Data.Matplot_opt.to_numpy(dtype='str')
    except:     optplt=''
    
    try :       Gain = Data.Gain.to_numpy()
    except :    Gain = 0;
    
    try :       Frequence = Data.Frequence.to_numpy()
    except :    Frequence = 0;
    
    try :       Puissance = Data.Puissance.to_numpy()
    except :    Puissance = 0;
    
    try :       Amplitude_Mod = Data.Amplitude_Mod.to_numpy()
    except :    Amplitude_Mod = 0;
    
    try :       Nb_scan = Data.Nb_scan.to_numpy()
    except :    Nb_scan = 0;
    
    try :       Hauteur = Data.Hauteur.to_numpy()
    except :    Hauteur = 0;
    
    try :       Masse = Data.Masse.to_numpy()
    except :    Masse = 0;
    
    try :       LimiteInf = Data.LimiteSup.to_list()
    except :    LimiteInf = -np.inf
    
    try :       LimiteSup = Data.LimiteInf.to_list()
    except :    LimiteSup = np.inf
        
    
    
    Limite = np.zeros((np.size(Liste), 2)) # Créaction de la matrice limite
    LimiteInf = mono2tab(LimiteInf, np.size(Liste))
    LimiteSup = mono2tab(LimiteSup, np.size(Liste))
    
    
    for i, lim in enumerate(LimiteInf): # on parcous la liste pour mettre les trucs au bon endroit

        if LimiteInf[i] == '':
            Limite[i, 0] = -np.inf
        else :
            try : 
                Limite[i, 0]=LimiteInf[i]
            except :
                Limite[i, 0]=LimiteInf[i][0]
        
        if LimiteSup[i] == '' :
            Limite[i, 1]= np.inf
        else :
            try :
                Limite[i, 1]=LimiteSup[i]
            except :
                Limite[i, 1]=LimiteSup[i][0]
    
    
    
    gain = Gain
    m = Masse
    h = Hauteur
    p = Puissance*0.0001
    am = Amplitude_Mod
    ns = Nb_scan
    Facteur_cavite = 1.5929 + 0.019438*h - 0.0094083*h*h+0.00067139-2.1875e-05*h**4+1.1815e-07*h**5
    valeurnorm = Facteur_cavite*gain*m*am*ns*np.sqrt(p)

    try : 
        gain = Gain
        m = Masse
        h = Hauteur
        p = Puissance*0.0001
        am = Amplitude_Mod
        ns = Nb_scan
        Facteur_cavite = 1.5929 + 0.019438*h - 0.0094083*h*h+0.00067139-2.1875e-05*h**4+1.1815e-07*h**5
        valeurnorm = Facteur_cavite*gain*m*am*ns*np.sqrt(p)
    except :
        valeurnorm=0
        
    try : # Récupération de la liste des fichier traité si applicable
        Liste_temp=copy(Liste);
        for i in np.arange(0, np.size(Liste), 1):
            if Correction[i] != 0:
                try:
                    (cheminfichier, nomfichier) = os.path.split(Liste[i])
                    DOSSIER = 'Data_corr'
                    Liste_temp[i]=DOSSIER+os.sep+nomfichier[:-4]+Corr2Str_RPE(Correction[i])
                except IndexError:
                    warnings.warn('\nLe dossier '+DOSSIER+' est vide ou il y a trop de fichier par rapport au fichier d\'input\n')
            else:
                    Liste_temp[i]=Liste[i]            
    
        Liste_corr=np.array(Liste_temp);
    except FileNotFoundError:
        pass
    
    
    
    return(Liste, Legende, Liste_ref, Correction, optplt, Frequence, valeurnorm, Liste_corr, Limite, TITRE)
    

#%% RIXS


def Corr2Str_RIXS(Correction_number):
    '''
    Cette fonction sert Ã  donner l'extension d'un fichier en fonction de la correction demandÃ©

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
    Correction_number = int(Correction_number)
    # print(Correction_number)
    if Correction_number == 0:
        Correction_NAME=''
    
    elif Correction_number == 1:
       Correction_NAME = '_baseline_polyline_normalise_area.csv'
    
    elif Correction_number == 11:
          Correction_NAME = '_baseline_alteration_polyline_normalise_area.csv'

    elif Correction_number == 2 :
          Correction_NAME = '_baseline_gcvspline_normalise_area.csv'
    
    #elif Correction_number == 10:
        #Correction_NAME = '_cor_I100_I0.csv'
    
    elif Correction_number == 6:
        Correction_NAME='_despiked.csv'
        
    elif Correction_number == 66:
        Correction_NAME = '_traitement_IMPMC.csv'
        
    elif Correction_number == 99 :
        Correction_NAME = "_averaged.csv"
        
    else :
        raise ValueError("Pb la correction n'existe pas")
    
    return Correction_NAME

def Corr2folder_RIXS(Correction_number):
    '''
    Cette fonction sert à donner le repertoire de sauvegarde
    en fonction de la correction

    Parameters
    ----------
    Correction_number : int
        numéro de la correction.
    Returns
    -------
        Repertoire de la correction

    '''

    
    if Correction_number == 1 :
        Dossier = './Data_corr_baseline'
    else :
        Dossier = './Data_corr'
       
        
    return(Dossier)
    

def Chemin2input_RIXS(TITRE, CLEFDETRIE='', CHEMIN='.', LISTEDOSSIER=['Data_trait'], mode='portable', correction='0'):
    '''
    Cette fonction sert à créer un ficher d'input contenant le chemin relatif
    vers les fichiers matchant avec la clef de trie et étant dans les dossier 
    dans la liste LISTEDOSSIER
    
    Parameters
    ----------
    TITRE : string
        Titre.
    CLEFDETRIE : string, optional
        Chaine de caractère pour le trie, si vide le trie s'effectuer via le titre. The default is ''.
    CHEMIN : TYPE, optional
        DESCRIPTION. The default is '.'.
    LISTEDOSSIER : tab de string, optional
        Tableau des différents dossier à visiter. The default is ['Data_trait'].
    mode : string, optional
        choix du mode, PERKIN, portable ou ID20. The default is 'portable'.
    correction : TYPE, optional
        valeur de remplissage de la case correction. The default is '0'.
    
    Returns
    -------
    None.
    
    '''

    if CLEFDETRIE== '' : 
        CLEFDETRIE = '*'+TITRE+'*'
        
    Liste       =[]
    Liste_ref   =[]
    Legende     =[]
    
    for DOSSIER in LISTEDOSSIER:
       

        
        if mode == 'RIXS_alpha' :
            Liste_temp       = [x for x in list_folder_filter(CLEFDETRIE, CHEMIN, DOSSIER) if not( 'long' in x)];
            Liste_ref_temp   = list_folder_filter(CLEFDETRIE+'*long*', CHEMIN, DOSSIER);
            Legende_temp     = [re.sub("_RIXS_alpha.pkl", ' Kalpha', x) for x in Liste_temp]
            Correction       = correction
            Seuil            = 'alpha'
        
        elif mode == 'RIXS_beta' :
            Liste_temp       = [x for x in list_folder_filter(CLEFDETRIE+'*RIXS*', CHEMIN, DOSSIER) if not( 'long' in x)];
            Liste_ref_temp   = list_folder_filter(CLEFDETRIE+'*Kb*', CHEMIN, DOSSIER);
            Legende_temp     = [re.sub("_RIXS_beta.pkl", ' Kbeta', x) for x in Liste_temp]
            Correction       = correction
            Seuil            = 'beta'
        
        
        else :
            Liste_temp       = list_folder_filter(CLEFDETRIE, CHEMIN, DOSSIER);
            Legende_temp     = [re.sub("RIXS_*.pkl", '', x) for x in Liste_temp]
            Liste_ref_temp   = mono2tab('', np.size(Liste_temp))
            Correction       = correction
            Seuil            = ''
            
        Liste     = Liste     + ['.'+os.sep + DOSSIER + os.sep + x for x in Liste_temp]
        Liste_ref = Liste_ref + ['.'+os.sep + DOSSIER + os.sep + x for x in Liste_ref_temp]
        Legende   = Legende   + Legende_temp

    if np.size(Liste_ref) == 0:
        Liste_ref = [None]
        
    taille=np.size(Liste)
    df = dict(        Nom_fichier = Liste,
                      Legende       = Legende,
                      Nom_ref       = Liste_ref,
                      Seuil         = mono2tab(Seuil, taille),
                      Correction    = np.ones(taille)*Correction,
                      Matplot_opt   = mono2tab('', taille),
                      PosCEE   = mono2tab('', taille),
                      PosCIE   = mono2tab('', taille),
                      PosCET   = mono2tab('', taille),
                      PosITE   = mono2tab('', taille),
                      PosIIE   = mono2tab('', taille),
                      Clim     = mono2tab('', taille))
                      
    df = pd.DataFrame(df)
    
    SAVE=df.to_csv(sep=';', index=False);
    
    print(SAVE)
    
    text_file = open('input_'+TITRE+'.csv', "w")

    text_file.write(SAVE)
    
    text_file.close()
    

def readinput_RIXS(INPUTNAME='input.csv'):
    '''
    Cette fonction sert à lire les inputs de traitements et d'affichage des spectres'

    Parameters
    ----------
    INPUTNAME : TYPE, optional
        Nom du fichier d'input. The default is 'input.csv'.
    mode : TYPE, optional
        numpy ou pandas. The default is 'numpy'.
    concentrationinput : TYPE, optional
       

    Returns
    -------
    None.

    '''
    
    TITRE = INPUTNAME[6:-4]
    
    POSEcoupe={}
    
    Data=pd.read_csv(INPUTNAME, delimiter=';', comment='#', header=0)
        
    Liste =     Data.Nom_fichier.tolist()
    Legende =   Data.Legende.tolist()
    Liste_ref = Data.Nom_ref.tolist() # Chemin de la référence, si vide le fichier à déja été traiter
    
    Tailletab = np.size(Liste)
    
    try:
        Correction = Data.Correction.tolist()  
    except AttributeError:
        Correction = np.zeros(np.size(Liste))
    
    try:        optplt = Data.Matplot_opt.to_numpy(dtype='str')
    except:     optplt=''
    
    try: Seuil = Data.Seuil.to_list()
    except : Seuil = ''
    


    try : PosCEE = Data.PosCEE.to_list()
    except : PosCEE = mono2tab(None, Tailletab)
    POSEcoupe['CEE'] = PosCEE

    try : PosCIE = Data.PosCIE.to_list()
    except : PosCIE = mono2tab(None, Tailletab)
    POSEcoupe['CIE'] = PosCIE
    
    try : PosCET = Data.PosCET.to_list()
    except : PosCET = mono2tab(None, Tailletab)
    POSEcoupe['CET'] = PosCET
    
    try :
        PosITE = Data.PosITE.to_list()
        PosITE = [ast.literal_eval(i) for i in PosITE]
    except : PosITE = mono2tab(None, Tailletab)
    POSEcoupe['ITE'] = PosITE
   
    try :
        PosIIE = Data.PosIIE.to_list()
        PosIIE = [ast.literal_eval(i) for i in PosIIE]
    except : PosIIE = mono2tab(None, Tailletab)
    POSEcoupe['IIE'] = PosIIE
    
    try :
        clim = Data.Clim.to_list()
        clim = [ast.literal_eval(i) for i in clim]
    except : clim = mono2tab(None, Tailletab)
    POSEcoupe['clim'] = clim

    try : # Récupération de la liste des fichier traité si applicable
        Liste_temp=copy(Liste);
        for i in np.arange(0, np.size(Liste), 1):
            if Correction[i] != 0:
                
                try:
                    (cheminfichier, nomfichier) = os.path.split(Liste[i])
                    DOSSIER = Corr2folder_RAMAN(Correction[i])
                    Liste_temp[i]=DOSSIER+os.sep+nomfichier[:-4]+Corr2Str_RAMAN(Correction[i])
                except IndexError:
                    warnings.warn('\nLe dossier '+DOSSIER+' est vide ou il y a trop de fichier par rapport au fichier d\'input\n')
            else:
                    Liste_temp[i]=Liste[i]            
    
        Liste_corr=np.array(Liste_temp);
    except FileNotFoundError:
        pass

    return(Liste, Legende, Liste_ref, Seuil, Correction, optplt, Liste_corr, TITRE, POSEcoupe)


#%% Boite de dialogue 

def askinput():
    '''
    Fonction servant à récupérer le nom du fichier d'input

    Returns
    -------
    None.

    '''
    tk.Tk().withdraw() # Pour créer un fenètre Tkinter fantôme
    filepath = askopenfilename(initialdir='.', filetypes= (("input","*.csv"), ("all files","*.*")))
    path, file = os.path.split(filepath)
    os.chdir(path)
    TITRE = file[6:-4] # pour supprimer le input_*.csv
    # tkdialog.destroy()
    print(TITRE)
    return(TITRE)
