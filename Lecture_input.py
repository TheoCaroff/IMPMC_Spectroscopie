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
from tkinter.filedialog import askdirectory
import tempfile
import fnmatch

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
       Correction_NAME = '_cor_I100_saut_filtre.csv'
       #Correction_NAME = '_cor_I100_saut.csv'
   
    elif Correction_number == 5:
       #Correction_NAME = '_Tr.csv'
       Correction_NAME = '.csv'
    
    elif Correction_number == 6:
        Correction_NAME = '_jointNIR.csv'
    
    elif Correction_number == 7:
        Correction_NAME = '_ABScm.csv'
    
    elif Correction_number == 8:
        #Correction_NAME = '_cor_DO.csv'
        Correction_NAME = 'corr.csv'
    elif Correction_number == 9:
        Correction_NAME = '_soustrait_norm.csv'

    elif Correction_number == 0:
        Correction_NAME=''
        
    else :
        raise ValueError("Pb la correction n'existe pas")
    
    return(Correction_NAME)

def list_folder_filter(CLEFDETRIE='', CHEMIN='.', DOSSIER='/Data_trait'):
    os.chdir(CHEMIN)
    CHEMIN = CHEMIN + os.sep + DOSSIER
    Liste = os.listdir(CHEMIN); #Récupère la liste des fichiers
    Liste = natural_sort(Liste)
    Liste = fnmatch.filter(Liste, CLEFDETRIE)
    return(Liste)

def Chemin2input(TITRE, CLEFDETRIE='', CHEMIN='.', DOSSIER='Data_trait', mode='portable', correction='0'):
    if CLEFDETRIE== '' : 
        CLEFDETRIE = '*'+TITRE+'*'
    
    if mode == 'portable':
        Liste       = list_folder_filter(CLEFDETRIE+'*VIS_*', CHEMIN, DOSSIER);
        Liste_ref   = list_folder_filter(CLEFDETRIE+'*NIR*', CHEMIN, DOSSIER);
        Legende     = [re.sub("_Transmission.*csv", '', x) for x in Liste]
        Correction  = 6;
        
    if mode == 'PERKIN':
        Liste       = list_folder_filter(CLEFDETRIE, CHEMIN, DOSSIER);
        Liste_ref   = mono2tab('', np.size(Liste))
        Legende     = [re.sub("Sample.*csv", '', x) for x in Liste]
        Correction  = 1
        
    Liste=['.'+os.sep + DOSSIER + os.sep + x for x in Liste]
    Liste_ref=['.'+os.sep + DOSSIER + os.sep + x for x in Liste_ref]
    
    taille=np.size(Liste)

    
    df = pd.DataFrame(dict(Nom_fichier=Liste,
                      Legende       = Legende,
                      Nom_ref       = Liste_ref,
                      Correction    = np.ones(taille)*Correction,
                      Matplot_opt   = mono2tab('', taille),
                      Marqueur_CIE  = mono2tab('', taille),
                      Epaisseur     = np.zeros(taille),
                      Concentration = np.ones(taille)*-1,
                      Tx_redox      = np.ones(taille),
                      Addition_Tr   = np.zeros(taille)))
    
    SAVE=df.to_csv(sep=';', index=False);
    
    print(SAVE)
    
    text_file = open('input_'+TITRE+'.csv', "w")

    text_file.write(SAVE)
    
    text_file.close()
    
def Chemin2Liste(TITRE, Recuperer_nom_dossier_temporaire=False, CLEFDETRIE='', CHEMIN='.', DOSSIER='Data_trait'):
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
    
    Liste = list_folder_filter(CLEFDETRIE, RACINE, DOSSIER);
    
    Legende=[x[:-4] for x in Liste]
    
    Liste=[DOSSIER + os.sep + x for x in Liste]
    Liste_ref=''
    Liste_corr=''
    MarqueurCIE=''
    Addition_Tr=0
    valeurnorm=1
    Correction=0
    optplt=''
    
    return(Liste, Legende, Liste_ref, Correction, optplt, MarqueurCIE, Addition_Tr, valeurnorm, Liste_corr, TITRE)


def readinput(INPUTNAME='input.csv', mode='numpy', concentrationinput='none'):
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
            Data = np.genfromtxt(INPUTNAME, comments='#', skip_header=1, delimiter=';', dtype='str'); 
        except UnicodeDecodeError:
            print("INPUT PAS UNICODE")
            Data = np.genfromtxt(INPUTNAME, comments='#', skip_header=1, delimiter=';', dtype='str', encoding='latin-1');
    # Récupération des données de base du fichier d'input
        #print(Data)
        Liste = (Data[:, 0].tolist())
        Legende = (Data[:, 1])
        Liste_ref = Data[:, 2]; # Chemin de la référence, si vide le fichier à déja été traiter
        try:
            Correction =  Data[:, 3].astype(np.float); 
        except ValueError:
            Correction = np.zeros(np.size(Liste))
        
        try:    optplt = Data[:, 4]
        except: optplt=''
        
        try :   MarqueurCIE = Data[:, 5]
        except: MarqueurCIE =''    
        
        
        try :       Addition_Tr = Data[:, 9].astype(np.float)
        except :    Addition_Tr = 0;
        
        try :       Epaisseur = Data[:, 6].astype(np.float)*1E-1 
        except :    Epaisseur = 0;
        
        try :       Concentration = Data[:, 7].astype(np.float)
        except :    Concentration=-1
        
        try :       Tx_redox = Data[:, 8].astype(np.float) # Le taux de rédox = [Fe2+]/[Fetot]
        except :    Tx_redox=1;
        
    elif mode == 'pandas' :
        Data=pd.read_csv(INPUTNAME, delimiter=';', comment='#', header=0)
            
        Liste =     Data.Nom_fichier.tolist()
        Legende =   Data.Legende.tolist()
        Liste_ref = Data.Nom_ref.tolist(); # Chemin de la référence, si vide le fichier à déja été traiter
        
        try:
            Correction = Data.Correction.tolist()  
        except AttributeError:
            Correction = np.zeros(np.size(Liste))
        
        try:        optplt = Data.Matplot_opt.tolist()
        except:     optplt=''
        
        try :       MarqueurCIE = Data.Marqueur_CIE.tolist()
        except:     MarqueurCIE ='';  
        
        
        try :       Addition_Tr = Data.Addition_Tr.tolist()
        except :    Addition_Tr = 0;

        try :       Epaisseur = Data.Epaisseur.tolist()
        except :    Epaisseur = 0;
        
        try :       Concentration = Data.Concentration.to_list()
        except :    Concentration=-1
        
        try :       Tx_redox = Data.Tx_redox.to_list()# Le taux de rédox = [Fe2+]/[Fetot]
        except :    Tx_redox=1;
        
    print(concentrationinput)
    
    if concentrationinput == 'epaisseur':
        valeurnorm = Epaisseur
    
    elif concentrationinput == 'molaire' : # Récupération des information suplémentaire du fichier d'input
        valeurnorm=Concentration*Epaisseur*Tx_redox
    
    elif concentrationinput == 'massique' : # Récupération des information suplémentaire du fichier d'input
        # Information suplémentaire    
        Tx_El = 1 # Taux de d'élément dans le composer (Ici fer dans FeO_2 : Fe/(Fe+2O))
        Masse_molaire = 58.69; # En g/mol pour le colorant (ici Ni)
        Densite = 2.488 #Densité du verre en kg/L (ici pour un verre à vitre)
        
        # On convertie les concentration en molaire par état d'oxydation
        Concentration_molaire = Concentration*Tx_El * Densite*1E3/(Masse_molaire)*1E-2
        # On passe d'une masse d'oxyde (ou d'élement brut) à un concentration molaire/volumique en mol/L,
 
        valeurnorm=Concentration_molaire*Epaisseur * Tx_redox
    else:
        valeurnorm = 1

    
    try : # Récupération de la liste des fichier traité si applicable
        Liste_temp=copy(Liste);
        for i in np.arange(0, np.size(Liste), 1):
            if Correction[i] != 0:
                try:
                    (cheminfichier, nomfichier) = os.path.split(Liste[i])
                    Liste_temp[i]='./Data_corr'+os.sep+nomfichier[:-4]+Corr2Str(Correction[i])
                except IndexError:
                    warnings.warn('\nLe dossier Data_corr est vide ou il y a trop de fichier par rapport au fichier d\'input\n')
            else:
                    Liste_temp[i]=Liste[i]            
    
        Liste_corr=np.array(Liste_temp);
    except FileNotFoundError:
        pass

    return(Liste, Legende, Liste_ref, Correction, optplt, MarqueurCIE, Addition_Tr, valeurnorm, Liste_corr, TITRE)



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


def Readspectre(Fichier, skip_header=2, delimiter=';'):
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
    Y = Data[:, 1]
    return([X, Y])

def natural_sort(l):#Fonction qui sert à trier dans l'ordre naturel (ou humain) les listes.
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

def nm2cm1(X):
    return(1/(X*1E-7))

def normYminmax(X, Y, COUPURENORMminmax=[400, 2500]):
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
    INDEX = np.logical_and(X>COUPURENORMminmax[0], X<COUPURENORMminmax[1])             
    Y = Y - np.min(Y[INDEX])
    Y = Y/np.max(Y[INDEX])
    return(Y)