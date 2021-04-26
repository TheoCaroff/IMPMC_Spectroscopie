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
from Nettoyage_spectre import Corr2Str
import pandas as pd

def readinput(INPUTNAME='input.csv', mode='numpy', concentrationinput='none'):
    
    TITRE = INPUTNAME[6:-4]

    if mode == 'numpy':
        try: # on ouvre en string pour récuperer le nom et le chemin des fichiers
            Data = np.genfromtxt(INPUTNAME, comments='#', skip_header=1, delimiter=';', dtype='str'); 
        except UnicodeDecodeError:
            print("INPUT PAS UNICODE")
            Data = np.genfromtxt(INPUTNAME, comments='#', skip_header=1, delimiter=';', dtype='str', encoding='latin-1');
    
    # Récupération des données de base du fichier d'input
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
        
        try :       Tx_redox = Data[:, 8].astype(np.float) # Le taux de rédox = [Fe2+]/[Fetot]
        except :    Tx_redox=1;
        
        
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