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
import panda as pd

def readinput(INPUTNAME='input.csv', mode='numpy', calculconcentration='none'):

    if mode == 'numpy':
        try: # on ouvre en string pour récuperer le nom et le chemin des fichiers
            Data = np.genfromtxt(INPUTNAME, comments='#', skip_header=1, delimiter=';', dtype='str'); 
        except UnicodeDecodeError:
            print("INPUT PAS UNICODE")
            Data = np.genfromtxt(INPUTNAME, comments='#', skip_header=1, delimiter=';', dtype='str', encoding='latin-1');
    
    # Récupération des données de base du fichier d'input
        Liste = (Data[:, 0])
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
        
        if calculconcentration == 'molaire' : # Récupération des information suplémentaire du fichier d'input
        
            Epaisseur = Data[:, 6].astype(np.float)*1E-1 # conversion en float, data d'origine en mm, conversion en cm
            Concentration_massique = Data[:, 7].astype(np.float) # 
            Tx_redox = Data[:, 8].astype(np.float) # Le taux de rédox = [Fe2+]/[Fetot]
        
            # Information suplémentaire    
            Tx_El = 1 # Taux de d'élément dans le composer (Ici fer dans FeO_2 : Fe/(Fe+2O))
            Masse_molaire = 58.69; # En g/mol pour le Nickel
            Densite = 2.488 #Densité du verre en kg/L
            
            # On convertie les concentration en molaire par état d'oxydation
            #ConFe2 = Concentration_massique*Tx_El*Tx_redox * Densite*1E3/(Masse_molaire)*1E-2
            ConFe2 = Concentration_massique*1E-2
            # On passe d'une masse d'oxyde (ou d'élement brut) à un concentration molaire/volumique en mol/L,
            #1E3 converti des g en kg, 1E-2 parceque pourcentage
            ConFe2 = Concentration_massique*1E-2 
            
            Concentration_molaire = Data[:, 7].astype(np.float) # 
            valeurnorm=Concentration_molaire*Epaisseur
        else:
            valeurnorm = 1
        
        
        # Récupération de la liste des fichier traité si applicable
    try :
        Liste_temp=copy(Liste.tolist());
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

    return(Liste, Legende, Liste_ref, Correction, optplt, MarqueurCIE, Addition_Tr, valeurnorm, Liste_corr)