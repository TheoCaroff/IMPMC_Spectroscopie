# -*- coding: utf-8 -*-
"""
@author: Théo Caroff
Se code sert au traitement des spectre optiques dans le forme nm; Tr avec 2 ligne de texte avant les donénes.
Placer le fichier input a la racine du repertoire contenant le dossier Data_brut contenant les spectres de la forme indiqué ci-dessus
Les fonction d'affichage sont dans Affichage_spectre.py
Les fonction de traitement dans Nettoyage_spectre.py
"""
#Chargement des modules et lecture du fichier contenant les paramètre exp
import numpy as np
import os
import warnings
from copy import copy
from tkinter.filedialog import askdirectory
import fnmatch
from matplotlib import pyplot as plt

from Affichage_spectre import Affichage_abs
from Affichage_spectre import AffichageCIE1931
from Affichage_spectre import plt_BeerLambert_xy

from Nettoyage_spectre import Nettoyage_spectre
from Nettoyage_spectre import Corr2Str


INPUTNAME = 'input.csv'
TITRE = INPUTNAME[6:-4]
TITRE='pouet'

CALCULEPSI=False

# folder = askdirectory()
# os.chdir(folder)

if True: # Ouverture fichier d'input
    try: # on ouvre en string pour récuperer le nom et le chemin des fichier
        Data = np.genfromtxt(INPUTNAME, comments='#', skip_header=1, delimiter=';', dtype='str'); 
    except UnicodeDecodeError:
        print("INPUT PAS UNICODE")
        Data = np.genfromtxt(INPUTNAME, comments='#', skip_header=1, delimiter=';', dtype='str', encoding='latin-1');

# Récupération des données de base du fichier d'input
    Liste = (Data[:, 0])
    Legende = (Data[:, 1])
    Liste_ref = Data[:, 2]; # Chemin de la référence, si vide le fichier à déja été traiter
    try:
        Correction =  Data[:, 3].astype(np.float); # On enlève la bande uv, 1=oui, 0 =non;
    except ValueError:
        Correction = np.zeros(np.size(Liste))
    optplt = Data[:, 4]
    MarqueurCIE = Data[:, 5]
    Addition_Tr= Data[:, 6].astype(np.float)
if CALCULEPSI: # Récupération des information suplémentaire du fichier d'input

    Epaisseur = Data[:, 4].astype(np.float)*1E-1 # conversion en float, data d'origine en mm, convesion en cm
    Concentration_massique = Data[:, 3].astype(np.float) # 
    Tx_redox = Data[:, 6].astype(np.float) # Le taux de rédox = [Fe2+]/[Fetot]


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


# Récupération de la liste des fichier traité si applicable
try :
    Liste_temp=copy(Liste.tolist());
    for i in np.arange(0, np.size(Liste), 1):
        if Correction[i] != 0:
            try:
                (cheminfichier, nomfichier) = os.path.split(Liste[i])
                Liste_temp[i]='./Data_corriger'+os.sep+nomfichier[:-4]+Corr2Str(Correction[i])
                CORRIGER=True
            except IndexError:
                warnings.warn('\nLe dossier Data_corriger est vide ou il y a trop de fichier par rapport au fichier d\'input\n')
        else:
                Liste_temp[i]=Liste[i]            

    Liste_corr=np.array(Liste_temp);
except FileNotFoundError:
    pass
#%% Partie absorbance
ABScm=True
Autoaxe = False;
CORRIGER=CORRIGER

nbonde_min = 4000;
nbonde_max = 25000;
Amin = -2;
Amax = 0.5;

AdditionTr=0.0000
RAJOUT=''
modecouleurs='auto'; # 'auto', 'bigdata', 'manuel'


if CORRIGER:
    RAJOUT = RAJOUT #+ '_corriger'

if not ABScm: # Si on se met en transmittance
    nbonde_min = 350;
    nbonde_max = 2500;
    Amin = 0;
    Amax = 1.2;

if not (AdditionTr == 0):
    RAJOUT =  RAJOUT + '_+tr=' + str(AdditionTr)[2:]
else:
    RAJOUT = RAJOUT

SecondAxe=True; #choix double échelle mettre false pour avoir uniquement en cm^-1
TITRE = TITRE

if CORRIGER:
    Affichage_abs(Liste_corr, Legende, Autoaxe, nbonde_min, nbonde_max, Amin, Amax, SecondAxe,
              TITRE + RAJOUT, AdditionTr, ABScm=ABScm, modecouleurs=modecouleurs, optionplot=optplt)
else:
    Affichage_abs(Liste, Legende, Autoaxe, nbonde_min, nbonde_max, Amin, Amax, SecondAxe,
              TITRE + RAJOUT, AdditionTr, ABScm=ABScm, modecouleurs=modecouleurs, optionplot=optplt)



#%% Traitement spectre

Liste_corr=Nettoyage_spectre(Liste, Legende, Liste_ref, Correction, Addition_Tr=Addition_Tr)

#%% Affichage CIE

x, y = AffichageCIE1931(Liste_corr, Legende, TITRE, Marqueur=MarqueurCIE,Fleche=True, xylim=[0, 0.25, 0, 0.2], show=True)

Fichierref='/media/veracrypt1/These_principal/Manip/Test_CIE/Effet_lc_3311XX/Beer_lambert_Bleu_3311.csv'

#plt_BeerLambert_xy(Fichierref, 'lc Co2+ pur', optionplot='', show='True')



#%% Affichage espi

nbonde_min=5000;
nbonde_max=20000;
epsimin=0;
epsimax=40;

SecondAxe=True; #choix double échelle mettre false pour avoir uniquement en cm^-1

Index=(Correction!=0); # On ne prend en compte que les spectres voulu

Affichage_epsi(nbonde_min, nbonde_max, epsimin, epsimax, SecondAxe, Liste_corr[Index],
               Legende[Index], 0, ConFe2[Index], Epaisseur[Index], '$Fe^{2+}$')
