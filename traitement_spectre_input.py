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
from Affichage_spectre import Affichage_Lab
from Affichage_spectre import mono2tab

from Nettoyage_spectre import Nettoyage_spectre
from Nettoyage_spectre import Corr2Str


INPUTNAME = 'input_Bleu_Musee.csv'
TITRE = INPUTNAME[6:-4]
#TITRE='pouet'

CALCULEPSI=False
CORRIGER=False;

# folder = askdirectory()
# os.chdir(folder)

if True: # Ouverture fichier d'input
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


    optplt = Data[:, 4]
    MarqueurCIE = Data[:, 5]
    #Addition_Tr= Data[:, 9].astype(np.float)

if CALCULEPSI: # Récupération des information suplémentaire du fichier d'input

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
#%% Partie absorbance
Modeaff='ABSnormep' # ABScm, Reflectance, Transmittance, Epsilon, ABSnormep

Autoaxe = False;
CORRIGER=CORRIGER
RAJOUT=''
TITRE=TITRE

X_min = 4000;
X_max = 30000;
Y_min = -0.1;
Y_max = 2;

AdditionTr=0#Addition_Tr

AdditionTr=mono2tab(AdditionTr, np.size(Liste))

modecouleurs='auto'; # 'auto', 'bigdata', 'manuel'
optplt=optplt

if not CORRIGER:
    RAJOUT = RAJOUT + '_brut'

if Modeaff == 'Transmittance': # Si on se met en transmittance
    X_min = 250;
    X_max = 2500;
    Y_min = 0;
    Y_max = 0.2;

if not (np.sum(AdditionTr)== 0):
    RAJOUT =  RAJOUT + '_+tr=' + str(AdditionTr)[2:]
else:
    RAJOUT = RAJOUT

SecondAxe=True; #choix double échelle mettre false pour avoir uniquement en cm^-1
TITRE = TITRE

if CORRIGER:
    Affichage_abs(Liste_corr, Legende, Autoaxe, [X_min, X_max], [Y_min,Y_max], SecondAxe,
              TITRE + RAJOUT, AdditionTr, valeurnorm=valeurnorm, Modeaff=Modeaff, modecouleurs=modecouleurs, optionplot=optplt, SHOW=True)
else:
    Affichage_abs(Liste, Legende, Autoaxe, [X_min, X_max], [Y_min,Y_max], SecondAxe,
              TITRE + RAJOUT, AdditionTr, Modeaff=Modeaff, modecouleurs=modecouleurs, optionplot=optplt)



#%% Traitement spectre

Liste_corr=Nettoyage_spectre(Liste, Legende, Liste_ref, Correction)
CORRIGER=True


#%% Affichage CIE

x, y = AffichageCIE1931(Liste_corr, Legende, TITRE, Marqueur=MarqueurCIE,Fleche=True,
                        xylim=[0, 0.33, 0, 0.33], show=False)

Fichierref='/media/veracrypt1/These_principal/Manip/Test_CIE/Effet_lc_3311XX/Beer_lambert_Bleu_3311.csv'

plt_BeerLambert_xy(Fichierref, 'lc Co2+ pur', optionplot='', show='True', TITRE=TITRE)

#%% Affichage Lab
Lab = Affichage_Lab(Liste_corr, Legende, TITRE, MarqueurCIE)


#%% Affichage espi

nbonde_min=5000;
nbonde_max=20000;
epsimin=0;
epsimax=40;

SecondAxe=True; #choix double échelle mettre false pour avoir uniquement en cm^-1

Index=(Correction!=0); # On ne prend en compte que les spectres voulu

Affichage_epsi(nbonde_min, nbonde_max, epsimin, epsimax, SecondAxe, Liste_corr[Index],
               Legende[Index], 0, ConFe2[Index], Epaisseur[Index], '$Fe^{2+}$')
