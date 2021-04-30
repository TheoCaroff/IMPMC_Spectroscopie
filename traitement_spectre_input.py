# -*- coding: utf-8 -*-
"""
@author: Théo Caroff
Se code sert au traitement des spectre optiques dans le forme nm; Tr avec 2 ligne de texte avant les donénes.
Placer le fichier input a la racine du repertoire contenant le dossier Data_brut contenant les spectres de la forme indiqué ci-dessus
Les fonction d'affichage sont dans Affichage_spectre.py
Les fonction de traitement dans Nettoyage_spectre.py
"""
#%%Chargement des modules et lecture du fichier contenant les paramètre exp
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

from Lecture_input import readinput
from Lecture_input import nm2cm1
from Lecture_input import Chemin2Liste

mode='chemin' # input ou chemin

TITRE='Fibre' # Pour trier les spectres à afficher

INPUTNAME = 'input_BK2.csv'

Recuperer_nom_dossier_temporaire=False;


# folder = askdirectory()
# os.chdir(folder)

if mode == 'input' : 
    (Liste, Legende, Liste_ref, Correction, optplt, MarqueurCIE, Addition_Tr,
             valeurnorm, Liste_corr, TITRE) = readinput(INPUTNAME, mode='numpy', concentrationinput='epaisseur')
else :
    (Liste, Legende, Liste_ref, Correction, optplt, MarqueurCIE, Addition_Tr,
             valeurnorm, Liste_corr, TITRE) = Chemin2Liste(TITRE, Recuperer_nom_dossier_temporaire)

Correction = mono2tab(3, np.size(Liste));

CORRIGER=False;

#%% Partie absorbance
Modeaff='Transmittance' # ABScm, Reflectance, Transmittance, Epsilon, ABSnormep

Autoaxe = True;
SecondAxe=False; #choix double échelle mettre false pour avoir uniquement en cm^-1

CORRIGER=False
RAJOUT=''
TITRE=TITRE

X_min = 5000;
X_max = 32000;
Y_min = -0;
Y_max = 1;

Addition_Tr=Addition_Tr

Addition_Tr=mono2tab(Addition_Tr, np.size(Liste))

modecouleurs='auto'; # 'auto', 'bigdata', 'manuel'
optplt=optplt

if not CORRIGER:
    RAJOUT = RAJOUT + '_brut'

if Modeaff == 'Transmittance': # Si on se met en transmittance
    X_min = 250;
    X_max = 2500;
    Y_min = 0;
    Y_max = 1.5;

if not (np.sum(Addition_Tr)== 0):
    RAJOUT =  RAJOUT + '_+tr_' + str(Addition_Tr[0])[2:]
else:
    RAJOUT = RAJOUT


TITRE = TITRE

if CORRIGER:
    Affichage_abs(Liste_corr, Legende, Autoaxe, [X_min, X_max], [Y_min,Y_max], SecondAxe,
              TITRE + RAJOUT, Addition_Tr, valeurnorm=valeurnorm, Modeaff=Modeaff, modecouleurs=modecouleurs, optionplot=optplt, SHOW=True)
else:
    Affichage_abs(Liste, Legende, Autoaxe, [X_min, X_max], [Y_min,Y_max], SecondAxe,
              TITRE + RAJOUT, Addition_Tr, valeurnorm=valeurnorm, Modeaff=Modeaff, modecouleurs=modecouleurs, optionplot=optplt, SHOW=True)



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

