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
from Lecture_input import Chemin2input

mode='input' # input, chemin, create_input

TITRE='POUET' # Pour trier les spectres à afficher, si input_XXX.csv mettre TITRE=XXX et mode = input

Recuperer_nom_dossier_temporaire=False;

# folder = askdirectory()
# os.chdir(folder)

if mode == 'input' :
    INPUTNAME='input_'+TITRE+'.csv';
    (Liste, Legende, Liste_ref, Correction, optplt, MarqueurCIE, Addition_Tr,
             valeurnorm, Liste_corr, TITRE) = readinput(INPUTNAME, mode='numpy',
                                                        concentrationinput='epaisseur')
elif mode == 'create_input':
    Chemin2input(TITRE, DOSSIER='Data_trait')
    
else :
    (Liste, Legende, Liste_ref, Correction, optplt, MarqueurCIE, Addition_Tr,
             valeurnorm, Liste_corr, TITRE) = Chemin2Liste(TITRE,
             Recuperer_nom_dossier_temporaire,DOSSIER='Data_trait')
    
CORRIGER=False;

#%% Partie absorbance
Modeaff='ABScm' # ABScm, ABSnm, ABSnorm_min_max, Reflectance, Transmittance, Epsilon, ABSnormep
modecouleurs='manuel'; # 'auto', 'bigdata', 'manuel'


Autoaxe     = True;
SecondAxe   = True; #choix double échelle mettre false pour avoir uniquement en cm^-1

CORRIGER    = CORRIGER
RAJOUT=''
TITRE=TITRE

X_min = nm2cm1(2500);
X_max = nm2cm1(300);
Y_min = -0.1;
Y_max = 1.5;

Addition_Tr=Addition_Tr
Addition_Tr=mono2tab(Addition_Tr, np.size(Liste))

if Modeaff == 'Transmittance' or Modeaff=='ABSnm': # Si on se met en nm
    X_min = 300;
    X_max = 1000;
    Y_min = 0;
    Y_max = 1;

if not (np.sum(Addition_Tr)== 0):
    RAJOUT =  RAJOUT + '_+tr_' + str(Addition_Tr[0])[2:]

if CORRIGER :
    Liste_aff=Liste_corr
else:
    Liste_aff=Liste
    RAJOUT = RAJOUT + '_brut'


Affichage_abs(Liste_aff, Legende, Autoaxe, [X_min, X_max], [Y_min,Y_max], SecondAxe,
              TITRE + RAJOUT, Addition_Tr, valeurnorm=valeurnorm, Modeaff=Modeaff,
              modecouleurs=modecouleurs, optionplot=optplt, SHOW=True, COUPURENORMminmax=[400, 800])



#%% Traitement spectre

Liste_corr=Nettoyage_spectre(Liste, Legende, Liste_ref, Correction)
CORRIGER=True


#%% Affichage CIE

x, y = AffichageCIE1931(Liste_corr, Legende, TITRE, Marqueur=MarqueurCIE,Fleche=False,
                        xylim=[0, 0.33, 0, 0.33], show=False)

Fichierref='/media/veracrypt1/These_principal/Manip/Test_CIE/Effet_lc_3311XX/Beer_lambert_Bleu_3311.csv'

plt_BeerLambert_xy(Fichierref, 'lc Co2+ pur', optionplot='', show='True', TITRE=TITRE)

#%% Affichage Lab
Lab = Affichage_Lab(Liste_corr, Legende, TITRE, MarqueurCIE)

