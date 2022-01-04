# -*- coding: utf-8 -*-
"""
@author: Théo Caroff
Se code sert au traitement des spectre optiques dans le forme nm; Tr avec 2 ligne de texte avant les données.
Placer le fichier input a la racine du repertoire contenant le dossier Data_trait contenant les spectres de la forme indiqué ci-dessus
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
from Affichage_spectre import Affichage_Lab2D
from Affichage_spectre import mono2tab
from Affichage_spectre import Courbe_BeerLambert_XY
from Affichage_spectre import Calcul_BeerLambert_XYZ
from Affichage_spectre import Courbe_BeerLambert_Lab3D
from Affichage_spectre import Affichage_Lab3D
from Affichage_spectre import Affichage_gauss

from Nettoyage_spectre import Nettoyage_spectre

from Lecture_input import readinput
from Lecture_input import nm2cm1
from Lecture_input import Chemin2Liste
from Lecture_input import Chemin2input

mode='input' # input, chemin, create_input

TITRE='test' # Pour trier les spectres à afficher, si input_XXX.csv mettre TITRE=XXX et mode = input
CLEFDETRIE='*'#'*ref[0-9]*.csv'

Recuperer_nom_dossier_temporaire=False;

# folder = askdirectory()
# os.chdir(folder)

DOSSIER='Data_trait'

if mode == 'create_input':
    # Si document dans plusieurs repertoir, indiquer le chemin dans un liste de string.
    DOSSIER=[DOSSIER]
    Chemin2input(TITRE, CLEFDETRIE=CLEFDETRIE, LISTEDOSSIER=DOSSIER,
                 mode='portable', correction=6)
    mode='input'
    
if mode == 'input' :
    INPUTNAME='input_'+TITRE+'.csv';
    (Liste, Legende, Liste_ref, Liste_refN, Correction, optplt, MarqueurCIE, Addition_Tr,
             valeurnorm, Liste_corr, TITRE) = readinput(INPUTNAME, mode='numpy',
                                                        concentrationinput='epaisseur')
else :
    (Liste, Legende, Liste_ref, Liste_refN, Correction, optplt, MarqueurCIE, Addition_Tr,
             valeurnorm, Liste_corr, TITRE) = Chemin2Liste(TITRE,
             Recuperer_nom_dossier_temporaire, DOSSIER=DOSSIER, CLEFDETRIE=CLEFDETRIE)

#Legende=mono2tab(TITRE, 1)
                                                           
CORRIGER=True;

#%% Partie absorbance
Modeaff='ABScm' # ABScm, ABSnm, ABSnorm_min_max, SubBaseline, Reflectance
#Transmittance, Epsilon, ABSnormep, Gradient
modecouleurs='auto'; # 'auto', 'bigdata', 'manuel'


Autoaxe     = True;
SecondAxe   = True; #choix double échelle mettre false pour avoir uniquement en cm^-1

CORRIGER    = CORRIGER
RAJOUT=''
TITRE=TITRE
optplt=''

X_min = nm2cm1(3000);
X_max = nm2cm1(370);
Y_min = -0.0;
Y_max = 3;


Addition_Tr=0
Addition_Tr=mono2tab(Addition_Tr, np.size(Liste))

if Modeaff == 'Transmittance' or Modeaff=='ABSnm': # Si on se met en nm
    X_min = 300;
    X_max = 1100;
    Y_min = 0;
    Y_max = 1.5;
    
elif Modeaff == 'ABSnorm_min_max' or Modeaff == 'SubBaseline':
    Autoaxe = False;
    X_min = nm2cm1(10000);
    X_max = nm2cm1(370);
    Y_min = -0.1;
    Y_max = 1.2;
COUPURENORMminmax=[1000, 2500]

if not (np.sum(Addition_Tr)== 0):
    RAJOUT =  RAJOUT + '_+tr_' + str(Addition_Tr[0])[2:]

if CORRIGER :
    Liste_aff=Liste_corr
else:
    Liste_aff=Liste
    RAJOUT = RAJOUT + '_brut'


Affichage_abs(Liste_aff, Legende, Autoaxe, [X_min, X_max], [Y_min,Y_max], SecondAxe,
              TITRE + RAJOUT, Addition_Tr, valeurnorm=valeurnorm, Modeaff=Modeaff,
              modecouleurs=modecouleurs, optionplot=optplt, SHOW=True, GRAPHOUT=2,
              COUPURENORMminmax=COUPURENORMminmax, newgraph=True)

#%%
#Affichage_gauss(Liste, Legende, TITRE+'_fitgauss', SHOW=True, optplt='--')

#%% Affichage CIE

x, y = AffichageCIE1931(Liste_corr, Legende, TITRE, Marqueur=MarqueurCIE,Fleche=False,
                        xylim=[0.3, 0.8, 0.3, 0.7], show=True)

# Courbe_BeerLambert_XY(Liste_corr[0], Legende[0], TITRE=TITRE, show=False)
# Courbe_BeerLambert_XY(Liste_corr[1], Legende[1], TITRE=TITRE, show=True)


#%% Affichage Lab

Lab = Affichage_Lab2D(Liste_corr, Legende, TITRE, MarqueurCIE)

Lab_2 = Affichage_Lab3D(Liste, Legende, TITRE, SHOW=True)

#Lab_3 = Courbe_BeerLambert_Lab3D(Liste_corr[0], TITRE=TITRE, newfig=False)

#%% Traitement spectre
#Correction = [2]
Liste_corr=Nettoyage_spectre(Liste, Legende, Liste_ref, Correction)
CORRIGER=True
