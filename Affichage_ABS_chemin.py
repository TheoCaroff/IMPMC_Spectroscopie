# -*- coding: utf-8 -*-
"""
@author: Théo Caroff
"""
#Ce script supperpose des fichier d'absorbance.
import numpy as np
import os
from matplotlib import pyplot as plt

from tkinter.filedialog import askdirectory
import tempfile
import fnmatch

from Affichage_spectre import Affichage_abs
from Affichage_spectre import AffichageCIE1931


from Nettoyage_spectre import natural_sort

Recuper_nom_dossier_temporaire=True;

TITRE='med' # Titre du graph
CLEFDETRIE = '*bleu*'

#DOSSIER = 'Data_corriger_appareil' # Dossier avec les data et UNIQUEMENT LES DATAs
DOSSIER = 'Data_trait'

Autoaxe = False;
SecondAxe = True;
RAJOUT=''

X_min = 4000;
X_max = 30000;
Y_min = -0.1;
Y_max = 2;

Modeaff='ABScm' # ABScm, Reflectance, Transmittance, Epsilon, ABSnormep
AdditionTr=0#Addition_Tr
modecouleurs='auto'; # 'auto', 'bigdata', 'manuel'

# Soustration_ref =   False; # Pour soustraire par rapport à la première coubre dans l'ordre naturelle de trie
# Courbenorm      =   0; # Courbe qui sert de référence, attention commence à 0.

if Modeaff == 'Transmittance': # Si on se met en transmittance
    X_min = 250;
    X_max = 2500;
    Y_min = 0;
    Y_max = 0.2;

if not (np.sum(AdditionTr)== 0):
    RAJOUT =  RAJOUT + '_+tr=' + str(AdditionTr)[2:]
else:
    RAJOUT = RAJOUT



if Recuper_nom_dossier_temporaire:
    repspectro = open(tempfile.gettempdir()+os.sep+'repspectro', 'r');
    folder = repspectro.read();
    repspectro.close();

else:
    folder = '.'
    #folder = askdirectory()


os.chdir(folder)
folder = folder + os.sep + DOSSIER
Liste = os.listdir(folder); #Récupère la liste des fichiers
Liste = natural_sort(Liste)
Liste = fnmatch.filter(Liste, CLEFDETRIE)

Legende=[x[:-4] for x in Liste]

Liste=[DOSSIER + os.sep + x for x in Liste]



Affichage_abs(Liste, Legende, Autoaxe, [X_min, X_max], [Y_min,Y_max], SecondAxe,
              TITRE + RAJOUT, AdditionTr, Modeaff=Modeaff, modecouleurs=modecouleurs, SHOW=True)
