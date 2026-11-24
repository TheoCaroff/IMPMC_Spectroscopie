#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Théo Caroff
Se code sert au traitement des spectre Raman dans le forme cm; intensité avec 2 ligne de texte avant les données.
Placer le fichier input a la racine du repertoire contenant le dossier Data_trait contenant les spectres de la forme indiqué ci-dessus
Les fonction d'affichage sont dans Affichage_spectre.py
Les fonction de traitement dans Nettoyage_spectre.py

"""
#%%Chargement des modules et lecture du fichier contenant les paramètre exp
import numpy as np
import os
import warnings
from copy import copy
import fnmatch
from matplotlib import pyplot as plt
import pickle
import tkinter as tk
from tkinter.filedialog import askdirectory, askopenfilename
from Affichage_spectre  import circle_plot, rectangle_plot, annotate_clic_text, click_print_arrow
plt.style.use('default')
#plt.style.use('ggplot')
#plt.rcParams.keys() pour afficher les paramètres modifiable
setfont = lambda color : plt.style.use({'figure.facecolor': color, # yellow, white or none
                                       'savefig.transparent' : False,
                                       'lines.markersize': 6.0})  
from Affichage_spectre import Affichage_spectre
from Affichage_spectre import mono2tab
from Affichage_spectre import Affichage_gauss
from Affichage_spectre import Sav_fig
from Affichage_spectre import annot_max
from Affichage_spectre import IpVsMaxnu
from Affichage_spectre import DecnuVsMaxnu
from Nettoyage_spectre import Nettoyage_spectre_Raman

from Lecture_input import create_raman_txt
from Lecture_input import readinput_RAMAN as readinput
from Lecture_input import Chemin2input_RAMAN as Chemin2input
from Lecture_input import nm2cm1
from Lecture_input import Chemin2Liste

from Lecture_input import Readspectre

mode='input' # input, chemin, create_input

# TITRE = 'CDV_Blanc_alt' # Pour trier les spectres à afficher, si input_XXX.csv mettre TITRE=XXX et mode = input
TITRE = 'Jaune' # Pour trier les spectres à afficher, si input_XXX.csv mettre TITRE=XXX et mode = input
# TITRE = 'CDV_Vert_blanc_these'
# TITRE='temp' # Pour trier les spectres à afficher, si input_XXX.csv mettre TITRE=XXX et mode = input
CLEFDETRIE='*Jaune*'


if mode == 'askinput':
    tk.Tk().withdraw()
    # tkdialog = Tk().withdraw()
    filepath = askopenfilename(initialdir='.', filetypes= (("input","*.csv"), ("all files","*.*")))
    path, file = os.path.split(filepath)
    os.chdir(path)
    TITRE = file[6:-4]
    # tkdialog.destroy()
    print(TITRE)
    mode= 'input'
    
# DOSSIER = ['Data_trait_J1', 'Data_trait_J2', 'Data_trait_J3']
# DOSSIER = ['Data_trait_J4']
# DOSSIER = ['REF']
DOSSIER = ['Data_brut_txt']

if mode == 'create_input':
    # Si document dans plusieurs repertoir, indiquer le chemin dans un liste de string.
    Chemin2input(TITRE, CLEFDETRIE=CLEFDETRIE, LISTEDOSSIER=DOSSIER,
                 mode='RAMAN_portable', correction=3) #RAMAN_IPG, RAMAN_portable, RAMAN_IMPMC
    mode='input'
    
if mode == 'input' :
    INPUTNAME='input_'+TITRE+'.csv';
    (Liste, Legende, Liste_ref, Correction, optplt,
                     Liste_corr, TITRE, Liste_RefN, Marker) = readinput(INPUTNAME, comment='#')

#Legende=mono2tab(TITRE, 1)
                                            
CORRIGER=False;
plt.style.use({'figure.facecolor': 'white', 'savefig.transparent' : False, 'lines.markersize': 6.0})   
TitreOri = TITRE 

#%%
# %matplotlib auto
#%
plt.style.use({'figure.facecolor': 'white', 'savefig.transparent' : False, 'lines.markersize': 6.0})   
# plt.style.use({'figure.facecolor': 'none', 'savefig.transparent' : False, 'lines.markersize': 6.0})   

CORRIGER=True;


#%% Partie absorbance

# plt.style.use({'legend.loc': 'lower left'})
#plt.style.use({'lines.linestyle' : ':'})
CHEMINSAVFIG = '/home/theo/These/01_Redaction/Manuscrit/Chapter2/Figs/'
CHEMINSAVFIG = r'D:\Documents\These\01_Redaction\Manuscrit\_Fig_partA'

ModeX='cm' # 'cm', 'nm'
ModeY='Int' # Int, SubRubber, Gradient, Tr
normmode='' # min_max
   
 
modecouleurs='manuel' # 'auto', 'bigdata', 'manuel'
colormap   = 'plasma' # Utiliser pour mode bigdata nipy ou plasma
# pouet = [matplotlib.colors.rgb2hex(i) for i in colorsbigdata]

plt.style.use({'legend.fontsize': '7'})

TITRE = TitreOri

SHOW     = True ; LEGEND = True

Autoaxe     = True
legendin    = True; #False pour mettre la legende à l'extérieur de la figure
GRAPHOUT    = 5;     legncol = 1;
LISSAGE     = False ; # True, False ou un int >=4
INVERTAXE   = False;
CHECKBUTTON = True
SUBPLOT_TRANSITION  = False; ListeEl_subplot = ['Sibending', 'SiStretch', 'Postretching']# 

Sousmin     = False; #Pour soustraire le minimum
etagage     = False;

CORRIGER    = CORRIGER


X_min = 530 ;
X_max = 1250;
Ylim = [-0.0001, 0.002]
Ylim = None # [-4E-4, 1E-3]

COUPURENORMminmax = [300, 450]


if ModeY == 'Tr' :
    Y_min = -0.1;
    Y_max = 1;

if normmode == 'min_max' :
    Autoaxe = False;
    Y_min = -0.1;
    Y_max = 1.1;
    Ylim = [Y_min, Y_max]
if CORRIGER :
    Liste_aff=Liste_corr
    RAJOUT = ''
else:
    Liste_aff=Liste
    RAJOUT = '_BRUT'

fig, ax, TITRE, otherpara = Affichage_spectre(Liste_aff, Legende, Autoaxe, Xlim=[X_min, X_max], Ylim=Ylim,
                              TITRE=TITRE+RAJOUT, ModeX=ModeX, ModeY=ModeY, normmode=normmode,
                              spectromode='RAMAN', modecouleurs=modecouleurs, optionplot=optplt,
                              GRAPHOUT=GRAPHOUT, legendin=legendin,
                              ax1=None, newgraph=True, SHOW=SHOW, markevery=1,
                              COUPURENORMminmax=COUPURENORMminmax, colormap=colormap,
                              Sousmin=Sousmin,  etagage=etagage, 
                              Langue = 'anglais', Lissage=LISSAGE, third_dim=False,
                              Z_dim=Liste_RefN, Z_Label='%mol de $K_2O$',
                              CHECKBUTTON=CHECKBUTTON, SAVE_HTML=False, legncol=legncol,
                              SUBPLOT_TRANSITION=SUBPLOT_TRANSITION, ListeElement = ListeEl_subplot)
ax.set_title('')

if INVERTAXE : ax.invert_xaxis()

ax.set_ylabel('Absorbance (u.a)')

# annotate_clic_text(fig, ax)
# click_print_arrow(fig,ax)

if not SHOW and True:
    plt.title('')
    # ax.set_ylabel('Absorbance linéaire ($cm^{-1}$)') 
    # ax.set_ylabel('') 
    Sav_fig(TITRE, Repertoire = CHEMINSAVFIG, legend=LEGEND, legendin=legendin,
                      legncol=legncol, save_format='.png', PDF=True)



# ax1.set_title('Spectre Raman aprés soustraction du fond de fluorescence et normalisation à l\'air', y=-0.2)
# Sav_fig(TITRE, legendin=legendin)

#%% inset ax
#Affichage_gauss(Liste, Legende, TITRE+'_fitgauss', SHOW=True, optplt='--')
# plt.legend(loc='upper left')
# plt.show()

axins2 = ax.inset_axes([0.35, 0.1, 0.50, 0.40])
# axins2 = ax1.inset_axes([800, 0.7, 1500, 800], transform=ax.transData)
Autoaxe=False
SHOW = False
X_min = 80
X_max = 400
Y_min = 2000
Y_max = 8700
Ylim=[Y_min, Y_max]

pouet, axinstemp, pouet2 = Affichage_spectre(Liste_aff, Legende, Autoaxe, Xlim=[X_min, X_max], Ylim=Ylim,
                              TITRE=TITRE+RAJOUT, ModeX=ModeX, ModeY=ModeY, normmode=normmode,
                              spectromode='RAMAN', modecouleurs=modecouleurs, optionplot=optplt,
                              GRAPHOUT=GRAPHOUT, legendin=legendin,
                              ax1=axins2, newgraph=False, SHOW=SHOW, markevery=1,
                              COUPURENORMminmax=COUPURENORMminmax, colormap=colormap,
                              Sousmin=Sousmin,  etagage=etagage, 
                              Langue = 'Francais')

ax.indicate_inset_zoom(axins2, edgecolor="black", alpha=1)
axins2.set_xlabel('')
axins2.set_ylabel('')

#axins2.set_xticklabels([X_min])
axins2.set_yticklabels([])

Sav_fig(TITRE+'zoom', legendin=legendin)
if False:
    Sav_fig(TITRE, Repertoire = CHEMINSAVFIG, legend=LEGEND, legendin=legendin,
                      legncol=legncol, save_format='.png', PDF=True)

#%% Ip vs Numax

# Mettre la concentration X dans colonne refN.
print(Liste_RefN)

Liste_RefN = np.array(Liste_RefN)

IpVsMaxnu(Liste_corr, Legende, optplt, Marker=None, legncol=1)



percentagemol = [float(i) for i in Liste_ref]
DecnuVsMaxnu(Liste_corr, Legende, optplt, Marker=None, molpercentage=Liste_RefN)

#%% Traitement spectre
%matplotlib inline
plt.style.use({'legend.loc': 'lower right'})
# !! sinon beug dans spyder, a n'utiliser qu'avec spyder

DOSSAVE=''

Liste_corr=Nettoyage_spectre_Raman(Liste, Legende, Liste_ref, Correction, DOSSAVE)
#Liste_corr=Nettoyage_spectre_Raman(Liste, Legende, Liste_ref, mono2tab(2, np.size(Liste)), DOSSAVE)

CORRIGER=True
plt.style.use({'figure.facecolor': 'white', 'savefig.transparent' : False, 'lines.markersize': 6.0})   

#%%
# def create_empty_txt(files):
#     """
#     Pour chaque chemin dans `files` :
#     - on retire l'extension d'origine
#     - on crée un fichier vide avec l'extension .txt
#     """
#     for f in files:
#         p = Path(f)
#         txt_file = p.with_suffix('.txt')
#         try:
#             txt_file.touch(exist_ok=False)  # crée le fichier, lève une erreur si existe déjà
#             print(f"Créé : {txt_file}")
#         except FileExistsError:
#             print(f"Existant, non écrasé : {txt_file}")
            
# create_empty_txt(Liste)
#%% Manipulation des fichier !!! copie des trucs

# Repsav = 'CDV_C3_Vert'
# try : os.mkdir(Repsav)
# except : pass

# for element in Liste_corr :
#     head, tail = os.path.split(element)
#     os.system('cp ' + element + ' ' + Repsav +os.sep + tail)

#%%
create_raman_txt('../Data_brut')