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
from tkinter.filedialog import askdirectory, askopenfilename

import tkinter as tk
import fnmatch
import pickle
from matplotlib import pyplot as plt
plt.style.use('default')
#plt.style.use('ggplot')
#plt.rcParams.keys() pour afficher les paramètres modifiable
setfont = lambda color : plt.style.use({'figure.facecolor': color, # yellow, white or none
                                       'savefig.transparent' : False,
                                       'lines.markersize': 6.0})  

from Affichage_spectre import Affichage_spectre
from Affichage_spectre import AffichageCIE1931
from Affichage_spectre import Affichage_Lab2D
from Affichage_spectre import mono2tab
from Affichage_spectre import Courbe_BeerLambert_XY
from Affichage_spectre import Calcul_BeerLambert_XYZ
from Affichage_spectre import Courbe_BeerLambert_Lab3D
from Affichage_spectre import Affichage_Lab3D
from Affichage_spectre import Affichage_gauss
from Affichage_spectre import Sav_fig
from Affichage_spectre import annot_max
from Affichage_spectre import Liste2dendrogram
from Affichage_spectre  import circle_plot, rectangle_plot, annotate_clic_text

from Nettoyage_spectre import Nettoyage_spectre_OAS as Nettoyage_spectre

from Lecture_input import readinput_OAS as readinput
from Lecture_input import nm2cm1
from Lecture_input import Chemin2Liste
from Lecture_input import Chemin2input_OAS as Chemin2input
from Lecture_input import Readspectre

close = lambda : plt.close('all')

mode='input' # input, chemin, create_input, askinput

TITRE = 'grisaille_blanche' # Pour trier les spectres à afficher, si input_XXX.csv mettre TITRE=XXX et mode = input

CLEFDETRIE='*'


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
    
Recuperer_nom_dossier_temporaire=False;

# folder = askdirectory()
# os.chdir(folder)

# DOSSIER = ['Data_corr']
DOSSIER = ['20240418_Mesure_sphere']
# DOSSIER = ['C3']
# DOSSIER = ['Data_brut']
# DOSSIER =['Ref']
# DOSSIER =['Data_trait_Mn_PT1', 'Data_trait_Mn_PT2', 'Data_trait_Mn_PT3', 'Data_trait_Mn_PT4']


if mode == 'create_input':
    # Si document dans plusieurs repertoir, indiquer le chemin dans un liste de string.
    Chemin2input(TITRE, CLEFDETRIE=CLEFDETRIE, LISTEDOSSIER=DOSSIER,
                 mode='PERKIN', correction=0) #PERKIN, portable, portable_trait
    mode='input'
    
if mode == 'input' :
    INPUTNAME='input_'+TITRE+'.csv';
    (Liste, Legende, Liste_ref, Liste_refN, Correction, optplt, MarqueurCIE, Addition_Tr,
             valeurnorm, Liste_corr, TITRE, annot, Limite) = readinput(INPUTNAME, mode='pandas',
                                                        concentrationinput='epaisseur') # epaisseur ou massique
else :
    (Liste, Legende, Liste_ref, Liste_refN, Correction, optplt, MarqueurCIE, Addition_Tr,
             valeurnorm, Liste_corr, TITRE, annot) = Chemin2Liste(TITRE,
             Recuperer_nom_dossier_temporaire, LISTEDOSSIER=DOSSIER, CLEFDETRIE=CLEFDETRIE)
             

CORRIGER=False;
plt.style.use('default')   
plt.style.use({'figure.facecolor': 'yellow', 'savefig.transparent' : False, 'lines.markersize': 6.0})  
TitreOri=TITRE 
#%
#Pour save env travaile :   dill.dump_session(TITRE + '.pkl')
# pour load :               dill.load_session(TITRE + '.pkl')

#%
# %matplotlib auto
#%
# plt.style.use({'figure.facecolor': 'white', 'savefig.transparent' : False, 'lines.markersize': 6.0})   
plt.style.use({'figure.facecolor': 'none', 'savefig.transparent' : False, 'lines.markersize': 6.0})   

CORRIGER=True;

#%% Partie absorbance

CHEMINSAVFIG = '/home/theo/These/01_Redaction/Manuscrit/Chapter6/Figs/'
CHEMINSAVFIG = r'D:\Documents\These\01_Redaction\Manuscrit\_Fig_partB'
# CHEMINSAVFIG = '/home/theo/These/01_Redaction/Manuscrit/Appendix2/Figs/'
plt.style.use('default')
# plt.style.use({'figure.facecolor': 'none', 'savefig.transparent' : False, 'lines.markersize': 6.0})   

# plt.style.use({'legend.loc': 'upper right'})
# plt.style.use({'legend.fontsize': '8'})

SHOW     = True ; SAVEthese = True; LEGEND = True

ModeX    = 'nm' # 'cm', 'nm', 'eV'
ModeY    = 'Tr' # 'ABS', 'Tr', SubRubber, Reflectance, Kubelka,  Epsilon, Gradient, ABSnormep, FondRubber, DoubleGradient
normmode = '' # min_max , ep


modecouleurs = 'manuel' # 'auto', 'bigdata', 'manuel'
colormap   = 'nipy' # Utiliser pour mode bigdata nipy ou plasma, ocean, GnBu
     
Autoaxe     = True
SecondAxe   = True; SpectreVIS  = True; 
legncol     = 1; legendin = True; #False pour mettre la legende à l'extérieur de la fiure
GRAPHOUT    = 6; #5 pour 2im col g/d, 6 pour im plein page large

Sousmin     = False; #Pour soustraire le minimum
etagage     = False;

CORRIGER    = CORRIGER
LWARN               = False
SUBPLOT_TRANSITION  = False ; ListeEl_subplot = ['Fe2+']#['Co2+', 'Cu2+', 'Fe2+', 'Fe3+']# ['Mn3+', 'Fe3+']
ANNOT               = False
LISAGE              = False
CHECKBUTTON         = False # Pour afficher la liste des courbes à afficher ou non

X_min = 4E3; 
X_max = 33E3
Y_min = 0
Y_max = 0.7
Ylim = [Y_min, Y_max]
# Ylim = None
 
COUPURENORMminmax = [5, 19000]

if ModeX == 'nm' : # Si on se met en nm
    SecondAxe = False;
    X_min = 300;
    X_max = 2500;
    COUPURENORMminmax = [X_min, X_max]

if normmode == 'min_max' :
    Autoaxe = True;
    # X_min = 16000;
    # X_max = nm2cm1(400);
    Y_min = 0; Y_max = 2.5;

if ModeY == 'Tr' :
    Y_min = -0.1; Y_max = 100; Ylim = [Y_min, Y_max]

if CORRIGER : Liste_aff=Liste_corr; RAJOUT = ''
else:
    Liste_aff=Liste
    if not int(Correction.sum()) == 0: RAJOUT = ''
    else : RAJOUT = ''

fig, ax1, TITRE, otherpara = Affichage_spectre(Liste_aff, Legende, Autoaxe, Xlim=[X_min, X_max], Ylim=Ylim,
                              TITRE=TITRE+RAJOUT, ModeX=ModeX, ModeY=ModeY, normmode=normmode,
                              spectromode='OAS', modecouleurs=modecouleurs, optionplot=optplt,
                              colormap=colormap, GRAPHOUT=GRAPHOUT, legendin=legendin,
                              ax1=None, newgraph=True, SHOW=SHOW, markevery=1,
                              COUPURENORMminmax=COUPURENORMminmax, 
                              Sousmin=Sousmin,  etagage=etagage,
                              Langue = 'francais', valeurnorm=valeurnorm, annotation=annot,
                              SecondAxe=SecondAxe, Dosage=False, Lissage=LISAGE, LWARN=LWARN, legncol=legncol,
                              SpectreVIS = SpectreVIS, SAVE_HTML=False, LimiteAff=Limite,
                              SUBPLOT_TRANSITION=SUBPLOT_TRANSITION, ListeElement=ListeEl_subplot, CHECKBUTTON=CHECKBUTTON)

# ax.set_title('', y=-0.2)
# Sav_fig(TITRE, legendin=legendin)
# annotate_clic_text(fig, ax1)
# lasso = rectangle_plot(ax1, 'red', 0.5)
# lasso = rectangle_plot(ax1, 'green', 0.5)
# ax1.legend()

    
if not SHOW and SAVEthese:
    plt.title('')
    # ax1.set_ylabel('Absorbance linéaire ($cm^{-1}$)') 
    ax1.set_ylabel('Absorbance (u.a)')
    # fig.adjust_subplot(top=0.875, bottom=0.015, left=0.125, right=0.925, hspace=0.1, wspace=0.0)
    # ax1.set_ylabel('$\\varepsilon (L.mol^{-1}.cm^{-1})$') 
    # ax1.set_ylabel('') 
    ax1.legend()
    Sav_fig(TITRE, Repertoire = CHEMINSAVFIG, legend=LEGEND, legendin=legendin,
                      legncol=legncol, save_format='.png', PDF=True, fig=fig)
    # file = open('Para_'+TITRE, 'wb')
    # pickle.dump((TITRE, X_min, X_max, Ylim, ModeX, ModeY, normmode, COUPURENORMminmax,
    #              SecondAxe, SpectreVIS, legncol, legendin, GRAPHOUT, Sousmin), file)
    # file.close()
TITRE = TitreOri



#%% Affichage CIE
plt.style.use('default')
SHOW        = True;
legendin    = True;
x, y = AffichageCIE1931(Liste_aff, Legende, TITRE, Marqueur=MarqueurCIE,Fleche=False,
                        xylim=[0, 1, 0, 1], show=SHOW,
                        legendin=legendin, GRAPHOUT=1) #  xylim = [0.05, 0.25, 0, 0.20] ou bleu ou pourpre

plt.style.use({'lines.linestyle' : ':', 'lines.linewidth' : 1.5})
plt.style.use({'legend.loc': 'upper left'})

# Ncourbe=0
# Courbe_BeerLambert_XY(Liste_corr[Ncourbe], Legende[Ncourbe],
#                     TITRE=TITRE, show=SHOW, legendin=legendin, style=':w')
# Ncourbe=0
# Courbe_BeerLambert_XY(Liste_corr[Ncourbe], Legende[Ncourbe],
#                       TITRE=TITRE, show=True, legendin=legendin)

if not SHOW : 
    plt.title('')
    plt.tight_layout()
    Sav_fig(TITRE+'_CIE1931', Repertoire = CHEMINSAVFIG,
                      legendin=legendin, legend=True, legncol=1, colorplot=True, PDF=True)

plt.style.use('default')
#%% Affichage Lab
plt.style.use({'figure.facecolor': 'none', 'savefig.transparent' : False})
# Llim = None
# alim = None
# blim = None

legendin=False;
SHOW = True
Data = Affichage_Lab2D(Liste_aff, Legende, TITRE, MarqueurCIE, legendin=legendin, SHOW=SHOW)
if not SHOW : Sav_fig(TITRE+'_LAB2D', Repertoire = CHEMINSAVFIG, legendin=legendin,
                      legncol=legncol, Tight=False,legend=False)

plt.style.use('default')

legendin=True;
Data2 = Affichage_Lab3D(Liste_aff, Legende, TITRE, MarqueurCIE, SHOW=True,
                legendin=legendin, GRAPHOUT=1)

# Courbe_BeerLambert_Lab3D(Liste_corr[0], TITRE=TITRE,newfig=False,
#                           legendin=legendin, show=False)
# Courbe_BeerLambert_Lab3D(Liste_corr[3], TITRE=TITRE,newfig=False,
#                           legendin=legendin)
#%% Traitement spectre  
%matplotlib inline
# !! sinon beug dans spyder, a n'utiliser qu'avec spyder

DOSSAVE=''

Borne=[11E3, 31E3]

Liste_corr=Nettoyage_spectre(Liste, Legende, Liste_ref, Correction,
                             Liste_refN, DOSSAVE, Borne=Borne,
                             valeurnorm=valeurnorm, DataExcel=False)
CORRIGER=True
plt.style.use({'figure.facecolor': 'white', 'savefig.transparent' : False, 'lines.markersize': 6.0})   

#%% Manipulation des fichier !!! copie des trucs

# Repsav = 'CDV_C3_Vert'
# try : os.mkdir(Repsav)
# except : pass

# for element in Liste_corr :
#     head, tail = os.path.split(element)
#     os.system('cp ' + element + ' ' + Repsav +os.sep + tail)    