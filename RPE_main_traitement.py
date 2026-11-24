#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 14:14:46 2023

@author: deajais
"""

#%%Chargement des modules et lecture du fichier contenant les paramètre exp
import numpy as np
import os
import warnings
from copy import copy
import tkinter as tk
from tkinter.filedialog import askdirectory, askopenfilename
import fnmatch
from matplotlib import pyplot as plt
from scipy.integrate import quad, dblquad


plt.style.use('default')
#plt.style.use('ggplot')
#plt.rcParams.keys() pour afficher les paramètres modifiable
setfont = lambda color : plt.style.use({'figure.facecolor': color, # yellow, white or none
                                       'savefig.transparent' : False,
                                       'lines.markersize': 6.0})  
close = lambda : plt.close('all')

from Affichage_spectre import Affichage_spectre


from Affichage_spectre import mono2tab

from Affichage_spectre import Sav_fig
from Affichage_spectre import annot_max
from Affichage_spectre import GRAPHOUT_mode
from Affichage_spectre  import circle_plot, rectangle_plot, annotate_clic_text

from Nettoyage_spectre import Nettoyage_spectre_RPE
from Nettoyage_spectre import Integrale

from Lecture_input import readinput_RPE as readinput
from Lecture_input import nm2cm1
from Lecture_input import Chemin2Liste
from Lecture_input import Chemin2input_RPE as Chemin2input
from Lecture_input import Readspectre
from Lecture_input import Releve_parametre

mode='askinput' # input, chemin, create_input, askinput

TITRE='RTvs100K' # Pour trier les spectres à afficher, si input_XXX.csv mettre TITRE=XXX et mode = input
# TITRE = 'KMn0_25comp100K-TA'
# TITRE = 'Vincent'
#TITRE = 'X1_1350C'
CLEFDETRIE='*0,2*'


if mode == 'askinput':
    tkdialog = tk.Tk().withdraw()
    filepath = askopenfilename(initialdir='.')
    path, file = os.path.split(filepath)
    os.chdir(path)
    TITRE = file[6:-4]
    # tkdialog.destroy()
    mode= 'input'


DOSSIER = ['Data_trait']
DOSSIER = ['100K', 'RT']

if mode == 'create_input':
    # Si document dans plusieurs repertoir, indiquer le chemin dans un liste de string.
    Chemin2input(TITRE, CLEFDETRIE=CLEFDETRIE, LISTEDOSSIER=DOSSIER,
                 mode='RPE', correction=31) # RPE ou fityk
    mode='input'
    
if mode == 'input' :
    INPUTNAME='input_'+TITRE+'.csv';
    (Liste, Legende, Liste_ref, Correction, optplt, Frequence,
         valeurnorm, Liste_corr, Limite, TITRE) = readinput(INPUTNAME)                               

CORRIGER=False;
plt.style.use({'figure.facecolor': 'white', 'savefig.transparent' : True, 'lines.markersize': 6.0})   

TITRE_ORI = TITRE

# Liste=Liste[2:4]
# Legende=Legende[2:4]
#%%
%matplotlib auto

#%%
plt.style.use({'figure.facecolor': 'white', 'savefig.transparent' : False, 'lines.markersize': 6.0})   
plt.style.use({'figure.facecolor': 'none', 'savefig.transparent' : True, 'lines.markersize': 6.0})   

CORRIGER=False;

#%% Tracé
# close()
CHEMINSAVFIG = r'C:\Users\theoc\Documents\These\Temporaire\Fe-Mn-Transfert_charge'

# plt.style.use({'legend.loc': 'upper left'})
#plt.style.use({'lines.linestyle' : ':'})

ModeX= 'mT_norm' # g , Champ(G), mT, mT_norm
ModeY='Int' # Int, Integrate, Gradient, DoubleGradient
normmode='' # min_max

modecouleurs = 'auto' # 'auto', 'bigdata', 'manuel'
colormap   = 'nipy' # Utiliser pour mode bigdata nipy ou plasma ou yellow


SHOW        = True; LEGEND = True ; legncol=1
Autoaxe     = True

legendin    = True; #False pour mettre la legende à l'extérieur de la figure
GRAPHOUT    = 5;

Sousmin     = False; #Pour soustraire le minimum
etagage     = False;

SAVE_HTML   = False
CORRIGER    = CORRIGER
CHECKBUTTON = True

X_min = 8;
X_max = 0;
Y_min = -11000;
Y_max = 19000;
Ylim = None
COUPURENORMminmax = [300, 400]


if ModeX == 'g' :
    Autoaxe     = False
    X_min       = 6;
    X_max       = 1;

if CORRIGER :
    Liste_aff=Liste_corr
    RAJOUT = ''
else:
    Liste_aff=Liste
    RAJOUT = '_BRUT'

fig, ax1, TITRE, otherpara = Affichage_spectre(Liste_aff, Legende, Autoaxe, Xlim=[X_min, X_max], Ylim=Ylim,
                              TITRE=TITRE+RAJOUT, ModeX=ModeX, ModeY=ModeY, normmode=normmode,
                              spectromode='RPE', modecouleurs=modecouleurs, optionplot=optplt,
                              colormap=colormap, valeurnorm=valeurnorm, valeurnormX=Frequence,
                              GRAPHOUT=GRAPHOUT, legendin=legendin,
                              ax1=None, newgraph=True, SHOW=SHOW, SAVE_HTML=SAVE_HTML,markevery=1,
                              COUPURENORMminmax=COUPURENORMminmax, 
                              Sousmin=Sousmin,  etagage=etagage, SecondAxe=True,
                              Langue = 'anglais', third_dim=False, CHECKBUTTON=CHECKBUTTON, LimiteAff=Limite)

rect = rectangle_plot(ax1)

if not SHOW and True:
    plt.title('')
    # ax1.set_ylabel('Signal RPE normalisé (u.a)')
    # ax1.set_ylabel('Signal RPE (u.a)')
    
    # ax1.yaxis.set_major_locator(ticker.NullLocator())
    if ModeY== 'Gradient':
        ax1.set_ylabel('Dérivée du signal RPE (u.a)')
    ax1.legend()
    Sav_fig(TITRE, Repertoire = CHEMINSAVFIG, legend=LEGEND, legendin=legendin,
                      legncol=legncol, save_format='.png', PDF=True, fig=fig)
    
TITRE = TITRE_ORI


#%% inset ax
#Affichage_gauss(Liste, Legende, TITRE+'_fitgauss', SHOW=True, optplt='--')
# plt.legend(loc='upper left')
# plt.show()

axins2 = ax1.inset_axes([0.20, 0.7, 0.30, 0.20])
Autoaxe=False
X_min = 25000
X_max = 27000
Y_min = 2
Y_max = 4

Affichage_spectre(Liste_aff, Legende, Autoaxe, Xlim=[X_min, X_max], Ylim=[Y_min,Y_max],
                  TITRE=TITRE+RAJOUT, ModeX=ModeX, ModeY=ModeY, normmode=normmode,
                  spectromode='OAS', modecouleurs=modecouleurs, optionplot=optplt,
                  GRAPHOUT=GRAPHOUT, legendin=legendin,
                  ax1=axins2, newgraph=False, SHOW=False, markevery=1,
                  COUPURENORMminmax=COUPURENORMminmax, 
                  Sousmin=Sousmin,  etagage=etagage, 
                  Langue = 'Francais', valeurnorm=valeurnorm,
                  SecondAxe=False, SpectreVIS = False)
ax.indicate_inset_zoom(axins2, edgecolor="black")
axins2.set_xlabel('')
axins2.set_ylabel('')

#axins2.set_xticklabels([X_min])
#axins2.set_yticklabels([])

Sav_fig(TITRE+'zoom', legendin=legendin)

#%% Calcul de l'intégrale
from Lecture_input import integrate2curve
Index_curvre  = 3# Index curve réduite
Index_curvatm = 2# Index curve atmosphère

COR = False
LimiteG = None#[2000, 8000]

if COR : 
    file1 = Liste_corr[Index_curvre]
    file2 = Liste_corr[Index_curvatm]
else :
    file1 = Liste[Index_curvre]
    file2 = Liste[Index_curvatm]

FIGSIZE, DPI =  GRAPHOUT_mode(5)

fig, axs = plt.subplots(2, 1, sharex=True, figsize=FIGSIZE, dpi = DPI)
X1, Y1 = Readspectre(file1, skip_header=3, delimiter='\t')
X2, Y2 = Readspectre(file2, skip_header=3, delimiter='\t')

Y1 = Y1 / valeurnorm[Index_curvre] # Noramlisation par les paramètre exp
Y2 = Y2 / valeurnorm[Index_curvatm]

if LimiteG is not None :
    INDEX1 = np.logical_and(X1 > np.min(LimiteG), X1 < np.max(LimiteG))
    X1 = X1[INDEX1]
    Y1 = Y1[INDEX1]
    
    INDEX2 = np.logical_and(X2 > np.min(LimiteG), X2 < np.max(LimiteG))
    X2 = X2[INDEX2]
    Y2 = Y2[INDEX2]
    
    
X1, Y1 = integrate2curve(X1, Y1)
nbspin_red = np.trapz(y=Y1, x=X1)
X1, Y1 = integrate2curve(X1, Y1)

X2, Y2 = integrate2curve(X2, Y2)
nbspin_atm = np.trapz(y=Y2, x=X2)
X2, Y2 = integrate2curve(X2, Y2)


axs[0].plot(X1/10,Y1, label = Legende[Index_curvre])
axs[0].plot(X2/10, Y2, label = Legende[Index_curvatm])
axs[0].legend()
axs[1].plot(X2/10, Y2/Y1,  label = Legende[Index_curvatm] + '/ ' +Legende[Index_curvre] )
plt.legend()
axs[1].set_ylim([1, 1.3])
# axs[1].set_ylim([0.8, 0.9])
# axs[1].set_ylim([0.9, 1])
axs[1].grid(True)

axs[1].set_xlabel('Champs magnétique (mT)')
axs[1].set_ylabel('Rapport redox')
axs[0].set_ylabel('$\iint$ signal RPE')

fig.tight_layout(pad=0.01)
print('Rapport redox = ' + str(nbspin_atm/nbspin_red))



#%%% Calcul de l'intégrale DEA

import os
    
x_min = 1000
x_max = 2500
N_spin_ala = 1.75e17
h_al = 3 #en mmm

if CORRIGER :
    Liste_calc=Liste_corr
else:
    Liste_calc=Liste

Liste_Integrale = Integrale(Liste_calc, Legende,x_min, x_max, Gain=Gain,Frequence =Frequence, 
                     Puissance =Puissance, Masse=Masse, Hauteur=Hauteur,Dossier='Data_Integ' ) 

Int_alanine = Integrale(Liste_ref, Legende, x_min=2500, x_max=4500, Gain=Gain, Frequence =Frequence, 
                     Puissance =Puissance, Masse=Masse, Hauteur=Hauteur)
print(Int_alanine)
#Gain_al,Frequence_al,Puissance_al,Amplitude_Mod_al = Releve_parametre(CLEFDETRIE='*alanine*',CHEMIN='.', DOSSIER = 'Data_parametre')
#gain_al,p_al,AM_al = Gain_al[0],Puissance_al[0],Amplitude_Mod_al[0]
#f_al =  1.5929 + 0.019438*h_al - 0.0094083*h_al*h_al+0.00067139-2.1875e-05*h_al**4+1.1815e-07*h_al**5

for i in np.arange(0, np.size(Liste_calc)) :
   # N_spin_ech = Liste_Integrale[i]*(N_spin_ala*gain_al*f_al*AM_al*g_al)/Double_Int_alanine #nb de spin dans 1g d'échantillon
  #probleme  
    try :
        text_file=open('Quantification.csv', "x")
        ENTETE='Nom fichier; Integrale; Nb_spin;\n'
        text_file.write(ENTETE)
    except FileExistsError :
        text_file=open('Quantification.csv', "a")
    
    DATA= Legende + ';' + str(Liste_Integrale[i]) + ';' + str(N_spin_ech) +'\n'
    text_file.write(DATA)
    
    text_file.close()

    

#%% Traitement spectre
# %matplotlib inline
# !! sinon beug dans spyder, a n'utiliser qu'avec spyder

DOSSAVE='./Data_corr'


Liste_corr=Nettoyage_spectre_RPE(Liste, Legende, Liste_ref,Correction, 
                                 Frequence = Frequence, valeurnorm=valeurnorm, Dossier=DOSSAVE)

CORRIGER=True
plt.style.use({'figure.facecolor': 'white', 'savefig.transparent' : False, 'lines.markersize': 6.0}) 

#%% Dosage : hateur pic à pic

# NaMn0.25 : 
hautred = 97.4+109.7
hautatm = 60.5+64.20

redox = hautatm/hautred
print(redox)


#%% 