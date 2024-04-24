#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 12:51:50 2022

@author: theo
"""
import pandas as pd
import numpy as np

ConvertOxEl={'Ions':["Fe", "Fe", "Cu", "Cu", "Mn", "Mn", "Co", "Co", "Cr", "Sn", "Sn", "Ni", "Mo", "Ce", "Ti", "Sb", "Sb", 'S'],
   'Facteurs':[2, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 1, 2, 2, 1]}
ConvertOxEl=pd.DataFrame(ConvertOxEl,
                         index=["Fe2O3", "FeO", "CuO", "Cu2O", "MnO2", "MnO", "CoO",
                                "Co3O4", "Cr2O3", "SnO", "SnO2", "NiO", "MoO2", "CeO2",
                                "TiO2", "Sb2O3", "Sb2O5", 'SO2'])

MasseMol ={
    "SiO2" 	 : 60.0843,
    "B2O3" 	 : 69.6202,
    "Al2O3"  : 101.9613,
    "Li2O" 	 : 29.8814,
    "Na2O" 	 : 61.9789,
    "K2O" 	 : 94.1954,
    "MgO" 	 : 40.3044,
    "CaO" 	 : 56.0774,
    "SrO" 	 : 103.6194,
    "BaO" 	 : 153.3394,
    "ZrO2" 	 : 123.2188,
    "ZnO" 	 : 81.3794,
    "PbO" 	 : 223.1894,
    "F2" 	 : 37.9968,
    "Bi2O3"  : 465.9590,
    "TiO2" 	 : 79.8788,
    "CeO2" 	 : 172.1188,
    "Ce2O3"  : 328.2382,
    "Cr2O3"  : 151.9904,
    "MnO2" 	 : 86.9368,
    "MnO" 	 : 70.9374,
    "Mn2O3"  : 157.8742,
    "Co3O4"  : 240.7972,
    "CoO" 	 : 74.9326,
    "SnO2" 	 : 150.7088,
    "As2O3"  : 197.8414,
    "Sb2O3"  : 291.4982,
    "SO3" 	 : 80.0642,
    "H2O" 	 : 18.0153,
    "Se" 	 : 78.9600,
    "Fe2O3"  : 159.6922,
    "FeO" 	 : 71.8464,
    "P2O5" 	 : 141.9445,
    "La2O3"  : 325.8182,
    "CuO" 	 : 79.5454,
    "Nd2O3"  : 336.4782,
    "V2O5" 	 : 181.8800,
    "Nb2O5"  : 265.8098,
    "WO3" 	 : 231.8482,
    "CdO" 	 : 128.3994,
    "MoO3" 	 : 143.9382,
    "WO3" 	 : 231.8482,
    "Ga2O3"  : 187.4442,
    "Er2O3"  : 382.5182,
    "Gd2O3"  : 362.4982,
    "GeO2" 	 : 104.6088,
    "Ta2O5"  : 441.8928,
    "Y2O3" 	 : 225.8099,
    "In2O3"  : 277.6382,
    "TeO2" 	 : 159.5988,
    "Au2O3"  : 441.9313,
    "Ag2O" 	 : 231.7358,
    "Tl2O" 	 : 424.7660,
    "NiO" 	 : 74.6894,
    "V2O3"   : 149.881,
    "As2O5"  : 229.839,
    "Rb2O"   : 186.935,
    "Sb2O5"  : 323.52,
    "HgO"    : 216.59,
    "H" : 1.007,
    "He" : 4.002,
    "Li" : 6.941,
    "Be" : 9.012,
    "B" : 10.811,
    "C" : 12.011,
    "N" : 14.006,
    "O" : 15.999,
    "F" : 18.998,
    "Ne" : 20.179,
    "Na" : 22.989,
    "Mg" : 24.305,
    "Al" : 26.981,
    "Si" : 28.085,
    "P" : 30.973,
    "S" : 32.066,
    "Cl" : 35.452,
    "Ar" : 39.948,
    "K" : 39.098,
    "Ca" : 40.078,
    "Sc" : 44.955,
    "Ti" : 47.88,
    "V" : 50.941,
    "Cr" : 51.996,
    "Mn" : 54.938,
    "Fe" : 55.847,
    "Co" : 58.933,
    "Ni" : 58.69,
    "Cu" : 63.546,
    "Zn" : 65.39,
    "Ga" : 69.723,
    "Ge" : 72.61,
    "As" : 74.921,
    "Se" : 78.96,
    "Br" : 79.904,
    "Kr" : 83.8,
    "Rb" : 85.467,
    "Sr" : 87.62,
    "Y" : 88.905,
    "Zr" : 91.224,
    "Nb" : 92.906,
    "Mo" : 95.94,
    "Tc" : 98.906,
    "Ru" : 101.07,
    "Rh" : 102.905,
    "Pd" : 106.42,
    "Ag" : 107.868,
    "Cd" : 112.411,
    "In" : 114.82,
    "Sn" : 118.71,
    "Sb" : 121.75,
    "Te" : 127.6,
    "I" : 126.904,
    "Xe" : 131.29,
    "Cs" : 132.905,
    "Ba" : 137.327,
    "La" : 138.905,
    "Ce" : 140.115,
    "Pr" : 140.907,
    "Nd" : 144.24,
    "Pm" : 146.915,
    "Sm" : 150.36,
    "Eu" : 151.965,
    "Gd" : 157.25,
    "Tb" : 158.925,
    "Dy" : 162.5,
    "Ho" : 164.930,
    "Er" : 167.26,
    "Tm" : 168.934,
    "Yb" : 173.04,
    "Lu" : 174.967,
    "Hf" : 178.49,
    "Ta" : 180.947,
    "W" : 183.85,
    "Re" : 186.207,
    "Os" : 190.2,
    "Ir" : 192.22,
    "Pt" : 195.08,
    "Au" : 196.966,
    "Hg" : 200.59,
    "Tl" : 204.383,
    "Pb" : 207.2,
    "Bi" : 208.980,
    "Po" : 208.982,
    "At" : 209.987,
    "Rn" : 222.017,
    "Fr" : 223.019,
    "Ra" : 226.025,
    "Ac" : 227.027,
    "Th" : 232.038,
    "Pa" : 231.035,
    "U" : 238.028,
    "Np" : 237.048,
    "Pu" : 244.064,
    "Am" : 243.061,
    "Cm" : 247.070,
    "Bk" : 247.070,
    "Cf" : 251.079,
    "Es" : 252.082,
    "Fm" : 257.095,
    "Md" : 258.098,
    "No" : 259.100,
    "Lr" : 260.105,
    "Rf" : 261.108,
    "Db" : 262.113,
    "Sg" : 263.118,
    "Bh" : 262.122,
    "Hs" : 265,
    "Mt" : 266,
    "Ds" : 269,
    "Rg" : 272,
    "Cn" : 277,
    "Nh" : 286,
    "Fl" : 289,
    "Mc" : 289,
    "Lv" : 293,
    "Ts" : 294,
    "Og" : 294}
MasseMol = pd.Series(MasseMol)


def Masse2mol(CompoMasse : pd.core.series.Series):
    '''
    Convertie rapport massique en rapport molaire, compo entré et
    sortie sous forme de dataframe

    Parameters
    ----------
    CompoMasse : pd.core.series.Series
        DESCRIPTION.

    Returns
    -------
    CompoMol : TYPE
        DESCRIPTION.

    '''
    global MasseMol
    norm = (CompoMasse/MasseMol[CompoMasse.index]).sum()
    CompoMol = CompoMasse/MasseMol[CompoMasse.index]*1/norm*100
    return CompoMol

def Mol2masse(CompoMol : pd.core.series.Series):
    '''
    Convertie rapport molaire en rapport massique

    Parameters
    ----------
    CompoMol : pd.core.series.Series
        DESCRIPTION.

    Returns
    -------
    CompoMasse : TYPE
        DESCRIPTION.

    '''
    global MasseMol
    norm = (CompoMol*MasseMol[CompoMol.index]).sum()
    CompoMasse = CompoMol*MasseMol[CompoMol.index]*1/norm*100
    return CompoMasse

def Pmass2concentrationMol(Oxyde, Pmass, Densite):
    '''
    Transforme un pourcentage massique en concentration molaire (mol.L⁻1)

    Parameters
    ----------
    Oxyde : TYPE
        DESCRIPTION.
    Pmass : TYPE
        DESCRIPTION.
    Densite : TYPE
        DESCRIPTION.

    Returns
    -------
    Concentration : TYPE
        DESCRIPTION.

    '''
    Concentration = ConvertOxEl['Facteurs'][Oxyde] * Pmass/100 * Densite*1000 * 1/MasseMol[Oxyde]
    return (Concentration)


def Formule2molEl(Formule : pd.core.series.Series) :
    '''
    Cette fonction convertie une formule chimique en pourcentage molaire élémentaire

    Parameters
    ----------
    Formule : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    
    ntot=(Formule).sum() # Calcul nombre de mol tot
    CompoMol=(Formule)/ntot
    return(CompoMol)

def Formule2massEl(Formule : pd.core.series.Series) :
    '''
    Cette fonction convertie une formule chimique en pourcentage massique élémentaire

    Parameters
    ----------
    Formule : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    ntot=(Formule*MasseMol[Formule.index]).sum() # Calcul nombre de mol tot
    CompoMol=(Formule*MasseMol[Formule.index])/ntot
    
    return(CompoMol)

def Formule2MolOx(Formule : pd.core.series.Series):
    '''
    Convertie en %molaire une formule chimique

    Parameters
    ----------
    Formule : pd.core.series.Series
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    ConvertEl2Ox={'Oxyde':["Na2O", "CaO", "SiO2", "Al2O3", "MgO", "K2O", "CuO", "SnO2", "Sb2O3",
                           "PbO", "Y2O3", "Mn2O3", "BaO", "H2O", "Fe2O3", "P2O5"],
                  'Facteurs':[2, 1, 1, 2, 1, 2, 1, 1, 2, 1, 2, 2, 1, 2, 2, 2]}
    ConvertEl2Ox=pd.DataFrame(ConvertEl2Ox,
                         index=["Na", "Ca", "Si", "Al", "Mg", "K ", "Cu", "Sn", "Sb",
                                "Pb", "Y", "Mn", "Ba", "H", "Fe", "P "])

    
    Formule=Formule.drop('O') # On supprime la valeur en oxygène
    noxy = 1/(ConvertEl2Ox['Facteurs'][Formule.index])*Formule
    #On passe en nombre de mol par oxyde.
    
    noxy.index = ConvertEl2Ox['Oxyde'][Formule.index] # Pour avoir les bons noms
    
    CompoMolOxide=noxy/noxy.sum() #Expression en %
    
    return(CompoMolOxide)
 
def ConvertMassOx (MassOxDep):
    """Pour convertir le pourcentage massique d'un oxyde en pourcentage massique d'un autre oxyde,
    sans passer par les pourcentages molaires.
    Ex : FeO to Fe2O3, MnO to Mn2O3"""

    MassOx={'Oxyde':["Fe2O3", "Mn2O3"],
       'Facteurs':[1+MasseMol['O']/(2*MasseMol['FeO']), 1+MasseMol['O']/(2*MasseMol['MnO'])]}
    MassOx=pd.DataFrame(MassOx, index=["FeO", "MnO"]) 

    MassOxFin = MassOx['Facteurs'][MassOxDep.columns]*MassOxDep
    
    return MassOxFin

#%%% Utilisation masse to mol
if False :
    Densite = 2.564
    Pmassique = 0.9901
    Oxyde_ori = 'MnO2'
    
    CompoMasse = {'SiO2': 57, 'Na2O' : 17.10, 'CaO': 20.90,'MnO2' : 5.00}
    
    CompoMol = CompoMasse
    CompoMasse = pd.Series(CompoMasse)
    
    Compovoulu = ['SiO2', 'Na2O', 'MnO']
    
    
    #%%% Utilisation mol to masse
    
    Compo = {'B2O3': 64.67, 'Li2O' : 32.33, 'NiO': 3}
    Compo = pd.Series(Compo)
    
    CompoMol = Compo
    
    CompoMass = Mol2masse(CompoMol)
    
    
    
    #%% Utilisation Formule2mol
    
    Formule = {'Na': 0.1, 'Mn' : 0.1, 'Si': 0.2,'O' : 0.6}
    Formule = {'Na': 0.69, 'Mg' : 0, 'Si': 0.65, 'O' : 56.26, 'Al' : 0.04, 'K' : 0.10, 'Ca' : 16.38, 'Cu' : 0.09, 'Sn' : 0.33, 'Sb' : 25.38, 'Pb' : 0 }
    
    Formule = {'Si' : 4, 'Na': 1.5*2, 'Ca' : 1.5, 'O' : 0}
    
    
    Formule = {'Ba' : 2, 'Mn': 5, 'H' : 4, 'O' : 12}
    
    Formule = {'Si' : 1, 'Fe': 2, 'O' : 4}
    
    Formule = {'Mn' : 2, 'Si': 1, 'O' : 4}
    
    Formule = pd.Series(Formule)
    
    print(Formule)
    CompoMol = Formule2MolOx(Formule)
    print("\nComposition molaire :")
    #print(CompoMol)
    print (round(CompoMol*100, 2))
    
    print("\nComposition massique :")
    print(Mol2masse(CompoMol))
    


