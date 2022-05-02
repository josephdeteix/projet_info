#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 15:42:55 2022

@author: josephdeteix
"""

''' librairies et constantes'''

import numpy as np

N = 10
L = 1*(10**(-7))

'''les fonctions'''

def test_mol_rebondit(coord_mol):
    ''' renvoie l indice de la coordonnee qui depasse si il y en a une ainsi que si c est trop grand (True) ou trop petit (False) et None sinon'''
    if coord_mol[2,0,0] > L :
        return 0,True
    elif coord_mol[2,0,0] < -L :
        return 0,False
    elif coord_mol[2,0,1] > L :
        return 1,True
    elif coord_mol[2,0,1] < -L :
        return 1,False
    elif coord_mol[2,0,2] > L :
        return 2,True
    elif coord_mol[2,0,2] < -L :
        return 2,False
    return None

def applique_rebond(coord_mol):
    '''applique le rebond aux molecules qui y sont eligibles'''
    if test_mol_rebondit(coord_mol) != None:
        indice = test_mol_rebondit(coord_mol)[0]
        sens = test_mol_rebondit(coord_mol)[1]
        if sens:
            coord_mol[2,0,indice] = L - (coord_mol[2,0,indice] - L)
            coord_mol[2,1,indice] = -coord_mol[2,1,indice]
        else: 
           coord_mol[2,0,indice] = -L + ( -coord_mol[2,0,indice] - L)
           coord_mol[2,1,indice] = -coord_mol[2,1,indice] 
    return coord_mol

def applique_rebond_tab(coord_tab):
    '''applique le rebond sur toutes les molecules du nouveau tableau'''
    for i in range(coord_tab.shape[0]):
        coord_tab[i,:,:,:] = applique_rebond(coord_tab[i,:,:,:])
    return coord_tab

'''programme principale'''
if __name__ == "__main__":
    taille_boite = 1*(10**(-7))