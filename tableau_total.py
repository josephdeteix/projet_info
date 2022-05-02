#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 00:21:31 2022

@author: josephdeteix
"""

''' librairies et constantes'''

import numpy as np
import integration
import initialisation_moi

nombre_iteration = 1000
N = 10
L = 1*(10**(-7))

'''les fonctions'''

def tous_les_tableau(nombre_iteration,N,L):
    tab_total = np.zeros((nombre_iteration,N,3,4,3))
    tab_total[0,:,:,:,:] = initialisation_moi.tableau_tot(N, L)
    for i in range(1,nombre_iteration):
        tab_total[i,:,:,:,:] = integration.nouveau_tableau(tab_total[i-1,:,:,:,:])
    return tab_total

'''programme principale'''

if __name__ == "__main__" :
    premier = tous_les_tableau(nombre_iteration,N,L)
    coords_G1 = premier[:,0,2,0,:]
    print(coords_G1)
