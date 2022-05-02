#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 12 08:02:40 2022

@author: josephdeteix
"""
"""
toutes les valeurs des constants sont fausses pour l'instant
"""

''' librairies et constantes'''

import numpy as np
from scipy.spatial.distance import euclidean as distance
from math import pi
from numpy.linalg import norm as norm
import itertools
from math import cos
from math import sin
import initialisation_moi

N = 10
L = 1*(10**(-9))
epsilon_0 = 	8.85418782*(10**(-3))
'''charge fictive plus'''
cp= 8.78*(10**(-21))
'''charge_fictive_moins'''
cm = -8.78*(10**(-21))
'''parametre du potentiel de lennard-jones'''
sigma = 2.85*(10**(-10))
e_min = 1.3895715709066758e-21
M = 18*(10**(-3)) # masse molaire en kg/mol
Na = 6.022*(10**23) 
m = M/Na 
'''la base canonique'''
e_x=np.array([1,0,0])
e_y=np.array([0,1,0])
e_z=np.array([0,0,1])
base_cano = np.array([e_x,e_y,e_z])

''' les fonctions'''
def vecteur_unitaire_norme(point_1,point_2):
    '''renvoie la norme et le vecteur unitaire selon la direction de la droite entre les deux points'''
    norme = distance(point_1,point_2)
    return (point_1 - point_2)/norme,norme

def decomposition_vecteur_base(vecteur,base):
    '''renvoie la norme et la décomposition d'un vecteur selon la base choisie'''
    norme = norm(vecteur)
    return norme,np.array([np.vdot(vecteur,base[0]),np.vdot(vecteur,base[1]),np.vdot(vecteur,base[2])])

def r_decomposition(point_1,point_2):
    '''renvoie norme et decomposition unitaire dans la base canonique entre les deux points'''
    vecteur = point_1 - point_2
    norme , vecteur_decompose = decomposition_vecteur_base(vecteur, base_cano)
    return norme,vecteur_decompose/norme

def calcul_force_coulomb(point_1,point_2 ,charge_1,charge_2):
    '''calcul la force de coulomb (sur les trois dim de l espace) entre point_1 et point_2'''
    r,decomposition = r_decomposition(point_1, point_2) 
    norme_force = (charge_1*charge_2)/(4*pi*epsilon_0*(r**2))
    return -decomposition * norme_force #la c est ce qui s'aplique sur le point 1

def calcul_force_lennard_jhones(point_1,point_2):
    '''calcul la force de lennard-jhones (sur les trois dim de l espace) entre deux centre de masses entre point_1 et point_2'''
    r,decomposition = r_decomposition(point_1, point_2)
    norme = -(24*e_min/r)*(2*((sigma/r)**12) - (sigma/r)**6)
    return norme*decomposition # la c est ausi ce qui s'aplique sur le point 1 

def forces_entre_deux_molecules(coord_mol_1,coord_mol_2):
    '''calcul toutes les forces qui s'appliquent entre deux molecules , renvoit deux tableau 3*3 avec en entrée deux 3*2*3'''
    tab_forces_1=np.zeros((3,3))
    tab_forces_2=np.zeros((3,3))
    '''on calcul les forces qui s'appliquent sur chaque points'''
    f_plus_2_sur_plus_1 = calcul_force_coulomb(coord_mol_1[1,0,:] , coord_mol_2[1,0,:],cp,cp) # comme les forces s apliquent sur A
    f_moins_2_sur_plus_1 = calcul_force_coulomb(coord_mol_1[1,0,:] , coord_mol_2[0,0,:],cp,cm)
    f_plus_2_sur_moins_1 = calcul_force_coulomb(coord_mol_1[0,0,:] , coord_mol_2[1,0,:],cm,cp)
    f_moins_2_sur_moins_1 = calcul_force_coulomb(coord_mol_1[0,0,:] ,coord_mol_2[0,0,:],cm,cm)
    f_cdm_sur_cdm = calcul_force_lennard_jhones(coord_mol_1[2,0,:] , coord_mol_2[2,0,:])
    tab_forces_1[1,:] = f_plus_2_sur_plus_1 + f_moins_2_sur_plus_1
    tab_forces_1[0,:] = f_plus_2_sur_moins_1 + f_moins_2_sur_moins_1
    tab_forces_1[2,:] = f_cdm_sur_cdm
    tab_forces_2[1,:] = - f_plus_2_sur_plus_1 - f_plus_2_sur_moins_1
    tab_forces_2[0,:] = - f_moins_2_sur_plus_1 - f_moins_2_sur_moins_1
    tab_forces_2[2,:] = - f_cdm_sur_cdm
    '''maintenant on projette les forces des charges plus selon l axe principale et les plans perpendiculaire a cet axe pour la molecule 1'''
    ''' on a les coordonnées euleriennes donc on choisit la base ex',ey',ez'  et on connait la matrice de passage'''
    psi = coord_mol_1[1,1,0] # la c est un 1 car on a fait attention a envoyer ce qui faut dans coord mol
    theta = coord_mol_1[1,1,1]
    phi = coord_mol_1[1,1,2]
    matrice_passage = np.array([[cos(psi)*cos(phi) - sin(psi)*cos(theta)*sin(phi) , -cos(psi)*sin(phi) - sin(psi)*cos(theta)*cos(phi) , sin(psi)*sin(theta)] , [sin(psi)*cos(phi) + cos(psi)*cos(theta)*sin(phi) , -sin(psi)*sin(phi) + cos(psi)*cos(theta)*cos(phi) , -cos(psi)*sin(theta)] , [sin(theta)*sin(phi) , sin(theta)*cos(phi) , cos(theta)]])
    matrice_passage_inv = np.linalg.inv(matrice_passage)
    tab_force_1_base_euler = matrice_passage_inv.dot(tab_forces_1) #ici revoir comment on fait le changement de base
    '''maintenant on fait les petits changement qui s'imposent'''
    tab_force_1_base_euler[2,0] =+ tab_force_1_base_euler[0,0] + tab_force_1_base_euler[1,0]
    tab_force_1_base_euler[0,0] , tab_force_1_base_euler[1,0] = 0 , 0
    '''puis on remet dans la base canonique pour G et A et B reste dans la base de la molecule'''
    tab_forces_G_remanie = matrice_passage.dot(tab_force_1_base_euler[2,:])
    tab_forces_1[2,:] , tab_forces_1[0:2:1,:] = tab_forces_1[2,:] + tab_forces_G_remanie , tab_force_1_base_euler[0:2:1,:]
    '''on fait la meme chose pour la molecule 2 '''
    psi = coord_mol_2[1,1,0]
    theta = coord_mol_2[1,1,1]
    phi = coord_mol_2[1,1,2]
    matrice_passage = np.array([[cos(psi)*cos(phi) - sin(psi)*cos(theta)*sin(phi) , -cos(psi)*sin(phi) - sin(psi)*cos(theta)*cos(phi) , sin(psi)*sin(theta)] , [sin(psi)*cos(phi) + cos(psi)*cos(theta)*sin(phi) , -sin(psi)*sin(phi) + cos(psi)*cos(theta)*cos(phi) , -cos(psi)*sin(theta)] , [sin(theta)*sin(phi) , sin(theta)*cos(phi) , cos(theta)]])
    matrice_passage_inv = np.linalg.inv(matrice_passage)
    tab_force_2_base_euler = matrice_passage_inv.dot(tab_forces_2)
    '''maintenant on fait les petits changement qui s'imposent'''
    tab_force_2_base_euler[2,0] =+ tab_force_2_base_euler[0,0] + tab_force_2_base_euler[1,0]
    tab_force_2_base_euler[0,0] , tab_force_2_base_euler[1,0] = 0 , 0 
    '''puis on remet dans la base canonique pour G et A et B reste dans la base de la molecule'''   
    tab_forces_G_remanie = matrice_passage.dot(tab_force_2_base_euler[2,:])
    tab_forces_2[2,:] , tab_forces_2[0:2:1,:] = tab_forces_2[2,:] + tab_forces_G_remanie , tab_force_2_base_euler[0:2:1,:]
    return tab_forces_1 , tab_forces_2

def calcul_force(tab):
    '''calcul le tableau des forces qui s exerce sur les trois dimensions de l espace sur les trois points de toutes les molecules renvoie un 10*3*3 avec en entrée un 10*3*4*3'''
    tab_forces = np.zeros((N,3,3))
    for couple in itertools.combinations(range(10),2) : 
        coord_mol_1 , coord_mol_2 = tab[couple[0],:,0:3:2,:] , tab[couple[1],:,0:3:2,:]
        force_mol_1 , force_mol_2 = forces_entre_deux_molecules(coord_mol_1 , coord_mol_2)
        tab_forces[couple[0],:,:] =+ force_mol_1
        tab_forces[couple[1],:,:] =+ force_mol_2
    return tab_forces

    '''programme principale'''

if __name__ == "__main__" :
    tab_init = initialisation_moi.tableau_tot(N, L)
    tab_deux = calcul_force(tab_init)
