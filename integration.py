#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 16:32:00 2022

@author: josephdeteix
"""

''' librairies et constantes'''
import calculs_forces
import numpy as np
from math import cos
from math import sin
from math import pi
import rebond
import initialisation_moi

N = 10
L = 1*(10**(-7))
dt = 1e-11#pas du temps d'integrations
d_AG = 2.57*10**(-11) # 0,1 microm par ex pour une boîte carrée
d_AB = 3.85*10**(-11) # distance pôle - (A) à pôle + (B)
M = 18*(10**(-3)) # masse molaire en kg/mol
Na = 6.022*(10**23) 
m = M/Na 
a = d_AG
b = d_AB - d_AG
I1 = 3.025017871798798e-47#celui de l'axe x et y
I2 = 2.6328029130757733e-47#celui de l'axe Z

'''
pour toutes les fonctions incluant uniquement le TMC on est deja dans le referentiel barycentrique
faire gaffe voir relire ce qui entre et ressort des fonctions pour etre sur de pas se tromper avec cart et euler
'''
    
'''les fonctions'''

def separation_force_PFD_TMC(tab_forces):
    '''prends le tableau des forces du tab de shape 10*3*3 et renvoie une tableau des forces pour le PFD (shape 10*1*3) et le TMC (shape 10*2*3)'''
    tab_forces_PFD = tab_forces[:,2,:]
    tab_forces_TMC = tab_forces[:,0:2:1,:]
    return tab_forces_PFD , tab_forces_TMC
       
def forces_en_moment(tab_coord_mol_cart,tab_forces_TMC_mol):
    '''calcul les moments imposé a la molecule a partir des forces et des positions des trois points en entree (2,2,3) et (2,3)'''
    GA = np.array([0,0, a ])
    GB = np.array([0,0, -b ])
    M1 = np.cross(GA,tab_forces_TMC_mol[1,:])
    M2 = np.cross(GB,tab_forces_TMC_mol[0,:])
    return M1 + M2

def matrice_rotation(psi,theta,phi):
    return np.array([[cos(psi)*cos(phi) - sin(psi)*cos(theta)*sin(phi) , -cos(psi)*sin(phi) - sin(psi)*cos(theta)*cos(phi) , sin(psi)*sin(theta)] , [sin(psi)*cos(phi) + cos(psi)*cos(theta)*sin(phi) , -sin(psi)*sin(phi) + cos(psi)*cos(theta)*cos(phi) , -cos(psi)*sin(theta)] , [sin(theta)*sin(phi) , sin(theta)*cos(phi) , cos(theta)]])

def coord_eulerienne_en_cartesienne_barycentrique(tab_angle_euler):
    '''fait juste un changement de coord dans le ref barycentrique avec en entree un 2*2*3 et ressort la meme chose'''
    tab_coord_mol_cart = np.zeros((2,2,3))
    psi = tab_angle_euler[1,0,0]
    theta = tab_angle_euler[1,0,1]
    phi = tab_angle_euler[1,0,2]
    position_A_rm = np.array([0,0,a])
    position_B_rm = np.array([0,0,b])
    tab_coord_mol_cart[0,0,:] = matrice_rotation(psi, theta ,phi).dot(position_B_rm)
    tab_coord_mol_cart[1,0,:] = matrice_rotation(psi, theta ,phi).dot(position_A_rm)
    return tab_coord_mol_cart
    
def mouvement_TMC_une_mol(tab_coord_mol,tab_forces_TMC_mol):
    '''calcule les rotations pour une mol dans les coordonnees psi theta et phi en entree un 3*4*3 et 2*3 sort un 2*4*3'''
    nouveau_tab_coord_mol = np.zeros((2,4,3))
    tab_moment = forces_en_moment(tab_coord_mol[:,0:2:1,:], tab_forces_TMC_mol) #on chope les coord cartesiennes
    m1 = tab_moment[0]
    m2 = tab_moment[1]
    #on attribue chaque valeurs du tableau a sa coordonnee
    psi = tab_coord_mol[1,2,0]
    theta = tab_coord_mol[1,2,1]
    phi = tab_coord_mol[1,2,2]
    psi_point = tab_coord_mol[1,3,0]
    theta_point = tab_coord_mol[1,3,1]
    phi_point = tab_coord_mol[1,3,2]
    # on calcule les d(quelquechose)
    Lz = (psi_point*cos(theta) + phi_point)*I2
    dtheta_point = dt*(m1 - psi_point*sin(theta)*Lz + I1*(psi_point**2)*sin(theta)*cos(theta))/I1
    dpsi = psi*dt
    dtheta = theta*dt
    dphi = phi*dt
    dpsi_point_sin_theta = dt*(m2 - I1*theta_point*psi_point*cos(theta) + Lz*theta_point)/I1
    #on integre
    psi_point_sin_theta = psi_point*sin(theta) + dpsi_point_sin_theta
    theta_point = theta_point + dtheta_point
    theta = theta +dtheta
    psi_point = psi_point_sin_theta / sin(theta)
    psi = psi +dpsi
    phi_point = Lz/I2 - psi_point*cos(theta)
    phi = phi + dphi
    nouveau_tab_coord_mol[:,2::1,:] = np.array([[[psi,theta,phi],[psi_point,theta_point,phi_point]],[[psi,theta,phi],[psi_point,theta_point,phi_point]]])
    nouveau_tab_coord_mol[:,0:2:1,:] = coord_eulerienne_en_cartesienne_barycentrique(nouveau_tab_coord_mol[:,2::1,:])
    return nouveau_tab_coord_mol

def mouvement_TMC_tot(tab_coord_TMC,tab_forces_TMC):
    '''calcule uniquement les rotations de toutes les molecules, en entree un 10*3*4*3 pour les coordonnees et un 10*2*3 pour les forces et sort un 10*2*4*3 pour A et B'''
    nouveau_tab_mouvement_TMC = np.zeros((10,2,4,3))
    for i in range(tab_coord_TMC.shape[0]):
        tab_coord_mol = tab_coord_TMC[i,:,:,:]
        tab_forces_TMC_mol = tab_forces_TMC[i,:,:]
        nouveau_tab_mouvement_TMC_mol = mouvement_TMC_une_mol(tab_coord_mol, tab_forces_TMC_mol)
        nouveau_tab_mouvement_TMC[i,:,:,:] = nouveau_tab_mouvement_TMC_mol
    return nouveau_tab_mouvement_TMC
        
def mouvement_PFD_une_mol(tab_coord_mol,tab_forces_PFD_mol):
    '''calcule le mouvement uniquement a cause du PFD d une seule molecule sur G en entree (4*3) et (3) et sort (1*4*3) du temps d'après'''
    tab_coord_mol_cart = tab_coord_mol[0:2:1,:]
    nouveau_tab_coord_mol = np.zeros((1,4,3)) #on remet un 1 en shape juste pour etre en accord avec la fonctions en dessous
    tab_dV = tab_forces_PFD_mol * dt / m
    tab_dR = tab_coord_mol_cart[1,:] * dt
    nouveau_tab_coord_mol[0,1,:] = tab_coord_mol_cart[1,:] + tab_dV
    nouveau_tab_coord_mol[0,0,:] = tab_coord_mol_cart[0,:] + tab_dR
    
    return nouveau_tab_coord_mol
    
def mouvement_PFD_tot(tab_coord_PFD, tab_forces_PFD):
    '''calcule le mouvement uniquement sur le centre d inertie de toutes les molecules, en entree un 10*1*4*3 pour les coordonnees et un 10*1*3 pourles forces, ressort un 10*1*4*3'''
    nouveau_tab_mouvement_PFD = np.zeros((10,1,4,3))
    for i in range(tab_coord_PFD.shape[0]):
        tab_coord_mol = tab_coord_PFD[i,:,:] 
        tab_forces_PFD_mol = tab_forces_PFD[i,:]
        nouveau_tab_mouvement_PFD_mol = mouvement_PFD_une_mol(tab_coord_mol, tab_forces_PFD_mol)
        nouveau_tab_mouvement_PFD[i,:,:,:] = nouveau_tab_mouvement_PFD_mol
    return nouveau_tab_mouvement_PFD
    
def nouveau_tableau(tab_coord):
    '''la fonctions finale qui donne le nouveau tableau des positions prends un tab 10*3*4*3 et renvoit un tableau de meme shape '''
    nouveau_tab = np.zeros((10,3,4,3))
    tab_forces = calculs_forces.calcul_force(tab_coord)
    tab_forces_PFD , tab_forces_TMC = separation_force_PFD_TMC(tab_forces)
    tab_coord_PFD , tab_coord_TMC = tab_coord[:,2,:,:] , tab_coord
    nouveau_tab_mouvement_PFD_tot = mouvement_PFD_tot(tab_coord_PFD , tab_forces_PFD)
    nouveau_tab_mouvement_TMC_tot = mouvement_TMC_tot(tab_coord_TMC , tab_forces_TMC)
    nouveau_tab_mouvement_AB_complet = nouveau_tab_mouvement_TMC_tot + nouveau_tab_mouvement_PFD_tot #on compose les mouvements
    nouveau_tab[:,0:2:1,:,:] = nouveau_tab_mouvement_AB_complet
    nouveau_tab[:,2,:,:] = nouveau_tab_mouvement_PFD_tot[:,0,:,:]
    nouveau_tab = rebond.applique_rebond_tab(nouveau_tab)
    return nouveau_tab

'''programme principale'''

if __name__ == "__main__" :
    tab_init = initialisation_moi.tableau_tot(N, L)
    tab_deux = nouveau_tableau(tab_init)
    print(tab_init[0,:,:,:])
    print(tab_deux[0,:,:,:])