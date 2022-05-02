#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 14:38:27 2022

@author: josephdeteix
"""

''' librairies et constantes'''

import numpy as np
from math import cos
from math import sin
from math import pi

N = 10
L = 1*(10**(-9))
Kb = 1.38*(10**(-23)) # en J.K-1
T = 293 # K
vmax = 500 # pas trop d'idée ordre de grandeur : à vérifier
M = 18*(10**(-3)) # masse molaire en kg/mol
Na = 6.022*(10**23) 
m = M/Na # masse d'une molécule 
d_AG = 2.57*10**(-11) # 0,1 microm par ex pour une boîte carrée
d_AB = 3.85*10**(-11) # distance pôle - (A) à pôle + (B)
d_BG = d_AB - d_AG

'''les fonctions'''

def matrice_passage(psi,theta,phi):
    mat = np.array([[cos(psi)*cos(phi) - sin(psi)*cos(theta)*sin(phi) , -cos(psi)*sin(phi) - sin(psi)*cos(theta)*cos(phi) , sin(psi)*sin(theta)] , [sin(psi)*cos(phi) + cos(psi)*cos(theta)*sin(phi) , -sin(psi)*sin(phi) + cos(psi)*cos(theta)*cos(phi) , -cos(psi)*sin(theta)] , [sin(theta)*sin(phi) , sin(theta)*cos(phi) , cos(theta)]])
    return mat
    
def positions(N,L):
    ''' renvoie un 10*3*3 des positions pour tous les points et toutes les mol ainsi les angles d euler '''
    positions_possible_G = np.linspace(-L,L,endpoint = False)[1::1]
    positions_G = np.random.choice(positions_possible_G,(N,3),replace = False)
    angles_euler = np.zeros((N,3,2,3))
    angles_euler[:,1,0,0] = abs(np.random.uniform(0,2*pi,size=(N,)))
    angles_euler[:,1,0,1] = abs(np.random.uniform(0,2*pi,size=(N,)))
    angles_euler[:,1,0,2] = abs(np.random.uniform(0,2*pi,size=(N,)))
    angles_euler[:,0,0,0] = angles_euler[:,1,0,0]
    angles_euler[:,0,0,1] = angles_euler[:,1,0,0] 
    angles_euler[:,0,0,2] = angles_euler[:,1,0,2]
    position_B_euler_selon_z = np.full(N, -d_BG)
    position_B_euler = np.zeros((N,3))
    position_B_euler[:,2] = position_B_euler_selon_z
    position_A_euler_selon_z = np.full(N, d_AG)
    position_A_euler = np.zeros((N,3))
    position_A_euler[:,2] = position_A_euler_selon_z
    positions_euler = np.zeros((N,3,1,3))
    positions_euler[:,0,0,:] = position_B_euler
    positions_euler[:,1,0,:] = position_A_euler
    positions_euler[:,2,0,:] = positions_G
    positions_cart = np.zeros((N,3,1,3))
    for i in range (N):
        psi = angles_euler[i,1,0,0]
        theta = angles_euler[i,1,0,1]
        phi = angles_euler[i,1,0,2]
        mat = matrice_passage(psi,theta,phi)
        position_cart_B_bary = mat.dot(positions_euler[i,0,0,:])
        position_cart_A_bary = mat.dot(positions_euler[i,1,0,:])
        positions_cart[i,0,0,:] = position_cart_B_bary + positions_G[i,:]
        positions_cart[i,1,0,:] = position_cart_A_bary + positions_G[i,:]
    positions_cart[:,2,0,:] = positions_G
    return positions_cart , angles_euler
    
def vitesse(N):
    '''donne le tableau (10*3*1*3) des vitesses'''
    sigma = np.sqrt((Kb*T)/m)
    mean = 500/2                # vmax = 500 m/s
    res = np.random.normal(mean, sigma, size = (N,3,1,3))
    return res


def vitesse_euler(N):
    vitesse_min = 40*pi
    vitesse_max = 80*pi
    tab_vitesse = np.random.uniform(vitesse_min,vitesse_max,size = (N,2,3))
    return tab_vitesse
    
def tableau_tot(N,L):
    '''renvoie le tableau 10*3*4*3 de l'initialisation'''
    tab_tot = np.zeros((N,3,4,3))
    positions_cart , angles_euler = positions(N,L)
    vitesse_cart = vitesse(N)
    tab_tot[:,:,0,:] = positions_cart[:,:,0,:]
    tab_tot[:,:,1,:] = vitesse_cart[:,:,0,:]
    tab_tot[:,:,2::1,:] = angles_euler
    tab_tot[:,0:2:1,3,:] = vitesse_euler(N)
    return tab_tot
    
'''programme principale'''

if __name__ == "__main__" :
    tab_tot = tableau_tot(N, L)
    print(tab_tot[0,:,:,:])
