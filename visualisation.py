#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  2 18:39:32 2022

@author: clarabollaert
"""

#https://stackoverflow.com/questions/41602588/matplotlib-3d-scatter-animations


from matplotlib import pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation

import tableau_total

# constantes
N = 10
L = 1*(10**(-9)) # pour le calcul

nbre_frames = 1000


# Définition des fonctions

def update_graph(num):
    Liste_G_x =[]
    Liste_G_y = []
    Liste_G_z = []
    Liste_G_x=Tableau[num,:,2,0,0]
    Liste_G_y=Tableau[num,:,2,0,1]
    Liste_G_z=Tableau[num,:,2,0,2]
    
    
    Liste_A_x =[]
    Liste_A_y = []
    Liste_A_z = []
    Liste_A_x=Tableau[num,:,1,0,0]
    Liste_A_y=Tableau[num,:,1,0,1]
    Liste_A_z=Tableau[num,:,1,0,2]
    
    Liste_B_x =[]
    Liste_B_y = []
    Liste_B_z = []
    Liste_B_x = Tableau[num,:,0,0,0]
    Liste_B_y = Tableau[num,:,0,0,1]
    Liste_B_z = Tableau[num,:,0,0,2]

    
    
    graphG._offsets3d= (Liste_G_x, Liste_G_y, Liste_G_z)
    graphA._offsets3d= (Liste_A_x, Liste_A_y, Liste_A_z)
    graphB._offsets3d= (Liste_B_x, Liste_B_y, Liste_B_z)
    
    #return(graph)


def augmentation_taille_mol(Tableau, facteur):
    '''Augmente la taille des molécules (norme du vecteur AB) d'un facteur pour la bonne visualisation'''
    x_a, y_a, z_a = Tableau[:,:,1,0,0], Tableau[:,:,1,0,1], Tableau[:,:,1,0,2]
    x_g, y_g, z_g = Tableau[:,:,2,0,0], Tableau[:,:,2,0,1], Tableau[:,:,2,0,2]
    x_b, y_b, z_b = Tableau[:,:,0,0,0], Tableau[:,:,0,0,1], Tableau[:,:,0,0,2]
    
    
    # Aux = np.sqrt((x_a-x_b)**2 + (y_a-y_b)**2 + (z_a-z_b)**2)
    
    # print(Aux) 
    
    # On augmente la taille du vecteur GA
    
    
    x_a= facteur*(x_a-x_g) +x_g
    y_a= facteur*(y_a-y_g) +y_g
    z_a= facteur*(z_a-z_g) +z_g
    Tableau[:,:,1,0,0] = x_a
    Tableau[:,:,1,0,1] = y_a
    Tableau[:,:,1,0,2] = z_a
    
    # Idem pour GB
    
    x_b= facteur*(x_b-x_g) +x_g
    y_b= facteur*(y_b-y_g) +y_g
    z_b= facteur*(z_b-z_g) +z_g
    Tableau[:,:,0,0,0] = x_b
    Tableau[:,:,0,0,1] = y_b
    Tableau[:,:,0,0,2] = z_b
    
    # Aux1 = np.sqrt((x_a-x_b)**2 + (y_a-y_b)**2 + (z_a-z_b)**2)
    # print(Aux1) # pour vérifier la modification de taille
    
    return Tableau



# Programme principal



if __name__ =="__main__":
    # Définition de la figure
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    # Propriétés des axes
    ax.set_xlim3d([-L, L])
    ax.set_xlabel('X')
    ax.set_ylim3d([-L, L])
    ax.set_ylabel('Y')
    ax.set_zlim3d([-L, L])
    ax.set_zlabel('Z')
    ax.set_title('3D Test')



    # on récupère le tableau des résultats
    Tableau = tableau_total.tous_les_tableau(nbre_frames, N, L)

    Tableau = augmentation_taille_mol(Tableau, 10)

    # Initialisation!

    Liste_G_x = Tableau[0,:,2,0,0]
    Liste_G_y = Tableau[0,:,2,0,1]
    Liste_G_z = Tableau[0,:,2,0,2]

    Liste_A_x = Tableau[0,:,1,0,0]
    Liste_A_y = Tableau[0,:,1,0,1]
    Liste_A_z = Tableau[0,:,1,0,2]

    Liste_B_x = Tableau[0,:,0,0,0]
    Liste_B_y = Tableau[0,:,0,0,1]
    Liste_B_z = Tableau[0,:,0,0,2]

    graphG = ax.scatter([],[], color = 'red')
    graphA = ax.scatter([],[], color = 'blue')
    graphB = ax.scatter([], [], color = 'green')


    # Creating the Animation object
    line_ani = animation.FuncAnimation(fig, update_graph, nbre_frames, 
                               interval=200, blit=False)
    plt.show()