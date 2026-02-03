# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 12:37:54 2023

@author: bduplessis
"""
import os
import pandas as pd
import numpy as np

def CreationRDC(dossier, Fichier_entrees, Donnees_SST, Fichier_Tubes):
    "Récupération des données réseau"
    os.chdir(dossier+'\\Entrées')
    dossier_entrees = os.getcwd()
    Fichier_entrees = 'Entrees_Albi23_short.xlsx'
    Donnees_SST = 'DonneesSST_Albi23.xlsx'
    Fichier_Tubes = 'DonneesReseau_Albi23.xlsx'
    os.chdir(dossier)
    
    entrees = pd.read_excel(dossier_entrees+'\\'+Fichier_entrees)
    df_SST=pd.read_excel(dossier_entrees+'\\'+Donnees_SST)
    df_Tubes=pd.read_excel(dossier_entrees+'\\'+Fichier_Tubes)
    
    Input = np.array(entrees) # Import des données d'entrée : Temperature extérieure, appels de puissance
    Infos_SST = np.array(df_SST)
    Carac_tubes = np.array(df_Tubes) #Import des caractéristiques du réseau
    
    Nb_SST = len(df_SST['ID_SST'])
    Nb_tubes = len(df_Tubes['ID_Tube'])
    
    temps_simu = len(Input)
    Branches = list(map(int,df_SST['Boucle'][0].split(",")))
    
    ###################################################Création du réseau"
    "Création des sous-stations"
    SST = [ 0 for _ in range(Nb_SST) ]
    for i in range(Nb_SST):
        SST[i]=DHC_Components.SousStation(i+1,df_SST,temps_simu)
    
    "Création des tubes"
    TUBES = [ 0 for _ in range(Nb_tubes) ]
    Nb_noeuds = 0
    for i in range(Nb_tubes):
        TUBES[i]=DHC_Components.Tuyau(i+1, df_Tubes)
        Nb_noeuds = max(Nb_noeudsoeuds,TUBES[i].NodeOut)
    Nb_noeuds = Nb_noeuds+1
    
    "Création de la matrice des noeuds"
    NOEUDS = [0 for _ in range(Nb_noeuds)]
    for i in range(Nb_noeuds):
        NOEUDS[i]=DHC_Components.Noeud(temps_simu)
        
    "Calcul du débit total du réseau"
    Mp_vector=np.zeros(Nb_SST+2)
    BB = np.zeros(Nb_noeuds) # Vecteur des débits imposés aux noeuds des SST
    for i in range(Nb_SST):
        Mp_vector[i+1] = SST[i].Mp_nom
        BB[SST[i].NodeIn] = SST[i].Mp_nom
        BB[SST[i].NodeOut] =- SST[i].Mp_nom
    Mp_vector[Nb_SST+1] = Mp_vector.sum()
    Mp_vector[0] = Mp_vector[Nb_SST+1]
    BB[0]=-Mp_vector[0]
    BB[-1] = Mp_vector[-1]
    NOEUDS[0].State[0] = (75,Mp_vector[0])
    
    Struct_ReseauT = np.zeros((Nb_tubes+Nb_SST,Nb_noeuds))
    for i in range(Nb_tubes):
        Struct_ReseauT[i,TUBES[i].NodeIn]=-1
        Struct_ReseauT[i,TUBES[i].NodeOut]=1
    for i in range(Nb_SST):
        Struct_ReseauT[Nb_tubes+i,SST[i].NodeIn]=-1
        Struct_ReseauT[Nb_tubes+i,SST[i].NodeOut]=1                       
    
    "Dimensionnement et initialisation des tubes"
    structT= Struct_ReseauT[:Nb_tubes,1:-1]
    Mp_tubes = np.linalg.solve(structT.T,BB[1:-1])  #Calcul des débits nominaux dans les tuyaux
    for i in range(Nb_tubes):
        TUBES[i].Dimensionnement(i+1, Mp_tubes, N, Delta_t, temps_simu)
        
    "Détermination des débits aux noeuds"
    Mp_reseau=np.concatenate((Mp_tubes,Mp_vector[:Nb_SST]))
    AA=np.array(Struct_ReseauT)
    AA[AA<0]=0                              # Matrice de mélange aux noeuds
    Mp_noeuds=abs(np.dot(AA.T,Mp_reseau))
    
    "Initialisation des températures et Débits aux noeuds"
    for k in range(Nb_noeuds-1):
        NOEUDS[k+1].State[0] = (50,Mp_noeuds[k+1])
 