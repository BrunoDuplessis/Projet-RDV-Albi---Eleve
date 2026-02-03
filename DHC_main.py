# -*- coding: utf-8 -*-
"""
Energy modelling of DHC

Authors: Bruno Duplessis, Antoine Fabre
Date: 29/03/24
"""


import numpy as np
import os
import sys
import matplotlib.pyplot as plt

dossier = os.getcwd()
sys.path.append(dossier+'\\Modèles')
import SimulationReseau

"Input files loading"
Fichier_entrees = 'Entrees_Albi_10min_30j.xlsx'  #Température extérieure et demande de chauffage des bâtiments
Donnees_SST = 'DonneesSST_Albi24.xlsx'           #Caractéristiques des sous-stations
Fichier_Tubes = 'DonneesReseau_Albi24.xlsx'      #Caractéristiques des canalisations 

"Main parameters"
T_geo_base = 65                          #Température de l'eau en sortie du doublet géothermal
T_sol = 10                               #Température du sol
A_chauffe = -1.85                        #coefficient A de la loi d'eau du circuit primaire
B_chauffe = 77                           #coefficient B de la loi d'eau du circuit primaire
Tprimaire_min = 60                       #Température minimale au départ du circuit primaire

"Output files"
Fichier_resultats = 'Results.xlsx'

"Core program"
parametres = (T_geo_base,T_sol,A_chauffe, B_chauffe, Tprimaire_min)

(P_geo, P_gaz, T_cons_prod, T_res_ret, Sorties_SST, Input, temps_simu) = SimulationReseau.Run(dossier, Fichier_entrees, Donnees_SST, Fichier_Tubes, parametres, Fichier_resultats)

"Figures"
ax_t = np.arange(6,temps_simu,1)

fig1, ax0 = plt.subplots(figsize=(10,6))
ax0.plot(ax_t,T_cons_prod[6:], label = "Tconsigne départ (°C)")
ax0.plot(ax_t,T_res_ret[6:], label = "T retour (°C)")
ax0.plot(ax_t,Input[6:,1], label ="Text (°C)")
ax0.set_title('Température chaufferie')
ax0.legend()
plt.show()

production = {
    'Géothermie': P_geo[7:],
    'Gaz': P_gaz[7:],
}

fig2,ax = plt.subplots(figsize=(10,6))
ax.stackplot(ax_t,production.values(), labels=production.keys())
ax.legend()

debitprimaire = np.sum(Sorties_SST[2],axis=1)

fig3, axes = plt.subplots(nrows=4, ncols=2, figsize=(10,15))
for i, ax in enumerate(axes.flat):
    ax.plot(ax_t,Sorties_SST[4][6:,i], label = "Ts out")
    ax.plot(ax_t,Sorties_SST[3][6:,i], label = "Tcons")
    ax.set_title('SST_%d' %(i+1))
    ax.legend()
