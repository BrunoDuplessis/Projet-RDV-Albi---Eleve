# -*- coding: utf-8 -*-
"""
Created on Tue Dec 27 10:00:27 2022

@author: Bruno Duplessis, Antoine Fabre
"""
import numpy as np
from math import exp, pi

"Constantes"
Cp= 4180 # Capacité calorifique eau (kg/m3)
Rho = 1000 # Masse volumique eau (kg/m3)
Nu = 0.000001 # Viscosité cinématique de l'eau

" Modèle de la sous-station avec régulation "
class Noeud:
    def __init__(self,temps_simulation):
        self.State = np.zeros((temps_simulation+1,2)) # Température et débits aux noeuds

class SousStation:
    def __init__(self, ID, df_SST,temps_simulation):
        self.ID = ID
                 
        " Input : Efficacité nominale, Puissance nominale,Pertes de charges nominales, Regime de temperature de l'emetteur "
        self.E_nom = df_SST['E_nominale'][ID-1]
        self.P_nom = df_SST['P_nominale'][ID-1]
        self.Delta_P_nom = df_SST['DeltaP_nom'][ID-1]
        self.Trad_in_nom = df_SST['Trad_in_nom'][ID-1]
        self.Trad_out_nom = df_SST['Trad_out_nom'][ID-1]
        self.Kp = df_SST['Kp_PID'][ID-1]
        self.NodeIn = int(df_SST['Noeud_arrivee'][ID-1])
        self.NodeOut = int(df_SST['Noeud_départ'][ID-1])
        self.Mp_nom = df_SST['Mp_nom'][ID-1]
        self.Ki = df_SST['Ki_PID'][ID-1]
               
        "Calcul des points nominaux"
        self.a = (self.Trad_in_nom-30) / (-7-20) # Calcul de la loi d'eau secondaire en fonction du type de radiateur (Régime Trad_in_nom/Trad_out_nom) température de non chauffage Text=20°C
        self.b = self.Trad_in_nom + 7*self.a # Calcul de la loi d'eau secondaire
        self.NUT_nom = self.E_nom / (1 - self.E_nom) # Calcul du NUT nominal
        self.UA_nom = Cp * self.Mp_nom * self.NUT_nom # Calcul du coefficient de transfert nominal
        self.hAp_nom = 2 * self.UA_nom # Calcul du coefficient de convection côté primaire nominal
        self.hAs_nom = self.hAp_nom # Calcul du coefficient de convection côté secondaire nominal
        self.k = self.Delta_P_nom / self.Mp_nom**2  # Coefficient linaire pertes de charge  
        
        "Etat de la sous-station : puissance fournie, Tconsigne secondaire, Tsortie secondaire, Tsortie primaire, débit primaire, Tentrée secondaire, écart à la consigne, intégrale de l'écart, efficacité, pertes de charge côté primaire, ouverture vanne primaire"
        self.State = np.zeros((temps_simulation+1,11))
        
        "Initialisation de l'état de la sous-station" # puissance fournie, Tconsigne secondaire, Tsortie secondaire, Tsortie primaire, débit primaire, Tentrée secondaire, écart à la consigne, intégrale de l'écart, efficacité, pertes de charge côté primaire, ouverture vanne primaire
        self.State[0] = [0.5*self.P_nom, 50, 50, 40, self.Mp_nom, 30, 0, 0, self.E_nom, self.Delta_P_nom, 0.5]
        
        
    def Run_2 (self, Text, P_appelée, Tp_in, temps,dt): #Inputs : Temperature extérieure, puissance appelée, matrice des noeuds, pas de temps
    
        "récupération ouverture de vanne"
        Youv = self.State[temps,10]
        IntegraleEcartT= self.State[temps,7]
 
        "Bilan échangeur"
        Mp = Youv * self.Mp_nom # Calcul du débit primaire en fonction de l'ouverture de vanne
        Delta_P = self.k * Mp**2 # Calcul pertes de charge
        hAp = self.hAp_nom * (Mp / self.Mp_nom)**(0.8) # Calcul coefficient de convection primaire
        
        UA = (1 / hAp + 1 / self.hAs_nom)**(-1) # Calcul du coefficient de transfert
        NUT = UA /(Cp * Mp) # Calcul du NUT
        R = Mp / self.Mp_nom # Calcul du rapport des débits
        if R==1:
            E = NUT/(1+NUT)
        else :
            E = (1 - exp(-NUT * (1 - R)))/(1- R * exp(-NUT * (1 - R)))  # Calcul de l'efficacité    
        Ts_in = Tp_in - P_appelée/(E*Mp*Cp) # Calcul température entrée côté secondaire
        
        "Contrôle puissance réellement fournie"
        if Ts_in <25.0 :
            Ts_in = 25.0 # on suppose que le pincement des émetteurs de chaleur est de 5K (et une température intérieure de 20 °C)
            P_fournie = E*Mp*Cp*(Tp_in-Ts_in)  #Calcul de la puissance réellement fournie
        else:
            P_fournie = P_appelée  #Calcul de la puissance réellement fournie
    
        Tp_out = Tp_in - P_fournie / (Mp * Cp) # Calcul de la température sortie primaire 
        Ts_out = Ts_in + P_fournie/(self.Mp_nom*Cp) # Calcul température de sortie côté secondaire        
            
        " Régulation "
        Tcons = self.a * Text + self.b # Loi d'eau secaondaire
        EcartT = Tcons-Ts_out
        IntegraleEcartT_next = IntegraleEcartT + EcartT*dt
        
        Youv_next = Youv + self.Kp*EcartT  + self.Ki*IntegraleEcartT  # Calcul de l'incrément d'ouverture de la vanne
        Youv_next = np.clip(Youv_next, 0.001, 1.0)# Calcul de l'ouverture de vanne finale comprise entre 0 et 1
       
        " Output : Etat de la sous-station au pas de temps t et ouverture de vanne au suivant"    
        self.State[temps]=[P_fournie, Tcons, Ts_out, Tp_out, Mp, Ts_in, EcartT, IntegraleEcartT, E, Delta_P, Youv]
        self.State[temps+1,7] = IntegraleEcartT_next
        self.State[temps+1,10] = Youv_next


class TuyauPF:
    " Modèle de canalisation "
    def __init__(self, ID, df_Tubes) :
        self.ID = ID
        self.Lambda_iso = df_Tubes['Lambda_isolant'][ID-1]  
        self.V_lim = df_Tubes['V_lim'][ID-1]  
        self.Long = df_Tubes['Longueur'][ID-1] 
        self.NodeIn=int(df_Tubes['Noeud_entree'][ID-1])
        self.NodeOut=int(df_Tubes['Noeud_sortie'][ID-1])
        self.Mtub_nom = 0     
                
    def Dimensionnement(self,ID,matrice_debits,temps_simulation):  
        self.Mtub_nom = matrice_debits[ID-1]
        self.D_tube = np.round( (4 * self.Mtub_nom / (Rho * self.V_lim *3.14))**(0.5), 2) # Calcul du diamètre des canalisation
        self.e_iso = 30.77 * self.D_tube**3 - 11.43 * self.D_tube**2 + 1.45 * self.D_tube
        self.Delta_p_nom = 8 * self.Lambda_iso * self.Long * self.Mtub_nom**2 /( 3.14**2 * self.D_tube**5 * Rho) # Calcul des pertes de charge nominale
              
        "Etat du tube : température en sortie, débit, pertes de charge "
        self.State = np.zeros((temps_simulation+1,3))
        self.State[0]=[50, self.Mtub_nom, self.Delta_p_nom]            #Initialisation du tube 
        
    def Run(self, T_in, M_tub, temps, Delta_temps,Tsol):
        " Input :  Température d'entrée, Débit, Temps(indice), pas de temps des entrées"
       
        "paramètres du schéma numérique de résolution du modèle"
        delta_t = 1 # pas de temps de calcul en minutes
        N=int(self.Long/10)      #Nombre initial de volumes dans le tuyau
        nb_dt = int(Delta_temps/delta_t)
        dt = 60*delta_t        #conversation pas de temps de calcul en secondes 
        
        "Initialisation du maillage du tuyau"
        L_vol = np.ones(N)*self.Long/N                                                          # longueur de chaque volume 
        T_vol=np.ones(N)*50                                                                     # Température de chaque volume
        M_vol=np.ones(N)*Rho*Cp*pi * self.D_tube**2 / 4 * L_vol                                 # Masse thermique de chaque volume 
        R_vol = np.log(( 2 * self.e_iso + self.D_tube) / self.D_tube)/(2 * pi *self.Lambda_iso * L_vol)   # Résistance thermique de chaque volume

        if M_tub==0:
            M_tub=0.0001
        
        for k in range(nb_dt):
            M_0 = np.array([Cp*M_tub*dt])
            T_0 = np.array([T_in])
            L_0 = np.array([M_tub*dt/(Rho*pi*(self.D_tube/2)**2)])            

            T_vol2 = np.concatenate((T_0,T_vol*(1-dt/(R_vol*M_vol)) + dt/(R_vol*M_vol)*Tsol))
            L_vol2 = np.concatenate((L_0,L_vol))
            M_vol2 = np.concatenate((M_0,M_vol))

            L_vol3 = np.concatenate((L_vol2[np.cumsum(L_vol2)<self.Long],np.array([self.Long-np.sum(L_vol2[np.cumsum(L_vol2)<self.Long])])))
            T_vol3 = np.resize(T_vol2,np.shape(L_vol3))
            M_vol3 = Rho*Cp*pi * self.D_tube**2 / 4 *L_vol3
            M_vol22 = np.copy(M_vol3)
            M_vol22.resize(np.shape(M_vol2),refcheck=False)

            T_vol = T_vol3
            L_vol = L_vol3
            M_vol = M_vol3
            R_vol = np.log(( 2 * self.e_iso + self.D_tube) / self.D_tube)/(2 * pi *self.Lambda_iso * L_vol)
    
        Tout = np.sum((M_vol2-M_vol22)*T_vol2)/np.sum(M_vol2-M_vol22)
        
        "Calcul des pertes de charges"
        Re = M_tub * self.D_tube**2 /( 3.14 * Nu ) # Calcul du nombre de Reynolds
        if Re <2000 : # Calcule du coefficient de pertes thermiques       
            Lambda = 64/Re
        else :     
            Lambda = 0.3164 * Re**(-0.25)

        Delta_p = 8 * Lambda * self.Long * M_tub**2 /( 3.14**2 * self.D_tube**5 * Rho) # Calcul des pertes de charge

        "outputs : température en sortie, débit, pertes de charge "
        self.State[temps+1] = [Tout, M_tub, Delta_p]
        
class Tuyau:
    " classe de tuyau avec modèle Plug Flow"
    def __init__(self, ID, df_Tubes):
        self.ID = ID
        self.Lambda_iso = df_Tubes['Lambda_isolant'][ID-1]  
        self.V_lim = df_Tubes['V_lim'][ID-1]  
        self.L = df_Tubes['Longueur'][ID-1] 
        self.NodeIn=int(df_Tubes['Noeud_entree'][ID-1])
        self.NodeOut=int(df_Tubes['Noeud_sortie'][ID-1])
        self.Mtub_nom = 0
    
    def Dimensionnement(self,ID,matrice_debits,N,Delta_temps,temps_simulation):
        " Output : Diamètre tube, Resistance thermique, Coeff beta "
        self.Mtub_nom = matrice_debits[ID-1]
        self.D_tube = np.round( (4 * self.Mtub_nom / (Rho * self.V_lim *3.14))**(0.5), 2) # Calcul du diamètre des canalisation
        self.e_iso = 30.77 * self.D_tube**3 - 11.43 * self.D_tube**2 + 1.45 * self.D_tube
        self.Vol = 3.14 * self.D_tube**2 / 4 * self.L / N # Volume de chaque tronçon de tuyau
        self.R = np.log(( 2 * self.e_iso + self.D_tube) / self.D_tube)/(2 * 3.14 * self.Lambda_iso * self.L / N) # Calcul de la résistance thermique de la canalisation 
        self.Beta = self.Vol * Rho * Cp /Delta_temps 
        self.Delta_p_nom = 8 * self.Lambda_iso * self.L * self.Mtub_nom**2 /( 3.14**2 * self.D_tube**5 * Rho) # Calcul des pertes de charge nominale
        
        "Etat du tube : débit, pertes de charge "
        self.State = np.zeros((temps_simulation+1,3))
        "Champs de température"
        self.T_champs = np.ones((temps_simulation+1, N))*70 
        
        "Initialisation du tube"
        self.State[0]=[self.T_champs[0,-1], self.Mtub_nom, self.Delta_p_nom]
        
    def Run(self, N, T_in, M_tub, temps,Delta_temps):
        
        " Input : Pas de discretisation, Température d'entrée, Débit, Temps(indice), Température au pas de temps précédent"
        " Output, Température du tuyau, Température au dernier pas de temps, Pertes de charge "             
        alpha = self.Beta + M_tub * Cp + 1/self.R
        Re = M_tub * self.D_tube**2 /( 3.14 * Nu ) # Calcul du nombre de Reynolds

        if Re <2000 : # Calcule du coefficient de pertes thermiques       
            Lambda = 64/Re
        else :     
            Lambda = 0.3164 * Re**(-0.25)

        Delta_p = 8 * Lambda * self.L * M_tub**2 /( 3.14**2 * self.D_tube**5 * Rho) # Calcul des pertes de charge

        T_end = round(600/Delta_temps ) # calcul du pas de temps final
        Ttube=np.zeros((T_end+1,N)) # création de la matrice
        Ttube[0,:]=self.T_champs[temps] # Récupération du champs de température au pas de temps précédent (n)
        Ttube[:,0]=T_in*np.ones((1,T_end+1),float) # On impose T_in en entrée du tuyau sur tout le pas de temps

        alpha = self.Beta + self.State[temps,0] * Cp + 1/self.R
        for j in range(1,T_end+1):
            for i in range(1,N-1):                       
                Ttube[j,i] = (+Tsol/self.R + self.Beta * Ttube[j-1,i] + self.State[temps,0] * Cp * Ttube[j,i-1]) * 1/alpha

        "Champs de température"
        for k in range (N-1):
            self.T_champs[temps+1] = Ttube[T_end,k]

        "Etat du tube : débit, pertes de charge "
        self.State[temps+1] = [self.T_champs[temps+1,-1], M_tub, Delta_p]

def Pompe(Pdc_list, M_tot, eta):
        
    "Input : Pertes de charges de chacun des tronçons, Débit total, Rendement pompe"
    "Output : Puissance pompe"
        
    Pdc_max = np.max(Pdc_list) # Calcul du chemin le plus défavorable
    P_pompe = M_tot/1000 * Pdc_max # Calcul puissance pompe
        
    return P_pompe

" Modèle production de chaleur ( Geothermie + Gaz )"

def Production_chaleur(T_res_ret, E_geo, M_res, eta_gaz, T_cons_prod, T_geo_base, Text):
    " Input : Temperature d'entrée du site de production et de consigne, Débit réseau, Efficacité echangeur, Rendement chaudière "
    " Output : Puissance geothermie, gaz et totale, Température de départ du réseau, température de sortie de l'échangeur geothermique"    
    if T_res_ret > T_cons_prod:
        T_geo_out = T_res_ret  
        T_res_dep = T_res_ret
        P_gaz = 0
        P_geo = 0
    else:
        if T_res_ret > T_geo_base:
            T_geo_out = T_res_ret
        else:
            T_geo_out = E_geo * ( T_geo_base - T_res_ret ) + T_res_ret # Calcul de la température de sortie réseau de l'echangeur geothermie
            if T_geo_out > T_cons_prod:
                T_geo_out = T_cons_prod                 
        P_geo = Cp * M_res * ( T_geo_out - T_res_ret ) # Calcul puissance geothermie
        P_gaz = M_res * Cp * ( T_cons_prod - T_geo_out ) # Calcul puissance chaudière gaz
        
    P_gaz_reel = P_gaz / eta_gaz # calcul puissance gaz réelle
    P_res_tot = P_geo + P_gaz_reel  # calcul puissance totale
    
    return P_geo, P_gaz_reel, P_res_tot, T_cons_prod, T_geo_out

   