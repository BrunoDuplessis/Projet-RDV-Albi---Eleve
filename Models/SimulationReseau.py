import numpy as np
import pandas as pd
import os
import DHC_Components
from tqdm import tqdm

def Run(dossier, Fichier_entrees, Donnees_SST, Fichier_Tubes, parametres, Fichier_resultats):
    
    "Récupération des données d'entrée"
    os.chdir(dossier+'\\Inputs')
    dossier_entrees = os.getcwd()
    os.chdir(dossier)

    entrees = pd.read_excel(dossier_entrees+'\\'+Fichier_entrees)
    df_SST=pd.read_excel(dossier_entrees+'\\'+Donnees_SST)
    df_Tubes=pd.read_excel(dossier_entrees+'\\'+Fichier_Tubes)
    temps_simu = len(np.array(entrees))
    Input = np.array(entrees) # Import des données d'entrée : Temperature extérieure, appels de puissance
    Infos_SST = np.array(df_SST)
    Carac_tubes = np.array(df_Tubes)
    
    "Paramètres"
    Nb_SST = len(df_SST['ID_SST'])
    Nb_tubes = len(df_Tubes['ID_Tube'])
    Delta_t = 10 # Pas de temps des données d'entrées en minutes
    Branches = list(map(int,df_SST['Boucle'][0].split(",")))
    (T_geo_base,Tsol,A_chauffe, B_chauffe, Tprimaire_min) = parametres
    
    "Constantes"
    Cp= 4180 # Capacité calorifique eau (kg/m3)
    
    "Création du réseau"
    "Sous-stations"
    SST = [ 0 for _ in range(Nb_SST) ]             
    for i in range(Nb_SST):
        SST[i]=DHC_Components.SousStation(i+1,df_SST,temps_simu)
    
    "Canalisations"
    TUBES = [ 0 for _ in range(Nb_tubes) ]          
    Nb_noeuds = 0
    for i in range(Nb_tubes):
        TUBES[i]=DHC_Components.TuyauPF(i+1, df_Tubes)
        Nb_noeuds = max(Nb_noeuds,TUBES[i].NodeOut)
    Nb_noeuds = Nb_noeuds+1
    
    "Noeuds (réseaux)"
    NOEUDS = [0 for _ in range(Nb_noeuds)]          
    for i in range(Nb_noeuds):
        NOEUDS[i]=DHC_Components.Noeud(temps_simu)
    
    "Calcul des grandeurs nominales"
    "Débit nominal"
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
    
    "Dimensionnement des canalisations"
    "Définition de la structure du réseau"
    Struct_ReseauT = np.zeros((Nb_tubes+Nb_SST,Nb_noeuds))
    for i in range(Nb_tubes):
        Struct_ReseauT[i,TUBES[i].NodeIn]=-1
        Struct_ReseauT[i,TUBES[i].NodeOut]=1
    for i in range(Nb_SST):
        Struct_ReseauT[Nb_tubes+i,SST[i].NodeIn]=-1
        Struct_ReseauT[Nb_tubes+i,SST[i].NodeOut]=1
    ##
    structT= Struct_ReseauT[:Nb_tubes,1:-1]
    Mp_tubes = np.linalg.solve(structT.T,BB[1:-1])  #Calcul des débits nominaux dans les tuyaux
    for i in range(Nb_tubes):
        TUBES[i].Dimensionnement(i+1, Mp_tubes,temps_simu)
    
    "Détermination des débits aux noeuds"
    Mp_reseau=np.concatenate((Mp_tubes,Mp_vector[:Nb_SST]))
    AA=np.array(Struct_ReseauT)
    AA[AA<0]=0                              # Matrice de mélange aux noeuds
    Mp_noeuds=abs(np.dot(AA.T,Mp_reseau))

    "Initialisation des températures et Débits aux noeuds"
    for k in range(Nb_noeuds-1):
        NOEUDS[k+1].State[0] = (50,Mp_noeuds[k+1])
    
    "Fonctionnement du réseau"
    "Initialisation"
    P_pompe = [0] # Création vecteur puissance pompe
    P_geo = [0] # Création vecteur puissance geothermie
    P_gaz = [0] # Création vecteur puissance chaudière gaz
    P_tot = [0] # Création vecteur puissance thermique totale
    T_geo_out = [0] # Création vecteur température sortie géothermie

    Tp_in_SST = np.zeros((temps_simu,Nb_SST)) #Initialisation des sorties à afficher
    Tp_out_SST = np.zeros((temps_simu,Nb_SST)) 
    Mp_SST = np.zeros((temps_simu,Nb_SST))
    Tcons_SST = np.zeros((temps_simu,Nb_SST))
    Ts_out_SST = np.zeros((temps_simu,Nb_SST))
    T_cons_prod = [] 
    P_out_SST = np.zeros((temps_simu,Nb_SST))
    Y_ouv_SST = np.ones((temps_simu,Nb_SST))*0.5
    EcartT = np.zeros((temps_simu,Nb_SST))
    T_res_ret = []
    
    for i in range (Nb_SST):
        Tp_in_SST[0,i] =NOEUDS[SST[i].NodeIn].State[0][0]
        Tp_out_SST[0,i]=SST[i].State[0][2]
        Ts_out_SST[0,i]=SST[i].State[0][3]
        Tcons_SST[0,i] = SST[i].a * Input[0][1] + SST[i].b
        Tcons_SST[-1,i] = SST[i].a * Input[-1][1] + SST[i].b
        P_out_SST[0,i] = Cp*SST[i].State[0][4]*(SST[i].State[0][3] - SST[i].State[0][1])
        Mp_SST[0,i] = SST[i].Mp_nom
        
    for n in range(temps_simu):

        T_cons_prod1 = A_chauffe * Input[n,1] + B_chauffe #loi de chauffe primaire (loi d'eau 80/60 (a,b) = ( -1.85, 77) ET min = 60 ; loi d'eau 60/40 (a,b) = (-1.11, 62) et min = 40
        if T_cons_prod1 < Tprimaire_min:
            T_cons_prod1 = Tprimaire_min
        else:
            T_cons_prod1 = A_chauffe * Input[n,1] + B_chauffe
        NOEUDS[0].State[n+1][0]=T_cons_prod1
        
        "Calcul des sous-stations"
        T_SST_out = np.zeros(Nb_SST)# Vecteur des températures de sorties des SST
        BB = np.zeros(Nb_noeuds) # Vecteur des débits imposés aux noeuds des SST
        for i in range(Nb_SST):
            SST[i].Run_2(Input[n,1], Input[n,i+2], NOEUDS[SST[i].NodeIn].State[n][0], n)
            NOEUDS[SST[i].NodeOut].State[n][0]=SST[i].State[n][2]                  #Température sortie primaire SST
            BB[SST[i].NodeIn] = SST[i].State[n][4]                               #Débit coté primaire            
            BB[SST[i].NodeOut] = -SST[i].State[n][4]                              #Débit coté primaire 
            Mp_vector[i+1] = SST[i].State[n][4]                              #Débit coté primaire 
            T_SST_out[i] = SST[i].State[n][2]                                #Température sortie primaire SST
            Tp_in_SST[n,i] = np.round(NOEUDS[SST[i].NodeIn].State[n][0],1)
            Tp_out_SST[n,i] = np.round(SST[i].State[n][2],1)                         #Température sortie primaire SST
            Mp_SST[n,i] = np.round(SST[i].State[n][4],2)                                 #Débit côté primaire SST
            Ts_out_SST[n,i] = np.round(SST[i].State[n][3],1)                            #Température sortie secondaire SST
            Tcons_SST[n,i] = np.round(SST[i].a * Input[n][1] + SST[i].b,1)
            P_out_SST[n,i] = np.round((Cp*SST[i].Mp_nom*(SST[i].State[n][3] - SST[i].State[n][1]))/1000,1)
            Y_ouv_SST[n,i] = np.round(100*SST[i].State[n][0],1)
            EcartT[n,i] = SST[i].State[n][7]
        
        Mp_vector[Nb_SST] = Mp_vector.sum() # calcul du débit total du réseau
        BB[0]=-Mp_vector[Nb_SST]
        BB[-1] = Mp_vector[Nb_SST] 
        NOEUDS[0].State[n][1] = Mp_vector[Nb_SST]
    
        "Calcul des débits dans les tuyaux"
        structT= Struct_ReseauT[:Nb_tubes,1:-1]
        Mp_tubes = np.linalg.solve(structT.T,BB[1:-1]) 
    
        "Calcul des tuyaux"
        T_tube_out = np.zeros(Nb_tubes)
        for i in range(Nb_tubes):
            TUBES[i].Run(NOEUDS[TUBES[i].NodeIn].State[n][0], Mp_tubes[i], n, Delta_t, Tsol)
            NOEUDS[TUBES[i].NodeOut].State[n+1][0]=TUBES[i].State[n+1][0]
            T_tube_out [i] = TUBES[i].State[n+1][0]

        "Détermination des débits aux noeuds"
        Mp_reseau = np.concatenate((Mp_tubes, Mp_vector[1:Nb_SST+1]))
        Mp_noeuds = np.dot(AA.T,Mp_reseau)
    
        "Calcul des températures aux noeuds"
        T_nodes=np.concatenate((T_tube_out,T_SST_out))
        M_Matrix=()
        for k in range(Nb_noeuds-1):
            M_Matrix=np.concatenate((M_Matrix,Mp_reseau))
        M_Matrix=M_Matrix.reshape(Nb_noeuds-1,len(Mp_reseau))
        structure_M=AA[:,1:].T*M_Matrix #résolution température (bilan masses)
        structure_M=(1/np.sum(structure_M,axis=1)*structure_M.T).T
        T_nodes2 = np.dot(structure_M,T_nodes)
    
        for k in range(Nb_noeuds-1):
            NOEUDS[k+1].State[n+1] = (T_nodes2[k],Mp_noeuds[k]) 
        
        " Calcul des pertes de charges sur tout les chemins "
        Delta_P_branches = np.zeros(Nb_SST)
        for k in range(Nb_SST):
            for i in Branches:
                Delta_P_branches[k] = SST[k].State[n+1][6] + TUBES[i-1].State[n+1][2]
            
        " Calcul puissance pompe "
        Pompe1 = DHC_Components.Pompe(Delta_P_branches, Mp_vector[Nb_SST], 0.9)
        P_pompe.append(Pompe1) # Calcul puissance pompe
    
        " Calcul site de production "
    
        Prod_chaleur = DHC_Components.Production_chaleur(NOEUDS[-1].State[n+1][0], 0.95, NOEUDS[-1].State[n+1][1], 0.85, T_cons_prod1,T_geo_base, Input[n,1])
        P_geo.append(np.round(Prod_chaleur[0]/1000,1)) # Calcul puissance geothermie
        P_gaz.append(np.round(Prod_chaleur[1]/1000,1)) # Calcul puissance chaudière gaz
        P_tot.append(np.round(Prod_chaleur[2]/1000,1)) # Calcul puissance totale
        T_cons_prod.append(np.round(Prod_chaleur[3],1)) # Calcul Température de consigne réseau
        T_geo_out.append(np.round(Prod_chaleur[4],1)) # Calcul Température sortie echangeur geothermie
        T_res_ret.append(np.round(NOEUDS[-1].State[n+1][0],1))
    
    "Sorties"
    Part_geo = np.sum(P_geo)
    Part_gaz = np.sum(P_gaz)
    # Inconfort_0=round(np.sum(inconfort[:,0])*10/60,1)
    # Inconfort_1=round(np.sum(inconfort[:,1])*10/60,1)
    # Inconfort_2=round(np.sum(inconfort[:,2])*10/60,1)
    # Inconfort_3=round(np.sum(inconfort[:,3])*10/60,1)
    # Inconfort_4=round(np.sum(inconfort[:,4])*10/60,1)
    # Inconfort_5=round(np.sum(inconfort[:,5])*10/60,1)
    # Inconfort_6=round(np.sum(inconfort[:,6])*10/60,1)
    # Inconfort_7=round(np.sum(inconfort[:,7])*10/60,1)
    
    Sorties_SST = (Tp_in_SST, Tp_out_SST, Mp_SST, Tcons_SST, Ts_out_SST, P_out_SST, Y_ouv_SST, EcartT)
    
    "Calcul des pertes en conduite"                     
    P_all_SST= np.sum(Cp*Sorties_SST[2]*(Sorties_SST[0]-Sorties_SST[1]),axis=1) #Somme des puissances délivrées au primaire des sous-stations
    P_primaire = Cp*np.sum(Sorties_SST[2],axis=1)*(np.array(T_cons_prod)-np.array(T_res_ret)) #Puissance fournie au départ primaire
    pertes = 100*(1-np.sum(P_all_SST)/np.sum(P_primaire))
    
    "Affichage des indicateurs de performance"
    print("La consommation de gaz s'élève à ", np.round(Part_gaz/6,1)," kWh")
    print("Le taux d'énergies renouvelables est de ",np.round(100*Part_geo/(Part_geo+Part_gaz),1)," %")
    print("Les pertes en conduite s'élèvent à ",np.round(pertes,1)," % de la puissance totale délivrée")
    # print("Le taux d'inconfort dans les SST est respectivement de ",Inconfort_0,",",Inconfort_1, ",", Inconfort_2,",",Inconfort_3,",",Inconfort_4,",",Inconfort_5,",",Inconfort_6,",",Inconfort_7,"heures")
        
    "Sauvegarde des résultats"
    os.chdir(dossier+'\\Outputs')
    outputs = pd.DataFrame({"Text (°C)": Input[:,1],"T° départ réseau (°C)": T_cons_prod[:], "T° retour réseau (°C)": T_res_ret[:], " Puissance geothermie (W)" : P_geo[1:], " Puissance gaz (W)" : P_gaz[1:], \
                        "Tp_in_SST1 (°C)" : Sorties_SST[0][:,0], "Tp_out_SST1 (°C)" : Sorties_SST[1][:,0], "Débit primaire_SST1 (kg/s)" : Sorties_SST[2][:,0], "Température de consigne_SST1 (°C)":Sorties_SST[3][:,0], "Ts_out_SST1 (°C)" : Sorties_SST[4][:,0], "Puissance_SST1 (kW)":Sorties_SST[5][:,0],"Ouverture_vanne_SST1 (%)":Sorties_SST[6][:,0], "Indicateur d'inconfort potentiel":Sorties_SST[7][:,0],\
                        "Tp_in_SST2 (°C)" : Sorties_SST[0][:,1], "Tp_out_SST2 (°C)" : Sorties_SST[1][:,1], "Débit primaire_SST2 (kg/s)" : Sorties_SST[2][:,1], "Température de consigne_SST2 (°C)":Sorties_SST[3][:,1], "Ts_out_SST2 (°C)" : Sorties_SST[4][:,1], "Puissance_SST2 (kW)":Sorties_SST[5][:,1],"Ouverture_vanne_SST2 (%)":Sorties_SST[6][:,1],"Indicateur d'inconfort potentiel":Sorties_SST[7][:,1],\
                        "Tp_in_SST3 (°C)" : Sorties_SST[0][:,2], "Tp_out_SST3 (°C)" : Sorties_SST[1][:,2], "Débit primaire_SST3 (kg/s)" : Sorties_SST[2][:,2], "Température de consigne_SST3 (°C)":Sorties_SST[3][:,2], "Ts_out_SST3 (°C)" : Sorties_SST[4][:,2], "Puissance_SST3 (kW)":Sorties_SST[5][:,2],"Ouverture_vanne_SST3 (%)":Sorties_SST[6][:,2],"Indicateur d'inconfort potentiel":Sorties_SST[7][:,2],\
                        "Tp_in_SST4 (°C)" : Sorties_SST[0][:,3], "Tp_out_SST4 (°C)" : Sorties_SST[1][:,3], "Débit primaire_SST4 (kg/s)" : Sorties_SST[2][:,3], "Température de consigne_SST4 (°C)":Sorties_SST[3][:,3], "Ts_out_SST4 (°C)" : Sorties_SST[4][:,3], "Puissance_SST4 (kW)":Sorties_SST[5][:,3],"Ouverture_vanne_SST4 (%)":Sorties_SST[6][:,3],"Indicateur d'inconfort potentiel":Sorties_SST[7][:,3],\
                        "Tp_in_SST5 (°C)" : Sorties_SST[0][:,4], "Tp_out_SST5 (°C)" : Sorties_SST[1][:,4], "Débit primaire_SST5 (kg/s)" : Sorties_SST[2][:,4], "Température de consigne_SST5 (°C)":Sorties_SST[3][:,4], "Ts_out_SST5 (°C)" : Sorties_SST[4][:,4], "Puissance_SST5 (kW)":Sorties_SST[5][:,4],"Ouverture_vanne_SST5 (%)":Sorties_SST[6][:,4],"Indicateur d'inconfort potentiel":Sorties_SST[7][:,4],\
                        "Tp_in_SST6 (°C)" : Sorties_SST[0][:,5], "Tp_out_SST6 (°C)" : Sorties_SST[1][:,5], "Débit primaire_SST6 (kg/s)" : Sorties_SST[2][:,5], "Température de consigne_SST6 (°C)":Sorties_SST[3][:,5], "Ts_out_SST6 (°C)" : Sorties_SST[4][:,5], "Puissance_SST6 (kW)":Sorties_SST[5][:,5],"Ouverture_vanne_SST6 (%)":Sorties_SST[6][:,5],"Indicateur d'inconfort potentiel":Sorties_SST[7][:,5],\
                        "Tp_in_SST7 (°C)" : Sorties_SST[0][:,6], "Tp_out_SST7 (°C)" : Sorties_SST[1][:,6], "Débit primaire_SST7 (kg/s)" : Sorties_SST[2][:,6], "Température de consigne_SST7 (°C)":Sorties_SST[3][:,6], "Ts_out_SST7 (°C)" : Sorties_SST[4][:,6], "Puissance_SST7 (kW)":Sorties_SST[5][:,6],"Ouverture_vanne_SST7 (%)":Sorties_SST[6][:,6],"Indicateur d'inconfort potentiel":Sorties_SST[7][:,6],\
                        "Tp_in_SST8 (°C)" : Sorties_SST[0][:,7], "Tp_out_SST8 (°C)" : Sorties_SST[1][:,7], "Débit primaire_SST8 (kg/s)" : Sorties_SST[2][:,7], "Température de consigne_SST8 (°C)":Sorties_SST[3][:,7], "Ts_out_SST8 (°C)" : Sorties_SST[4][:,7], "Puissance_SST8 (kW)":Sorties_SST[5][:,7],"Ouverture_vanne_SST8 (%)":Sorties_SST[6][:,7],"Indicateur d'inconfort potentiel":Sorties_SST[7][:,7],})

    outputs.to_excel(Fichier_resultats)
    os.chdir(dossier)
    
    return P_geo, P_gaz, T_cons_prod, T_res_ret, Sorties_SST, Input, temps_simu


def Run_discret(dossier, Fichier_entrees, Donnees_SST, Fichier_Tubes, parametres, Fichier_resultats):
    
    "Récupération des données d'entrée"
    print("Lecture des données d'entrée")
    os.chdir(dossier+'\\Inputs')
    dossier_entrees = os.getcwd()
    os.chdir(dossier)

    entrees = pd.read_excel(dossier_entrees+'\\'+Fichier_entrees)
    df_SST=pd.read_excel(dossier_entrees+'\\'+Donnees_SST)
    df_Tubes=pd.read_excel(dossier_entrees+'\\'+Fichier_Tubes)
    temps_simu = len(np.array(entrees))
    Input = np.array(entrees) # Import des données d'entrée : Temperature extérieure, appels de puissance
    Infos_SST = np.array(df_SST)
    Carac_tubes = np.array(df_Tubes)

    "Paramètres"
    Nb_SST = len(df_SST['ID_SST'])
    Nb_tubes = len(df_Tubes['ID_Tube'])
    Delta_t = 10 # Pas de temps des données d'entrées en minutes
    Branches = list(map(int,df_SST['Boucle'][0].split(",")))
    (T_geo_base,Tsol,A_chauffe, B_chauffe, Tprimaire_min) = parametres
    
    "Paramètres discretisation temporelle"
    dt = 0.2 # facteur de discretisation temporel
    temps_simu_dt = int(temps_simu/dt)

    "Constantes"
    Cp= 4180 # Capacité calorifique eau (kg/m3)
    
    "Création du réseau"
    print("Création du réseau")
    "Sous-stations"
    SST = [ 0 for _ in range(Nb_SST) ]             
    for i in range(Nb_SST):
        SST[i]=DHC_Components.SousStation(i+1,df_SST,temps_simu_dt)
    
    "Canalisations"
    TUBES = [ 0 for _ in range(Nb_tubes) ]          
    Nb_noeuds = 0
    for i in range(Nb_tubes):
        TUBES[i]=DHC_Components.TuyauPF(i+1, df_Tubes)
        Nb_noeuds = max(Nb_noeuds,TUBES[i].NodeOut)
    Nb_noeuds = Nb_noeuds+1
    
    "Noeuds (réseaux)"
    NOEUDS = [0 for _ in range(Nb_noeuds)]          
    for i in range(Nb_noeuds):
        NOEUDS[i]=DHC_Components.Noeud(temps_simu_dt)
    
    "Calcul des grandeurs nominales"
    "Débit nominal"
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
    
    "Dimensionnement des canalisations"
    "Définition de la structure du réseau"
    Struct_ReseauT = np.zeros((Nb_tubes+Nb_SST,Nb_noeuds))
    for i in range(Nb_tubes):
        Struct_ReseauT[i,TUBES[i].NodeIn]=-1
        Struct_ReseauT[i,TUBES[i].NodeOut]=1
    for i in range(Nb_SST):
        Struct_ReseauT[Nb_tubes+i,SST[i].NodeIn]=-1
        Struct_ReseauT[Nb_tubes+i,SST[i].NodeOut]=1
    
    structT= Struct_ReseauT[:Nb_tubes,1:-1]
    Mp_tubes = np.linalg.solve(structT.T,BB[1:-1])  #Calcul des débits nominaux dans les tuyaux
    for i in range(Nb_tubes):
        TUBES[i].Dimensionnement(i+1, Mp_tubes,temps_simu_dt)
    
    "Détermination des débits aux noeuds"
    Mp_reseau=np.concatenate((Mp_tubes,Mp_vector[:Nb_SST]))
    AA=np.array(Struct_ReseauT)
    AA[AA<0]=0                              # Matrice de mélange aux noeuds
    Mp_noeuds=abs(np.dot(AA.T,Mp_reseau))

    "Initialisation des températures et Débits aux noeuds"
    print("Initialisation")
    for k in range(Nb_noeuds-1):
        NOEUDS[k+1].State[0] = (50,Mp_noeuds[k+1])
    
    "Fonctionnement du réseau"
    "Initialisation"
    P_pompe = [0] # Création vecteur puissance pompe
    P_geo = [0] # Création vecteur puissance geothermie
    P_gaz = [0] # Création vecteur puissance chaudière gaz
    P_tot = [0] # Création vecteur puissance thermique totale
    T_geo_out = [0] # Création vecteur température sortie géothermie

    Tp_in_SST = np.zeros((temps_simu_dt,Nb_SST)) #Initialisation des variables de calcul
    Tp_out_SST = np.zeros((temps_simu_dt,Nb_SST)) 
    Mp_SST = np.zeros((temps_simu_dt,Nb_SST))
    Tcons_SST = np.zeros((temps_simu_dt,Nb_SST))
    Ts_out_SST = np.zeros((temps_simu_dt,Nb_SST))
    T_cons_prod = [] 
    P_out_SST = np.zeros((temps_simu_dt,Nb_SST))
    P_appelee = np.zeros((temps_simu_dt,Nb_SST))
    Y_ouv_SST = np.ones((temps_simu_dt,Nb_SST))*0.5
    EcartT = np.zeros((temps_simu_dt,Nb_SST))
    T_res_ret = []
    P_non_fournie = np.zeros((temps_simu_dt,Nb_SST))

    Text = np.interp((np.arange(0,temps_simu,dt)),np.linspace(0,temps_simu-1,temps_simu),Input[:,1])   

    for i in range (Nb_SST):
        Tp_in_SST[0,i] =NOEUDS[SST[i].NodeIn].State[0][0]
        Tp_out_SST[0,i]=SST[i].State[0][2]
        Ts_out_SST[0,i]=SST[i].State[0][3]
        Tcons_SST[0,i] = SST[i].a * Input[0][1] + SST[i].b
        Tcons_SST[-1,i] = SST[i].a * Input[-1][1] + SST[i].b
        P_out_SST[0,i] = Cp*SST[i].State[0][4]*(SST[i].State[0][3] - SST[i].State[0][1])
        Mp_SST[0,i] = SST[i].Mp_nom
        P_appelee[:,i] = np.interp((np.arange(0,temps_simu,dt)),np.linspace(0,temps_simu-1,temps_simu),Input[:,i+2])

    print("Simulation en cours")    
    for n in tqdm(range(temps_simu_dt)):
        T_cons_prod1 = A_chauffe * Text[n] + B_chauffe #loi de chauffe primaire (loi d'eau 80/60 (a,b) = ( -1.85, 77) ET min = 60 ; loi d'eau 60/40 (a,b) = (-1.11, 62) et min = 40
        if T_cons_prod1 < Tprimaire_min:
            T_cons_prod1 = Tprimaire_min
        else:
            T_cons_prod1 = A_chauffe * Text[n] + B_chauffe
        NOEUDS[0].State[n+1][0]=T_cons_prod1
        
        "Calcul des sous-stations"
        T_SST_out = np.zeros(Nb_SST)# Vecteur des températures de sorties des SST
        BB = np.zeros(Nb_noeuds) # Vecteur des débits imposés aux noeuds des SST
        for i in range(Nb_SST):
            SST[i].Run_2(Text[n], P_appelee[n,i], NOEUDS[SST[i].NodeIn].State[n][0], n,dt)
            NOEUDS[SST[i].NodeOut].State[n][0]=SST[i].State[n][3]                  #Température sortie primaire SST
            BB[SST[i].NodeIn] = SST[i].State[n][4]                               #Débit coté primaire            
            BB[SST[i].NodeOut] = -SST[i].State[n][4]                              #Débit coté primaire 
            Mp_vector[i+1] = SST[i].State[n][4]                              #Débit coté primaire 
            T_SST_out[i] = SST[i].State[n][3]                                #Température sortie primaire SST
            Tp_in_SST[n,i] = np.round(NOEUDS[SST[i].NodeIn].State[n][0],decimals=1)
            Tp_out_SST[n,i] = np.round(SST[i].State[n][3],decimals=1)                         #Température sortie primaire SST
            Mp_SST[n,i] = np.round(SST[i].State[n][4],decimals=2)                                 #Débit côté primaire SST
            Ts_out_SST[n,i] = np.round(SST[i].State[n][2],decimals=1)                            #Température sortie secondaire SST
            Tcons_SST[n,i] = np.round(SST[i].State[n][1],decimals=1)
            P_out_SST[n,i] = np.round(SST[i].State[n][0],decimals=11)
            Y_ouv_SST[n,i] = np.round(SST[i].State[n][10],decimals=3)
            EcartT[n,i] = np.round(SST[i].State[n][6],decimals=3)
            P_non_fournie[n,i] = np.clip(0, P_appelee[n,i], P_appelee[n,i] - P_out_SST[n,i])

        Mp_vector[Nb_SST] = Mp_vector.sum() # calcul du débit total du réseau
        BB[0]=-Mp_vector[Nb_SST]
        BB[-1] = Mp_vector[Nb_SST] 
        NOEUDS[0].State[n][1] = Mp_vector[Nb_SST]
    
        "Calcul des débits dans les tuyaux"
        structT= Struct_ReseauT[:Nb_tubes,1:-1]
        Mp_tubes = np.linalg.solve(structT.T,BB[1:-1]) 
    
        "Calcul des tuyaux"
        T_tube_out = np.zeros(Nb_tubes)
        for i in range(Nb_tubes):
            TUBES[i].Run(NOEUDS[TUBES[i].NodeIn].State[n][0], Mp_tubes[i], n, Delta_t, Tsol)
            NOEUDS[TUBES[i].NodeOut].State[n+1][0]=TUBES[i].State[n+1][0]
            T_tube_out [i] = TUBES[i].State[n+1][0]

        "Détermination des débits aux noeuds"
        Mp_reseau = np.concatenate((Mp_tubes, Mp_vector[1:Nb_SST+1]))
        Mp_noeuds = np.dot(AA.T,Mp_reseau)
    
        "Calcul des températures aux noeuds"
        T_nodes=np.concatenate((T_tube_out,T_SST_out))
        M_Matrix=()
        for k in range(Nb_noeuds-1):
            M_Matrix=np.concatenate((M_Matrix,Mp_reseau))
        M_Matrix=M_Matrix.reshape(Nb_noeuds-1,len(Mp_reseau))
        structure_M=AA[:,1:].T*M_Matrix #résolution température (bilan masses)
        structure_M=(1/np.sum(structure_M,axis=1)*structure_M.T).T
        T_nodes2 = np.dot(structure_M,T_nodes)
    
        for k in range(Nb_noeuds-1):
            NOEUDS[k+1].State[n+1] = (T_nodes2[k],Mp_noeuds[k]) 
        
        " Calcul des pertes de charges sur tout les chemins "
        Delta_P_branches = np.zeros(Nb_SST)
        for k in range(Nb_SST):
            for i in Branches:
                Delta_P_branches[k] = SST[k].State[n+1][6] + TUBES[i-1].State[n+1][2]
            
        " Calcul puissance pompe "
        Pompe1 = DHC_Components.Pompe(Delta_P_branches, Mp_vector[Nb_SST], 0.9)
        P_pompe.append(Pompe1) # Calcul puissance pompe
    
        " Calcul site de production "
    
        Prod_chaleur = DHC_Components.Production_chaleur(NOEUDS[-1].State[n+1][0], 0.95, NOEUDS[-1].State[n+1][1], 0.85, T_cons_prod1,T_geo_base, Text[n])
        P_geo.append(np.round(Prod_chaleur[0]/1000,1)) # Calcul puissance geothermie
        P_gaz.append(np.round(Prod_chaleur[1]/1000,1)) # Calcul puissance chaudière gaz
        P_tot.append(np.round(Prod_chaleur[2]/1000,1)) # Calcul puissance totale
        T_cons_prod.append(np.round(Prod_chaleur[3],1)) # Calcul Température de consigne réseau
        T_geo_out.append(np.round(Prod_chaleur[4],1)) # Calcul Température sortie echangeur geothermie
        T_res_ret.append(np.round(NOEUDS[-1].State[n+1][0],1))
  
    "Sorties"
    print("Post-processing et sauvegarde")
    Part_geo = np.sum(P_geo)
    Part_gaz = np.sum(P_gaz)
    Tp_in_SST_red = np.zeros((temps_simu,Nb_SST))
    Tp_out_SST_red =np.zeros((temps_simu,Nb_SST))
    Mp_SST_red = np.zeros((temps_simu,Nb_SST))
    Tcons_SST_red = np.zeros((temps_simu,Nb_SST))
    Ts_out_SST_red = np.zeros((temps_simu,Nb_SST))
    P_out_SST_red = np.zeros((temps_simu,Nb_SST))
    Y_ouv_SST_red = np.zeros((temps_simu,Nb_SST))
    EcartT_red = np.zeros((temps_simu,Nb_SST))
    E_non_fournie = np.zeros((Nb_SST))

    for i in range(Nb_SST):
        Tp_in_SST_red[:,i] = np.interp(np.linspace(0, temps_simu_dt,temps_simu), np.arange(0, temps_simu_dt,1),Tp_in_SST[:,i])
        Tp_out_SST_red[:,i] = np.interp(np.linspace(0, temps_simu_dt,temps_simu), np.arange(0, temps_simu_dt,1),Tp_out_SST[:,i])
        Mp_SST_red[:,i] = np.interp(np.linspace(0, temps_simu_dt,temps_simu), np.arange(0, temps_simu_dt,1),Mp_SST[:,i])
        Tcons_SST_red[:,i] = np.interp(np.linspace(0, temps_simu_dt,temps_simu), np.arange(0, temps_simu_dt,1),Tcons_SST[:,i])
        Ts_out_SST_red[:,i] = np.interp(np.linspace(0, temps_simu_dt,temps_simu), np.arange(0, temps_simu_dt,1),Ts_out_SST[:,i])
        P_out_SST_red[:,i] = np.interp(np.linspace(0, temps_simu_dt,temps_simu), np.arange(0, temps_simu_dt,1),P_out_SST[:,i])
        Y_ouv_SST_red[:,i] = np.interp(np.linspace(0, temps_simu_dt,temps_simu), np.arange(0, temps_simu_dt,1),Y_ouv_SST[:,i])
        EcartT_red[:,i]= np.interp(np.linspace(0, temps_simu_dt,temps_simu), np.arange(0, temps_simu_dt,1),EcartT[:,i])
        E_non_fournie[i] = np.round(np.sum(P_non_fournie[:,i])*60/(10*dt),1)

    Sorties_SST = (Tp_in_SST_red, Tp_out_SST_red, Mp_SST_red, Tcons_SST_red, Ts_out_SST_red, P_out_SST_red, Y_ouv_SST_red, EcartT_red, )
    T_cons_prod_red = np.interp(np.linspace(0, temps_simu_dt,temps_simu), np.arange(0, temps_simu_dt,1),T_cons_prod)
    P_geo_red = np.interp(np.linspace(0, temps_simu_dt,temps_simu), np.arange(0, temps_simu_dt+1,1),P_geo)
    P_gaz_red = np.interp(np.linspace(0, temps_simu_dt,temps_simu), np.arange(0, temps_simu_dt+1,1),P_gaz)
    T_res_ret_red = np.interp(np.linspace(0, temps_simu_dt,temps_simu), np.arange(0, temps_simu_dt,1),T_res_ret)
    
    

    "Calcul des pertes en conduite"                     
    P_all_SST= np.sum(Cp*Sorties_SST[2]*(Sorties_SST[0]-Sorties_SST[1]),axis=1) #Somme des puissances délivrées au primaire des sous-stations
    P_primaire = Cp*np.sum(Sorties_SST[2],axis=1)*(np.array(T_cons_prod_red)-np.array(T_res_ret_red)) #Puissance fournie au départ primaire
    pertes = 100*(1-np.sum(P_all_SST)/np.sum(P_primaire))
    
    

    print("Sauvegarde des résultats")

    "Sauvegarde des résultats"
    os.chdir(dossier+'\\Outputs')
    outputs = pd.DataFrame({"Text (°C)": Input[:,1],"T° départ réseau (°C)": T_cons_prod_red[:], "T° retour réseau (°C)": T_res_ret_red[:], " Puissance geothermie (W)" : P_geo_red[:], " Puissance gaz (W)" : P_gaz_red[:], \
                        "Tp_in_SST1 (°C)" : Sorties_SST[0][:,0], "Tp_out_SST1 (°C)" : Sorties_SST[1][:,0], "Débit primaire_SST1 (kg/s)" : Sorties_SST[2][:,0], "Température de consigne_SST1 (°C)":Sorties_SST[3][:,0], "Ts_out_SST1 (°C)" : Sorties_SST[4][:,0], "Puissance_SST1 (kW)":Sorties_SST[5][:,0],"Ouverture_vanne_SST1 (%)":Sorties_SST[6][:,0], "Indicateur d'inconfort potentiel":Sorties_SST[7][:,0],\
                        "Tp_in_SST2 (°C)" : Sorties_SST[0][:,1], "Tp_out_SST2 (°C)" : Sorties_SST[1][:,1], "Débit primaire_SST2 (kg/s)" : Sorties_SST[2][:,1], "Température de consigne_SST2 (°C)":Sorties_SST[3][:,1], "Ts_out_SST2 (°C)" : Sorties_SST[4][:,1], "Puissance_SST2 (kW)":Sorties_SST[5][:,1],"Ouverture_vanne_SST2 (%)":Sorties_SST[6][:,1],"Indicateur d'inconfort potentiel":Sorties_SST[7][:,1],\
                        "Tp_in_SST3 (°C)" : Sorties_SST[0][:,2], "Tp_out_SST3 (°C)" : Sorties_SST[1][:,2], "Débit primaire_SST3 (kg/s)" : Sorties_SST[2][:,2], "Température de consigne_SST3 (°C)":Sorties_SST[3][:,2], "Ts_out_SST3 (°C)" : Sorties_SST[4][:,2], "Puissance_SST3 (kW)":Sorties_SST[5][:,2],"Ouverture_vanne_SST3 (%)":Sorties_SST[6][:,2],"Indicateur d'inconfort potentiel":Sorties_SST[7][:,2],\
                        "Tp_in_SST4 (°C)" : Sorties_SST[0][:,3], "Tp_out_SST4 (°C)" : Sorties_SST[1][:,3], "Débit primaire_SST4 (kg/s)" : Sorties_SST[2][:,3], "Température de consigne_SST4 (°C)":Sorties_SST[3][:,3], "Ts_out_SST4 (°C)" : Sorties_SST[4][:,3], "Puissance_SST4 (kW)":Sorties_SST[5][:,3],"Ouverture_vanne_SST4 (%)":Sorties_SST[6][:,3],"Indicateur d'inconfort potentiel":Sorties_SST[7][:,3],\
                        "Tp_in_SST5 (°C)" : Sorties_SST[0][:,4], "Tp_out_SST5 (°C)" : Sorties_SST[1][:,4], "Débit primaire_SST5 (kg/s)" : Sorties_SST[2][:,4], "Température de consigne_SST5 (°C)":Sorties_SST[3][:,4], "Ts_out_SST5 (°C)" : Sorties_SST[4][:,4], "Puissance_SST5 (kW)":Sorties_SST[5][:,4],"Ouverture_vanne_SST5 (%)":Sorties_SST[6][:,4],"Indicateur d'inconfort potentiel":Sorties_SST[7][:,4],\
                        "Tp_in_SST6 (°C)" : Sorties_SST[0][:,5], "Tp_out_SST6 (°C)" : Sorties_SST[1][:,5], "Débit primaire_SST6 (kg/s)" : Sorties_SST[2][:,5], "Température de consigne_SST6 (°C)":Sorties_SST[3][:,5], "Ts_out_SST6 (°C)" : Sorties_SST[4][:,5], "Puissance_SST6 (kW)":Sorties_SST[5][:,5],"Ouverture_vanne_SST6 (%)":Sorties_SST[6][:,5],"Indicateur d'inconfort potentiel":Sorties_SST[7][:,5],\
                        "Tp_in_SST7 (°C)" : Sorties_SST[0][:,6], "Tp_out_SST7 (°C)" : Sorties_SST[1][:,6], "Débit primaire_SST7 (kg/s)" : Sorties_SST[2][:,6], "Température de consigne_SST7 (°C)":Sorties_SST[3][:,6], "Ts_out_SST7 (°C)" : Sorties_SST[4][:,6], "Puissance_SST7 (kW)":Sorties_SST[5][:,6],"Ouverture_vanne_SST7 (%)":Sorties_SST[6][:,6],"Indicateur d'inconfort potentiel":Sorties_SST[7][:,6],\
                        "Tp_in_SST8 (°C)" : Sorties_SST[0][:,7], "Tp_out_SST8 (°C)" : Sorties_SST[1][:,7], "Débit primaire_SST8 (kg/s)" : Sorties_SST[2][:,7], "Température de consigne_SST8 (°C)":Sorties_SST[3][:,7], "Ts_out_SST8 (°C)" : Sorties_SST[4][:,7], "Puissance_SST8 (kW)":Sorties_SST[5][:,7],"Ouverture_vanne_SST8 (%)":Sorties_SST[6][:,7],"Indicateur d'inconfort potentiel":Sorties_SST[7][:,7],})

    outputs.to_excel(Fichier_resultats)
    os.chdir(dossier)
    print("Simulation terminée")

    "Affichage des indicateurs de performance"
    print("La consommation de gaz s'élève à ", np.round(Part_gaz/6000,1)," MWh")
    print("L'énergie dévivrée est de ", np.round((Part_gaz+Part_geo)/6000,1), "MWh")
    print("Le taux d'énergies renouvelables est de ", np.round(100*Part_geo/(Part_geo+Part_gaz),1)," %")
    print("Les pertes en conduite s'élèvent à ",np.round(pertes,1)," % de la puissance totale délivrée")
    print("L'énergie non-fournie par les SST est respectivement de ",E_non_fournie[0],"kWh,",E_non_fournie[1], "kWh,",E_non_fournie[2],"kWh,",E_non_fournie[3],"kWh,",E_non_fournie[4],"kWh,",E_non_fournie[5],"kWh,",E_non_fournie[6],"kWh,",E_non_fournie[7],"kWh")

    return P_geo_red, P_gaz_red, T_cons_prod_red, T_res_ret_red, Sorties_SST, Input, temps_simu
