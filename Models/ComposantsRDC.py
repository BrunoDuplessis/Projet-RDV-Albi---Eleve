from math import exp, pi
import numpy as np

"constantes"
Cp = 4180
Rho = 1000

def SST_nom(Kp,Ki,Kd,Trad_in_nom,Trad_out_nom,E_nom,P_nom): # Détermination des caractéristiques nominales d'une sous-station"

    "entrées : gain régulateur, température nominales au secondaire (sortie et entrée), efficacité nominale, puissance nominale"

    "sorties : vecteur des caractéristiques nominales (Débit, NUT, UA, coefficients de convection côté primaire et secondaire, paramètres de la loi de chauffe, gain régulateur)"

    "Calcul des caractéristiques nominales"
    Mp_nom = P_nom/(Cp*(Trad_in_nom-Trad_out_nom))
    NUT_nom = E_nom/(1-E_nom)
    UA_nom = Cp*Mp_nom*NUT_nom
    hAp_nom = 2*UA_nom
    hAs_nom = hAp_nom
    
    "Calcul de la loi de chauffe"
    a = (Trad_in_nom-30)/(-7-20)
    b = Trad_in_nom-a*(-7)
    
    return Mp_nom, NUT_nom, UA_nom, hAp_nom, hAs_nom, a, b, Kp, Ki, Kd, Trad_in_nom, Trad_out_nom, E_nom, P_nom 

def SST_fonct(Nom, Youv, Tp_in, P_appelée, Text, EcartT_1, IntegraleEcartT, dt): # Calcul du point de fonctionnement d'une sous-station"
    
    "Entrées : vecteur des caractéristiques nominales de la SST, T° entrée primaire, ouverture de vanne (du pas de temps précédent !), Puissance appelée, T° extérieure, temps (indice)"
    "Sorties : puissance fournie, T° de consigne, T° sortie secondaire, T° sortie primaire, débit, T° entrée secondaire, écart à la consigne, intégrale de l'écart à la consigne, ouverture de vanne côté primaire" 
    
    "Bilan échangeur"
    Mp_nom = Nom[0]
    NUT_nom = Nom[1]
    UA_nom = Nom[2]
    hAp_nom = Nom[3]
    hAs_nom = Nom[4]
    a = Nom[5]
    b = Nom[6]
    Kp = Nom[7]
    Ki = Nom[8]
    Kd = Nom [9]
    Trad_in_nom = Nom[10]
    Trad_out_nom = Nom[11]
    E_nom = Nom[12]
    P_nom = Nom[13]

    Mp = Youv * Mp_nom # Calcul du débit primaire en fonction de l'ouverture de vanne
    hAp = hAp_nom * (Mp / Mp_nom)**(0.8) # Calcul coefficient de convection primaire
    UA = (1 / hAp + 1 / hAs_nom)**(-1) # Calcul du coefficient de transfert
    NUT = UA /(Cp * Mp) # Calcul du NUT
    R = Mp / Mp_nom # Calcul du rapport des débits
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
    Ts_out = Ts_in + P_fournie/(Mp_nom*Cp) # Calcul température de sortie côté secondaire        

    " Régulation "
    Tcons = a * Text + b # Loi d'eau secondaire
    EcartT = Tcons-Ts_out
    IntegraleEcartT = IntegraleEcartT + EcartT*dt
    derivee = (EcartT-EcartT_1)/dt 
    Youv = Youv+ Kp*EcartT + Kd*derivee + Ki*IntegraleEcartT  # Calcul de l'incrément d'ouverture de la vanne
    Youv = np.clip(Youv, 0.001, 1.0)# Calcul de l'ouverture de vanne finale comprise entre 0 et 1

    return (P_fournie, Tcons, Ts_out, Tp_out, Mp, Ts_in, EcartT, IntegraleEcartT, Youv)

def Tuyau_dim(Lambda_iso, Longueur, Mtube_nominal) : # Dimensionnement d'une canalisation"
    
    "Entrées : conductivité isolant (W/(m.K)), longueur de la canalisation (m), débit nominal (kg/s)"
    
    "Sortie : vecteur des caractéristiques de la canalisation (Longueur, débit nominal, diamètre, épaisseur d'isolant"
    
    V_lim = 2
    Long = Longueur
    Mtub_nom = Mtube_nominal 
    D_tube = np.round( (4 * Mtub_nom / (Rho * V_lim *3.14))**(0.5), 2) # Calcul du diamètre des canalisation
    e_iso = 30.77 * D_tube**3 - 11.43 * D_tube**2 + 1.45 * D_tube #épaisseur d'isolant
    
    return Long, Mtub_nom, D_tube, e_iso, Lambda_iso

def Tuyau_fonct(carac_tube, T_in, M_tub): # Calcul de l'évolution de la température dans une canalisation"
    
    "Entrées :  vecteur des caractéristiques du tuyau, Température en entrée, Débit"
    
    "Sorties : température en sortie"
    
    "Initialisation"
    Delta_temps = 10 
    Long = carac_tube[0]
    Mtub_nom = carac_tube[1]
    D_tube =carac_tube[2]
    e_iso =carac_tube[3]
    Lambda_iso = carac_tube[4]
    T_sol = 10
    
    "paramètres du schéma numérique de résolution du modèle"
    delta_t = 1 # pas de temps de calcul en minutes
    N=int(Long/10)      #Nombre initial de volumes dans le tuyau
    nb_dt = int(Delta_temps/delta_t)
    dt = 60*delta_t        #conversation pas de temps de calcul en secondes 
        
    "Initialisation du maillage du tuyau"
    L_vol = np.ones(N)*Long/N                                                          # longueur de chaque volume 
    T_vol=np.ones(N)*50                                                                     # Température de chaque volume
    M_vol=np.ones(N)*Rho*Cp*pi * D_tube**2 / 4 * L_vol                                 # Masse thermique de chaque volume 
    R_vol = np.log(( 2 * e_iso + D_tube) / D_tube)/(2 * pi *Lambda_iso * L_vol)   # Résistance thermique de chaque volume

    if M_tub==0:
        M_tub=0.0001

    for k in range(nb_dt):
        M_0 = np.array([Cp*M_tub*dt])
        T_0 = np.array([T_in])
        L_0 = np.array([M_tub*dt/(Rho*pi*(D_tube/2)**2)])            

        T_vol2 = np.concatenate((T_0,T_vol*(1-dt/(R_vol*M_vol)) + dt/(R_vol*M_vol)*T_sol))
        L_vol2 = np.concatenate((L_0,L_vol))
        M_vol2 = np.concatenate((M_0,M_vol))

        L_vol3 = np.concatenate((L_vol2[np.cumsum(L_vol2)<Long],np.array([Long-np.sum(L_vol2[np.cumsum(L_vol2)<Long])])))
        T_vol3 = np.resize(T_vol2,np.shape(L_vol3))
        M_vol3 = Rho*Cp*pi * D_tube**2 / 4 *L_vol3
        M_vol22 = np.copy(M_vol3)
        M_vol22.resize(np.shape(M_vol2),refcheck=False)

        T_vol = T_vol3
        L_vol = L_vol3
        M_vol = M_vol3
        R_vol = np.log(( 2 * e_iso + D_tube) / D_tube)/(2 * pi *Lambda_iso * L_vol)

    Tout = np.sum((M_vol2-M_vol22)*T_vol2)/np.sum(M_vol2-M_vol22)
            
    return Tout

def Production_chaleur(T_in_prod, M_res, Text): # Calcul des puissances et consommations d'énergie géothermique et de gaz
    "Entrées : Temperature d'entrée du site de production, T° extérieure, Débit réseau"
    
    "Sorties : Puissance geothermie, gaz et totale, Température de départ du réseau, température de sortie de l'échangeur geothermique"
    
    E_geo = 0.9
    eta_gaz = 0.85
    Delta_temps = 10

    
    "Détermination de la loi de chauffe primaire"
    T_cons_prod = -1.85 *Text + 77 #loi de chauffe primaire
    if T_cons_prod < 60:
        T_cons_prod = 60
    else:
        T_cons_prod = -1.85 * Text +77
            
    T_geo_out = E_geo * (70 - T_in_prod) + T_in_prod# Calcul de la température de sortie réseau de l'echangeur geothermie
    if T_geo_out > T_cons_prod:
        T_geo_out = T_cons_prod
        P_geo = Cp * M_res * (T_geo_out - T_in_prod) # Calcul puissance geothermie
    else:
        P_geo = Cp * M_res * (T_geo_out - T_in_prod) # Calcul puissance geothermie
    
    P_gaz = M_res * Cp * (T_cons_prod - T_geo_out) # Calcul puissance chaudière gaz
    P_gaz_reel = P_gaz / eta_gaz # calcul puissance gaz consommée
    
    P_res_tot = P_geo + P_gaz  # calcul puissance totale fournie au réseau
    Conso_gaz = P_gaz_reel*Delta_temps/60 # Calcul consommation de gaz
    Conso_geo = P_geo*Delta_temps/60 # Calcul consommation géothermie
    
    return P_geo, P_gaz, P_res_tot, Conso_gaz, Conso_geo
   
def Production_chaleur(T_res_ret, E_geo, M_res, eta_gaz, T_cons_prod, Text):
    " Input : Temperature d'entrée du site de production (retour réseau) et de consigne, Débit réseau, Efficacité echangeur géothermal, Rendement chaudière"
    " Output : Puissance geothermie, gaz et totale, Température de départ du réseau, température de sortie de l'échangeur geothermique"
        
    T_geo_out = E_geo * ( 65 - T_res_ret ) +T_res_ret # Calcul de la température de sortie réseau de l'echangeur geothermie
    if T_geo_out > T_cons_prod:
        T_geo_out = T_cons_prod
        P_geo = Cp * M_res * ( T_geo_out - T_res_ret ) # Calcul puissance geothermie
    else:
         P_geo = Cp * M_res * ( T_geo_out - T_res_ret ) # Calcul puissance geothermie
    

    P_gaz = M_res * Cp * ( T_cons_prod - T_geo_out ) # Calcul puissance chaudière gaz
    P_gaz_reel = P_gaz / eta_gaz # calcul puissance gaz réelle
    
    P_res_tot = P_geo + P_gaz  # calcul puissance totale fournie au réseau
    
    return P_geo, P_gaz_reel, P_res_tot, T_cons_prod, T_geo_out, T_cons_prod
   
