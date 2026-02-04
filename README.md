## Contenu

Cet ensemble de notebooks sert de support au projet "Analyse de la performance d'un réseau de chaleur" proposé dans le cadre du cours "Réseaux de chaleur" donné dans l'option Energies du cycle ingénieur de l'IMT Albi-Carmaux.

Il est composé de 5 notebooks à réaliser dans l'ordre suivant :

- 1_SousStation

- 2_Analyse_SST

- 3_Analyse_Canalisation

- 4_Analyse_Production

- 5_Analyse_réseau

Les répertoires /Images et /Models ne doivent pas être modifiés.

Toutes les sauvegardes de simulations ainsi que les graphiques proposés dans les notebooks sont enregistrés dans le répertoire /Ouputs au fur et à mesure de leur production.

Dans certain cas, en particulier pour le notebook 5_Analyse_réseau, vous serez invités à modifier des paramètres des modèles qui se trouvent dans le répertoire /Inputs.

## Téléchargement du dépôt

Pour télécharger l'ensemble du dossier, deux options :

- utiliser "Create a new fork" pour travailler directement sur votre propre dépôt Git

- ou télécharger le [dossier compressé](https://github.com/BrunoDuplessis/Projet-RDV-Albi---Eleve/archive/refs/heads/Projet-2026.zip) pour un usage en local sur votre PC

NE PAS CLONER directement le dépôt !

## Création d'un environnement virtuel

Avant utilisation des notebooks, il convient de créer puis d'utiliser un environnement virtuel. Pour cela, on suivra la procédure suivante. Dans le terminal, lancer les commandes :

- sous Windows

```bash
    python -m venv venv
    venv\Scripts\activate
    pip install -r requirements.txt
```

- sous macOS / Linux

```bash
    python -m venv venv
    source venv/bin/activate
    pip install -r requirements.txt
```
