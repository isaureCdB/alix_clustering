==========================================================================================

rmsdXXX : Calcul de la matrice R
Crée un array numpy contenant la matrice des RMSD 2 à 2 d'un trinucléotide XXX donné. 

Utilisation : 
python rmsdXXX.py XXXr.npy

XXXr est le fichier des conformations pseudo-atome du trinucléotide de séquence XXX de votre choix. Le fichier XXXr.npy doit se trouver dans le même dossier que le script. L'array numpy créé est enregistré dans le même dossier sous le nom 'rmsdXXX.npy'
ATTENTION : L'exécution de ce script est longue. Par exemple, sur 4771 conformations de AAA, le script tourne en 23 minutes.

==========================================================================================

distribCI : Distribution des coordonnées internes
Calcule les coordonnées internes des conformations, enregistre les graphiques de leur distribution

Utilisation : 
python distribCI.py XXXr.npy 

XXr.npy est le fichier des conformations pseudo-atomes du trinucléotide XXX de votre choix. Le fichier XXXr.npy doit se trouver dans le même dossier que le script. 

==========================================================================================

repCI : Représentation en coordonnées internes
Crée un array numpy contenant les représentations en coordonnées internes d'un trinucléotide donné. 


Utilisation : 
python repCI.py XXXr.npy

XXXr.npy est le fichier des conformations pseudo-atome du trinucléotide XXX de votre choix. Le fichier XXXr.npy doit se trouver dans le même dossier que le script. L'array numpy créé est enregistré dans le même dossier sous le nom 'XXX_trigo.npy'

==========================================================================================

RMSDetCI : Étude de la relation entre similarité de coordonnées internes et RMSD
Renvoie le dictionnaire des "seuils" de similarité de coordonnées internes. Au dessus de ces seuils, on ne peut pas avoir un RMSD inférieur à 1Å.

Utilisation : 
python RMSDetCI.py ou lancement depuis un IDE. 

Le script peut encore être amélioré et étudié pour produire de meilleurs seuils. Pour l'instant, il utilise les données du trinucléotide AAA. Il enregistre le dictionnaire des seuils sous le nom "thresholds_CI.npy" dans le dossier courant.

==========================================================================================

rmsdXXX_fast.py : Calcul de la matrice R accéléré grâce à l'utilisation des seuils de similarité de coordonnées internes


Utilisation : 
python rmsdXXX_fast.py XXXr.npy XXX_trigo.npy thresholds_CI.npy

XXXr.npy est le fichier des conformations pseudo-atome du trinucléotide XXX.
XXX_trigo.npy est le fichier des conformations en coordonnées internes du trinucléotide XXX.
thresholds_CI.npy est le fichier contenant le dictionnaire qui définit les seuils de similarité maximaux à ne pas dépasser pour avoir un RMSD < 1Å. Il peut être calculé à l'aide du script RMSDetCI.py


==========================================================================================

CAH_withRMSD : Classification hiérarchique ascendante à partir de la matrice des RMSD 2à2
Effectue la classification hiérarchique avec le seuil de distance 1Å à partir de la matrice des similarités RMSD. Affiche le dendogramme et le nombre de clusters. 

Utilisation : 
python CAH_withRMSD rmsdXXX.npy

rmsdXXX.npy est le fichier contenant la matrice des RMSD 2 à 2 calculés pour les conformations de la séquence XXX. On peut aussi utiliser les matrices rmsdXXX_fast.npy calculés par le script précédent.

==========================================================================================
