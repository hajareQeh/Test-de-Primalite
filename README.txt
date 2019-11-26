## Projet Cryptographie - 2018
Dans un premier temps, pour avoir accès au projet, il est nécessaire de 
décompresser l'archive projet_se_shade.tar.gz

## Compiler et executer le projet
Le projet se compile et s'execute en utilisant les commandes suivant :

	gcc projetfin.c -o projet -Wall -lgmp
	./projet
Lors de l'exécution, ce dernier attendra que l'utilisateur saisisse un nombre d'iteration, pour le nombre a tester il faut l'ecrire sur le fichier "test.txt".
attention faut pas avoire des retoure a la ligne dans le nombre écrit dans le fichier "teste.txt".

Aucune installation préalable n'est nécessaire.

notre projet :
- une methode pgcd qui calcule le pgcd de deux grand nombres
-square_and_multiply
-mode: utiliser dans les cas ou l'un des parametres est négative
-facteur :fais la composition d'un nombre [ex : a=2^d * c] nous permet de trouver le d et c 
-reciprocite_quad : Loi de r´eciprocit´e quadratique Si pgcd(m, n) = 1 
-jacobi1
-solovay_strassen
	
