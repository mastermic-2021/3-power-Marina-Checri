encodegln(s,n)={
  my(v);
  v=[if(x==32,0,x-96)|x<-Vec(Vecsmall(s))];
  if(#v>n^2,warning("string truncated to length ",n^2));
  v = Vec(v,n^2);
  matrix(n,n,i,j,v[(i-1)*n+j]);
}

decodegln(M)={
  my(v);
  \\on concatène les entiers les uns à la suite des autres.
  v= Vec( concat(Vec(M~))~, (matsize(M)[1])*(matsize(M)[1]) - 1);
  v;
}
\\On traduire en lettres les entiers.
decode(v) = {
  v= Strchr([ if(c==0,32,c+96) | c <- v]);
  v;
}


/* ****************************************************************************** */
/* ******************************** Fonctions : ********************************* */
/* ****************************************************************************** */

\\(fonctions écrites précédemment, pour l'exponentiation rapide d'une matrice)
ExpMat(M,e) = {
	if(e == 0, return (matid(matsize(M)[1])),
	    if(e == 1, return (M),
	        if(e%2 == 0,return (ExpMat(M^2,e/2)), return (M*ExpMat(M^2,(e-1)/2)) )) );
}

IterN(M,n)={
	my(f);
	f=factor(n);
	for(i=1,matsize(f)[1], M=ExpMat(M,f[i,1]^f[i,2]));
	M;
}


/*
 * On cherche à trouver l'exposant k pour lequel la matrice M,
 * élevée à la puissance k est la matrice identité.
*/
ord(M)={
  my(Id);
  Id=matid(matsize(M)[1]);
  my(k);
  k=1;
  my(puissance);
  puissance=Mod(M,27);
  while(puissance!=Id,puissance=puissance*M;
	k++);
  k;
}


/*
 * Une idée que je n'ai pas encore réussi à mettre en pratique est de
 * trigonaliser la matrice avant de chercher son ordre, puisque dans ce cas,
 * on aurait M = PTP^(-1), où T est triangulaire, et
 * M^k = P T^k P^(-1), avec des calculs bien moins longs pour T^k.
 * *Mais il faudrait déjà s'assurer que la matrice est trigonalisable,*
 * *ce, à l'aide de ses vecteurs/valeurs propres (utiliser mateigen(x) ... ?), ...*

 * Remarque -5 est un générateur de (Z/27Z)*, puisque d'ordre 18 ...
*/


/* ****************************************************************************** */
/* **************************** Posons nos données : **************************** */
/* ****************************************************************************** */

n=65537; \\ n = 2^16 + 1, n premier
\\on récupère la matrice donnée dans le input.
M  = readstr("input.txt")[1];
M  = encodegln(M,12);


/* ****************************************************************************** */
/* ******************************* Application : ******************************** */
/* ****************************************************************************** */

/*
 * Le principe est le suivant :
 * On cherche la puissance ord pour laquelle M^ord = Id.
 * On trouve l'inverse (n') de n modulo ord,
   afin de pouvoir ensuite n'avoir qu'à élever M à cette puissance
   (en effet, on a M = message^65537, donc (message^n)^n'=message^(n*n')=message).
 * Nota Bene : c'est une simple relation de Bézout qu'il faut trouver.
 * Ensuite, il suffira d'élever notre matrice à la puissance n'
 * Il ne reste qu'à décoder !
*/

\\La matrice est dans Z/27Z :
M = Mod(M,27);

\\ On récupère l'ordre de la matrice
ord=ord(M);

\\ On résout n'*n=n'*65537=1[ord]
\\ (c'est une relation de Bézout : n*u = 1 + v*ord)
inverse_n=gcdext(n,ord)[1];

M=IterN(M,inverse_n);

\\on ramène M dans Z pour pouvoir décoder :
M=lift(M);

\\Puis, on décode !
message=decodegln(M);
message=decode(message);

print(message);
