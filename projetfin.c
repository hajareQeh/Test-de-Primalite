#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

int jacobi1(mpz_t a,mpz_t n);
void mode(mpz_t res,mpz_t a,mpz_t n);


//~ //square and multiply
void square_and_multiply(mpz_t r,const mpz_t a,const mpz_t n,char * h,int t){
		int i;
		mpz_set(r,a);
			for(i=1;i<t;i++){
				mpz_mul(r,r,r);
				mpz_mod(r,r,n);
					if(h[i] == '1'){
						mpz_mul(r,r,a);
						mpz_mod(r,r,n);
					}
			}
	}

//calcul pgcd

void pgcd(mpz_t res,const mpz_t a,const mpz_t b){
	mpz_t y;
	mpz_init(y);
	if( mpz_cmp_ui(b,0) == 0) 	mpz_set(res,a);
	else{ 
								mpz_mod(y,a,b);
								pgcd(res,b,y);
	}
	mpz_clear(y);
}
//Jacobi
	//etape 1
void mode(mpz_t res,mpz_t a,mpz_t n){
		mpz_t b;
		mpz_init(b);
		mpz_set(b,a);
			if(mpz_cmp_si(b,0)<0){
				mpz_mul_si(b,b,-1);
				mpz_mod(res,b,n);
				mpz_mul_si(res,res,-1);
			}else mpz_mod(res,a,n);
			
		mpz_clear(b);	
	}
	
	//~ etape2: rassemble la propriété 3 et 5
void facteur(mpz_t o,const mpz_t a,const mpz_t n){
	mpz_t k, res,d;
	int nb = 0;
	int j;
	mpz_init(k);
	mpz_init(d);
	mpz_init(res);
	mpz_set(k,a);
	if(mpz_cmp_si(k,0)<0)
			mpz_mul_si(k,k,-1);
	
	if(mpz_cmp_ui(k,1) > 0){
			mpz_mod_ui(res,k,2);
			j = mpz_cmp_ui(res,0);
		while(j==0){
			mpz_div_ui(k,k,2);
			nb++;
			mpz_mod_ui(res,k,2);
			j = mpz_cmp_ui(res,0);
		}
	}
	//pro 5	
	if(nb!=0){
			if(mpz_cmp_si(a,0)<0)
						mpz_mul_si(k,k,-1);
			mpz_mod_ui(d,n,8);
			if(mpz_cmp_si(d,1)==0 || mpz_cmp_si(d,7)==0){
						mpz_set(o,k);	
			}else{
					if(mpz_cmp_si(d,3)==0 || mpz_cmp_si(d,5)==0){
							if(nb%2==0) j=1; 
							else j=-1;
							
							mpz_mul_si(k,k,j);
							mpz_set(o,k);	
					}
			}
	}else
		mpz_set(o,a);
		mpz_clear(k);
		mpz_clear(res);
		mpz_clear(d);
}

		
		
	
//etape 3
int reciprocite_quad(mpz_t a,mpz_t n){
	mpz_t d;
	mpz_t c;
	mpz_t x;
	mpz_init(x);
	mpz_init(d);
	mpz_init(c);
	mpz_set(x,a);
			mpz_mod_ui(d,n,4);
			if (mpz_cmp_si(a,0)<0)
			mpz_mul_si(x,a,-1);
			mpz_mod_ui(c,x,4);
			if(mpz_cmp_si(a,0)<0) mpz_mul_si(c,c,-1);
			if(mpz_cmp_si(d,1)==0 || mpz_cmp_si(c,1)==0){ 
					if(mpz_cmp_si(a,0)<0){
						mpz_mul_si(n,n,-1);
						mpz_mul_si(a,a,-1);
					}
						mpz_clear(d);
						mpz_clear(c);
						return jacobi1(n,a);//(n/a);
			
			}else{ 
				if(mpz_cmp_si(c,3)==0 && mpz_cmp_si(d,3)==0){ 
					if(mpz_cmp_si(a,0)<0)
						mpz_mul_si(a,a,-1);
					else mpz_mul_si(n,n,-1);
						mpz_clear(d);
						mpz_clear(c);
						return jacobi1(n,a);//(-n/a);
				}else 
				{
						mpz_div(c,a,n);
						return mpz_get_si(c) ;
						
				}
		}
	}


//deput de methode de jacobi appel de etape1 apres 2 ensuit appliqué etape 3 puis appler fonction de l'etape 4
int jacobi1(mpz_t a,mpz_t n){
	mpz_t k;
	mpz_init(k);
	mode(a,a,n);
	facteur (a,a,n);
	if (mpz_cmp_si(a,0)==0) return 0;
	if (mpz_cmp_si(a,1)==0) return 1;
	else{
			if(mpz_cmp_si(a,-1)==0) return -1;
			else{
					pgcd(k,a,n);
					if(mpz_cmp_si(k,1) !=0){ 
							mpz_clear(k);
							return 0;
					}else{
							mpz_clear(k);
							return reciprocite_quad(a,n);
					}
			}		
		}
	}




// Solovay-Strassen
int solovay_strassen(mpz_t n,int k){
		mpz_t a,tmp,d,H;
		int r,i,t;
		gmp_randstate_t gmpRandState;
 		mpz_init(a);
 		mpz_init(d);
		mpz_init(tmp);
		mpz_init(H);
		gmp_randinit_default(gmpRandState);
		mpz_sub_ui(H,n,1);  //H=n-1
		mpz_div_ui(H,H,2);	//H=H/2
		t=mpz_sizeinbase(H,2);	//taille de H en binaire
		char *h= malloc(sizeof(char) * t);
		mpz_get_str(h,2,H);
		
		for(i=1;i<=k;i++){
					mpz_set(tmp,n);
				do{
		 			mpz_urandomm(a,gmpRandState,tmp);
				}while(mpz_cmp_ui(a,2)< 0);
				gmp_printf("a=%Zd\n",a);
					square_and_multiply(d,a,tmp,h,t);
					mpz_set(tmp,n);
					r=jacobi1(a,tmp);
					mode(d,d,n);
					mpz_sub(tmp,d,n);
					if((mpz_cmp_si(tmp,r)!=0 && mpz_cmp_si(d,r)!=0) || r==0) {
					mpz_clear(a);
					mpz_clear(tmp);
					free(h);
					mpz_clear(H);
					mpz_clear(d);
					return 0;
				}
		}
		mpz_clear(tmp);
		mpz_clear(a);
		mpz_clear(d);
		mpz_clear(H);
		free(h);
		return 1;
	}

//main
int main(){
    mpz_t temp;
	mpz_init(temp);
	int k,j;
	printf("entre le nombre d'itération :");
	scanf("%d",&k);
	FILE *fp;
	fp = fopen("test.txt", "r");
	mpz_inp_str(temp, fp,0);
	j=solovay_strassen(temp,k);
	if(j==0)	gmp_printf("composé\n");
	else 		gmp_printf("premier\n");
	mpz_clear(temp);
	return 0;
}
