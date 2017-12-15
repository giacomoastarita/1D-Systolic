/*
 =====================================================================================================================================
 Name        : MatrixMultiplication3S.c
 Author      : Antonio Tino & Giacomo Astarita
 Version     : 2.0
 Description : Bisogna effettuare il calcolo del prodotto C = A x B con A[m][h], B[h][n] e C[m][n].
			   Le matrici A, B e C sono suddivise in blocchi di colonne e bisogna assegnare al processo Pi, con i = 0, ..., p-1,
			   i blocchi Ai, Bi e Ci, rispettivamente di dimensione m x (h/p), h x (n/p) e m x (n/p).

			   Strategia: ANELLO
 =====================================================================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#define STAMPA 0

/* Assegna lo spazio di memoria alla matrice */
int* inizializza_spazio_matrice(int rows, int cols){
	int* Matrix = calloc(rows * cols, sizeof(int));
	return Matrix;
}

/* Creazione della matrice contenente numeri casuali */
void create_matrix(int *matrix, int row, int col){
	int i, j;

	for(i=0;i<row;i++){
		for(j=0;j<col;j++)
			matrix[col*i+j] = rand() % 10;
	}
}

/* Stampa a video la matrice */
void stampa_matrice(int* matrice, int rows, int cols, char* nome){
	if (STAMPA==0)
		return;

	printf("\nStampa matrice %s\n",nome);
	fflush(stdout);

	for(int i=0;i<rows;i++){
		for(int j=0;j<cols;j++){
			printf("%d\t",matrice[i*cols+j]);
			fflush(stdout);
		}
		printf("\n");
		fflush(stdout);
	}
}

/* Calcola il lavoro da assegnare ad ogni processore */
void calculateWork(int *work, int nproc, int colonne){
	int w, i;
	
	if(colonne%nproc==0)
		for(i=0;i<nproc;i++)
			work[i]=colonne/nproc;
	else{
		int resto=colonne-(colonne/nproc)*nproc;
		for(i=0;i<nproc;i++){
			if (resto>0){
				work[i]=colonne/nproc+1;
				resto--;
			}else 
				work[i]=colonne/nproc;
		}
	}
}

/* Effettua il prodotto tra matrici */
void prodotto_matrici(int* A, int* B, int* C, int rowA, int colA, int rowB, int colB, int righeBfatte){

	for(int i=0;i<rowA;i++){
		for(int j=0;j<colB;j++){
			for(int k=0;k<colA;k++){
				//colA==rowB altrimenti non si può effettuare il prodotto
				C[i*colB+j]+=A[i*colA+k]*B[(righeBfatte+k)*colB+j];
			}
		}
	}
}

/* Copia le varie sottomatrici da inviare */
void copia_sottomatrici(int* originale, int* copia, int row, int col, int *colW, int idproc){
	int colonna_iniziale, colonna_finale, offset=0,i;
	int riga=0, colonna=0;
	
	for(i=0;i<idproc;i++)
		offset+=colW[i];

	colonna_iniziale=offset;
	colonna_finale=colonna_iniziale+colW[idproc];

	for(i=colonna_iniziale;i<colonna_finale;i++){
		for(int j=0;j<row;j++){
			copia[riga*colW[idproc]+colonna]=originale[j*col+i];
			riga++;
		}
		colonna++;
		riga=0;
	}
}

/* Ripristina le varie sottomatrici */
void ripristina_sottomatrici(int* risultato, int* parziale, int rowR, int colR, int *colW, int idproc, int colP){
	int i=0, j=0, offset=0;
	
	for(i=0;i<idproc;i++)
		offset+=colW[i];

	for(i=0;i<rowR;i++)
		for(j=0;j<colP;j++)
			risultato[i*colR+j+offset]=parziale[i*colP+j];
}

int main(int argc, char *argv[]){
	
	/* Id del processo e numero di processi  */
	int nproc, myid;

	/* Puntatori alle matrici */
	int *A, *B, *C;
	
	/* Puntatori */
	int *Aloc, *Bloc, *Cloc, *Mtemp, *colAp, *colBp;
	
	/* Dimensioni delle matrici A e B */
	int rowA, rowB, colA, colB;

	/* Id del processo a cui inviare e id del processo da cui ricevere*/
	int idSend, idRecv;

	/* Variabili per il calcolo del tempo */
	double start, finish;

	/* Indici */
	int  i, j, k;

	/* Restituisce lo status*/
	MPI_Status status;

	/* Start MPI */
	MPI_Init(&argc, &argv);

	/* Trova l'id del processo */
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	/* Trova il numero di processi */
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	
	/* Controllo che l'utente abbia inserito correttamente i dati, in caso negativo utilizzo quelli di default */
	if(argc < 5){
		rowA=5;
		colA=4;
		rowB=4;
		colB=3;

		if(myid==0){
			printf("PARAMETRI DI DEFAULT\n");
			printf("La matrice A è di dimensione %d x %d\n", rowA, colA);
			printf("La matrice B è di dimensione %d x %d\n", rowB, colB);
			fflush(stdout);
		}
		
	}else{
		rowA=atoi(argv[1]);
		colA=atoi(argv[2]);
		rowB=atoi(argv[3]);
		colB=atoi(argv[4]);

		if(myid==0){
			printf("PARAMETRI INSERITI DALL'UTENTE\n");
			printf("La matrice A è di dimensione %d x %d\n", rowA, colA);
			printf("La matrice B è di dimensione %d x %d\n", rowB, colB);
			fflush(stdout);
		}
	}

	/* Controllo che il numero di colonne della matrice A siano uguali al numero di righe della matrice B */
	if(colA!=rowB) {
		if(myid==0){
			printf("Impossibile effettuare la moltiplicazione (colA != rowB)\n");
			fflush(stdout);
		}
		exit(0);
	}
	
	/* Controllo che il numero di processori non sia maggiore del numero di colonne di A */
	if(nproc>colA){
		if(myid==0){
			printf("Numero di processi %d è inferiore alle colonne di A %d. Impossibile distribuire il lavoro\n", nproc, colA);
			fflush(stdout);
		}
		exit(0);
	}

	/* Vettore in cui mi salvo il numero di colonne delle matrici per ogni processore*/
	colAp=inizializza_spazio_matrice(1, nproc);
	colBp=inizializza_spazio_matrice(1, nproc);

	/* Calcola il lavoro da distribuire ad ogni processore */
	calculateWork(colAp, nproc, colA);
	calculateWork(colBp, nproc, colB);

	/* Allocazione di memoria per le matrici A, B e C locali */
	Aloc=inizializza_spazio_matrice(rowA, colAp[myid]);
	Bloc=inizializza_spazio_matrice(rowB, colBp[myid]);
	Cloc=inizializza_spazio_matrice(rowA, colBp[myid]);

	if(myid == 0){
		/* Inizio tempo di calcolo */
		printf("\nInizio tempo di calcolo...\n\n");
		fflush(stdout);
		start = MPI_Wtime();

		/* Allocazione di memoria per le matrici A e B*/
		A=inizializza_spazio_matrice(rowA, colA);
		B=inizializza_spazio_matrice(rowB, colB);

		srand(13);
		
		/*Creazione delle matrici A e B*/
		create_matrix(A, rowA, colA);
		create_matrix(B, rowB, colB);

		/*Se il flag STAMPA è settato ad 1, si effettua la stampa delle matrici A e B*/
		stampa_matrice(A, rowA, colA, "A");
		stampa_matrice(B, rowB, colB, "B");

		/* Assegno i primi elementi delle matrici A e B al processo P0 */
		copia_sottomatrici(A, Aloc, rowA, colA, colAp, myid);
		copia_sottomatrici(B, Bloc, rowB, colB, colBp, myid);

		/* Il processo P0 assegna il lavoro agli altri processori */
		for(i=1;i<nproc;i++){
			/* Allocazione di memoria per la matrice A temporanea*/
			Mtemp=inizializza_spazio_matrice(rowA, colAp[i]);

			/* Copia le sottomatrici da inviare */
			copia_sottomatrici(A, Mtemp, rowA, colA, colAp, i);

			MPI_Send(Mtemp, rowA*colAp[i], MPI_INT, i, 0, MPI_COMM_WORLD);

			/* Libera la memoria utilizzata per la matrice temporanea*/
			free(Mtemp);

			/* Allocazione di memoria per la matrice B temporanea*/
			Mtemp=inizializza_spazio_matrice(rowB, colBp[i]);

			/* Copia le sottomatrici da inviare */
			copia_sottomatrici(B, Mtemp, rowB, colB, colBp, i);

			MPI_Send(Mtemp, rowB*colBp[i], MPI_INT, i, 0, MPI_COMM_WORLD);

			/* Libera la memoria utilizzata per la matrice temporanea*/
			free(Mtemp);
		}

		/* Libera la memoria utilizzata per la matrice A  e la matrice B*/
		free(A);
		free(B);
	}
	else {
		MPI_Recv(Aloc, rowA*colAp[myid], MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(Bloc, rowB*colBp[myid], MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	}

	/* Computazione */
	int passo=0, righeBfatte=0;

	/* Ogni processo calcola l'identificativo del processo a cui inviare */
	idSend=(myid+1)%nproc;

	/* Ogni processo calcola l'identificativo del processo da cui ricevere */
	if(myid-1<0)
		idRecv=nproc-1;
	else
		idRecv=myid-1;

	int realIdRecv=myid;
	int tempid;

	/* Fase di calcolo */
	for(i=0;i<myid;i++)
		righeBfatte+=colAp[i];
	
	prodotto_matrici(Aloc, Bloc, Cloc, rowA, colAp[realIdRecv], rowB, colBp[myid], righeBfatte);
	
	int *Aloctemp=inizializza_spazio_matrice(rowA, colAp[myid]);
	
	for(i=0;i<nproc-1;i++){
		passo++;
		
		if(myid-passo<0)
			tempid=nproc-passo;
		else 
			tempid=myid-passo;

		if(colAp[realIdRecv]!=colAp[tempid]){
			Aloctemp=inizializza_spazio_matrice(rowA, colAp[tempid]);
		}
		
		/* Si utilizza la MPI_Sendrecv per effettuare le send e le recv in modo parallelo, quindi evitando problemi di buffer*/
		MPI_Sendrecv(Aloc, rowA*colAp[realIdRecv], MPI_INT, idSend, 0, Aloctemp, rowA*colAp[tempid], MPI_INT, idRecv, 0 ,MPI_COMM_WORLD, &status);
		
		realIdRecv=tempid;
		Aloc=Aloctemp;
		
		if(righeBfatte-colAp[realIdRecv]<0)
			righeBfatte=rowB-righeBfatte-colAp[realIdRecv];
		else
			righeBfatte=(righeBfatte-colAp[realIdRecv])%rowB;

		/* Fase di calcolo */
		prodotto_matrici(Aloc, Bloc, Cloc, rowA, colAp[realIdRecv], rowB, colBp[myid], righeBfatte);
	}

	/* Ogni processo libera la memoria che ha utilizzato */
	free(Aloc);
	free(Bloc);

	if(myid==0){

		/* Allocazione di memoria per la matrice C */
		C=inizializza_spazio_matrice(rowA, colB);

		/* Ripristina la propria sottomatrice */
		ripristina_sottomatrici(C, Cloc, rowA, colB, colBp, myid, colBp[myid]);

		for(i=1;i<nproc;i++){
			MPI_Recv(Cloc, rowA*colBp[i], MPI_INT, i, 0, MPI_COMM_WORLD, &status);

			/* Ripristina le sottomatrici degli altri processi*/
			ripristina_sottomatrici(C, Cloc, rowA, colB, colBp, i, colBp[i]);
		}

		/* Fine tempo di calcolo */
		printf("\n...Fine tempo di calcolo\n");
		fflush(stdout);
		finish = MPI_Wtime();

		stampa_matrice(C,rowA,colB,"C");

		printf("\nTempo di computazione: %f\n", finish - start);
		fflush(stdout);

		free(C);
	}else{
		/* Ogni processo invia a P0 la sua parte della matrice C e libera la memoria utilizzata*/
		MPI_Send(Cloc, rowA*colBp[myid], MPI_INT, 0, 0, MPI_COMM_WORLD);
		free(Cloc);
	}
	
	/* Shut down MPI */
	MPI_Finalize();

	return 0;
}