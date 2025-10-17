/* NEIGHBORHOOD.c (29/09/2025) */

/* ***************************************************** */

#include "libhdr"
#define NN 90001  /* max number of NIND */
#define GG 10001  /* max number of GEN */
#define MM 701  /* max number of NCRO */

int BW, SW, ran_i, ran_h, fil, col, Sfil, Scol, f[NN], c[NN], Vk;
int marker, NIND, NCRO, NLOCI, SNP, numLOCI[GG];
int gen, generations, i, k, l;
int RM[31], p1, p2;
int gm[NN][MM][2], sm[NN][MM][2];
int slocus[NN][MM][2], locus[NN][MM][2];
int crom, chrom[MM][31];
unsigned long long int pos[MM][31];

double AA, Aa, aa, q[MM][31];
double fibd[NN], fLH1[NN];
double EHom, EHet, Hom[NN];
double NeFibd, NeHet, NeHetb, NeSk2;
double FISmax, FSTmax, FITmax;
double FISmin, FSTmin, FITmin;
double qNS1[MM][31], qNS2[MM][31], qNS1NS2[MM][31];
double qNS3[MM][31], qNS4[MM][31], qNS3NS4[MM][31];
double EHetNS1, EHetNS2, EHetNS1NS2;
double EHetNS3, EHetNS4, EHetNS3NS4;
double HT, HS, HI, HetNS1, HetNS2, HetNS1NS2, lociNS1NS2;
double HetNS3, HetNS4, HetNS3NS4, lociNS3NS4;

struct acc Het[GG], Fibd[GG], FLH1[GG], AVE_FLH1, GEN, SK2;
struct covacc GENHET;

FILE *fptr, *frep, *fgen, *ffmap, *ffped, *ffpedR;

/* ***************************************************** */

main()
{
	fptr = fopen ("dfilename.dat","w");
	fgen = fopen ("genfile.dat","w");

	getinputs ();
	recombination_masks ();

	initial_generation ();

	for (gen=0; gen<=generations; gen++)
	{
		frep = fopen ("repfile.dat","a");
		if (gen%1000 == 0)	fprintf(frep,"GEN = %d\n", gen);
		fclose(frep);

		frequency_genes ();
		mutation_neutral ();
		if (gen == generations)	plink_files();
		if (gen == generations-15)	initial_multiallelic();
		if (gen >= generations-15)	estimates_of_F();
		if (gen == generations)	estimates_of_Fst_between_opposite_NS();
		if (gen == generations)	estimates_of_Fst_between_close_NS();

		mating ();

//		if (tracelevel!=0) dumpoffspring();
	}

 	printout();
}

/* ***************************************************** */

getinputs()
{
	tracestart();
	getseed();

	getintandskip("NIND (max 90000):",&NIND,2,90000);
	getintandskip("Vk (0 - 1):",&Vk,0,1);
	getintandskip("BW (odd; 99: RM):",&BW,0,1000);
	getintandskip("SW (odd):",&SW,0,1000);
	NCRO=680;
	NLOCI=30;
	getintandskip("Number of generations :",&generations,1,10000);
}

/* **************************************************** */

recombination_masks ()
{
	for (l=0; l<NLOCI; l++)   RM[l]=pow(2.0,(double)l);
//	if (tracelevel!=0)  for (l=0; l<NLOCI; l++)   fprintf (fptr,"RM[%d] = %d\n", l, RM[l]);

}

/* ***************************************************** */

initial_generation ()
{
	/***** Files and columns of individuals *****/

	i=(-1);
	for (fil=0; fil<sqrt(NIND); fil++)
	for (col=0; col<sqrt(NIND); col++)
	{
		i ++;
		f[i] = fil;
		c[i] = col;
//		if (tracelevel!=0)	fprintf(fptr,"i = %d    f = %d     c = %d\n", i, f[i], c[i]);
	}

	/* ***** Initialise the first generation ***** */

	//NOT SEGREGATING
	for (i=0; i<NIND; i++)
	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	{
	    	gm[i][k][0]=0;
	    	gm[i][k][1]=0;
	}

/*	if (tracelevel!=0)
	{
		fprintf (fptr,"gm0 = %d\n", gm[1][1][0]);

		for (k=0; k<NCRO; k++)
		{
			for (l=0; l<NLOCI; l++)
			{
				if((gm[0][k][0]&RM[l])==RM[l])		fprintf (fptr,"1 ");
				else				fprintf (fptr,"0 ");
			}
			fprintf (fptr,"\n");
		}
	}
*/
}

/* ***************************************************** */

frequency_genes ()
{
	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	{
		AA=0.0; Aa=0.0; aa=0.0;

		for (i=0; i<NIND; i++)
		{
			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))	aa+=1.0;
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	AA+=1.0;
		     	else	Aa+=1.0;
		}

		q[k][l] = (aa/(double)NIND)+(Aa/(2.0*(double)NIND));
//		if (tracelevel!=0)    fprintf(fptr,"\n q[%d][%d]=%f", k, l, q[k][l]);
	}
}

/* ***************************************************** */

mutation_neutral ()
{
	//if (tracelevel!=0)    fprintf(fptr,"\n New neutral mutants\n");

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	{
		if (q[k][l] == 0.0)
		{
			ran_i = (int)(uniform()*NIND);
			ran_h = (int)(uniform()*2.0);
			gm[ran_i][k][ran_h]=(gm[ran_i][k][ran_h] | RM[l]);
		}
		if (q[k][l] == 1.0)
		{
			ran_i = (int)(uniform()*NIND);
			ran_h = (int)(uniform()*2.0);
			gm[ran_i][k][ran_h]=(gm[ran_i][k][ran_h] & (~RM[l]));
		}
	}
}

/* **************************************************************** */

plink_files()
{
	ffped = fopen ("data.ped","w");
	ffpedR = fopen ("dataR.ped","w");
	ffmap = fopen ("data.map","w");

	// data.map

	SNP = 0;
	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	if ((q[k][l] != 0.0)&&(q[k][l] != 1.0))
	{
		SNP ++;
		for (crom=0; crom<20; crom++)
		{
			if ((k >= crom*34) && (k < (crom+1)*34))
			{
				chrom[k][l] = crom +  1;
				pos[k][l] = ( (l + NLOCI*k) - (1020*(chrom[k][l]-1)) ) * 100000;
			}
		}

		fprintf(ffmap,"%d SNP%d %f %llu\n", chrom[k][l], SNP, pos[k][l]/1000000.0, pos[k][l]);
	}

	for (Sfil=((sqrt(NIND)/2)-(SW/2)); Sfil<=(sqrt(NIND)/2)+(SW/2); Sfil++)
	for (Scol=((sqrt(NIND)/2)-(SW/2)); Scol<=(sqrt(NIND)/2)+(SW/2); Scol++)
	if ((Sfil >= 0) && (Scol >= 0) && (Sfil < sqrt(NIND)) && (Scol < sqrt(NIND)))
	{
		i = (Sfil*sqrt(NIND)) + Scol;

		// data.ped

		fprintf(ffped,"1 IND%d 0 0 1 -9 ", i);
		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		if ((q[k][l] != 0.0)&&(q[k][l] != 1.0))
		{
			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))		fprintf(ffped,"T T ");
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	fprintf(ffped,"A A ");
	    		else if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	fprintf(ffped,"T A ");
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))	fprintf(ffped,"A T ");
		}
		fprintf(ffped,"\n");
	}

	/* **************** ped file with five dispersed SW - dataR.ped ***************** */

	// ****** Upper left corner ******* */

	for (Sfil=0; Sfil<SW; Sfil++)
	for (Scol=0; Scol<SW; Scol++)
	if ((Sfil >= 0) && (Scol >= 0) && (Sfil < sqrt(NIND)) && (Scol < sqrt(NIND)))
	{
		i = (Sfil*sqrt(NIND)) + Scol;

		// dataR.ped

		fprintf(ffpedR,"1 INDUL%d 0 0 1 -9 ", i);
		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		if ((q[k][l] != 0.0)&&(q[k][l] != 1.0))
		{
			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))		fprintf(ffpedR,"T T ");
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	fprintf(ffpedR,"A A ");
	    		else if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	fprintf(ffpedR,"T A ");
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))	fprintf(ffpedR,"A T ");
		}
		fprintf(ffpedR,"\n");
	}

	// ****** Lower left corner ******* */

	for (Sfil=0; Sfil<SW; Sfil++)
	for (Scol=(sqrt(NIND)-SW); Scol<sqrt(NIND); Scol++)
	if ((Sfil >= 0) && (Scol >= 0) && (Sfil < sqrt(NIND)) && (Scol < sqrt(NIND)))
	{
		i = (Sfil*sqrt(NIND)) + Scol;

		// dataR.ped

		fprintf(ffpedR,"1 INDLL%d 0 0 1 -9 ", i);
		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		if ((q[k][l] != 0.0)&&(q[k][l] != 1.0))
		{
			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))		fprintf(ffpedR,"T T ");
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	fprintf(ffpedR,"A A ");
	    		else if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	fprintf(ffpedR,"T A ");
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))	fprintf(ffpedR,"A T ");
		}
		fprintf(ffpedR,"\n");
	}

	// ****** Upper right corner ******* */

	for (Sfil=(sqrt(NIND)-SW); Sfil<sqrt(NIND); Sfil++)
	for (Scol=0; Scol<SW; Scol++)
	if ((Sfil >= 0) && (Scol >= 0) && (Sfil < sqrt(NIND)) && (Scol < sqrt(NIND)))
	{
		i = (Sfil*sqrt(NIND)) + Scol;

		// dataR.ped

		fprintf(ffpedR,"1 INDUR%d 0 0 1 -9 ", i);
		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		if ((q[k][l] != 0.0)&&(q[k][l] != 1.0))
		{
			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))		fprintf(ffpedR,"T T ");
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	fprintf(ffpedR,"A A ");
	    		else if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	fprintf(ffpedR,"T A ");
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))	fprintf(ffpedR,"A T ");
		}
		fprintf(ffpedR,"\n");
	}

	// ****** Lower right corner ******* */

	for (Sfil=(sqrt(NIND)-SW); Sfil<sqrt(NIND); Sfil++)
	for (Scol=(sqrt(NIND)-SW); Scol<sqrt(NIND); Scol++)
	if ((Sfil >= 0) && (Scol >= 0) && (Sfil < sqrt(NIND)) && (Scol < sqrt(NIND)))
	{
		i = (Sfil*sqrt(NIND)) + Scol;

		// dataR.ped

		fprintf(ffpedR,"1 INDLR%d 0 0 1 -9 ", i);
		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		if ((q[k][l] != 0.0)&&(q[k][l] != 1.0))
		{
			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))		fprintf(ffpedR,"T T ");
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	fprintf(ffpedR,"A A ");
	    		else if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	fprintf(ffpedR,"T A ");
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))	fprintf(ffpedR,"A T ");
		}
		fprintf(ffpedR,"\n");
	}

	// *********** Centre ************ */

	for (Sfil=((sqrt(NIND)/2)-(SW/2)); Sfil<=(sqrt(NIND)/2)+(SW/2); Sfil++)
	for (Scol=((sqrt(NIND)/2)-(SW/2)); Scol<=(sqrt(NIND)/2)+(SW/2); Scol++)
	if ((Sfil >= 0) && (Scol >= 0) && (Sfil < sqrt(NIND)) && (Scol < sqrt(NIND)))
	{
		i = (Sfil*sqrt(NIND)) + Scol;

		// dataR.ped

		fprintf(ffpedR,"1 INDC%d 0 0 1 -9 ", i);
		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		if ((q[k][l] != 0.0)&&(q[k][l] != 1.0))
		{
			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))		fprintf(ffpedR,"T T ");
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	fprintf(ffpedR,"A A ");
	    		else if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	fprintf(ffpedR,"T A ");
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))	fprintf(ffpedR,"A T ");
		}
		fprintf(ffpedR,"\n");
	}

}

/* ************************************************************ */

initial_multiallelic ()
{
	/* ***** MULTIALLELIC GENES ***** */

	for (k=0; k<NCRO; k++)
	for (i=0; i<NIND; i++)
	{
		locus[i][k][0] = i;
		locus[i][k][1] = i+100000;

//		if (tracelevel!=0)	if (k==0) fprintf (fptr,"locus=%d ind=%d   %d %d\n", k, i, locus[i][k][0], locus[i][k][1]);
	}
}

/* ***************************************************** */

estimates_of_F ()
{
	numLOCI[gen] = 0.0;
	EHom = 0.0;
	EHet = 0.0;

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	if ((q[k][l] != 0.0)&&(q[k][l] != 1.0))
	{
		numLOCI[gen] ++;
		EHom += (1.0 - 2.0*q[k][l]*(1.0-q[k][l]));
		EHet += 2.0*q[k][l]*(1.0-q[k][l]);
	}
	accum(&Het[gen], EHet/numLOCI[gen]);

	for (i=0; i<NIND; i++)
	{
		// MULTIALLELIC GENES

		fibd[i] = 0.0;
		for (k=0; k<NCRO; k++)
		{
			if (locus[i][k][0] == locus[i][k][1])
			{
				fibd[i] ++;
//				if (tracelevel!=0)	if (k==0) fprintf(fptr," k=%d    locus[i][k][0]=%d locus[i][k][1]=%d \n", k, locus[i][k][0], locus[i][k][1]);
			}
		}
		fibd[i] = fibd[i] / NCRO;
//		if (tracelevel!=0)	fprintf(fptr," ******* fibd[i] = %f ******** \n", fibd[i]);

		Hom[i] = 0.0;

		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		if ((q[k][l] != 0.0)&&(q[k][l] != 1.0))
		{
//			if (tracelevel!=0) if((k==0)&&(l<5))   fprintf(fptr,"\n q[%d][%d]=%f   EHom=%f", k, l, q[k][l], EHom);

		    	if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))	Hom[i] ++;
			else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	Hom[i] ++;
		}

		fLH1[i] = (Hom[i]-EHom) / (numLOCI[gen]-EHom);

		accum(&Fibd[gen], fibd[i]);
		accum(&FLH1[gen], fLH1[i]);
	}
}

/* ***************************************************** */

estimates_of_Fst_between_opposite_NS ()
{
	/********* Frequencies of two extreme neighborhoods NS1 and NS2 **********/;

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	{
		/***** NS1_NS2 *****/;

		AA=0.0; Aa=0.0; aa=0.0;

		for (fil=0; fil<BW; fil++)
		for (col=0; col<BW; col++)
		{
			i = fil*sqrt(NIND) + col;

			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))	aa+=1.0;
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	AA+=1.0;
		     	else	Aa+=1.0;
		}

		for (fil=(sqrt(NIND)-BW); fil<sqrt(NIND); fil++)
		for (col=(sqrt(NIND)-BW); col<sqrt(NIND); col++)
		{
			i = fil*sqrt(NIND) + col;

			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))	aa+=1.0;
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	AA+=1.0;
		     	else	Aa+=1.0;
		}

		qNS1NS2[k][l] = (aa/(2*BW*BW))+(Aa/(2.0*(2*BW*BW)));

		/***** NS1 *****/;

		AA=0.0; Aa=0.0; aa=0.0;

		for (fil=0; fil<BW; fil++)
		for (col=0; col<BW; col++)
		{
			i = fil*sqrt(NIND) + col;

			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))	aa+=1.0;
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	AA+=1.0;
		     	else	Aa+=1.0;
		}

		qNS1[k][l] = (aa/(BW*BW))+(Aa/(2.0*(BW*BW)));

		/***** NS2 *****/;

		AA=0.0; Aa=0.0; aa=0.0;

		for (fil=(sqrt(NIND)-BW); fil<sqrt(NIND); fil++)
		for (col=(sqrt(NIND)-BW); col<sqrt(NIND); col++)
		{
			i = fil*sqrt(NIND) + col;

			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))	aa+=1.0;
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	AA+=1.0;
		     	else	Aa+=1.0;
		}

		qNS2[k][l] = (aa/(BW*BW))+(Aa/(2.0*(BW*BW)));
		if (tracelevel!=0) if (k < 100)   fprintf(fptr,"\nk=%d   l=%d    qNS1=%f    qNS2=%f    qNS1NS2=%f", k, l, qNS1[k][l], qNS2[k][l], qNS1NS2[k][l]);
	}

	/********* Frequencies of heterozygotes **********/;

	/***** NS1_NS2 *****/;

	EHetNS1NS2 = 0.0;
	lociNS1NS2 = 0.0;

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	if ((qNS1NS2[k][l] != 0.0)&&(qNS1NS2[k][l] != 1.0))
	{
		lociNS1NS2 ++;
		EHetNS1NS2 += 2.0*qNS1NS2[k][l]*(1.0-qNS1NS2[k][l]);
	}

	HetNS1NS2 = 0.0;

	for (fil=0; fil<BW; fil++)
	for (col=0; col<BW; col++)
	{
		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		if ((qNS1NS2[k][l] != 0.0)&&(qNS1NS2[k][l] != 1.0))
		{
		    	if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	HetNS1NS2 ++;
			else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))	HetNS1NS2 ++;
		}
	}

	for (fil=(sqrt(NIND)-BW); fil<sqrt(NIND); fil++)
	for (col=(sqrt(NIND)-BW); col<sqrt(NIND); col++)
	{
		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		if ((qNS1NS2[k][l] != 0.0)&&(qNS1NS2[k][l] != 1.0))
		{
		    	if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	HetNS1NS2 ++;
			else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))	HetNS1NS2 ++;
		}

//		if (tracelevel!=0) fprintf(fptr,"\nHetNS1NS2=%f\n", HetNS1NS2);
	}

	/***** NS1 *****/;

	EHetNS1 = 0.0;

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	if ((qNS1NS2[k][l] != 0.0)&&(qNS1NS2[k][l] != 1.0))
	{
		EHetNS1 += 2.0*qNS1[k][l]*(1.0-qNS1[k][l]);
	}

	HetNS1 = 0.0;

	for (fil=0; fil<BW; fil++)
	for (col=0; col<BW; col++)
	{
		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		if ((qNS1NS2[k][l] != 0.0)&&(qNS1NS2[k][l] != 1.0))
		{
		    	if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	HetNS1 ++;
			else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))	HetNS1 ++;
		}
//		if (tracelevel!=0) fprintf(fptr,"\nHetNS1=%f\n", HetNS1);
	}

	/***** NS2 *****/;

	EHetNS2 = 0.0;

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	if ((qNS1NS2[k][l] != 0.0)&&(qNS1NS2[k][l] != 1.0))
	{
		EHetNS2 += 2.0*qNS2[k][l]*(1.0-qNS2[k][l]);
	}

	HetNS2 = 0.0;

	for (fil=(sqrt(NIND)-BW); fil<sqrt(NIND); fil++)
	for (col=(sqrt(NIND)-BW); col<sqrt(NIND); col++)
	{
		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		if ((qNS1NS2[k][l] != 0.0)&&(qNS1NS2[k][l] != 1.0))
		{
		    	if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	HetNS2 ++;
			else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))	HetNS2 ++;
		}
		if (tracelevel!=0) fprintf(fptr,"\nHetNS2=%f\n", HetNS2);
	}

	HT = EHetNS1NS2 / lociNS1NS2;
	HI = ((HetNS1/(BW*BW)) + (HetNS2/(BW*BW))) / (2.0*lociNS1NS2);
	HS = (EHetNS1+EHetNS2)/(2.0*lociNS1NS2);
	FSTmax = ( HT - HS ) / HT;
	FISmax = ( HS - HI ) / HS;
	FITmax = ( HT - HI ) / HT;
}

/* ***************************************************** */

estimates_of_Fst_between_close_NS ()
{
	/********* Frequencies of two closed neighborhoods NS3 and NS4 **********/;

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	{
		/***** NS3_NS4 *****/;

		AA=0.0; Aa=0.0; aa=0.0;

		for (fil=((sqrt(NIND)/2)-(BW/2)); fil<=(sqrt(NIND)/2)+(BW/2); fil++)
		for (col=((sqrt(NIND)/2)-(BW/2)); col<=(sqrt(NIND)/2)+(BW/2); col++)
		{
			i = fil*sqrt(NIND) + col;

			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))	aa+=1.0;
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	AA+=1.0;
		     	else	Aa+=1.0;
		}

		for (fil=((sqrt(NIND)/2)+BW-(BW/2)); fil<=(sqrt(NIND)/2)+BW+(BW/2); fil++)
		for (col=((sqrt(NIND)/2)+BW-(BW/2)); col<=(sqrt(NIND)/2)+BW+(BW/2); col++)
		{
			i = fil*sqrt(NIND) + col;

			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))	aa+=1.0;
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	AA+=1.0;
		     	else	Aa+=1.0;
		}

		qNS3NS4[k][l] = (aa/(2*BW*BW))+(Aa/(2.0*(2*BW*BW)));

		/***** NS3 *****/;

		AA=0.0; Aa=0.0; aa=0.0;

		for (fil=((sqrt(NIND)/2)-(BW/2)); fil<=(sqrt(NIND)/2)+(BW/2); fil++)
		for (col=((sqrt(NIND)/2)-(BW/2)); col<=(sqrt(NIND)/2)+(BW/2); col++)
		{
			i = fil*sqrt(NIND) + col;

			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))	aa+=1.0;
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	AA+=1.0;
		     	else	Aa+=1.0;
		}

		qNS3[k][l] = (aa/(BW*BW))+(Aa/(2.0*(BW*BW)));

		/***** NS4 *****/;

		AA=0.0; Aa=0.0; aa=0.0;

		for (fil=((sqrt(NIND)/2)+BW-(BW/2)); fil<=(sqrt(NIND)/2)+BW+(BW/2); fil++)
		for (col=((sqrt(NIND)/2)+BW-(BW/2)); col<=(sqrt(NIND)/2)+BW+(BW/2); col++)
		{
			i = fil*sqrt(NIND) + col;

			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))	aa+=1.0;
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	AA+=1.0;
		     	else	Aa+=1.0;
		}

		qNS4[k][l] = (aa/(BW*BW))+(Aa/(2.0*(BW*BW)));
		if (tracelevel!=0) if (k < 100)   fprintf(fptr,"\nk=%d   l=%d    qNS3=%f    qNS4=%f    qNS3NS4=%f", k, l, qNS3[k][l], qNS4[k][l], qNS3NS4[k][l]);
	}

	/********* Frequencies of heterozygotes **********/;

	/***** NS3_NS4 *****/;

	EHetNS3NS4 = 0.0;
	lociNS3NS4 = 0.0;

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	if ((qNS3NS4[k][l] != 0.0)&&(qNS3NS4[k][l] != 1.0))
	{
		lociNS3NS4 ++;
		EHetNS3NS4 += 2.0*qNS3NS4[k][l]*(1.0-qNS3NS4[k][l]);
	}

	HetNS3NS4 = 0.0;

	for (fil=((sqrt(NIND)/2)-(BW/2)); fil<=(sqrt(NIND)/2)+(BW/2); fil++)
	for (col=((sqrt(NIND)/2)-(BW/2)); col<=(sqrt(NIND)/2)+(BW/2); col++)
	{
		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		if ((qNS3NS4[k][l] != 0.0)&&(qNS3NS4[k][l] != 1.0))
		{
		    	if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	HetNS3NS4 ++;
			else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))	HetNS3NS4 ++;
		}
	}

	for (fil=((sqrt(NIND)/2)+BW-(BW/2)); fil<=(sqrt(NIND)/2)+BW+(BW/2); fil++)
	for (col=((sqrt(NIND)/2)+BW-(BW/2)); col<=(sqrt(NIND)/2)+BW+(BW/2); col++)
	{
		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		if ((qNS3NS4[k][l] != 0.0)&&(qNS3NS4[k][l] != 1.0))
		{
		    	if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	HetNS3NS4 ++;
			else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))	HetNS3NS4 ++;
		}

//		if (tracelevel!=0) fprintf(fptr,"\nHetNS3NS4=%f\n", HetNS3NS4);
	}

	/***** NS3 *****/;

	EHetNS3 = 0.0;

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	if ((qNS3NS4[k][l] != 0.0)&&(qNS3NS4[k][l] != 1.0))
	{
		EHetNS3 += 2.0*qNS3[k][l]*(1.0-qNS3[k][l]);
	}

	HetNS3 = 0.0;

	for (fil=((sqrt(NIND)/2)-(BW/2)); fil<=(sqrt(NIND)/2)+(BW/2); fil++)
	for (col=((sqrt(NIND)/2)-(BW/2)); col<=(sqrt(NIND)/2)+(BW/2); col++)
	{
		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		if ((qNS3NS4[k][l] != 0.0)&&(qNS3NS4[k][l] != 1.0))
		{
		    	if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	HetNS3 ++;
			else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))	HetNS3 ++;
		}
//		if (tracelevel!=0) fprintf(fptr,"\nHetNS3=%f\n", HetNS3);
	}

	/***** NS4 *****/;

	EHetNS4 = 0.0;

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	if ((qNS3NS4[k][l] != 0.0)&&(qNS3NS4[k][l] != 1.0))
	{
		EHetNS4 += 2.0*qNS4[k][l]*(1.0-qNS4[k][l]);
	}

	HetNS4 = 0.0;

	for (fil=((sqrt(NIND)/2)+BW-(BW/2)); fil<=(sqrt(NIND)/2)+BW+(BW/2); fil++)
	for (col=((sqrt(NIND)/2)+BW-(BW/2)); col<=(sqrt(NIND)/2)+BW+(BW/2); col++)
	{
		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		if ((qNS3NS4[k][l] != 0.0)&&(qNS3NS4[k][l] != 1.0))
		{
		    	if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	HetNS4 ++;
			else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))	HetNS4 ++;
		}
		if (tracelevel!=0) fprintf(fptr,"\nHetNS4=%f\n", HetNS4);
	}

	HT = EHetNS3NS4 / lociNS3NS4;
	HI = ((HetNS3/(BW*BW)) + (HetNS4/(BW*BW))) / (2.0*lociNS3NS4);
	HS = (EHetNS3+EHetNS4)/(2.0*lociNS3NS4);
	FSTmin = ( HT - HS ) / HT;
	FISmin = ( HS - HI ) / HS;
	FITmin = ( HT - HI ) / HT;
}

/* ***************************************************** */

mating ()
{
	int a, p1, p2, EE[MM], FF[MM], family[NN];
	int pointrec[MM][31], ncrorec[MM], rndk, rndl;
	int rndfil, rndcol;
	double rnd, family_sum2=0.0;
	double fit2, fit3;

	for (i=0; i<NIND; i++)
	for (k=0; k<NCRO; k++)
	{
		sm[i][k][0]=gm[i][k][0];
		sm[i][k][1]=gm[i][k][1];

		// MULTIALLELIC GENES
		slocus[i][k][0] = locus[i][k][0];
		slocus[i][k][1] = locus[i][k][1];
	}

//	if (tracelevel!=0)	fprintf(fptr,"\n Parents \n");

	/***** Fitness (mortality) of parents *****/

	if (Vk == 0) 	{ fit2 = 0.0; fit3 = 0.0; }
	else 		{ fit2 = 0.666; fit3 = 0.9; }
	
	/***** Generation of offspring *****/

	for (i=0; i<NIND;i++)   family[i]=0;

	for (i=0; i<NIND; i++)
	{
		if (BW == 99)
		{
			/***** parents full RM *****/

			labelp1bw99: /* Repeat p1 */;
			p1 = (int)(uniform()*NIND);
			rnd = uniform();
			if ((p1%2==0)&&(rnd<fit2))	goto labelp1bw99;
			if ((p1%3==0)&&(rnd<fit3))	goto labelp1bw99;

			labelp2bw99: /* Repeat p2 */;
			do { p2 = (int)(uniform()*NIND); }
			while (p2 == p1);
			rnd = uniform();
			if ((p2%2==0)&&(rnd<fit2))	goto labelp2bw99;
			if ((p2%3==0)&&(rnd<fit3))	goto labelp2bw99;
		}
		else
		{
			/***** parents from BW *****/

			do
			{
				labelp1: /* Repeat p1 */;

				do { rndfil = (int)(uniform()*sqrt(NIND)); }
				while( (rndfil < (f[i] - (BW/2))) || (rndfil > (f[i] + (BW/2))) || (rndfil < 0) || (rndfil >= sqrt(NIND)) );

				do { rndcol = (int)(uniform()*sqrt(NIND)); }
				while( (rndcol < (c[i] - (BW/2))) || (rndcol > (c[i] + (BW/2))) || (rndcol < 0) || (rndcol >= sqrt(NIND)) );

				p1 = (rndfil *sqrt(NIND)) + rndcol;

				rnd = uniform();
				if ((p1%2==0)&&(rnd<fit2))	goto labelp1;
				if ((p1%3==0)&&(rnd<fit3))	goto labelp1;
			}
			while (p1 == i);
//			if (tracelevel!=0)	fprintf(fptr,"i = %d    rndfil = %d    rnd col = %d   p1 = %d \n", i, rndfil, rndcol, p1);

			do
			{
				labelp2: /* Repeat p2 */;

				do { rndfil = (int)(uniform()*sqrt(NIND)); }
				while( (rndfil < (f[i] - (BW/2))) || (rndfil > (f[i] + (BW/2))) || (rndfil < 0) || (rndfil >= sqrt(NIND)) );

				do { rndcol = (int)(uniform()*sqrt(NIND)); }
				while( (rndcol < (c[i] - (BW/2))) || (rndcol > (c[i] + (BW/2))) || (rndcol < 0) || (rndcol >= sqrt(NIND)) );

				p2 = (rndfil *sqrt(NIND)) + rndcol;

				rnd = uniform();
				if ((p2%2==0)&&(rnd<fit2))	goto labelp2;
				if ((p2%3==0)&&(rnd<fit3))	goto labelp2;
			}
			while ((p2 == i) || (p2 == p1));				
//			if (tracelevel!=0)	fprintf(fptr,"i = %d    rndfil = %d    rnd col = %d   p1 = %d \n", i, rndfil, rndcol, p1);
		}

//		if (tracelevel!=0)   fprintf (fptr,"%d\t%d\n", p1, p2);

	    	family[p1]++;    family[p2]++;

		/* ************** Restricted recombination ***************** */

		/* ****** Chromosome from father ****** */

		for (k=0; k<NCRO; k++)
		{
			ncrorec[k] = 0;
			for (l=0; l<NLOCI; l++)  pointrec[k][l] = 0;
		}

		// SEGREGATION CHROMOSOMES
 
//		if ((tracelevel!=0)&&(i==31))	fprintf (fptr,"SEGREGATION OF CHROMOSOMES\n");

		for (crom=0; crom<20; crom++)
		{
			ncrorec[(crom)*34] = 1;
			pointrec[(crom)*34][0] = 1;
		}

		// ONE CROSSINGOVER PER CHROMOSOME

		for (crom=0; crom<20; crom++)
		{
			rndk = (int)(uniform()*34);
			do { rndl = (int)(uniform()*NLOCI); }
			while (rndl == 0);

			ncrorec[(crom*34)+rndk] = 1;
			pointrec[(crom*34)+rndk][rndl] = 1;

//			if (tracelevel!=0)	fprintf (fptr,"Chrom %d rndk=%d  rndl=%d\n", crom, rndk, rndl);
		}

/*		if (tracelevel!=0)
		if (gen==0)
		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		{
			if ((ncrorec[k] == 1)&&(pointrec[k][l] == 1))	fprintf (fptr,"k=%d  l=%d\n", k, l);	
		}
*/
		marker = 1;

		for (k=0; k<NCRO; k++)
		{
			EE[k]=0;
			if (ncrorec[k] == 0)
			{
				if (marker==(-1))
				{
					EE[k] = ~EE[k];
				}
			}
			else
			{
				for (l=0; l<NLOCI; l++)
		      	{
					if (pointrec[k][l] == 0)
					{
						if (marker==(-1))  EE[k] = EE[k] | RM[l];
					}
					else
					{
						if (marker==1)
						{
							EE[k] = EE[k] | RM[l];
							marker = marker * (-1);
						}
						else
						{
							marker = marker * (-1);
						}
					}
				}
			}
		}

		rnd = uniform();
		for (k=0; k<NCRO; k++)
		{
			if (rnd < 0.5)	EE[k] = ~EE[k];
			FF[k] = ~EE[k];
			gm[i][k][0]=((EE[k]&sm[p1][k][0])|(FF[k]&sm[p1][k][1]));

			// MULTIALLELIC GENES
			if ((EE[k] & RM[0]) == RM[0])	locus[i][k][0] = slocus[p1][k][0];
			else						locus[i][k][0] = slocus[p1][k][1];

//			if (tracelevel!=0)	if (k<=3) fprintf(fptr," k=%d    i=%d    p1=%d    slocus[p1][k][0]=%d slocus[p1][k][1]=%d locus[i][k][0]=%d \n", k, i, p1, slocus[p1][k][0], slocus[p1][k][1], locus[i][k][0]);
		}
/*
		if (tracelevel!=0)
		if (i == 0)
		{
			fprintf (fptr,"EE\n");
			for (k=0; k<NCRO; k++)
			{
				for (l=0; l<NLOCI; l++)
				{
					if((EE[k]&RM[l])==RM[l])  	fprintf (fptr,"1 ");
					else			   	fprintf (fptr,"0 ");
				}
				fprintf (fptr,"\n");
			}
			fprintf (fptr,"FF\n");
			for (k=0; k<NCRO; k++)
			{
				for (l=0; l<NLOCI; l++)
				{
					if((FF[k]&RM[l])==RM[l])  	fprintf (fptr,"1 ");
					else			   	fprintf (fptr,"0 ");
				}
				fprintf (fptr,"\n");
			}
		}
*/
		/* ****** Chromosome from mother ****** */

		for (k=0; k<NCRO; k++)
		{
			ncrorec[k] = 0;
			for (l=0; l<NLOCI; l++)  pointrec[k][l] = 0;
		}

		for (crom=0; crom<20; crom++)
		{
			ncrorec[(crom)*34] = 1;
			pointrec[(crom)*34][0] = 1;
		}

		// ONE CROSSINGOVER PER CHROMOSOME

		for (crom=0; crom<20; crom++)
		{
			rndk = (int)(uniform()*34);
			do { rndl = (int)(uniform()*NLOCI); }
			while (rndl == 0);

			ncrorec[(crom*34)+rndk] = 1;
			pointrec[(crom*34)+rndk][rndl] = 1;

//			if (tracelevel!=0)	fprintf (fptr,"Chrom %d rndk=%d  rndl=%d\n", crom, rndk, rndl);
		}

		marker = 1;

		for (k=0; k<NCRO; k++)
		{
			EE[k]=0;
			if (ncrorec[k] == 0)
			{
				if (marker==(-1))
				{
					EE[k] = ~EE[k];
				}
			}
			else
			{
				for (l=0; l<NLOCI; l++)
		      	{
					if (pointrec[k][l] == 0)
					{
						if (marker==(-1))  EE[k] = EE[k] | RM[l];
					}
					else
					{
						if (marker==1)
						{
							EE[k] = EE[k] | RM[l];
							marker = marker * (-1);
						}
						else
						{
							marker = marker * (-1);
						}
					}
				}
			}
		}

		rnd = uniform();
		for (k=0; k<NCRO; k++)
		{
			if (rnd < 0.5)	EE[k] = ~EE[k];
			FF[k] = ~EE[k];
			gm[i][k][1]=((EE[k]&sm[p2][k][0])|(FF[k]&sm[p2][k][1]));

			// MULTIALLELIC GENES
			if ((EE[k] & RM[0]) == RM[0])	locus[i][k][1] = slocus[p2][k][0];
			else				locus[i][k][1] = slocus[p2][k][1];

//			if (tracelevel!=0)	if (k<=3) fprintf(fptr," k=%d    i=%d    p2=%d    slocus[p2][k][0]=%d slocus[p2][k][1]=%d locus[i][k][1]=%d \n", k, i, p2, slocus[p2][k][0], slocus[p2][k][1], locus[i][k][1]);
//			if (tracelevel!=0)	if (k<=3) if ((EE[k] & RM[0]) == RM[0]) fprintf(fptr,"EE[%d]=1\n", k);
		}
	}

	/* ***** VARIANCE OF FAMILY SIZE ***** */

	for (i=0; i<NIND; i++)	family_sum2 += (family[i] * family[i]);
	accum(&SK2, ((family_sum2-(4.0*NIND))/(NIND-1.0)) );
	//if (tracelevel!=0)	fprintf (fptr,"SK2 = %5.3f\n", ((family_sum2-(4.0*NIND))/(NIND-1.0)));
}

/* ***************************************************** */

printout()
{
	int g;
	double b;

	NeFibd = (10.0) / ( 2.0*(accmean(&Fibd[generations])-accmean(&Fibd[generations-10]))/(1.0-accmean(&Fibd[generations-10])) );
	NeHet = 1.0 / ( 1.0 - exp(10.0 * log(accmean(&Het[generations-10]) / accmean(&Het[generations]))) );

	fprintf(fgen, "\ngen     LOCI     Het            Fibd       FLH1\n");
	for(gen=generations-10; gen<=generations; gen++)
	{
		fprintf(fgen, "%3d  %d  %f  %6.4f  %6.4f\n", gen, numLOCI[gen], accmean(&Het[gen]), accmean(&Fibd[gen]), accmean(&FLH1[gen]));
		accum (&GEN, gen);
		covaccum (&GENHET, (double)gen, accmean(&Het[gen]));
		accum(&AVE_FLH1, accmean(&FLH1[gen]));
//		fprintf(fgen, "\ngen = %d  Het = %f  \n", gen, accmean(&Het[gen]));
	}
	b = covariance(&GENHET) / variance(&GEN);
	NeHetb = 1.0/(-2.0*b);
	NeSk2 = (4.0*NIND) / (2.0 + accmean(&SK2));

	fprintf(fgen,"\nN=%d  N.LOCI=%d  gens=%d\n", NIND, NCRO*NLOCI, generations);
	fprintf(fgen, "\nNeFibd = %f  NeHetb = %f  Sk2 = %f  NeSk2 = %f  FLH1 = %f\n", NeFibd, NeHetb, accmean(&SK2), NeSk2, accmean(&AVE_FLH1));
	fprintf(fgen, "\nlociNS1NS2 = %f  FISmax = %f  FSTmax = %f  FITmax = %f\n", lociNS1NS2, FISmax, FSTmax, FITmax);
	fprintf(fgen, "\nlociNS3NS4 = %f  FISmin = %f  FSTmin = %f  FITmin = %f\n", lociNS3NS4, FISmin, FSTmin, FITmin);
}

/* ***************************************************** */

