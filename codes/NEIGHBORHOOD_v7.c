/* NEIGHBORHOOD_v7.c (08/09/2026) print every 50000 gens and get 2000 progeny */

/* ***************************************************** */

#include "libhdr"
#include <stdio.h>
#include <stdlib.h> 
#include <string.h> 
#include <math.h>

#define NN 90001  /* max number of NIND */
#define GG 1000001  /* max number of GEN */
#define MM 701  /* max number of NCRO */
#define NSAMPLE 2000   /* número de hijos adicionales para el nuevo PED */

int BW, SW, SWR, ran_i, ran_h, fil, col, Sfil, Scol, Vk;
int marker, NIND, NCRO, NLOCI, SNP;
static int numLOCI[GG];
int gen, generations, genmulti, i, k, l;
int RM[31], p1, p2;

int (*gm)[MM][2] = NULL;
int (*sm)[MM][2] = NULL;
int (*locus)[MM][2] = NULL;
int (*slocus)[MM][2] = NULL;
/* Genotipos de los 2000 hijos adicionales */
int (*sample_gm)[MM][2] = NULL;

int crom;
static int chrom[MM][31];
static unsigned long long int pos[MM][31];

double AA, Aa, aa;
static double q[MM][31];
static double fibd[NN], fLH1[NN];
static double EHom, EHet;
static double Hom[NN];
double NeFibd, NeHet, NeHetb, NeSk2;
double FISmax, FSTmax, FITmax;
double FISmin, FSTmin, FITmin;

static double qNS1[MM][31], qNS2[MM][31], qNS1NS2[MM][31];
static double qNS3[MM][31], qNS4[MM][31], qNS3NS4[MM][31];

double EHetNS1, EHetNS2, EHetNS1NS2;
double EHetNS3, EHetNS4, EHetNS3NS4;
double HT, HS, HI, HetNS1, HetNS2, HetNS1NS2, lociNS1NS2;
double HetNS3, HetNS4, HetNS3NS4, lociNS3NS4;

int all, num_a;
static int qallele[2 * NN];
double SUM_q, SUM_q2;  

static int f[NN], c[NN];

static double qGroups[13][MM][31];
static double countGroups[13];

int generate_2000_sample(void);

double FST_matrix[13][13];
double FIS_matrix[13][13];
double FIT_matrix[13][13];
double FST_global_13;
double FIS_global_13;
double FIT_global_13;

struct acc Num_alleles[GG], Mean_q[GG], Var_q[GG];
struct acc Het[GG], Fibd[GG], FLH1[GG], AVE_FLH1, GEN, SK2;
struct covacc GENHET;

FILE *fptr, *frep, *fgen, *ffmap, *ffmapR, *ffped, *ffpedR;

/* ***************************************************** */

main()
{
	gm = malloc(sizeof(*gm) * NN);
	sm = malloc(sizeof(*sm) * NN);
	locus = malloc(sizeof(*locus) * NN);
	slocus = malloc(sizeof(*slocus) * NN);

	sample_gm = malloc(sizeof(*sample_gm) * NSAMPLE);

	if (gm == NULL || sm == NULL || locus == NULL || slocus == NULL ||
    		sample_gm == NULL)
	{
	    printf("No hay memoria suficiente.\n");
	    return 1;
	}

	memset(gm, 0, sizeof(*gm) * NN);
	memset(sm, 0, sizeof(*sm) * NN);
	memset(locus, 0, sizeof(*locus) * NN);
	memset(slocus, 0, sizeof(*slocus) * NN);
	memset(sample_gm, 0, sizeof(*sample_gm) * NSAMPLE);

	fptr = fopen ("dfilename.dat","w");
	fgen = fopen ("genfile.dat","w");
	frep = fopen ("repfile.dat","w");

	getinputs ();
	recombination_masks ();

	initial_generation ();

	for (gen=0; gen<=generations; gen++)
	{
		if (gen % 1000 == 0)
		{
			fprintf(frep, "GEN = %d\n", gen);
			fflush(frep);
		}

		frequency_genes ();
		mutation_neutral ();
		
		if (gen > 0 && gen % 50000 == 0)
		{
			plink_files();
		}
		else if (gen == generations && (gen % 50000 != 0))
		{
			plink_files();
		}

		if (gen == generations-genmulti)   initial_multiallelic();
		if (gen >= generations-genmulti)   estimates_of_F();
		if (gen == generations)        estimates_of_Fst_between_opposite_NS();
		if (gen == generations)        estimates_of_Fst_between_close_NS();
		if (gen == generations)        estimates_of_Fst_between_13_groups();

		mating ();
	}

	fclose(frep);
 	printout();
	
	fclose(fptr);
	fclose(fgen);

	free(gm);
	free(sm);
	free(locus);
	free(slocus);
	free(sample_gm);
}

/* ***************************************************** */

getinputs()
{
	tracestart();
	getseed();

	getintandskip("NIND (900, 10000, ... max 90000):", &NIND, 2, 90000);
	
	Vk = 0;

	getintandskip("BW (odd; 99: RM):", &BW, 0, 1000);
	
	if (BW != 99) 
	{
		int plano_lado = (int)(sqrt(NIND) + 0.5);
		if (BW >= plano_lado) 
		{
			printf("\nBW (%d) cannot be larger or equal ", BW);
			printf("to plano_lado (%d).\n", plano_lado);
			printf("Adjusting BW to (%d).\n\n", plano_lado - 1);
			BW = plano_lado - 1;
		}
	}

	getintandskip("SW (odd):", &SW, 0, 1000);
	getintandskip("SWR (odd):", &SWR, 0, 1000);
	
	NCRO = 680;
	NLOCI = 30;
	
	getintandskip("Number of generations :", &generations, 1, 1000000);
	getintandskip("Number of generations MULTI ALLELE LOCUS:", &genmulti, 1, 100000);
}

/* **************************************************** */

recombination_masks ()
{
	for (l=0; l<NLOCI; l++) 
	{
		RM[l] = 1 << l;
	}
}

/* ***************************************************** */

void initial_generation(void)
{
    int side = (int)(sqrt((double)NIND) + 0.5);

    /* Coordenadas espaciales */
    i = -1;
    for (fil = 0; fil < side; fil++)
    {
        for (col = 0; col < side; col++)
        {
            i++;
            if (i < NN) {
                f[i] = fil;
                c[i] = col;
            }
        }
    }

    /* Limpiar matrices genoti­picas antes de inicializar */
    memset(gm, 0, sizeof(*gm) * NN);

    /* Calculo de la CDF de la distribucion 1/x */
    int total_copies = 2 * NIND;
    int max_copies = total_copies - 1;

    double H_n = 0.0;
    for (int x = 1; x <= max_copies; x++) {
        H_n += 1.0 / (double)x;
    }

    double *cdf = (double *)malloc(sizeof(double) * (max_copies + 1));
    if (cdf == NULL) {
        printf("Memory allocation failed for CDF.\n");
        exit(1);
    }

    double cumulative = 0.0;
    cdf[0] = 0.0;
    for (int x = 1; x <= max_copies; x++) {
        cumulative += (1.0 / (double)x) / H_n;
        cdf[x] = cumulative;
    }

    /* Busqueda binaria */
    for (k = 0; k < NCRO; k++)
    {
        for (l = 0; l < NLOCI; l++)
        {
            double r = uniform();

            /* Buqueda binaria O(log N) */
            int x_copies = sample_x_copies(r, cdf, max_copies);

            /* Ubicacion aleatoria de las copias en la poblacion */
            int placed = 0;
            while (placed < x_copies)
            {
                int ind = (int)(uniform() * NIND);
                int hap = (int)(uniform() * 2.0);

                if ((gm[ind][k][hap] & RM[l]) == 0) {
                    gm[ind][k][hap] |= RM[l];
                    placed++;
                }
            }
        }
    }

    free(cdf);
}

/* ***************************************************** */

/* Busqueda binaria para la Transformada Inversa sobre el array CDF */
int sample_x_copies(double r, const double *cdf, int max_copies)
{
    int low = 1;
    int high = max_copies;
    int ans = max_copies;

    while (low <= high)
    {
        int mid = low + (high - low) / 2;

        if (cdf[mid] >= r)
        {
            ans = mid;
            high = mid - 1; /* Buscar un valor mas pequeño hacia la izquierda */
        }
        else
        {
            low = mid + 1;  /* Buscar en la mitad derecha */
        }
    }

    return ans;
}

/* ***************************************************** */

void frequency_genes(void)
{
    /* Limpieza inicial de la matriz de frecuencias */
    for (k = 0; k < NCRO; k++) {
        for (l = 0; l < NLOCI; l++) {
            q[k][l] = 0.0;
        }
    }

    /* Conteo de alelos mutados en la poblacion */
    for (i = 0; i < NIND; i++)
    {
        for (k = 0; k < NCRO; k++)
        {
            int gm0 = gm[i][k][0];
            int gm1 = gm[i][k][1];

            for (l = 0; l < NLOCI; l++)
            {
                int mask = RM[l];
                if (gm0 & mask) q[k][l] += 1.0;
                if (gm1 & mask) q[k][l] += 1.0;
            }
        }
    }

    /* Normalizacion por el total de haplotipos (2 * NIND) */
    double divisor = 2.0 * (double)NIND;
    for (k = 0; k < NCRO; k++) {
        for (l = 0; l < NLOCI; l++) {
            q[k][l] /= divisor;
        }
    }

    /* ************************************************************************
     * GUARDAR FRECUENCIAS: Generacion Inicial (0) y Final (generations)
     * *********************************************************************** */
    if (gen == 0 || gen == generations) 
    {
        FILE *fsnp = fopen("snp_frequencies.dat", (gen == 0) ? "w" : "a");
        
        if (fsnp != NULL) 
        {
            if (gen == 0) {
                fprintf(fsnp, "GEN\tCHROM\tLOCUS\tFREQ_ALT(p)\tFREQ_REF(1-p)\n");
            }

            for (k = 0; k < NCRO; k++) 
            {
                for (l = 0; l < NLOCI; l++) 
                {
                    double freq_alt = q[k][l];
                    double freq_ref = 1.0 - freq_alt;
                    fprintf(fsnp, "%d\t%d\t%d\t%.6f\t%.6f\n", gen, k + 1, l + 1, freq_alt, freq_ref);
                }
            }
            fclose(fsnp);
        }
        else
        {
            printf("No se pudo abrir el archivo snp_frequencies.dat\n");
        }
    }
}

/* ***************************************************** */

mutation_neutral ()
{
	double eps = 1e-9; 

	for (k=0; k<NCRO; k++)
	{
		for (l=0; l<NLOCI; l++)
		{
			double freq = q[k][l];

			if (freq <= eps)
			{
				ran_i = (int)(uniform() * NIND);
				ran_h = (int)(uniform() * 2.0);

				gm[ran_i][k][ran_h] |= RM[l];
			}
			else if (freq >= (1.0 - eps))
			{
				ran_i = (int)(uniform() * NIND);
				ran_h = (int)(uniform() * 2.0);

				gm[ran_i][k][ran_h] &= (~RM[l]);
			}
		}
	}
}

/* **************************************************************** */

plink_files()
{
	char filename_ped[100], filename_pedR[100], filename_map[100], filename_mapR[100];

	if (gen == generations) {
		sprintf(filename_ped, "data.ped");
		sprintf(filename_pedR, "dataR.ped");
		sprintf(filename_map, "data.map");
		sprintf(filename_mapR, "dataR.map");
	} else {
		sprintf(filename_ped, "data_gen%d.ped", gen);
		sprintf(filename_pedR, "dataR_gen%d.ped", gen);
		sprintf(filename_map, "data_gen%d.map", gen);
		sprintf(filename_mapR, "dataR_gen%d.map", gen);
	}

	ffped = fopen(filename_ped, "w");
	ffpedR = fopen(filename_pedR, "w");
	ffmap = fopen(filename_map, "w");
	ffmapR = fopen(filename_mapR, "w");

	if (ffped == NULL || ffpedR == NULL || ffmap == NULL || ffmapR == NULL)
	{
		printf("No se pudo abrir el fichero PLINK output files.\n");

		if (ffped != NULL) fclose(ffped);
		if (ffpedR != NULL) fclose(ffpedR);
		if (ffmap != NULL) fclose(ffmap);
		if (ffmapR != NULL) fclose(ffmapR);

		return 1;
	}

	int side = (int)(sqrt(NIND) + 0.5);
	int mid_side = side / 2;
	int quarter_side = side / 4;
	int three_quarters_side = (3 * side) / 4;
	int sw_half = SW / 2;
	int swr_half = SWR / 2;

	int *valid_k = (int *)malloc(NCRO * NLOCI * sizeof(int));
	int *valid_l = (int *)malloc(NCRO * NLOCI * sizeof(int));

	if (valid_k == NULL || valid_l == NULL)
	{
		printf("No hay valid SNP arrays.\n");

		if (valid_k != NULL) free(valid_k);
		if (valid_l != NULL) free(valid_l);

		fclose(ffped);
		fclose(ffpedR);
		fclose(ffmap);
		fclose(ffmapR);

		return 1;
	}

	int valid_count = 0;

	const char *genotipos_str[] =
	{
		"A A ",
		"T A ",
		"A T ",
		"T T "
	};

	/* **************** data.map y data.ped ******************* */

	SNP = 0;

	for (k = 0; k < NCRO; k++)
	{
		int crom = k / 34;
		int chrom_val = crom + 1;
		unsigned long long base_pos_sub = 1020 * crom;

		for (l = 0; l < NLOCI; l++)
		{
			if ((q[k][l] != 0.0) && (q[k][l] != 1.0))
			{
				SNP++;

				/*
				 * Guardamos el SNP válido.
				 */
				valid_k[valid_count] = k;
				valid_l[valid_count] = l;
				valid_count++;

				if (crom < 20)
				{
					chrom[k][l] = chrom_val;

					pos[k][l] =
						((l + NLOCI * k) - base_pos_sub)
						* 100000ULL;
				}

				double pos_mb =
					(double)pos[k][l] / 1000000.0;

				fprintf(ffmap,
						"%d SNP%d %f %llu\n",
						chrom[k][l],
						SNP,
						pos_mb,
						pos[k][l]);

				fprintf(ffmapR,
						"%d SNP%d %f %llu\n",
						chrom[k][l],
						SNP,
						pos_mb,
						pos[k][l]);
			}
		}
	}

	/* **************** data.ped *************** /*

	for (Sfil = (mid_side - sw_half);
		 Sfil <= (mid_side + sw_half);
		 Sfil++)
	{
		if (Sfil < 0 || Sfil >= side)
			continue;

		for (Scol = (mid_side - sw_half);
			 Scol <= (mid_side + sw_half);
			 Scol++)
		{
			if (Scol < 0 || Scol >= side)
				continue;

			i = (Sfil * side) + Scol;

			fprintf(ffped,
					"1 IND%d 0 0 1 -9 ",
					i);

			for (int idx = 0; idx < valid_count; idx++)
			{
				int vk = valid_k[idx];
				int vl = valid_l[idx];
				int mask = RM[vl];

				int bit0 =
					((gm[i][vk][0] & mask) == mask);

				int bit1 =
					((gm[i][vk][1] & mask) == mask);

				int lookup_idx =
					bit0 | (bit1 << 1);

				fputs(genotipos_str[lookup_idx], ffped);
			}

			fputc('\n', ffped);
		}
	}

	/* ******************** daraR.ped ****************** */

	for (int region = 1; region <= 13; region++)
	{
		int start_fil, end_fil, start_col, end_col;
		const char *prefix;

		switch(region)
		{
			case 1:
				start_fil = 0;
				end_fil = SWR - 1;
				start_col = 0;
				end_col = SWR - 1;
				prefix = "INDUL";
				break;

			case 2:
				start_fil = 0;
				end_fil = SWR - 1;
				start_col = mid_side - swr_half;
				end_col = mid_side + swr_half;
				prefix = "INDML";
				break;

			case 3:
				start_fil = 0;
				end_fil = SWR - 1;
				start_col = side - SWR;
				end_col = side - 1;
				prefix = "INDLL";
				break;

			case 4:
				start_fil = mid_side - swr_half;
				end_fil = mid_side + swr_half;
				start_col = 0;
				end_col = SWR - 1;
				prefix = "INDUC";
				break;

			case 5:
				start_fil = side - SWR;
				end_fil = side - 1;
				start_col = 0;
				end_col = SWR - 1;
				prefix = "INDUR";
				break;

			case 6:
				start_fil = side - SWR;
				end_fil = side - 1;
				start_col = mid_side - swr_half;
				end_col = mid_side + swr_half;
				prefix = "INDMR";
				break;

			case 7:
				start_fil = side - SWR;
				end_fil = side - 1;
				start_col = side - SWR;
				end_col = side - 1;
				prefix = "INDLR";
				break;

			case 8:
				start_fil = mid_side - swr_half;
				end_fil = mid_side + swr_half;
				start_col = side - SWR;
				end_col = side - 1;
				prefix = "INDLC";
				break;

			case 9:
				start_fil = quarter_side - swr_half;
				end_fil = quarter_side + swr_half;
				start_col = quarter_side - swr_half;
				end_col = quarter_side + swr_half;
				prefix = "INDMUL";
				break;

			case 10:
				start_fil = quarter_side - swr_half;
				end_fil = quarter_side + swr_half;
				start_col = three_quarters_side - swr_half;
				end_col = three_quarters_side + swr_half;
				prefix = "INDMUR";
				break;

			case 11:
				start_fil = three_quarters_side - swr_half;
				end_fil = three_quarters_side + swr_half;
				start_col = quarter_side - swr_half;
				end_col = quarter_side + swr_half;
				prefix = "INDMLL";
				break;

			case 12:
				start_fil = three_quarters_side - swr_half;
				end_fil = three_quarters_side + swr_half;
				start_col = three_quarters_side - swr_half;
				end_col = three_quarters_side + swr_half;
				prefix = "INDMLR";
				break;

			case 13:
				start_fil = mid_side - swr_half;
				end_fil = mid_side + swr_half;
				start_col = mid_side - swr_half;
				end_col = mid_side + swr_half;
				prefix = "INDC";
				break;

			default:
				continue;
		}

		for (Sfil = start_fil;
			 Sfil <= end_fil;
			 Sfil++)
		{
			if (Sfil < 0 || Sfil >= side)
				continue;

			for (Scol = start_col;
				 Scol <= end_col;
				 Scol++)
			{
				if (Scol < 0 || Scol >= side)
					continue;

				i = (Sfil * side) + Scol;

				fprintf(ffpedR,
						"1 %s%d 0 0 1 -9 ",
						prefix,
						i);

				for (int idx = 0;
					 idx < valid_count;
					 idx++)
				{
					int vk = valid_k[idx];
					int vl = valid_l[idx];
					int mask = RM[vl];

					int bit0 =
						((gm[i][vk][0] & mask) == mask);

					int bit1 =
						((gm[i][vk][1] & mask) == mask);

					int lookup_idx =
						bit0 | (bit1 << 1);

					fputs(genotipos_str[lookup_idx],
						  ffpedR);
				}

				fputc('\n', ffpedR);
			}
		}
	}

	/* ******** data 2000; SOLO en la última generación ***** */

	if (gen == generations)
	{
		if (generate_2000_sample() != 0)
		{
			printf("No se pudieron general los 2000 descendientes\n");

			free(valid_k);
			free(valid_l);

			fclose(ffped);
			fclose(ffpedR);
			fclose(ffmap);
			fclose(ffmapR);

			return 1;
		}

		FILE *ffped2000 = fopen("data2000.ped", "w");
		FILE *ffmap2000 = fopen("data2000.map", "w");

		if (ffped2000 == NULL || ffmap2000 == NULL)
		{
			printf("No se pudo abrir data2000.ped/map.\n");

			if (ffped2000 != NULL)
				fclose(ffped2000);

			if (ffmap2000 != NULL)
				fclose(ffmap2000);

			free(valid_k);
			free(valid_l);

			fclose(ffped);
			fclose(ffpedR);
			fclose(ffmap);
			fclose(ffmapR);

			return 1;
		}

		/* ******* data2000.map **************** */

		SNP = 0;

		for (int idx = 0;
			 idx < valid_count;
			 idx++)
		{
			int vk = valid_k[idx];
			int vl = valid_l[idx];

			SNP++;

			double pos_mb =
				(double)pos[vk][vl] / 1000000.0;

			fprintf(ffmap2000,
					"%d SNP%d %f %llu\n",
					chrom[vk][vl],
					SNP,
					pos_mb,
					pos[vk][vl]);
		}

		/* ******* data2000.ped **************** */

		for (int child = 0;
			 child < NSAMPLE;
			 child++)
		{
			fprintf(ffped2000,
					"1 SAMPLE2000_%d 0 0 1 -9 ",
					child + 1);

			for (int idx = 0;
				 idx < valid_count;
				 idx++)
			{
				int vk = valid_k[idx];
				int vl = valid_l[idx];

				int mask = RM[vl];

				int bit0 =
					((sample_gm[child][vk][0] & mask) == mask);

				int bit1 =
					((sample_gm[child][vk][1] & mask) == mask);

				int lookup_idx =
					bit0 | (bit1 << 1);

				fputs(genotipos_str[lookup_idx],
					  ffped2000);
			}

			fputc('\n', ffped2000);
		}

		fclose(ffped2000);
		fclose(ffmap2000);

		printf("data2000.ped and data2000.map generated: %d offspring.\n",
			   NSAMPLE);
	}

	free(valid_k);
	free(valid_l);

	fclose(ffped);
	fclose(ffpedR);
	fclose(ffmap);
	fclose(ffmapR);

	return 0;
}

/* ************************************************************ */

initial_multiallelic ()
{
	for (i=0; i<NIND; i++)
	{
		int id_father = i;
		int id_mother = i + NIND;

		for (k=0; k<NCRO; k++)
		{
			locus[i][k][0] = id_father;
			locus[i][k][1] = id_mother;
		}
	}
}

/* ***************************************************** */

estimates_of_F ()
{
	numLOCI[gen] = 0;
	EHom = 0.0;
	EHet = 0.0;

	for (k=0; k<NCRO; k++)
	{
		for (l=0; l<NLOCI; l++)
		{
			if (q[k][l] != 0.0 && q[k][l] != 1.0)
			{
				numLOCI[gen]++;
				EHom += (1.0 - 2.0 * q[k][l] * (1.0 - q[k][l]));
				EHet += 2.0 * q[k][l] * (1.0 - q[k][l]);
			}
		}
	}
	accum(&Het[gen], EHet / (double)numLOCI[gen]);

	for (i=0; i<NIND; i++)
	{
		fibd[i] = 0.0;
		for (k=0; k<NCRO; k++)
		{
			if (locus[i][k][0] == locus[i][k][1])
			{
				fibd[i] += 1.0;
			}
		}
		fibd[i] /= (double)NCRO;

		int hom_count = 0;
		for (k=0; k<NCRO; k++)
		{
			int gm0 = gm[i][k][0];
			int gm1 = gm[i][k][1];

			for (l=0; l<NLOCI; l++)
			{
				if (q[k][l] != 0.0 && q[k][l] != 1.0)
				{
					int mask = RM[l];
					if (((gm0 ^ gm1) & mask) == 0)
					{
						hom_count++;
					}
				}
			}
		}
		Hom[i] = (double)hom_count;

		fLH1[i] = (Hom[i] - EHom) / ((double)numLOCI[gen] - EHom);

		accum(&Fibd[gen], fibd[i]);
		accum(&FLH1[gen], fLH1[i]);
	}

	int total_alleles_limit = 2 * NIND;
	size_t clear_size = (total_alleles_limit + 1) * sizeof(int);

	for (k=0; k<NCRO; k++)
	{
		memset(qallele, 0, clear_size);

		for (i=0; i<NIND; i++)
		{
			qallele[locus[i][k][0]]++;
			qallele[locus[i][k][1]]++;
		}

		num_a = 0;
		for (all=0; all<total_alleles_limit; all++)
		{
			if (qallele[all] > 0)
			{
				num_a++;
			}
		}
		accum(&Num_alleles[gen], (double)num_a);
	}
}

/* ***************************************************** */

estimates_of_Fst_between_opposite_NS ()
{
    int side = (int)(sqrt(NIND) + 0.5); 
    
    for (k=0; k<NCRO; k++) {
        for (l=0; l<NLOCI; l++) {
            qNS1[k][l] = 0.0;
            qNS2[k][l] = 0.0;
        }
    }

    for (fil=0; fil<BW; fil++)
    {
        for (col=0; col<BW; col++)
        {
            i = fil * side + col; 
            for (k=0; k<NCRO; k++)
            {
                int gm0 = gm[i][k][0];
                int gm1 = gm[i][k][1];
                for (l=0; l<NLOCI; l++)
                {
                    int mask = RM[l];
                    int bit0 = ((gm0 & mask) == mask);
                    int bit1 = ((gm1 & mask) == mask);
                    if (bit0 && bit1)      qNS1[k][l] += 2.0;
                    else if (bit0 || bit1) qNS1[k][l] += 1.0;
                }
            }
        }
    }

    for (fil=(side-BW); fil<side; fil++)
    {
        for (col=(side-BW); col<side; col++)
        {
            i = fil * side + col;
            for (k=0; k<NCRO; k++)
            {
                int gm0 = gm[i][k][0];
                int gm1 = gm[i][k][1];
                for (l=0; l<NLOCI; l++)
                {
                    int mask = RM[l];
                    int bit0 = ((gm0 & mask) == mask);
                    int bit1 = ((gm1 & mask) == mask);
                    if (bit0 && bit1)      qNS2[k][l] += 2.0;
                    else if (bit0 || bit1) qNS2[k][l] += 1.0;
                }
            }
        }
    }

    double bw_size = (double)(BW * BW);
    for (k=0; k<NCRO; k++)
    {
        for (l=0; l<NLOCI; l++)
        {
            qNS1[k][l] /= (2.0 * bw_size);
            qNS2[k][l] /= (2.0 * bw_size);
            qNS1NS2[k][l] = (qNS1[k][l] + qNS2[k][l]) / 2.0;
        }
    }

    EHetNS1NS2 = 0.0; lociNS1NS2 = 0.0; EHetNS1 = 0.0; EHetNS2 = 0.0;
    for (k=0; k<NCRO; k++)
    {
        for (l=0; l<NLOCI; l++)
        {
            if (qNS1NS2[k][l] != 0.0 && qNS1NS2[k][l] != 1.0)
            {
                lociNS1NS2 += 1.0;
                EHetNS1NS2 += 2.0 * qNS1NS2[k][l] * (1.0 - qNS1NS2[k][l]);
                EHetNS1    += 2.0 * qNS1[k][l]    * (1.0 - qNS1[k][l]);
                EHetNS2    += 2.0 * qNS2[k][l]    * (1.0 - qNS2[k][l]);
            }
        }
    }

    HetNS1NS2 = 0.0; HetNS1 = 0.0; HetNS2 = 0.0;

    for (fil=0; fil<BW; fil++)
    {
        for (col=0; col<BW; col++)
        {
            i = fil * side + col; 
            for (k=0; k<NCRO; k++)
            {
                int gm0 = gm[i][k][0];
                int gm1 = gm[i][k][1];
                for (l=0; l<NLOCI; l++)
                {
                    if (qNS1NS2[k][l] != 0.0 && qNS1NS2[k][l] != 1.0)
                    {
                        int mask = RM[l];
                        if (((gm0 ^ gm1) & mask) == mask) { HetNS1NS2 += 1.0; HetNS1 += 1.0; }
                    }
                }
            }
        }
    }

    for (fil=(side-BW); fil<side; fil++)
    {
        for (col=(side-BW); col<side; col++)
        {
            i = fil * side + col;
            for (k=0; k<NCRO; k++)
            {
                int gm0 = gm[i][k][0];
                int gm1 = gm[i][k][1];
                for (l=0; l<NLOCI; l++)
                {
                    if (qNS1NS2[k][l] != 0.0 && qNS1NS2[k][l] != 1.0)
                    {
                        int mask = RM[l];
                        if (((gm0 ^ gm1) & mask) == mask) { HetNS1NS2 += 1.0; HetNS2 += 1.0; }
                    }
                }
            }
        }
    }

    HT = EHetNS1NS2 / lociNS1NS2;
    HI = ((HetNS1 / bw_size) + (HetNS2 / bw_size)) / (2.0 * lociNS1NS2);
    HS = (EHetNS1 + EHetNS2) / (2.0 * lociNS1NS2);
    FSTmax = (HT - HS) / HT;
    FISmax = (HS - HI) / HS;
    FITmax = (HT - HI) / HT;
}

/* ***************************************************** */

estimates_of_Fst_between_close_NS ()
{
    int side = (int)(sqrt(NIND) + 0.5); 
    int half_side = side / 2;
    int half_bw = BW / 2;

    int ns3_fil_start = half_side - half_bw;
    int ns3_fil_end   = half_side + half_bw;
    int ns3_col_start = half_side - half_bw;
    int ns3_col_end   = half_side + half_bw;

    int ns4_fil_start = half_side + BW - half_bw;
    int ns4_fil_end   = half_side + BW + half_bw;
    int ns4_col_start = half_side + BW - half_bw;
    int ns4_col_end   = half_side + BW + half_bw;

    if (ns3_fil_start < 0) ns3_fil_start = 0; if (ns3_fil_end >= side) ns3_fil_end = side - 1;
    if (ns3_col_start < 0) ns3_col_start = 0; if (ns3_col_end >= side) ns3_col_end = side - 1;
    if (ns4_fil_start < 0) ns4_fil_start = 0; if (ns4_fil_end >= side) ns4_fil_end = side - 1;
    if (ns4_col_start < 0) ns4_col_start = 0; if (ns4_col_end >= side) ns4_col_end = side - 1;

    double count_ns3 = (double)((ns3_fil_end - ns3_fil_start + 1) * (ns3_col_end - ns3_col_start + 1));
    double count_ns4 = (double)((ns4_fil_end - ns4_fil_start + 1) * (ns4_col_end - ns4_col_start + 1));

    for (k=0; k<NCRO; k++) {
        for (l=0; l<NLOCI; l++) {
            qNS3[k][l] = 0.0;
            qNS4[k][l] = 0.0;
        }
    }

    for (fil = ns3_fil_start; fil <= ns3_fil_end; fil++)
    {
        for (col = ns3_col_start; col <= ns3_col_end; col++)
        {
            i = fil * side + col;
            for (k=0; k<NCRO; k++)
            {
                int gm0 = gm[i][k][0];
                int gm1 = gm[i][k][1];
                for (l=0; l<NLOCI; l++)
                {
                    int mask = RM[l];
                    int bit0 = ((gm0 & mask) == mask);
                    int bit1 = ((gm1 & mask) == mask);
                    if (bit0 && bit1)      qNS3[k][l] += 2.0;
                    else if (bit0 || bit1) qNS3[k][l] += 1.0;
                }
            }
        }
    }

    for (fil = ns4_fil_start; fil <= ns4_fil_end; fil++)
    {
        for (col = ns4_col_start; col <= ns4_col_end; col++)
        {
            i = fil * side + col;
            for (k=0; k<NCRO; k++)
            {
                int gm0 = gm[i][k][0];
                int gm1 = gm[i][k][1];
                for (l=0; l<NLOCI; l++)
                {
                    int mask = RM[l];
                    int bit0 = ((gm0 & mask) == mask);
                    int bit1 = ((gm1 & mask) == mask);
                    if (bit0 && bit1)      qNS4[k][l] += 2.0;
                    else if (bit0 || bit1) qNS4[k][l] += 1.0;
                }
            }
        }
    }

    for (k=0; k<NCRO; k++)
    {
        for (l=0; l<NLOCI; l++)
        {
            qNS3[k][l] /= (2.0 * count_ns3);
            qNS4[k][l] /= (2.0 * count_ns4);
            qNS3NS4[k][l] = (qNS3[k][l] + qNS4[k][l]) / 2.0;
        }
    }

    EHetNS3NS4 = 0.0; lociNS3NS4 = 0.0; EHetNS3 = 0.0; EHetNS4 = 0.0;
    for (k=0; k<NCRO; k++)
    {
        for (l=0; l<NLOCI; l++)
        {
            if (qNS3NS4[k][l] != 0.0 && qNS3NS4[k][l] != 1.0)
            {
                lociNS3NS4 += 1.0;
                EHetNS3NS4 += 2.0 * qNS3NS4[k][l] * (1.0 - qNS3NS4[k][l]);
                EHetNS3    += 2.0 * qNS3[k][l]    * (1.0 - qNS3[k][l]);
                EHetNS4    += 2.0 * qNS4[k][l]    * (1.0 - qNS4[k][l]);
            }
        }
    }

    HetNS3NS4 = 0.0; HetNS3 = 0.0; HetNS4 = 0.0;

    for (fil = ns3_fil_start; fil <= ns3_fil_end; fil++)
    {
        for (col = ns3_col_start; col <= ns3_col_end; col++)
        {
            i = fil * side + col;
            for (k=0; k<NCRO; k++)
            {
                int gm0 = gm[i][k][0]; int gm1 = gm[i][k][1];
                for (l=0; l<NLOCI; l++)
                {
                    if (qNS3NS4[k][l] != 0.0 && qNS3NS4[k][l] != 1.0)
                    {
                        int mask = RM[l];
                        if (((gm0 ^ gm1) & mask) == mask) { HetNS3NS4 += 1.0; HetNS3 += 1.0; }
                    }
                }
            }
        }
    }

    for (fil = ns4_fil_start; fil <= ns4_fil_end; fil++)
    {
        for (col = ns4_col_start; col <= ns4_col_end; col++)
        {
            i = fil * side + col;
            for (k=0; k<NCRO; k++)
            {
                int gm0 = gm[i][k][0]; int gm1 = gm[i][k][1];
                for (l=0; l<NLOCI; l++)
                {
                    if (qNS3NS4[k][l] != 0.0 && qNS3NS4[k][l] != 1.0)
                    {
                        int mask = RM[l];
                        if (((gm0 ^ gm1) & mask) == mask) { HetNS3NS4 += 1.0; HetNS4 += 1.0; }
                    }
                }
            }
        }
    }

    HT = EHetNS3NS4 / lociNS3NS4;
    HI = ((HetNS3 / count_ns3) + (HetNS4 / count_ns4)) / (2.0 * lociNS3NS4);
    HS = (EHetNS3 + EHetNS4) / (2.0 * lociNS3NS4);
    FSTmin = (HT - HS) / HT;
    FISmin = (HS - HI) / HS;
    FITmin = (HT - HI) / HT;
}

/* ***************************************************** */

estimates_of_Fst_between_13_groups ()
{
    int side = (int)(sqrt(NIND) + 0.5); /* 300 fijo para NIND = 90000 */
    int half_side = side / 2;
    int half_bw = BW / 2;
    int g, g1, g2;
    int quarter_side = half_side / 2;

    struct GroupBounds {
        int r_start, r_end;
        int c_start, c_end;
    } bounds[13];

    
    // 0. Upper Left Corner (Esquina Superior Izquierda)
    bounds[0].r_start = 0;           bounds[0].r_end = BW - 1;
    bounds[0].c_start = 0;           bounds[0].c_end = BW - 1;
    
    // 1. Upper Right Corner (Esquina Superior Derecha)
    bounds[1].r_start = 0;           bounds[1].r_end = BW - 1;
    bounds[1].c_start = side - BW;   bounds[1].c_end = side - 1;
    
    // 2. Lower Left Corner (Esquina Inferior Izquierda)
    bounds[2].r_start = side - BW;   bounds[2].r_end = side - 1;
    bounds[2].c_start = 0;           bounds[2].c_end = BW - 1;
    
    // 3. Lower Right Corner (Esquina Inferior Derecha)
    bounds[3].r_start = side - BW;   bounds[3].r_end = side - 1;
    bounds[3].c_start = side - BW;   bounds[3].c_end = side - 1;
    
    // 4. Upper Border (Borde Central Norte)
    bounds[4].r_start = 0;           bounds[4].r_end = BW - 1;
    bounds[4].c_start = half_side - half_bw; bounds[4].c_end = half_side + half_bw;
    
    // 5. Lower Border (Borde Central Sur)
    bounds[5].r_start = side - BW;   bounds[5].r_end = side - 1;
    bounds[5].c_start = half_side - half_bw; bounds[5].c_end = half_side + half_bw;
    
    // 6. Left Border (Borde Central Oeste)
    bounds[6].r_start = half_side - half_bw; bounds[6].r_end = half_side + half_bw;
    bounds[6].c_start = 0;           bounds[6].c_end = BW - 1;
    
    // 7. Right Border (Borde Central Este)
    bounds[7].r_start = half_side - half_bw; bounds[7].r_end = half_side + half_bw;
    bounds[7].c_start = side - BW;   bounds[7].c_end = side - 1;

    // 8. Central Neighborhood (Centro Puro - NS3)
    bounds[8].r_start = half_side - half_bw; bounds[8].r_end = half_side + half_bw;
    bounds[8].c_start = half_side - half_bw; bounds[8].c_end = half_side + half_bw;

    // 9. Quadrant 1 (Intermedio Superior Izquierdo)
    bounds[9].r_start = quarter_side;           bounds[9].r_end = quarter_side + BW - 1;
    bounds[9].c_start = quarter_side;           bounds[9].c_end = quarter_side + BW - 1;
    
    // 10. Quadrant 2 (Intermedio Superior Derecho)
    bounds[10].r_start = quarter_side;          bounds[10].r_end = quarter_side + BW - 1;
    bounds[10].c_start = half_side + quarter_side; bounds[10].c_end = half_side + quarter_side + BW - 1;
    
    // 11. Quadrant 3 (Intermedio Inferior Izquierdo)
    bounds[11].r_start = half_side + quarter_side; bounds[11].r_end = half_side + quarter_side + BW - 1;
    bounds[11].c_start = quarter_side;          bounds[11].c_end = quarter_side + BW - 1;
    
    // 12. Quadrant 4 (Intermedio Inferior Derecho)
    bounds[12].r_start = half_side + quarter_side; bounds[12].r_end = half_side + quarter_side + BW - 1;
    bounds[12].c_start = half_side + quarter_side; bounds[12].c_end = half_side + quarter_side + BW - 1;

    for (g = 0; g < 13; g++) {
        if (bounds[g].r_start < 0) bounds[g].r_start = 0; if (bounds[g].r_end >= side) bounds[g].r_end = side - 1;
        if (bounds[g].c_start < 0) bounds[g].c_start = 0; if (bounds[g].c_end >= side) bounds[g].c_end = side - 1;
        countGroups[g] = (double)((bounds[g].r_end - bounds[g].r_start + 1) * (bounds[g].c_end - bounds[g].c_start + 1));
    }

    memset(qGroups, 0, sizeof(qGroups));

    for (g = 0; g < 13; g++)
    {
        for (fil = bounds[g].r_start; fil <= bounds[g].r_end; fil++)
        {
            for (col = bounds[g].c_start; col <= bounds[g].c_end; col++)
            {
                i = fil * side + col;
                for (k = 0; k < NCRO; k++)
                {
                    int gm0 = gm[i][k][0];
                    int gm1 = gm[i][k][1];
                    for (l = 0; l < NLOCI; l++)
                    {
                        int mask = RM[l];
                        int bit0 = ((gm0 & mask) == mask);
                        int bit1 = ((gm1 & mask) == mask);
                        if (bit0 && bit1)      qGroups[g][k][l] += 2.0;
                        else if (bit0 || bit1) qGroups[g][k][l] += 1.0;
                    }
                }
            }
        }
        for (k = 0; k < NCRO; k++) {
            for (l = 0; l < NLOCI; l++) {
                qGroups[g][k][l] /= (2.0 * countGroups[g]);
            }
        }
    }

    for (g1 = 0; g1 < 13; g1++)
    {
        for (g2 = g1; g2 < 13; g2++)
        {
            if (g1 == g2) { FST_matrix[g1][g2] = 0.0; FIS_matrix[g1][g2] = 0.0; FIT_matrix[g1][g2] = 0.0; continue; }

            double EHetPair = 0.0, lociPair = 0.0;
            double EHetG1 = 0.0, EHetG2 = 0.0;
            double HetG1 = 0.0, HetG2 = 0.0;

            for (k = 0; k < NCRO; k++) {
                for (l = 0; l < NLOCI; l++) {
                    double qBar = (qGroups[g1][k][l] + qGroups[g2][k][l]) / 2.0;
                    if (qBar != 0.0 && qBar != 1.0) {
                        lociPair += 1.0;
                        EHetPair += 2.0 * qBar * (1.0 - qBar);
                        EHetG1   += 2.0 * qGroups[g1][k][l] * (1.0 - qGroups[g1][k][l]);
                        EHetG2   += 2.0 * qGroups[g2][k][l] * (1.0 - qGroups[g2][k][l]);
                    }
                }
            }

            for (fil = bounds[g1].r_start; fil <= bounds[g1].r_end; fil++) {
                for (col = bounds[g1].c_start; col <= bounds[g1].c_end; col++) {
                    i = fil * side + col;
                    for (k = 0; k < NCRO; k++) {
                        int gm0 = gm[i][k][0]; int gm1 = gm[i][k][1];
                        for (l = 0; l < NLOCI; l++) {
                            double qBar = (qGroups[g1][k][l] + qGroups[g2][k][l]) / 2.0;
                            if (qBar != 0.0 && qBar != 1.0) {
                                int mask = RM[l];
                                if (((gm0 ^ gm1) & mask) == mask) HetG1 += 1.0;
                            }
                        }
                    }
                }
            }

            for (fil = bounds[g2].r_start; fil <= bounds[g2].r_end; fil++) {
                for (col = bounds[g2].c_start; col <= bounds[g2].c_end; col++) {
                    i = fil * side + col;
                    for (k = 0; k < NCRO; k++) {
                        int gm0 = gm[i][k][0]; int gm1 = gm[i][k][1];
                        for (l = 0; l < NLOCI; l++) {
                            double qBar = (qGroups[g2][k][l] + qGroups[g1][k][l]) / 2.0;
                            if (qBar != 0.0 && qBar != 1.0) {
                                int mask = RM[l];
                                if (((gm0 ^ gm1) & mask) == mask) HetG2 += 1.0;
                            }
                        }
                    }
                }
            }

            double HT_pair = EHetPair / lociPair;
            double HI_pair = ((HetG1 / countGroups[g1]) + (HetG2 / countGroups[g2])) / (2.0 * lociPair);
            double HS_pair = (EHetG1 + EHetG2) / (2.0 * lociPair);

            FST_matrix[g1][g2] = FST_matrix[g2][g1] = (HT_pair - HS_pair) / HT_pair;
            FIS_matrix[g1][g2] = FIS_matrix[g2][g1] = (HS_pair - HI_pair) / HS_pair;
            FIT_matrix[g1][g2] = FIT_matrix[g2][g1] = (HT_pair - HI_pair) / HT_pair;
        }
    }

    double global_EHetHT = 0.0, global_EHetHS = 0.0, global_HetHI = 0.0;
    double global_loci = 0.0;
    double total_individuals = 0.0;

    for (g = 0; g < 13; g++) total_individuals += countGroups[g];

    for (k = 0; k < NCRO; k++)
    {
        for (l = 0; l < NLOCI; l++)
        {
            double qMeanAllGroups = 0.0;
            double sumEHetSubpopulations = 0.0;

            for (g = 0; g < 13; g++) {
                qMeanAllGroups += qGroups[g][k][l] * (countGroups[g] / total_individuals);
                sumEHetSubpopulations += 2.0 * qGroups[g][k][l] * (1.0 - qGroups[g][k][l]);
            }

            if (qMeanAllGroups != 0.0 && qMeanAllGroups != 1.0)
            {
                global_loci += 1.0;
                global_EHetHT += 2.0 * qMeanAllGroups * (1.0 - qMeanAllGroups); 
                global_EHetHS += (sumEHetSubpopulations / 13.0);                
            }
        }
    }

    double sumObservedHeterozygosity = 0.0;
    for (g = 0; g < 13; g++)
    {
        double sub_het = 0.0;
        for (fil = bounds[g].r_start; fil <= bounds[g].r_end; fil++) {
            for (col = bounds[g].c_start; col <= bounds[g].c_end; col++) {
                i = fil * side + col;
                for (k = 0; k < NCRO; k++) {
                    int gm0 = gm[i][k][0]; int gm1 = gm[i][k][1];
                    for (l = 0; l < NLOCI; l++) {
                        if (qGroups[g][k][l] != 0.0 && qGroups[g][k][l] != 1.0) {
                            int mask = RM[l];
                            if (((gm0 ^ gm1) & mask) == mask) sub_het += 1.0;
                        }
                    }
                }
            }
        }
        sumObservedHeterozygosity += (sub_het / countGroups[g]);
    }
    global_HetHI = (sumObservedHeterozygosity / 13.0) / global_loci; 

    double HT_glob = global_EHetHT / global_loci;
    double HS_glob = global_EHetHS / global_loci;

    FST_global_13 = (HT_glob - HS_glob) / HT_glob;
    FIS_global_13 = (HS_glob - global_HetHI) / HS_glob;
    FIT_global_13 = (HT_glob - global_HetHI) / HT_glob;
}

/* ***************************************************** */

mating ()
{
    int a, p1, p2, EE[MM], FF[MM], family[NN];
    int pointrec[MM], ncrorec[MM], rndk, rndl;
    int rndfil, rndcol;
    double rnd, family_sum2=0.0;
    double fit2, fit3;

    int side = (int)(sqrt(NIND) + 0.5); 
    int half_bw = BW / 2;

    for (i=0; i<NIND; i++)
    {
        for (k=0; k<NCRO; k++)
        {
            sm[i][k][0] = gm[i][k][0];
            sm[i][k][1] = gm[i][k][1];
            slocus[i][k][0] = locus[i][k][0];
            slocus[i][k][1] = locus[i][k][1];
        }
    }

    fit2 = 0.0; fit3 = 0.0;
    for (i=0; i<NIND; i++) family[i] = 0;

    for (i=0; i<NIND; i++)
    {
        if (BW == 99)
        {
            p1 = (int)(uniform() * NIND);
            do { p2 = (int)(uniform() * NIND); } while (p2 == p1);
        }
        else
        {
            do {
                int min_fil = f[i] - half_bw; if (min_fil < 0) min_fil = 0;
                int max_fil = f[i] + half_bw; if (max_fil >= side) max_fil = side - 1;
                rndfil = min_fil + (int)(uniform() * (max_fil - min_fil + 1));

                int min_col = c[i] - half_bw; if (min_col < 0) min_col = 0;
                int max_col = c[i] + half_bw; if (max_col >= side) max_col = side - 1;
                rndcol = min_col + (int)(uniform() * (max_col - min_col + 1));

                p1 = rndfil * side + rndcol;
            } while (p1 == i);

            do {
                int min_fil = f[i] - half_bw; if (min_fil < 0) min_fil = 0;
                int max_fil = f[i] + half_bw; if (max_fil >= side) max_fil = side - 1;
                rndfil = min_fil + (int)(uniform() * (max_fil - min_fil + 1));

                int min_col = c[i] - half_bw; if (min_col < 0) min_col = 0;
                int max_col = c[i] + half_bw; if (max_col >= side) max_col = side - 1;
                rndcol = min_col + (int)(uniform() * (max_col - min_col + 1));

                p2 = rndfil * side + rndcol;
            } while (p2 == i || p2 == p1);
        }

        family[p1]++;    
        family[p2]++;

        memset(ncrorec, 0, sizeof(ncrorec));
        memset(pointrec, 0, sizeof(pointrec));

        for (crom=0; crom<20; crom++)
        {
            ncrorec[crom * 34] = 1;
            pointrec[crom * 34] = 0;
        }

        for (crom=0; crom<20; crom++)
        {
            rndk = (int)(uniform() * 34);
            do { rndl = (int)(uniform() * NLOCI); } while (rndl == 0);

            ncrorec[(crom * 34) + rndk] = 1;
            pointrec[(crom * 34) + rndk] = rndl;
        }

        /* Recombination father */
        marker = 1;
        for (k=0; k<NCRO; k++)
        {
            EE[k] = 0;
            if (ncrorec[k] == 0)
            {
                if (marker == (-1)) EE[k] = ~EE[k];
            }
            else
            {
                for (l=0; l<NLOCI; l++)
                {
                    if (l != pointrec[k])
                    {
                        if (marker == (-1)) EE[k] |= RM[l];
                    }
                    else
                    {
                        if (marker == 1)
                        {
                            EE[k] |= RM[l];
                            marker = -marker;
                        }
                        else
                        {
                            marker = -marker;
                        }
                    }
                }
            }
        }

        rnd = uniform();
        int bit_mask_RM0 = RM[0]; 
        for (k=0; k<NCRO; k++)
        {
            if (rnd < 0.5) EE[k] = ~EE[k];
            FF[k] = ~EE[k];
            gm[i][k][0] = ((EE[k] & sm[p1][k][0]) | (FF[k] & sm[p1][k][1]));

            if ((EE[k] & bit_mask_RM0) == bit_mask_RM0) 
                locus[i][k][0] = slocus[p1][k][0];
            else 
                locus[i][k][0] = slocus[p1][k][1];
        }

        /* Recombination mother */
        memset(ncrorec, 0, sizeof(ncrorec));
        memset(pointrec, 0, sizeof(pointrec));

        for (crom=0; crom<20; crom++)
        {
            ncrorec[crom * 34] = 1;
            pointrec[crom * 34] = 0;
        }

        for (crom=0; crom<20; crom++)
        {
            rndk = (int)(uniform() * 34);
            do { rndl = (int)(uniform() * NLOCI); } while (rndl == 0);

            ncrorec[(crom * 34) + rndk] = 1;
            pointrec[(crom * 34) + rndk] = rndl;
        }

        marker = 1;
        for (k=0; k<NCRO; k++)
        {
            EE[k] = 0;
            if (ncrorec[k] == 0)
            {
                if (marker == (-1)) EE[k] = ~EE[k];
            }
            else
            {
                for (l=0; l<NLOCI; l++)
                {
                    if (l != pointrec[k])
                    {
                        if (marker == (-1)) EE[k] |= RM[l];
                    }
                    else
                    {
                        if (marker == 1)
                        {
                            EE[k] |= RM[l];
                            marker = -marker;
                        }
                        else
                        {
                            marker = -marker;
                        }
                    }
                }
            }
        }

        rnd = uniform();
        for (k=0; k<NCRO; k++)
        {
            if (rnd < 0.5) EE[k] = ~EE[k];
            FF[k] = ~EE[k];
            gm[i][k][1] = ((EE[k] & sm[p2][k][0]) | (FF[k] & sm[p2][k][1]));

            if ((EE[k] & bit_mask_RM0) == bit_mask_RM0) 
                locus[i][k][1] = slocus[p2][k][0];
            else 
                locus[i][k][1] = slocus[p2][k][1];
        }
    }

    for (i=0; i<NIND; i++) family_sum2 += ((double)family[i] * (double)family[i]);
    accum(&SK2, ((family_sum2 - (4.0 * NIND)) / (NIND - 1.0)));
}

/* ********************* Genera 2000 hijos adicionales para data2000.ped *************************************** */

int generate_2000_sample(void)
{
    int side = (int)(sqrt((double)NIND) + 0.5);
    int mid_side = side / 2;
    int sw_half = SW / 2;
    int half_bw = BW / 2;

    int sw_start_fil = mid_side - sw_half;
    int sw_end_fil   = mid_side + sw_half;
    int sw_start_col = mid_side - sw_half;
    int sw_end_col   = mid_side + sw_half;

    int nsw = 0;

    for (int fil = sw_start_fil; fil <= sw_end_fil; fil++)
    {
        if (fil < 0 || fil >= side)
            continue;

        for (int col = sw_start_col; col <= sw_end_col; col++)
        {
            if (col < 0 || col >= side)
                continue;

            int id = fil * side + col;

            if (id >= 0 && id < NIND)
                nsw++;
        }
    }

    if (nsw <= 0)
    {
        printf("SW no tiene individuos \n");
        return 1;
    }

    int *sw_individuals =
        (int *)malloc((size_t)nsw * sizeof(int));

    if (sw_individuals == NULL)
    {
        printf("[critical error]: No memory for SW individuals.\n");
        return 1;
    }

    int nsw_check = 0;

    for (int fil = sw_start_fil; fil <= sw_end_fil; fil++)
    {
        if (fil < 0 || fil >= side)
            continue;

        for (int col = sw_start_col; col <= sw_end_col; col++)
        {
            if (col < 0 || col >= side)
                continue;

            int id = fil * side + col;

            if (id >= 0 && id < NIND)
            {
                sw_individuals[nsw_check] = id;
                nsw_check++;
            }
        }
    }

    if (nsw_check != nsw)
    {
        printf("[critical error]: SW individual count mismatch.\n");
        free(sw_individuals);
        return 1;
    }

    int EE[MM];
    int FF[MM];
    int pointrec[MM];
    int ncrorec[MM];

    int bit_mask_RM0 = RM[0];

    /* ********* Generación de los NSAMPLE hijos ************* */

    for (int child = 0; child < NSAMPLE; child++)
    {
        int i_target =
            sw_individuals[(int)(uniform() * nsw)];

        int p1, p2;

        if (BW == 99)
        {
            p1 = (int)(uniform() * NIND);

            do
            {
                p2 = (int)(uniform() * NIND);
            }
            while (p2 == p1);
        }
        else
        {
            do
            {
                int min_fil = f[i_target] - half_bw;

                if (min_fil < 0)
                    min_fil = 0;

                int max_fil = f[i_target] + half_bw;

                if (max_fil >= side)
                    max_fil = side - 1;

                int rndfil =
                    min_fil +
                    (int)(uniform() *
                          (max_fil - min_fil + 1));

                int min_col = c[i_target] - half_bw;

                if (min_col < 0)
                    min_col = 0;

                int max_col = c[i_target] + half_bw;

                if (max_col >= side)
                    max_col = side - 1;

                int rndcol =
                    min_col +
                    (int)(uniform() *
                          (max_col - min_col + 1));

                p1 = rndfil * side + rndcol;
            }
            while (p1 == i_target);

            do
            {
                int min_fil = f[i_target] - half_bw;

                if (min_fil < 0)
                    min_fil = 0;

                int max_fil = f[i_target] + half_bw;

                if (max_fil >= side)
                    max_fil = side - 1;

                int rndfil =
                    min_fil +
                    (int)(uniform() *
                          (max_fil - min_fil + 1));

                int min_col = c[i_target] - half_bw;

                if (min_col < 0)
                    min_col = 0;

                int max_col = c[i_target] + half_bw;

                if (max_col >= side)
                    max_col = side - 1;

                int rndcol =
                    min_col +
                    (int)(uniform() *
                          (max_col - min_col + 1));

                p2 = rndfil * side + rndcol;
            }
            while (p2 == i_target || p2 == p1);
        }

         /* *********** gameto del padre **************** */

        memset(ncrorec, 0, sizeof(ncrorec));
        memset(pointrec, 0, sizeof(pointrec));

        for (crom = 0; crom < 20; crom++)
        {
            ncrorec[crom * 34] = 1;
            pointrec[crom * 34] = 0;
        }

        for (crom = 0; crom < 20; crom++)
        {
            int rndk =
                (int)(uniform() * 34);

            int rndl;

            do
            {
                rndl =
                    (int)(uniform() * NLOCI);
            }
            while (rndl == 0);

            ncrorec[(crom * 34) + rndk] = 1;
            pointrec[(crom * 34) + rndk] = rndl;
        }

        marker = 1;

        for (k = 0; k < NCRO; k++)
        {
            EE[k] = 0;

            if (ncrorec[k] == 0)
            {
                if (marker == (-1))
                    EE[k] = ~EE[k];
            }
            else
            {
                for (l = 0; l < NLOCI; l++)
                {
                    if (l != pointrec[k])
                    {
                        if (marker == (-1))
                            EE[k] |= RM[l];
                    }
                    else
                    {
                        if (marker == 1)
                        {
                            EE[k] |= RM[l];
                            marker = -marker;
                        }
                        else
                        {
                            marker = -marker;
                        }
                    }
                }
            }
        }

        double rnd = uniform();

        for (k = 0; k < NCRO; k++)
        {
            if (rnd < 0.5)
                EE[k] = ~EE[k];

            FF[k] = ~EE[k];

            sample_gm[child][k][0] =
                ((EE[k] & gm[p1][k][0]) |
                 (FF[k] & gm[p1][k][1]));
        }

        /* *********** gameto de la madre **************** */

        memset(ncrorec, 0, sizeof(ncrorec));
        memset(pointrec, 0, sizeof(pointrec));

        for (crom = 0; crom < 20; crom++)
        {
            ncrorec[crom * 34] = 1;
            pointrec[crom * 34] = 0;
        }

        for (crom = 0; crom < 20; crom++)
        {
            int rndk =
                (int)(uniform() * 34);

            int rndl;

            do
            {
                rndl =
                    (int)(uniform() * NLOCI);
            }
            while (rndl == 0);

            ncrorec[(crom * 34) + rndk] = 1;
            pointrec[(crom * 34) + rndk] = rndl;
        }

        marker = 1;

        for (k = 0; k < NCRO; k++)
        {
            EE[k] = 0;

            if (ncrorec[k] == 0)
            {
                if (marker == (-1))
                    EE[k] = ~EE[k];
            }
            else
            {
                for (l = 0; l < NLOCI; l++)
                {
                    if (l != pointrec[k])
                    {
                        if (marker == (-1))
                            EE[k] |= RM[l];
                    }
                    else
                    {
                        if (marker == 1)
                        {
                            EE[k] |= RM[l];
                            marker = -marker;
                        }
                        else
                        {
                            marker = -marker;
                        }
                    }
                }
            }
        }

        rnd = uniform();

        for (k = 0; k < NCRO; k++)
        {
            if (rnd < 0.5)
                EE[k] = ~EE[k];

            FF[k] = ~EE[k];

            sample_gm[child][k][1] =
                ((EE[k] & gm[p2][k][0]) |
                 (FF[k] & gm[p2][k][1]));
        }
    }

    free(sw_individuals);

    printf("\nGenerated %d additional offspring \n", NSAMPLE);

    return 0;
}

/* ***************************************************** */

printout()
{
	int g;
	double b;
	int g1, g2; 

	NeFibd = (100.0) / ( 2.0*(accmean(&Fibd[generations])-accmean(&Fibd[generations-100]))/(1.0-accmean(&Fibd[generations-100])) );
	NeHet = 1.0 / ( 1.0 - exp(100.0 * log(accmean(&Het[generations-100]) / accmean(&Het[generations]))) );

	fprintf(fgen, "gen     LOCI     Het            Fibd           FLH1           NumA\n");
	for(gen=generations-genmulti; gen<=generations; gen++)
	{
		fprintf(fgen, "%3d  %d  %f  %f  %f  %f\n", gen, numLOCI[gen], accmean(&Het[gen]), accmean(&Fibd[gen]), accmean(&FLH1[gen]), accmean(&Num_alleles[gen]));
		
		if (gen >= generations-100) {
			accum (&GEN, (double)gen);
			covaccum (&GENHET, (double)gen, 1.0-accmean(&Fibd[gen]));
		}
		accum(&AVE_FLH1, accmean(&FLH1[gen]));
	}
	
	b = covariance(&GENHET) / variance(&GEN);
	NeHetb = 1.0/(-2.0*b);
	NeSk2 = (4.0*NIND) / (2.0 + accmean(&SK2));

	fprintf(fgen,"\nN=%d  N.LOCI=%d  gens=%d\n", NIND, NCRO*NLOCI, generations);
	fprintf(fgen, "\nNeFibd = %f  NeHetb = %f  Sk2 = %f  NeSk2 = %f  FLH1 = %f\n", NeFibd, NeHetb, accmean(&SK2), NeSk2, accmean(&AVE_FLH1));
	fprintf(fgen, "\nlociNS1NS2 = %f  FISmax = %f  FSTmax = %f  FITmax = %f\n", lociNS1NS2, FISmax, FSTmax, FITmax);
	fprintf(fgen, "\nlociNS3NS4 = %f  FISmin = %f  FSTmin = %f  FITmin = %f\n", lociNS3NS4, FISmin, FSTmin, FITmin);

	const char *group_names[13] = {
		"UL_Corn", "UR_Corn", "LL_Corn", "LR_Corn", "U_Bord ", 
		"L_Bord ", "L_Side ", "R_Side ", "Center ", 
		"Quad_1 ", "Quad_2 ", "Quad_3 ", "Quad_4 " 
	};

	fprintf(fgen, "\n=================================================================================\n");
	fprintf(fgen, "      F-STATISTICS (PAIRWISE & GLOBAL) FOR THE 13 GROUPS\n");
	fprintf(fgen, "=================================================================================\n\n");
	
	/* Imprime los estadÃ­sticos globales corregidos */
	fprintf(fgen, ">>> 13-GROUPS GLOBAL STATS:  FST_global = %7.5f  |  FIS_global = %7.5f  |  FIT_global = %7.5f\n\n", 
			FST_global_13, FIS_global_13, FIT_global_13);

	fprintf(fgen, "PAIRWISE FST:\n          ");
	for (g2 = 0; g2 < 13; g2++) {  
		fprintf(fgen, "%s ", group_names[g2]);
	}
	fprintf(fgen, "\n");

	for (g1 = 0; g1 < 13; g1++) {  
		fprintf(fgen, "%s ", group_names[g1]);
		for (g2 = 0; g2 < 13; g2++) {
			fprintf(fgen, "%7.4f ", FST_matrix[g1][g2]);
		}
		fprintf(fgen, "\n");
	}

	FILE *ffst13 = fopen("fst_matrix_13.dat", "w");
	if (ffst13 != NULL) {
		fprintf(ffst13, "# Global_13_Stats: FST=%7.5f\tFIS=%7.5f\tFIT=%7.5f\n", FST_global_13, FIS_global_13, FIT_global_13);
		for (g1 = 0; g1 < 13; g1++) {
			for (g2 = 0; g2 < 13; g2++) {
				fprintf(ffst13, "%7.5f", FST_matrix[g1][g2]);
				if (g2 < 12) fprintf(ffst13, "\t");
			}
			fprintf(ffst13, "\n");
		}
		fclose(ffst13);
	}
}

/* ***************************************************** */

