#define _GNU_SOURCE

#include <stdio.h>
#include <sys/time.h>
#include <zlib.h>
#include <pthread.h>
#include "kseq.h" // FASTA/Q parser
#include "khashl.h" // hash table
#include "logaddexp.h"

#include "kaci.h"
#include "operations128.h"
#include "constants.h"
#include "utils.h"

char* version = "0.0.1";

uint64_t MASK_OWNS_MEMORY = 1ULL << 63;
uint64_t MASK_FAMILY_NUMBER = ~(1ULL << 63);
typedef struct {
	uint64_t familyNumber;
	uint16_t* probs;
} ref_val;

typedef struct {
	uint64_t codon;
	int codonPos;
	uint128* cdnKmer;
	int seqNumber;
	uint64_t seqLength;
	uint64_t seqStart;
	int direction;
} query_detail;

#define MAX_N_QUERIES 10
typedef struct {
	uint64_t n;
	query_detail** details;
} query_val;

typedef struct {
	char* inputFileName;
	int* outputFileName;
} thr_args;


static kh_inline khint_t kh_nohash_uint64(khint64_t key) {
	return (khint_t)key;
}

KSEQ_INIT(gzFile, gzread)
KHASHL_MAP_INIT(, kc_famref_t, kc_famref, uint64_t, uint64_t, kh_hash_uint64, kh_eq_generic)
KHASHL_MAP_INIT(, kc_ref_t, kc_ref, uint64_t, ref_val*, kh_hash_uint64, kh_eq_generic)
KHASHL_MAP_INIT(, kc_store_t, kc_store, uint64_t, uint64_t, kh_nohash_uint64, kh_eq_generic)
KHASHL_MAP_INIT(, kc_query_t, kc_query, uint64_t, query_val*, kh_hash_uint64, kh_eq_generic)

// global variables
kc_famref_t *hFamRef = NULL;			// family hash table
kc_ref_t *hRef = NULL;					// reference hash table
kc_store_t *hStore = NULL;				// hash table used for saving the reference
kc_query_t *hQueryGlobal = NULL;		// query hash table

pthread_mutex_t mutexReport = PTHREAD_MUTEX_INITIALIZER;

// global variables
int nThreads = 1;
int binBlockSize = 4096;					// size of data block for writing into binary ref file
int verbose = 0;							// set to 1 for more output 
int maxQueryCount = 10;						// max number of query encodings for a given aa kmer
int reportQueryDetails = 0;					// set to 1 to report all found query kmers
int klen = 10;          					// kmer length
char* translationTable=NULL;				// current translation table
//int refMinAAperKmer = 1;					// min number of different amino acids per family kmer
int refMaxCountForOneAAperKmer = 100;		// max count per one amino acid in family kmer
double minFamrefKmerCount = 0.0; 	
int minProteinsPerFamily = 0;				// family must have > minProteinsPerFamily to be counted
int minLinkNumber = 1;						// min number of kmers in one aa substitution group
uint64_t nRemovedUncertainKmers = 0;		// number of kmers removed by _href_remove_uncertain
char* fileNameRefOrigin=NULL;				// file name used for reference building or loading

#define N_STRLEN 256						// string length for filename manipulations
char inputDir[N_STRLEN];					// input directory for appending to file name in batch mode
char outputDir[N_STRLEN];					// output directory for writing results

uint64_t nConsensusGlobal[CODON_ALPHABET_SIZE];	// number of codon occurrences in matched query kmers
double totsumGlobal[AA_ALPHABET_SIZE];			// internal use
double likelihoodsGlobal[CODON_ALPHABET_SIZE*(AA_ALPHABET_SIZE+1)];
double decodingProbsGlobal[CODON_ALPHABET_SIZE*(AA_ALPHABET_SIZE+1)];

#define MAX_N_THREADS 256
pthread_t thrArray[MAX_N_THREADS] = {0};
char fnArray[MAX_N_THREADS][2*N_STRLEN];

/* prints the family reference hastable value 'val' */
void _print_famref_hashtable_value(FILE* fOut, void* val) {
	uint64_t* p = (uint64_t*)val;
	fprintf(fOut,"\t%lu\n",*p);
}

/* prints the reference hastable value 'val';
showExp_P flag determines how reference probabilities are going to be displayed
(-log(prob) as stored or prob) */
void _print_ref_hashtable_value(FILE* fOut, void* val, int showExp_P) {
	ref_val* v = val;
	fprintf(fOut, "\nowns_memory=%lu\n",(v->familyNumber & MASK_OWNS_MEMORY)>>63);
	fprintf(fOut, "familyNumber=%lu\n",v->familyNumber & MASK_FAMILY_NUMBER);
	if (showExp_P) {
		double pv;
		int jj;
		for (jj=0; jj<AA_ALPHABET_SIZE-1; jj++) {
			pv = exp(-(v->probs[jj]/65535.0)*9.99999);
			fprintf(fOut, "%f,",pv);
		}
		pv = exp(-(v->probs[jj]/65535.0)*9.99999);
		fprintf(fOut, "%f\n",pv);
	} else {
		int jj;
		for (jj=0; jj<AA_ALPHABET_SIZE-1; jj++) {
			fprintf(fOut, "%d,",v->probs[jj]);
		}
		fprintf(fOut, "%d\n",v->probs[jj]);
	}
}

/* prints query hastable value 'val'.
detailN is the number of the detail to print
(takes effect only when reportQueryDetails is True). */
void _print_query_hashtable_value(FILE* fOut, void* val, int detailN) {
	query_val* v = val;
	if (reportQueryDetails) {
		query_detail** details = v->details;
		fprintf(fOut,"%u,",details[detailN]->seqNumber);
		fprintf(fOut,"%lu,",details[detailN]->seqLength);
		fprintf(fOut,"%lu,",details[detailN]->seqStart);
		fprintf(fOut,"%u,",details[detailN]->direction);
		fprintf(fOut,"%s,",codon_table[details[detailN]->codon]);
		fprintf(fOut,"%u,",details[detailN]->codonPos);
		_print_cdn_kmer(fOut,*(details[detailN]->cdnKmer),klen);
	}
}

/* prints to the terminal the current family reference hash table */
void _print_famref(int showExp_P) {
	khint_t el;
	uint64_t ak;
	for (el=0; el<kh_end(hFamRef); el++) {
		if (kh_exist(hFamRef,el)) {
			ak = kh_key(hFamRef, el);
			_print_aa_kmer(stdout, ak, klen);
			_print_famref_hashtable_value(stdout, &kh_val(hFamRef,el));
		}
	}
}

/* prints to the terminal the current reference hash table */
void tst_print_ref_kmers(int showExp_P) {
	khint_t el;
	uint64_t ak;
	for (el=0; el<kh_end(hRef); el++) {
		if (kh_exist(hRef,el)) {
			ak = kh_key(hRef, el);
			_print_aa_kmer(stdout, ak, klen);
			_print_ref_hashtable_value(stdout, kh_val(hRef,el), showExp_P);
		}
	}
}

/* prints to the terminal the current query hash table */
void tst_print_query_kmers() {
	if (hQueryGlobal==NULL) return;
	khint_t el;
	uint64_t ak;
	query_val* v;
	for (el=0; el<kh_end(hQueryGlobal); el++) {
		if (kh_exist(hQueryGlobal,el)) {
			ak = kh_key(hQueryGlobal, el);
			v = kh_val(hQueryGlobal,el);
			_print_aa_kmer(stdout, ak, klen);
			printf(",");
			fprintf(stdout,"%lu,",v->n);
			for (int dd=0; dd<v->n; dd++) {
				_print_query_hashtable_value(stdout,v,dd);
			}
			printf("\n");
		}
	}
}

int _write_query_kmers(kc_query_t* hq, char* filename) {
	FILE* fOut;
	if (filename[0]=='\0') {
		fOut = stdout;
	} else {
		if ((fOut = fopen(filename, "w")) == 0) return 1;	
	}
	fprintf(fOut,"aa_kmer,family_number,seq_number,seq_length,seq_position,is_fwd,codon,codon_pos,nt_kmer\n");

	if (hq==NULL || !reportQueryDetails) {
		fclose(fOut);
		return 0;
	}

	khint_t elRef;
	uint64_t ak;
	ref_val* refv;
	khint_t elQuery;
	query_val* queryv;
	for (elQuery=0; elQuery<kh_end(hq); elQuery++) {
		if (kh_exist(hq,elQuery)) {
			ak = kh_key(hq, elQuery);
			queryv = kh_val(hq,elQuery);
			elRef = kc_ref_get(hRef,ak);
			refv = kh_val(hRef,elRef);
			for (int dd=0; dd<queryv->n; dd++) {
				_print_aa_kmer(fOut, ak, klen);
				fprintf(fOut,",%lu,",refv->familyNumber & MASK_FAMILY_NUMBER);
				_print_query_hashtable_value(fOut, kh_val(hq,elQuery),dd);
				fprintf(fOut,"\n");
			}
		}
	}
	if (fOut!=stdout) fclose(fOut);
	return 0;
}

/* puts all amino acid kmers that are found in sequence 'seq' 
and are diverse into the current family reference hash table */
void _hfamref_put_seq(const int k, const int len, const char *seq) {
	int ii, kcur;
	uint64_t ak = 0; // current kmer
    uint64_t mask = (1ULL<<k*BITS_PER_AA) - 1;
	for (ii = kcur = 0; ii < len; ii++) {
		int absent, aa = aa_char2n[(uint8_t)seq[ii]];
		if (aa < AA_ALPHABET_SIZE) { // this is an amino acid
			ak = (ak << BITS_PER_AA | aa) & mask;
			if (++kcur >= k) { // we find a k-mer
				if (_is_aa_kmer_diverse(ak,refMaxCountForOneAAperKmer,klen)) {
					khint_t itr;
					itr = kc_famref_put(hFamRef, ak, &absent);
					if (absent) {
						kh_val(hFamRef, itr) = 1;
					} else {
						kh_val(hFamRef, itr)++;
					}
				}
			}
		} else kcur = 0, ak = 0; // if there is a character outside of the aa alphabet, restart
	}
}

/* removes any kmers from the family reference that have count below 'lowestCountToKeep'*/
void _hfamref_remove_below_count(double lowestCountToKeep) {
	if (kh_size(hFamRef)==0) return;
	if (lowestCountToKeep==0) return;
	khint_t itr = 0, itrPrev = 0;
	if (verbose) printf(" kmers initially %'d,", kh_size(hFamRef));
	while (itr<kh_end(hFamRef)) {
		for (itr=itrPrev; itr<kh_end(hFamRef); ++itr) {
			if (kh_exist(hFamRef,itr)) {
				itrPrev = itr;
				if (kh_val(hFamRef,itr)<lowestCountToKeep) {
					kc_famref_del(hFamRef,itr);
					break;
				}
			}
		}
	}
	if (verbose) printf(" finally %'d. ", kh_size(hFamRef));
}

/* puts kmers from the family reference hash table 
into the reference hash table. At the same time the probability arrays
for single amino acid substitutions are calculated. Kmers that are linked to fewer than 
minLinkNumber other kmers are discarded. If the kmer is already present in the reference 
(and therefore does not have an unambiguous family assignment) its family id is set 
to 0 (uncertain) */
uint64_t _href_put_hfamref(uint64_t famNumber) {
	khint_t itrCurrent, itrMutated; // iterators for the family reference table
	khint_t itr;					// iterator for the reference tables
	uint64_t ak, ak_mut;
	int absent;
	uint64_t mask = (1ULL<<klen*BITS_PER_AA) - 1;

	uint64_t nAdded=0;
	// for each kmer in famref hash table
	for (itrCurrent=0; itrCurrent<kh_end(hFamRef); itrCurrent++) {
		if (kh_exist(hFamRef,itrCurrent)) {
			ak = kh_key(hFamRef, itrCurrent);
			// for each amino acid position in the kmer
			for (int ii=0; ii<klen; ii++) {
				uint16_t* probArray = malloc(AA_ALPHABET_SIZE*sizeof(uint16_t));
				for (int jj=0; jj<AA_ALPHABET_SIZE; jj++) probArray[jj] = 1;
				double iniSum = 0;
				for (int jj=0; jj<AA_ALPHABET_SIZE; jj++) iniSum += probArray[jj];
				char currentAA = aa_table[((ak << ii*BITS_PER_AA) & mask) >> (klen-1)*BITS_PER_AA];
				// for all possible substitutions
				for (int jj=0; jj<AA_ALPHABET_SIZE; jj++) {
					if (aa_table[jj]!=currentAA) {
						ak_mut = _mutate_aa_kmer(ak, ii+1, aa_table[jj], klen);
						itrMutated = kc_famref_get(hFamRef,ak_mut); // look up this kmer in the family reference table
						if (itrMutated<kh_end(hFamRef)) {
							probArray[jj] += kh_val(hFamRef, itrMutated);
						}
					} else {
						probArray[jj] += kh_val(hFamRef, itrCurrent);
					}
				}
				// normalize the probability array
				double sum=0;
				double p;
				double maxP = 9.99999;
				for (int jj=0; jj<AA_ALPHABET_SIZE; jj++) sum += probArray[jj];
				if ((sum-iniSum)<minLinkNumber) {
					free(probArray);
					continue;
				}
				for (int jj=0; jj<AA_ALPHABET_SIZE; jj++) {
					p = -log(probArray[jj]/sum);
					if (p>maxP) p = maxP;
					p = round((p/maxP)*65535);
					probArray[jj] = (uint16_t)p;
				}
				// change the current amino acid to ?
				ak_mut = _mutate_aa_kmer(ak, ii+1, aa_table[AA_ALPHABET_SIZE], klen);
				// and put this kmer in the ref hash table
				itr = kc_ref_put(hRef,ak_mut,&absent);
				if (absent){
					ref_val* v = malloc(sizeof *v);
					v->familyNumber = famNumber | MASK_OWNS_MEMORY;
					v->probs = probArray;
					kh_val(hRef, itr) = v;
					nAdded++;
				} else {
					// set family number to 0 if this kmers is present in more than one family
					if ((kh_val(hRef,itr)->familyNumber & MASK_FAMILY_NUMBER) != famNumber) {
						kh_val(hRef,itr)->familyNumber = MASK_OWNS_MEMORY;
					}
					free(probArray);
				}
			}
		}
	}
	if (verbose) printf("Added %ld kmers to the reference table. ", nAdded);
	return nAdded;
}

/* removes uncertain kmers from the reference hash table 
(kmers with family id equal to 0) */
void _href_remove_uncertain() {
	nRemovedUncertainKmers = 0;
	if (hRef==NULL) return; 
	if (kh_size(hRef)==0) return;
	nRemovedUncertainKmers = kh_size(hRef);
	khint_t itr = 0, itrPrev = 0;
	if (verbose) printf("Removing uncertain kmers from reference: kmers initially %'d,", kh_size(hRef));
	while (itr<kh_end(hRef)) {
		for (itr=itrPrev; itr<kh_end(hRef); ++itr) {
			if (kh_exist(hRef,itr)) {
				itrPrev = itr;
				if ((kh_val(hRef,itr)->familyNumber & MASK_FAMILY_NUMBER)==0) {
					if ((kh_val(hRef,itr)->familyNumber & MASK_OWNS_MEMORY)!=0) free(kh_val(hRef,itr)->probs);
					free(kh_val(hRef,itr));
					kc_ref_del(hRef,itr);
					break;
				}
			}
		}
	}
	nRemovedUncertainKmers -= kh_size(hRef);
	if (verbose) printf(" finally %'d.\n", kh_size(hRef));
}

/* puts a kmer into the specified query hash table.
If the kmer is already present, its count is incremented.
If the count exceeds maxQueryCount, the function returns 
without making any changes. */
int _hquery_put_ak(kc_query_t* hq, uint64_t ak, query_detail* qd) {	
	int absent;
	khint_t query_itr = kc_query_put(hq, ak, &absent);
	if (absent) {
		query_val* v = malloc(sizeof *v);
		if (reportQueryDetails) {
			v->details = malloc(sizeof(query_detail**)*maxQueryCount);
			v->details[0] = qd;
		}
		v->n = 1;
		kh_val(hq, query_itr) = v;
		return 1;
	} else {
		if (kh_val(hq, query_itr)->n<maxQueryCount) {
			query_val* v = kh_val(hq, query_itr);
			if (reportQueryDetails) {
				query_detail** details = v->details;
				details[v->n] = qd;
			}
			v->n += 1;
			return 1;
		} else {
			return 0;
		}
	}
}

/* adds the substitution probabilities to the likelihood for the specified codon */
void _add_to_likelihoods(uint64_t cdn, uint16_t* probs, 
						 double* likelihoods, double* totsum, uint64_t* nConsensus) {
	for (int jj=0; jj<AA_ALPHABET_SIZE; jj++) { // calculations for all possible aa substitutions
		double p = (probs[jj]/65535.0)*9.99999;
		//printf("At %s for %c prob = %f\n", codon_table[cdn], aa_table[jj], p);
		likelihoods[cdn*(AA_ALPHABET_SIZE+1)+jj] += p;	// aka P(A|C_z_i)
		totsum[jj] = logaddexp(-p,totsum[jj]); 			// aka P(A)
	}
	nConsensus[cdn]++;
}

/* puts any amino acid kmers found in the query sequebce 'ntseq'
(and present in the reference hash table)
into the query hash table 'hQuery'. Also calculates the likelihoods 
for the encoding codons based on probabilities from the reference. */
void _hquery_put_seq(const int seqStart, const int seqEnd, const char *ntseq, uint32_t seqN, 
               	  const char* tr_table, const int isFwd, 
				  kc_query_t* hQuery, double* likelihoods, double* totsum, uint64_t* nConsensus) 
{
	uint64_t ak = 0;	// bit-encoded amino acid kmer 
	uint64_t ak_mut;	// mutated kmer
	uint128 ck = {0,0};	// bit-encoded codon kmer
	uint128 ck_copy = {0,0};	// bit-encoded codon kmer for reporting
	uint8_t nt1 = 0;	// first letter of the current codon
	uint8_t nt2 = 0;	// second letter of the current codon
	uint8_t nt3 = 0;	// third letter of the current codon
	int kcur = 0;  		// current length of the amino acid kmer
	uint64_t maskA = (1ULL<<klen*BITS_PER_AA) - 1;
	uint128 maskC;
	uint64_t maskLastCodon = (1ULL<<BITS_PER_CODON) - 1;
	if (klen>64/BITS_PER_CODON) {
		maskC.low = UINT64_MAX;
		maskC.high = (1ULL<<(klen*BITS_PER_CODON-64)) - 1;
	} else {
		maskC.low = (1ULL<<(klen*BITS_PER_CODON)) - 1;
		maskC.high = 0ULL;
	}
	uint128 codonXXX = {CODON_ALPHABET_SIZE,0};
	uint64_t shiftA = (klen-1)*BITS_PER_AA;
	uint64_t shiftC = (klen-1)*BITS_PER_CODON;
	int ii;
	for (ii = seqStart; ii < seqEnd; ii+=3) {
		nt1 = (uint8_t)ntseq[ii];
		nt2 = (uint8_t)ntseq[ii+1];
		nt3 = (uint8_t)ntseq[ii+2];
        uint64_t codon;
        if (isFwd) {
		    codon = (nt_char2n[nt1] << 4) + (nt_char2n[nt2] << 2) + nt_char2n[nt3];
        } else {
			codon = (_rc(nt_char2n[nt3]) << 4) + (_rc(nt_char2n[nt2]) << 2) + _rc(nt_char2n[nt1]);
        }
		if (codon < CODON_ALPHABET_SIZE) { // this is a codon
			uint64_t aa = aa_char2n[(uint8_t)tr_table[codon]];
			if (aa < AA_ALPHABET_SIZE) {
				uint128 codon128 = {codon,0};
                if (isFwd) {
					ck = _and_128(_or_128(_shift_left_128(ck,BITS_PER_CODON), codon128), maskC);
					if (reportQueryDetails) ck_copy = _and_128(_or_128(_shift_left_128(ck_copy,BITS_PER_CODON), codon128), maskC);
					//ck = (ck << BITS_PER_CODON | codon) & maskC;
					ak = (ak << BITS_PER_AA | aa) & maskA;
                } else {
					ck = _and_128(_or_128(_shift_right_128(ck,BITS_PER_CODON), _shift_left_128(codon128,shiftC)), maskC);
					if (reportQueryDetails) ck_copy = _and_128(_or_128(_shift_right_128(ck_copy,BITS_PER_CODON), _shift_left_128(codon128,shiftC)), maskC);
					//ck = ((ck >> BITS_PER_CODON) | (codon << shiftC)) & maskC;
					ak = ((ak >> BITS_PER_AA) | (aa << shiftA)) & maskA;
                }
				if (++kcur >= klen) { // we find a k-mer
					// for each possible uncertainty
					for (int jj=0; jj<klen; jj++) {
						ak_mut = _mutate_aa_kmer(ak,jj+1,aa_table[AA_ALPHABET_SIZE], klen);
						khint_t el = kc_ref_get(hRef,ak_mut); // look up this kmer in the reference table
						if (el<kh_end(hRef)) { // if present
							uint64_t uncertain_codon = _shift_right_128(ck,(klen-jj-1)*BITS_PER_CODON).low & maskLastCodon;
							if (uncertain_codon<CODON_ALPHABET_SIZE) { // skip for masked codons
								query_detail* qd = 0;
								if (reportQueryDetails) {
									qd = malloc(sizeof *qd);
									qd->codon = _shift_right_128(ck_copy,(klen-jj-1)*BITS_PER_CODON).low & maskLastCodon;
									qd->cdnKmer = malloc(sizeof(uint128));
									qd->codonPos = jj*3+1;
									qd->direction = isFwd;
									qd->cdnKmer->high = ck_copy.high;
									qd->cdnKmer->low = ck_copy.low;
									qd->seqNumber = seqN;
									qd->seqLength = seqEnd;
									qd->seqStart = ii-(klen-1)*3+1;
								}
								if (_hquery_put_ak(hQuery, ak_mut, qd)) { // query hash table keeps track of how many times this kmer has popped up
									uint16_t* probs = kh_val(hRef,el)->probs;
									// this is where the calculation happens
									_add_to_likelihoods(uncertain_codon, probs, likelihoods, totsum, nConsensus);
									// mask the uncertain codon by setting it to XXX so that it cannot be counted again
									ck = _mutate_cdn_kmer(ck,jj+1,codonXXX,klen);
								} else{
									if (reportQueryDetails) {
										free(qd->cdnKmer);
										free(qd);
									}
								}
							}
						}
					}
				}
			} else kcur = 0, ak = 0, memset(&ck, 0, sizeof(ck)); // if there is a character outside of the aa alphabet, restart
		} else kcur = 0, ak = 0, memset(&ck, 0, sizeof(ck)); // if there is a character outside of the nt alphabet, restart
	}
}

/* calls _hquery_put_seq above with a specified shift 'seqStart (0, 1, or 2)
and for fwd or rev sequence */
int _hquery_put_seq_uni_orf(char* fn, const int seqStart, const int is_fwd, 
							kc_query_t *hQuery, double* likelihoods, double* totsum, uint64_t* nConsensus) {
	gzFile fp;
	kseq_t *ks;
	if ((fp = gzopen(fn, "r")) == 0) return 1;
	ks = kseq_init(fp);

	if (translationTable==NULL) return 2;
	if (hRef==NULL) return 3;

	uint32_t seqNumber=1;
	while (kseq_read(ks) >= 0) {
		_hquery_put_seq(seqStart, ((ks->seq.l-seqStart)/3)*3, ks->seq.s, seqNumber, 
								   translationTable, is_fwd, 
								   hQuery, likelihoods, totsum, nConsensus);
		seqNumber++;
	}

	kseq_destroy(ks);
	gzclose(fp);
	return 0;
}

/* calculates decoding probabilitiesfor based on the provided likelihoods, totsum, and nConsensus */
void _calculate_decoding_probs(double* likelihoods, double* totsum, uint64_t* nConsensus, double* decodingProbs) {
	// normalize likelihoods
	double m;
	double denominator[AA_ALPHABET_SIZE];
	memset(denominator,0,sizeof(double)*AA_ALPHABET_SIZE);
	double sum_nConsensus=0;
	for (int cc=0; cc<CODON_ALPHABET_SIZE; cc++) sum_nConsensus += nConsensus[cc];
	for (int cc=0; cc<CODON_ALPHABET_SIZE; cc++) {
		m = nConsensus[cc];
		for (int jj=0; jj<AA_ALPHABET_SIZE+1; jj++) {
        	if (m==0) {
            	likelihoods[cc*(AA_ALPHABET_SIZE+1)+jj] = log(1.0/21.0);
			} else {
				if (jj!=AA_ALPHABET_SIZE) {
					//  for the 20 aa models
                	denominator[jj] = m * totsum[jj];
	                likelihoods[cc*(AA_ALPHABET_SIZE+1)+jj] = log(1.0/21.0) - likelihoods[cc*(AA_ALPHABET_SIZE+1)+jj] - denominator[jj];
				} else {
					// for the nonspecific model
                	likelihoods[cc*(AA_ALPHABET_SIZE+1)+jj] = log(1.0/21.0) - m * log(sum_nConsensus);
				}  
			}
		}
	}

    // compute decoding probabilities
	double post_denoms[CODON_ALPHABET_SIZE];
	for (int cc=0; cc<CODON_ALPHABET_SIZE; cc++) post_denoms[cc] = log(0);
	for (int cc=0; cc<CODON_ALPHABET_SIZE; cc++) {
		for (int jj=0; jj<AA_ALPHABET_SIZE+1; jj++) {
			post_denoms[cc] = logaddexp(post_denoms[cc],likelihoods[cc*(AA_ALPHABET_SIZE+1)+jj]);
		}
	}

	for (int cc=0; cc<CODON_ALPHABET_SIZE; cc++) {
		for (int jj=0; jj<AA_ALPHABET_SIZE+1; jj++) {
			//dcp[cc*(AA_ALPHABET_SIZE+1)+jj] = lkl[cc*(AA_ALPHABET_SIZE+1)+jj];
			decodingProbs[cc*(AA_ALPHABET_SIZE+1)+jj] = exp(likelihoods[cc*(AA_ALPHABET_SIZE+1)+jj] - post_denoms[cc]);
		}
	}

}

void tst_calculate_decoding_probs() {
	_calculate_decoding_probs(likelihoodsGlobal, totsumGlobal, nConsensusGlobal, decodingProbsGlobal);
}

/* write provided decoding probabilities to a json file */
int _write_decoding_probs(char* fnOut, double* decodingProbs, uint64_t* nConsensus) {
	FILE* fOut;
	if ((fOut = fopen(fnOut, "w")) == 0) return 1;
	struct timespec ts;
    timespec_get(&ts, TIME_UTC);

	fprintf(fOut, "{\n");
	fprintf(fOut,"    \"version\": \"%s\",\n", version);
	fprintf(fOut,"    \"timestamp\": %.6f,\n", ts.tv_sec + ts.tv_nsec/1e9);
	fprintf(fOut,"    \"referenceFile\": \"%s\",\n", fileNameRefOrigin);
	fprintf(fOut,"    \"inputTranslationTable\": \"%s\",\n", translationTable);
	fprintf(fOut,"    \"maxQueryCount\": \"%d\",\n", maxQueryCount);
	for (int cc=0; cc<CODON_ALPHABET_SIZE; cc++) {
		fprintf(fOut, "    \"%s\": {\n", codon_table[cc]);
		for (int jj=0; jj<AA_ALPHABET_SIZE+1; jj++) {
			fprintf(fOut, "        \"%c\": %.10f,\n", aa_table[jj], decodingProbs[cc*(AA_ALPHABET_SIZE+1)+jj]);
		}
		fprintf(fOut, "        \"n_\": %ld\n", nConsensus[cc]);
		if (cc!=CODON_ALPHABET_SIZE-1) {
			fprintf(fOut, "    },\n");
		} else {
			fprintf(fOut, "    }\n");
		}
	}
	fprintf(fOut, "}\n");

	fclose(fOut);

	return 0;
}

/* write the current global decoding probabilities to a json file */
int tst_write_decoding_probs(char* fnOut) {
	return _write_decoding_probs(fnOut, decodingProbsGlobal, nConsensusGlobal);
}

/* frees resources associated with the family reference hash table */
void _hfamref_destroy() {
	if (hFamRef==NULL) return;
	kc_famref_destroy(hFamRef);
	hFamRef = NULL;
}

/* frees resources associated with the reference hash table */
void _href_destroy() {
	if (hRef==NULL) return; 
	khint_t el;
	for (el=0; el<kh_end(hRef); el++) {
		if (kh_exist(hRef,el)) {
			if ((kh_val(hRef,el)->familyNumber & MASK_OWNS_MEMORY) != 0) free(kh_val(hRef,el)->probs);
			free(kh_val(hRef,el));
		}
	}
	kc_ref_destroy(hRef);
	hRef = NULL;
}

/* frees resources associated with the store hash table
(it is used while writing the reference to the hard drive) */
void _hstore_destroy() {
	if (hStore==NULL) return; 
	kc_store_destroy(hStore);
	hStore = NULL;
}

/* frees resources associated with the query hash table */
void _hquery_destroy(kc_query_t* hq) {
	if (hq==NULL) return; 
	khint_t el;
	for (el=0; el<kh_end(hq); el++) {
		if (kh_exist(hq,el)) {
			query_val* v = kh_val(hq,el);
			uint64_t nDetails = v->n;
			query_detail** details = v->details;
			if (reportQueryDetails) {
				for (int dd=0; dd<nDetails; dd++) {
					query_detail* qd = details[dd];
					free(qd->cdnKmer);
					free(qd);
				}
				free(details);
			}
			free(kh_val(hq,el));
		}
	}
	kc_query_destroy(hq);
	hq = NULL;
}

/* reset everything involved in the query calculation */
void tst_reset_query() {
	_hquery_destroy(hQueryGlobal);
	hQueryGlobal = kc_query_init();
	memset(nConsensusGlobal,0,sizeof(uint64_t)*CODON_ALPHABET_SIZE);
	for (int jj=0; jj<AA_ALPHABET_SIZE; jj++) totsumGlobal[jj] = log(0);
	memset(likelihoodsGlobal,0,sizeof(double)*CODON_ALPHABET_SIZE*(AA_ALPHABET_SIZE+1));
	memset(decodingProbsGlobal,0,sizeof(double)*CODON_ALPHABET_SIZE*(AA_ALPHABET_SIZE+1));
}

/* calls _hquery_put_seq_uni_orf on the forward orf (no shift) */
int tst_process_query_orf1f(char* filename) {
    return _hquery_put_seq_uni_orf(filename, 0, 1, hQueryGlobal, likelihoodsGlobal, totsumGlobal, nConsensusGlobal);
}

/* calls _hquery_put_seq_uni_orf on the reverse orf (no shift) */
int tst_process_query_orf1r(char* filename) {
    return _hquery_put_seq_uni_orf(filename, 0, 0, hQueryGlobal, likelihoodsGlobal, totsumGlobal, nConsensusGlobal);
}

/* calls _hquery_put_seq_uni_orf on the forward orf (1nt shift) */
int tst_process_query_orf2f(char* filename) {
    return _hquery_put_seq_uni_orf(filename, 1, 1, hQueryGlobal, likelihoodsGlobal, totsumGlobal, nConsensusGlobal);
}

/* calls _hquery_put_seq_uni_orf on the reverse orf (1nt shift) */
int tst_process_query_orf2r(char* filename) {
    return _hquery_put_seq_uni_orf(filename, 1, 0, hQueryGlobal, likelihoodsGlobal, totsumGlobal, nConsensusGlobal);
}

/* calls _hquery_put_seq_uni_orf on the forward orf (2nt shift) */
int tst_process_query_orf3f(char* filename) {
    return _hquery_put_seq_uni_orf(filename, 2, 1, hQueryGlobal, likelihoodsGlobal, totsumGlobal, nConsensusGlobal);
}

/* calls _hquery_put_seq_uni_orf on the reverse orf (2nt shift) */
int tst_process_query_orf3r(char* filename) {
    return _hquery_put_seq_uni_orf(filename, 2, 0, hQueryGlobal, likelihoodsGlobal, totsumGlobal, nConsensusGlobal);
}

void tst_get_n_consensus(uint64_t* ns) {
	memcpy(ns, nConsensusGlobal, sizeof(uint64_t)*CODON_ALPHABET_SIZE);
}

/* copies decoding probabilities into 'probs'  */
void tst_get_decoding_probs(double* probs) {
	memcpy(probs, decodingProbsGlobal, sizeof(double)*CODON_ALPHABET_SIZE*(AA_ALPHABET_SIZE+1));
}

int set_kmer_size(int klength) {
    if (klength<3 || klength>12) {
        return 1;
    } else {
        klen = klength;
        return 0;
    }
}

int set_max_query_count(int v) {
	if (v<1) {
		return 1;
	} else {
		maxQueryCount = v;
		return 0;
	}
}

void set_report_query_details(int value) {
	reportQueryDetails = value;
}

/*int set_ref_min_aa_per_kmer(int v) {
	if (v<1) {
		return 1;
	} else {
		refMinAAperKmer = v;
		return 0;
	}
}*/

int set_max_familyref_count_for_one_aa_per_kmer(int v) {
	if (v<2) {
		return 1;
	} else {
		refMaxCountForOneAAperKmer = v;
		return 0;
	}
}

int set_min_familyref_count_for_kmer(int v) {
	if (v<0) {
		return 1;
	} else {
		minFamrefKmerCount = v;
		return 0;
	}
}

int set_min_proteins_per_family(int value) {
	if (value<0) {
		return 1;
	} else {
		minProteinsPerFamily = value;
		return 0;
	}
}

int set_min_link_number(int value) {
	if (value<0) {
		return 1;
	} else {
		minLinkNumber = value;
		return 0;
	}
}

void set_verbose(int value) {
	verbose = value;
}

int set_bin_block_size(int value) {
	if (value<=0) {
		return 1;
	} else {
		binBlockSize = value;
		return 0;
	}
}

void set_translation_table(const char* translation_table) {
	if (translationTable!=NULL) free(translationTable);
    translationTable = strdup(translation_table);
}

int set_input_dir(char* dirPath) {
	if (strlen(dirPath)>sizeof(inputDir)-2) return 1;
	_trim_string(dirPath,inputDir);
	if (inputDir[0]=='\0') return 0;
	if (inputDir[strlen(inputDir)-1]!='/') {
		if (strlen(inputDir)>sizeof(inputDir)-2) return 1;
		strcat(inputDir,"/");
	}
	return 0;
}

int set_output_dir(char* dirPath) {
	if (strlen(dirPath)>sizeof(outputDir)-2) return 1;
	_trim_string(dirPath,outputDir);
	if (outputDir[0]=='\0') return 0;
	if (outputDir[strlen(outputDir)-1]!='/') {
		if (strlen(outputDir)>sizeof(outputDir)-2) return 1;
		strcat(outputDir,"/");
	}
	return 0;
}

int set_n_threads(int nThr) {
	if (nThr<=0 || nThr>=MAX_N_THREADS) {
		return 1;
	} else {
		nThreads = nThr;
		return 0;
	}
}

unsigned long get_ref_size() {
	return kh_size(hRef);
}

int build_ref(char* filename, int saveFamilyList) {
	FILE* fFamList=0;
	if (saveFamilyList) {
		char nameFamList[N_STRLEN+15];
		strcpy(nameFamList, outputDir);
		strcat(nameFamList, "families.txt");
		if ((fFamList = fopen(nameFamList, "w")) == 0 ) return 2;
	}

	// copy file name
	if (fileNameRefOrigin!=NULL) free(fileNameRefOrigin);
    fileNameRefOrigin = strdup(filename);
	// family reference hash table
	_hfamref_destroy();
	hFamRef = kc_famref_init();
	// query hash table
	tst_reset_query();
	// reference hash table
	_href_destroy();
	hRef = kc_ref_init();
	gzFile fp;
	kseq_t *ks;
	if ((fp = gzopen(filename, "r")) == 0) return 1;
	ks = kseq_init(fp);
	char* fullDescription = malloc(sizeof(char)*128);
	char* prevFamName = malloc(sizeof(char)*128);
	memset(prevFamName, 0, sizeof(char)*128);
	char* famName = malloc(sizeof(char)*128);
	uint64_t nn=1ULL;
	int nProtSeqs = 0;
	uint64_t nAddedKmers = 0;
	while (kseq_read(ks) >= 0) {
		strcpy(fullDescription, ks->name.s);	
		// in case the full description got truncated by a separator
		if (ks->comment.s!=NULL) strcat(fullDescription, ks->comment.s);	
		_parse_family_name(fullDescription, famName);
		//if (famName[0]!='L') continue;
		if (strcmp(famName,prevFamName)!=0) {
			if (kh_size(hFamRef)!=0) {
				if (verbose) printf("Processing family %s ... ", prevFamName);
				if (nProtSeqs>=minProteinsPerFamily) {
					_hfamref_remove_below_count(minFamrefKmerCount);
					nAddedKmers = _href_put_hfamref(nn);
					if (verbose) printf("Done.\n");
					if (saveFamilyList) fprintf(fFamList, "%s,%ld,%ld\n", prevFamName, nn, nAddedKmers);
				} else {
					if (verbose) printf(" skipping (too few proteins in this family).\n");
					if (saveFamilyList) fprintf(fFamList, "%s,%d,%d\n", prevFamName, -1, 0);
				}
				nn++;
				nProtSeqs = 0;
				_hfamref_destroy();
				hFamRef = kc_famref_init();
			}
			strcpy(prevFamName, famName);
		}
		_hfamref_put_seq(klen, ks->seq.l, ks->seq.s);
		nProtSeqs++;
	}
	if (verbose) printf("Processing family %s ... ", prevFamName);
	if (nProtSeqs>=minProteinsPerFamily) {
		_hfamref_remove_below_count(minFamrefKmerCount);
		nAddedKmers = _href_put_hfamref(nn);
		if (verbose) printf("Done.\n");
		if (saveFamilyList) fprintf(fFamList, "%s,%ld,%ld\n", famName, nn, nAddedKmers);
	} else {
		if (verbose) printf(" skipping (too few proteins in this family).\n");
		if (saveFamilyList) fprintf(fFamList, "%s,%d,%d\n", famName, -1, 0);
	}

	free(fullDescription);
	free(famName);
	free(prevFamName);
	kseq_destroy(ks);
	gzclose(fp);
	if (saveFamilyList) fclose(fFamList);

	_href_remove_uncertain();

	return 0;
}

int save_parameters_ref(char* fnOut) {
	FILE* fOut;
	if ((fOut = fopen(fnOut, "w")) == 0) return 1;

	fprintf(fOut, "{\n");
	fprintf(fOut,"    \"version\": \"%s\",\n", version);
	fprintf(fOut,"    \"inputFile\": \"%s\",\n", fileNameRefOrigin);
	fprintf(fOut,"    \"kmerLength\": \"%d\",\n", klen);
	fprintf(fOut,"    \"maxCountForAAperKmer\": \"%d\",\n", refMaxCountForOneAAperKmer);
	fprintf(fOut,"    \"minFamilyrefKmerCount\": \"%f\",\n", minFamrefKmerCount);
	fprintf(fOut,"    \"minProteinsPerFamily\": \"%d\",\n", minProteinsPerFamily);
	fprintf(fOut,"    \"minLinkNumber\": \"%d\",\n", minLinkNumber);
	fprintf(fOut,"    \"nRemovedUncertainKmers\": \"%ld\",\n", nRemovedUncertainKmers);
	fprintf(fOut,"    \"finalReferenceSize\": \"%d\"\n", kh_size(hRef));
	fprintf(fOut, "}\n");

	fclose(fOut);

	return 0;
}

int save_ref_txt(char* fnOut) {
	if (verbose) printf("Saving reference... ");
	gzFile fpOut;
	if (hRef==NULL) return 2;
	if ((fpOut = gzopen(fnOut, "w")) == 0) return 1;

	// create another hash table to shuffle the reference
	uint64_t ak;
	int absent;
	khint_t itrRef, itrStore;
	_hstore_destroy();
	hStore = kc_store_init();
	for (itrRef=0; itrRef<kh_end(hRef); itrRef++) {
		if (kh_exist(hRef,itrRef)) {
			ak = kh_key(hRef,itrRef);
			itrStore = kc_store_put(hStore, ak, &absent);
			if (absent) kh_val(hStore, itrStore) = ak;
		}
	}
	if (kh_size(hRef)!=kh_size(hStore)) {
		printf("Something went wrong!\n");
		return 3;
	}

	char* aa_kmer = malloc((klen+1)*sizeof(char));
	for (itrStore=0; itrStore<kh_end(hStore); itrStore++) {
		if (kh_exist(hStore,itrStore)) {
			ak = kh_val(hStore,itrStore);
			itrRef = kc_ref_get(hRef,ak);
			// write kmer description
			char* descr;
			if (asprintf(&descr, ">%ld;%hd,%hd,%hd,%hd,%hd,%hd,%hd,%hd,%hd,%hd,%hd,%hd,%hd,%hd,%hd,%hd,%hd,%hd,%hd,%hd\n", \
			 					 kh_val(hRef,itrRef)->familyNumber & MASK_FAMILY_NUMBER,
								 kh_val(hRef,itrRef)->probs[0],
								 kh_val(hRef,itrRef)->probs[1],
								 kh_val(hRef,itrRef)->probs[2],
								 kh_val(hRef,itrRef)->probs[3],
								 kh_val(hRef,itrRef)->probs[4],
								 kh_val(hRef,itrRef)->probs[5],
								 kh_val(hRef,itrRef)->probs[6],
								 kh_val(hRef,itrRef)->probs[7],
								 kh_val(hRef,itrRef)->probs[8],
								 kh_val(hRef,itrRef)->probs[9],
								 kh_val(hRef,itrRef)->probs[10],
								 kh_val(hRef,itrRef)->probs[11],
								 kh_val(hRef,itrRef)->probs[12],
								 kh_val(hRef,itrRef)->probs[13],
								 kh_val(hRef,itrRef)->probs[14],
								 kh_val(hRef,itrRef)->probs[15],
								 kh_val(hRef,itrRef)->probs[16],
								 kh_val(hRef,itrRef)->probs[17],
								 kh_val(hRef,itrRef)->probs[18],
								 kh_val(hRef,itrRef)->probs[19])<=0) return 2;
			gzwrite(fpOut,descr,strlen(descr));
			free(descr);
			// write kmer
			uint64_t mask = (1ULL<<klen*BITS_PER_AA) - 1;
			ak = kh_key(hRef,itrRef);
			for (int ii=0; ii<klen; ii++) {
				aa_kmer[ii] = aa_table[(ak & mask) >> (klen-1)*BITS_PER_AA];
				ak = ak << BITS_PER_AA;
			}
			aa_kmer[klen] = '\n';
			gzwrite(fpOut, aa_kmer, klen+1);
		}
	}

	_hstore_destroy();
	free(aa_kmer);
	gzclose(fpOut);
	if (verbose) printf("Done.\n");
	return 0;
}

int save_ref_bin(char* fnOut) {
	FILE* fpOut;
	if (hRef==NULL) return 2;
	if ((fpOut = fopen(fnOut, "wb")) == 0) return 1;
	fwrite(&klen, sizeof(klen),1,fpOut); // write the kmer size
	unsigned long refSize = get_ref_size();

	// create another hash table to shuffle the reference
	uint64_t ak;
	int absent;
	khint_t itrRef, itrStore;
	_hstore_destroy();
	hStore = kc_store_init();
	for (itrRef=0; itrRef<kh_end(hRef); itrRef++) {
		if (kh_exist(hRef,itrRef)) {
			ak = kh_key(hRef,itrRef);
			itrStore = kc_store_put(hStore, ak, &absent);
			if (absent) kh_val(hStore, itrStore) = ak;
		}
	}
	if (kh_size(hRef)!=kh_size(hStore)) {
		printf("Something went wrong!\n");
		return 3;
	}

	int blockSize = binBlockSize;

	uint64_t* buffAk = malloc(sizeof(uint64_t)*blockSize);
	uint64_t* buffFamN = malloc(sizeof(uint64_t)*blockSize);
	uint16_t* buffProbsPointers[blockSize];

	itrStore=0;
	int ii=0;
	if (blockSize>refSize) blockSize = refSize;
	while (refSize>0) {
		if (kh_exist(hStore,itrStore)) {
			ak = kh_val(hStore,itrStore);
			itrRef = kc_ref_get(hRef,ak);
			buffAk[ii] = ak;
			buffFamN[ii] = kh_val(hRef,itrRef)->familyNumber & MASK_FAMILY_NUMBER;
			buffProbsPointers[ii] = kh_val(hRef,itrRef)->probs;
			ii++;
		}
		itrStore++;
		if (itrStore>=kh_end(hStore)) blockSize = ii;
		if (ii==blockSize) {
			// write this block
			fwrite(&blockSize, sizeof(blockSize),1,fpOut);
			fwrite(buffAk, sizeof(uint64_t), blockSize, fpOut);
			fwrite(buffFamN, sizeof(uint64_t), blockSize, fpOut);
			for (int jj=0; jj<blockSize; jj++) {
				fwrite(buffProbsPointers[jj], sizeof(uint16_t)*AA_ALPHABET_SIZE, 1, fpOut);
			}
			refSize -= blockSize;
			ii = 0;
		}
	}

	_hstore_destroy();
	free(buffAk);
	free(buffFamN);
	fclose(fpOut);
	return 0;
}

int load_ref_txt(char* fnIn) {
	// reset hash tables
	tst_reset_query();
	if (hRef!=NULL) _href_destroy();
	hRef = kc_ref_init();

	gzFile fp;
	kseq_t *ks;
	if ((fp = gzopen(fnIn, "r")) == 0) return 1;
	ks = kseq_init(fp);
	char* fullDescription = malloc(sizeof(char)*1024);
	int famNumber;

	int ii, kcur;
	uint64_t ak = 0; // current kmer
    uint64_t mask = (1ULL<<klen*BITS_PER_AA) - 1;
	while (kseq_read(ks) >= 0) {
		strcpy(fullDescription, ks->name.s);
		// in case the full description got truncated by a separator
		if (ks->comment.s!=NULL) strcat(fullDescription, ks->comment.s);
		uint16_t* probArray = malloc(AA_ALPHABET_SIZE*sizeof(uint16_t));	
		_parse_ref_entry(fullDescription, &famNumber, probArray);
		// look for a kmer
		for (ii = kcur = 0; ii < ks->seq.l; ii++) {
			int absent, aa = aa_char2n[(uint8_t)ks->seq.s[ii]];
			if (aa <= AA_ALPHABET_SIZE) { // this is an amino acid or '?'
				ak = (ak << BITS_PER_AA | aa) & mask;
				if (++kcur >= klen) { // we find a k-mer
					khint_t itr;
					itr = kc_ref_put(hRef, ak, &absent);
					if (absent) {
						ref_val* v = malloc(sizeof *v);
						v->familyNumber = famNumber | MASK_OWNS_MEMORY;
						v->probs = probArray;
						kh_val(hRef, itr) = v;
					}
				}
			} else kcur = 0, ak = 0;
		}

	}

	free(fullDescription);
	kseq_destroy(ks);
	gzclose(fp);

	// copy file name
	if (fileNameRefOrigin!=NULL) free(fileNameRefOrigin);
    fileNameRefOrigin = strdup(fnIn);

	return 0;
}

int load_ref_bin(char* fnIn) {
	FILE* fpIn;
	if ((fpIn = fopen(fnIn, "rb")) == 0) return 1;
	int kmerSize;
	int ret = fread(&kmerSize, sizeof(kmerSize),1,fpIn); // read the kmer size
	if (ret!=1) return 2;
	if (kmerSize!=klen) return 3;

	int blockSize = binBlockSize;
	ret = fread(&blockSize, sizeof(blockSize),1,fpIn); // read the size of the next data block
	_href_destroy();
	hRef = kc_ref_init();
	
	uint64_t ak;
	uint64_t famN;

	uint64_t* buffAk = malloc(blockSize*sizeof(ak));
	uint64_t* buffFamN = malloc(blockSize*sizeof(famN));
	while(ret!=0) {
		// read in data
		ret = fread(buffAk, sizeof(ak), blockSize, fpIn);
		ret = fread(buffFamN, sizeof(famN), blockSize, fpIn);
		// allocate memory for probability arrays
		uint16_t* buffProbs = malloc(sizeof(uint16_t)*AA_ALPHABET_SIZE*blockSize);
		ret = fread(buffProbs, sizeof(uint16_t), AA_ALPHABET_SIZE*blockSize, fpIn);
		// put data in the hash table
		//for (int ii=0; ii<blockSize; ii++) {
		for (int ii=0; ii<blockSize; ii++) {
			ak = buffAk[ii];
			famN = buffFamN[ii];
			if (ii==0) famN += MASK_OWNS_MEMORY;

			int absent;
			khint_t itr;
			itr = kc_ref_put(hRef, ak, &absent);
			ref_val* v = malloc(sizeof *v);
			v->familyNumber = famN;
			v->probs = buffProbs;
			kh_val(hRef, itr) = v;

			buffProbs += AA_ALPHABET_SIZE;
		}
		// read the size of the next data block
		ret = fread(&blockSize, sizeof(blockSize),1,fpIn);
	}
	// free buffers that are not used any more
	free(buffAk);
	free(buffFamN);
	fclose(fpIn);

	// copy file name
	if (fileNameRefOrigin!=NULL) free(fileNameRefOrigin);
    fileNameRefOrigin = strdup(fnIn);

	return 0;
}

int tst_compare_ref_to_file_txt(char* fnIn) {
	// error codes:
	// 0 - success
	// 1 - unable to open file
	// 4 - mismatch
	// 5 - preloaded reference contains kmers not present in the file
	int ii, kcur;
	uint64_t ak = 0; // current kmer
    uint64_t mask = (1ULL<<klen*BITS_PER_AA) - 1;

	gzFile fp;
	kseq_t *ks;
	if ((fp = gzopen(fnIn, "r")) == 0) return 1;
	ks = kseq_init(fp);
	char* fullDescription = malloc(sizeof(char)*1024);
	int famNumber;
	uint16_t* probArray = malloc(AA_ALPHABET_SIZE*sizeof(uint16_t));	
	while (kseq_read(ks) >= 0) {
		strcpy(fullDescription, ks->name.s);
		// in case the full description got truncated by a separator
		if (ks->comment.s!=NULL) strcat(fullDescription, ks->comment.s);
		_parse_ref_entry(fullDescription, &famNumber, probArray);

		for (ii = kcur = 0; ii < ks->seq.l; ii++) {
			int aa = aa_char2n[(uint8_t)ks->seq.s[ii]];
			if (aa <= AA_ALPHABET_SIZE) { // this is an amino acid or ?
				ak = (ak << BITS_PER_AA | aa) & mask;
				if (++kcur >= klen) { // we find a k-mer
					khint_t itr = kc_ref_get(hRef,ak);
					if ( !(itr<kh_end(hRef)) ) {
						printf("Kmer ");
						_print_aa_kmer(stdout,ak, klen);
						printf(" from the file is absent from the preloaded reference table.\n");
						return 4;
					}
					// allocate memory for the next reference value	from the file	
					ref_val* vF = malloc(sizeof *vF);   // value from the file
					vF->familyNumber = famNumber;
					vF->probs = probArray;	
				
					ref_val* vR = kh_val(hRef,itr); // value from the current ref table

					// check that all numbers match
					if ((vR->familyNumber & MASK_FAMILY_NUMBER) != vF->familyNumber) {
						printf("Mismatch in family number at kmer ");
						_print_aa_kmer(stdout, ak, klen);
						printf("\n");
						free(vF);
						return 4;
					}
					for (int ii=0; ii<klen; ii++) {
						for (int jj=0; jj<AA_ALPHABET_SIZE; jj++) {
							if (fabs(vR->probs[jj]-vF->probs[jj]) > 1e-6) {
								printf("Mismatch in probability array at kmer ");
								_print_aa_kmer(stdout, ak, klen);
								printf(": built=%hd, loaded=%hd", vR->probs[jj], vF->probs[jj]);
								printf("\n");
								free(vF);
								return 4;
							}
						} 
					}
					free(vF);
				}
			} else kcur = 0, ak = 0; // if there is a character outside of the aa alphabet, restart
		}

	}
	free(probArray);

	free(fullDescription);
	kseq_destroy(ks);
	gzclose(fp);

	return 0;
}

int tst_compare_ref_to_file_bin(char* fnIn) {
	// error codes:
	// 0 - success
	// 1 - unable to open file
	// 2 - error while parsing the file
	// 3 - kmer size in the file does not match the current size
	// 4 - mismatch
	// 5 - preloaded reference contains kmers not present in the file
	FILE* fpIn;
	if ((fpIn = fopen(fnIn, "rb")) == 0) return 1; // unable to open file
	int kmerSize;
	int ret = fread(&kmerSize, sizeof(kmerSize),1,fpIn); // read the kmer size
	if (ret!=1) return 2; // error parsing the file
	if (kmerSize!=klen) return 3; // unexpected kmer length

	int blockSize = binBlockSize;
	ret = fread(&blockSize, sizeof(blockSize),1,fpIn); // read the size of the next data block
	khint_t itr;
	uint64_t ak;
	uint64_t famN;
	uint64_t fileRefSize = 0;

	uint64_t* buffAk = malloc(blockSize*sizeof(ak));
	uint64_t* buffFamN = malloc(blockSize*sizeof(famN));
	while(ret!=0) {
		// read in data
		ret = fread(buffAk, sizeof(ak), blockSize, fpIn);
		ret = fread(buffFamN, sizeof(famN), blockSize, fpIn);
		// allocate memory
		uint16_t* buffProbs = malloc(sizeof(uint16_t)*AA_ALPHABET_SIZE*blockSize);
		uint16_t* buffProbsInit = buffProbs;
		ret = fread(buffProbs, sizeof(uint16_t), AA_ALPHABET_SIZE*blockSize, fpIn);
		// compare the data to the hash table
		for (int ii=0; ii<blockSize; ii++) {
			// check ak
			ak = buffAk[ii];
			itr = kc_ref_get(hRef,ak);
			if ( !(itr<kh_end(hRef)) ) {
				printf("Kmer ");
				_print_aa_kmer(stdout, ak, klen);
				printf(" from the file is absent from the preloaded reference table.\n");
				free(buffAk);
				free(buffFamN);
				free(buffProbs);
				return 4;
			}

			// allocate memory for the next reference value	from the file	
			ref_val* vF = malloc(sizeof *vF);   // value from the file
			vF->familyNumber = buffFamN[ii];
			vF->probs = buffProbs;	
			buffProbs += AA_ALPHABET_SIZE;
		
			ref_val* vR = kh_val(hRef,itr); // value from the current ref table

			// check that all numbers match
			if ((vR->familyNumber & MASK_FAMILY_NUMBER) != vF->familyNumber) {
				printf("Mismatch in family number at kmer ");
				_print_aa_kmer(stdout, ak, klen);
				printf("\n");
				free(vF);
				free(buffAk);
				free(buffFamN);
				free(buffProbsInit);
				return 4;
			}
			for (int jj=0; jj<AA_ALPHABET_SIZE; jj++) {
				if (vR->probs[jj] != vF->probs[jj]) {
					printf("Mismatch in probability array at kmer ");
					_print_aa_kmer(stdout, ak, klen);
					printf("\n");
					free(vF);
					free(buffAk);
					free(buffFamN);
					free(buffProbsInit);
					return 4;
				} 
			}
			free(vF);
		}
		free(buffProbsInit);
		fileRefSize += blockSize;
		// read the size of the next data block
		ret = fread(&blockSize, sizeof(blockSize),1,fpIn);
	}

	free(buffAk);
	free(buffFamN);
	fclose(fpIn);
	if (fileRefSize!=get_ref_size()) return 5;
	return 0;
}

/* for the given kmer from the reference table
copies probabilities into 'probs'. */
int tst_get_ref_probs_for_aa_kmer(char* kmer, double* probs) {
	if (hRef==NULL) return 1;
	if (strlen(kmer)!=klen) return 2;
	int aa;
	uint64_t ak = 0;	// bit-encoded amino acid kmer 
	uint64_t maskA = (1ULL<<klen*BITS_PER_AA) - 1;
	for (int ii=0; ii<strlen(kmer); ii++) {
		aa = aa_char2n[(uint8_t)kmer[ii]];
		if (aa < AA_ALPHABET_SIZE+1) {
			ak = (ak << BITS_PER_AA | aa) & maskA;
		} else {
			return 3;
		}
	}
	khint_t el = kc_ref_get(hRef,ak); // look up this kmer in the reference table
	if (el<kh_end(hRef)) { // if present
		uint16_t* probsInt = kh_val(hRef,el)->probs;
		for (int jj=0; jj<AA_ALPHABET_SIZE; jj++) {
			probs[jj] = exp(-(probsInt[jj]/65535.0)*9.99999);
		}
	} else {
		return 4;
	}
	return 0;		
}

/* reports an error returned by a thread in a batch processing mode */
void _report_thread_error(int retVal, char* fileName) {
	pthread_mutex_lock(&mutexReport);
	if (retVal==0) {
		// not an error
	} else if (retVal==1) {
		printf("Unable to process \"%s\": input file not found.\n", fileName);
	} else if (retVal==2) {
		printf("Unable to process %s: translation table has not been set.\n", fileName);
	} else if (retVal==3) {
		printf("Unable to process %s: reference has not been set.\n", fileName);
	} else if (retVal==11) {
		printf("Unable to process %s: file name is too long.\n", fileName);
	} else if (retVal==21) {
		printf("Unable to process a query: couldn't save to %s.\n", fileName);
	} else if (retVal==31) {
		printf("Unable to process a query: couldn't save to %s.\n", fileName);
	} else {
		printf("Unknown error (%d) while processing query %s.\n", retVal, fileName);
	}
	pthread_mutex_unlock(&mutexReport);
}

void* process_query(void* filename) {

	char* inputFileNamePrelim = (char*)filename;
	if (strlen(inputFileNamePrelim)>N_STRLEN-1-6) { _report_thread_error(11,inputFileNamePrelim); return 0; }

	char inputFileName[N_STRLEN];
	_trim_string(inputFileNamePrelim,inputFileName);

	char outputFileName[2*N_STRLEN+4];
	memcpy(outputFileName,outputDir,N_STRLEN);
	char baseFileName[N_STRLEN];
	_get_filename_base(inputFileName, baseFileName);
	if (baseFileName[0]=='\0') return 0;
	strcat(outputFileName,baseFileName);
	strcat(outputFileName,".json.tmp");

	char outputFileNameFin[2*N_STRLEN];
	memcpy(outputFileNameFin,outputDir,N_STRLEN);
	strcat(outputFileNameFin,baseFileName);
	strcat(outputFileNameFin,".json");

	char outputFileNameDetails[2*N_STRLEN];
	memcpy(outputFileNameDetails,outputDir,N_STRLEN);
	strcat(outputFileNameDetails,baseFileName);
	strcat(outputFileNameDetails,"_kmers.csv");

	kc_query_t *hQueryLocal;
	hQueryLocal = kc_query_init();

	double likelihoodsLocal[CODON_ALPHABET_SIZE*(AA_ALPHABET_SIZE+1)];
	memset(likelihoodsLocal,0,sizeof(double)*CODON_ALPHABET_SIZE*(AA_ALPHABET_SIZE+1));

	double decodingProbsLocal[CODON_ALPHABET_SIZE*(AA_ALPHABET_SIZE+1)];
	memset(decodingProbsLocal,0,sizeof(double)*CODON_ALPHABET_SIZE*(AA_ALPHABET_SIZE+1));

	uint64_t nConsensusLocal[CODON_ALPHABET_SIZE];
	memset(nConsensusLocal,0,sizeof(uint64_t)*CODON_ALPHABET_SIZE);

	double totsumLocal[AA_ALPHABET_SIZE];
	for (int jj=0; jj<AA_ALPHABET_SIZE; jj++) totsumLocal[jj] = log(0);

	uint64_t ret;
	ret = _hquery_put_seq_uni_orf(inputFileName, 0, 1, hQueryLocal, likelihoodsLocal, totsumLocal, nConsensusLocal);
	if (ret!=0) { _report_thread_error(ret,inputFileName); return 0; }
	ret = _hquery_put_seq_uni_orf(inputFileName, 1, 1, hQueryLocal, likelihoodsLocal, totsumLocal, nConsensusLocal);
	if (ret!=0) { _report_thread_error(ret,inputFileName); return 0; }
	ret = _hquery_put_seq_uni_orf(inputFileName, 2, 1, hQueryLocal, likelihoodsLocal, totsumLocal, nConsensusLocal);
	if (ret!=0) { _report_thread_error(ret,inputFileName); return 0; }
	ret = _hquery_put_seq_uni_orf(inputFileName, 0, 0, hQueryLocal, likelihoodsLocal, totsumLocal, nConsensusLocal);
	if (ret!=0) { _report_thread_error(ret,inputFileName); return 0; }
	ret = _hquery_put_seq_uni_orf(inputFileName, 1, 0, hQueryLocal, likelihoodsLocal, totsumLocal, nConsensusLocal);
	if (ret!=0) { _report_thread_error(ret,inputFileName); return 0; }
	ret = _hquery_put_seq_uni_orf(inputFileName, 2, 0, hQueryLocal, likelihoodsLocal, totsumLocal, nConsensusLocal);
	if (ret!=0) { _report_thread_error(ret,inputFileName); return 0; }

	_calculate_decoding_probs(likelihoodsLocal, totsumLocal, nConsensusLocal, decodingProbsLocal);

	printf("Writing a file to %s \n", outputFileNameFin);
	ret = _write_decoding_probs(outputFileName, decodingProbsLocal, nConsensusLocal);
	if (ret!=0) { _report_thread_error(20+ret,outputFileName); return 0; }
	if (reportQueryDetails) {
		ret = _write_query_kmers(hQueryLocal,outputFileNameDetails);
		if (ret!=0) { _report_thread_error(30+ret,outputFileNameDetails); return 0; }
	}

	ret = rename(outputFileName,outputFileNameFin);
	if (ret!=0) { _report_thread_error(20+ret,outputFileName); return 0; }

	_hquery_destroy(hQueryLocal);

	return NULL;
}

int batch_process_query(char* filenameBatch) {

	FILE* fBatch; 
	if ((fBatch = fopen(filenameBatch, "r")) == 0) return 1;
	if (translationTable==NULL) return 2;
	if (hRef==NULL) return 3;

	int tt;
	memset(thrArray,0,sizeof(pthread_t)*MAX_N_THREADS);

	char line[N_STRLEN];
	char fn[N_STRLEN];
	int haveNewLine;
	if(fgets(line, sizeof(line), fBatch)!=NULL) {
		haveNewLine=1;
		_trim_string(line,fn);
		if (fn[0]=='\0') haveNewLine = 0;
	} else {
		haveNewLine = 0;
	}
	while (haveNewLine) {
		// find a thread that's not busy
		for (tt=0; tt<nThreads; tt++) {
			if (thrArray[tt]==0) break;
			if (pthread_tryjoin_np(thrArray[tt],NULL)==0) break;
		}
		// assign and read a new line
		if (tt<nThreads) {
			strcpy(fnArray[tt],inputDir);
			strcat(fnArray[tt],fn);
			pthread_create(&thrArray[tt], NULL, process_query, (void*)fnArray[tt]); 
			if(fgets(line, sizeof(line), fBatch)!=NULL) {
				haveNewLine=1;
				_trim_string(line,fn);
				if (fn[0]=='\0') haveNewLine = 0;
			} else {
				haveNewLine = 0;
			}
		}
	}

	for (tt=0; tt<nThreads; tt++) {
		if (thrArray[tt]==0) continue;
		pthread_join(thrArray[tt],NULL);
	}

	fclose(fBatch);
	
	return 0;
}

int thread_process_query(char* filename) {
	int tt;
	// find a thread that's not busy
	for (tt=0; tt<MAX_N_THREADS; tt++) {
		if (thrArray[tt]==0) break;
		if (pthread_tryjoin_np(thrArray[tt],NULL)==0) break;
	}
	if (tt<MAX_N_THREADS) {
		strcpy(fnArray[tt],filename);
		pthread_create(&thrArray[tt], NULL, process_query, (void*)fnArray[tt]); 
		return tt;
	} else {
		return -1;
	}
}

int is_thread_alive(int thrNumber) {
	if (thrNumber<0 || thrNumber>=MAX_N_THREADS) {
		return -1;
	}
	if (thrArray[thrNumber]==0) return 0;
	if (pthread_tryjoin_np(thrArray[thrNumber],NULL)==0) {
		thrArray[thrNumber]=0;
		return 0;
	} else {
		return 1;
	}
}

void free_resources() {
	if (translationTable!=NULL) free(translationTable);
	translationTable = NULL;
	if (fileNameRefOrigin!=NULL) free(fileNameRefOrigin);
	fileNameRefOrigin = NULL;
	_hquery_destroy(hQueryGlobal);
    _hfamref_destroy();
	_href_destroy();
}

int main() {

	set_kmer_size(11);
	set_translation_table("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG");
	set_verbose(1);
	set_report_query_details(1);
	set_input_dir("../test/genomes");
	set_output_dir("../examples/results");
	set_max_query_count(2);
	set_n_threads(2);
	load_ref_bin("../examples/ref/ref.bin");

	/*double probs[AA_ALPHABET_SIZE];
	int ret = tst_get_ref_probs_for_aa_kmer("KVQAARD?FYM",probs);
	printf("ret = %u\n", ret);
	for (int ii=0; ii<AA_ALPHABET_SIZE; ii++) {
		printf("%f, ", probs[ii]);
	}
	printf("\n");
	return 0;*/

	//process_query("../test/genomes/GCA_000146045.2_R64_genomic.fna.gz");
	batch_process_query("../test/genomes/genome_list.txt");

	batch_process_query("../test/genomes/genome_list.txt");

	free_resources();
	
	return 0;

}
