#include <stdint.h>

/* sets the size of the kmer for building 
reference and query hash tables 
3 <= klength <= 12
returns 0 if successful, 1 otherwise */
int set_kmer_size(int klength);

/* sets the limit on the count of aa kmer
in the query hash table
returns 0 if successful, 1 otherwise */
int set_max_query_count(int v);

/* set the flag that controls detailed output of query details */
void set_report_query_details(int v);

/* set the minimum number of different amino acids that have to be in the kmer
to be entered into the family hash table.
returns 0 is successful, 1 otherwise */
//int set_ref_min_aa_per_kmer(int aaPerKmer);

/* set the maximum number of times same amino acid can be present in the kmer
to be entered into the family hash table.
returns 0 is successful, 1 otherwise */
int set_max_familyref_count_for_one_aa_per_kmer(int v);

/* set the minimum count for kmers in one family to be kept.
For example, 2 keeps kmers that are present 
in at least two proteins in the family */
int set_min_familyref_count_for_kmer(int v);

/* set the minimum number of proteins 
that must be present in the family 
to be entered into the reference */
int set_min_proteins_per_family(int value);

/* when entering family reference kmers 
into the global reference hash table, only those
that are "linked" in a group of minLinkNumber will be kept
Linked kmers are those that are derived by a single  aa substitution. */
int set_min_link_number(int value);

/* set the flag that controls optional comments
sent to the terminal */
void set_verbose(int value);

/* set the size of processing block for bin ref files */
int set_bin_block_size(int value);

/* sets the number of threads that will be used 
to process queries in the batch mode. */
int set_n_threads(int value);

/* sets the input directory for 
 - query filesprocessed in the batch mode 
 */
int set_input_dir(char* dirPath);

/* sets the output directory for 
 - files created in the batch query processing mode 
 */
int set_output_dir(char* dirPath);

/* sets the translation table used for processing the query.
translation_table is a string of 64 capital amino acids + '*' for termination codons
The standard translation table for the codon_table from kaci.c is
FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
 */
void set_translation_table(const char* translation_table);

/* builds the reference hash table from a set of fasta sequences 
found in a file specified by filename. The sequences must have descriptions formatted as
>protein_name;family_name;protein_description
All sequences with the same family_name will be processed together
building a family hash table of kmers which is then used to calculate
probabilities stored in the reference hash table.
Overwrites any previously loaded references and resets the query.
If saveFamilyList==true writes a file with 
familyName, familyNumber, nContributedKmers to outputDir.
Behavior of this function is affected by 
set_kmer_size, set_max_ref_count_for_one_aa_per_kmer, set_min_familyref_count_for_kmer
returns 0 is successful, 1 if the file was not found */
int build_ref(char* filename, int saveFamilyList);

/* save current values of all parameters 
that affect refence construction */
int save_parameters_ref(char* fnOut);

/* saves the reference hash table built with build_ref as text
so that it can be later loaded with load_ref_txt */
int save_ref_txt(char* filename);

/* saves the reference hash table built with build_ref
as a custom format binary file 
so that it can be later loaded with load_ref_bin */
int save_ref_bin(char* filename);

/* loads a reference hash table previously saved as a text file 
using save_ref_txt() */
int load_ref_txt(char* filename);

/* loads a reference hash table previously saved as a binary file 
using save_ref_bin() */
int load_ref_bin(char* filename);

/* compares the reference built with build_ref
to a reference stored as a text file at filename.
Returns 0 if successful or a code in case of an error:
0 - success, 1 - unable to open file, 4 - mismatch,
5 - preloaded reference contains kmers not present in the file */
int tst_compare_ref_to_file_txt(char* filename);

/* for the given kmer from the reference table
copies probabilities into 'probs'.
Returns 
0 in case of success, 
1 if reference is not available,
2 if the length of kmer does not match klen
3 if the provided kmer contains characters that are not in the alphabet
4 if kmer is absent from the refrence  */
int tst_get_ref_probs_for_aa_kmer(char* kmer, double* probs);

/* compares the reference built with build_ref
to a reference stored as a binary file at filename.
0 - success, 1 - unable to open file, 2 - error while parsing the file,
3 - kmer size in the file does not match the current size, 4 - mismatch,
5 - preloaded reference contains kmers not present in the file */
int tst_compare_ref_to_file_bin(char* filename);

/* resets the query hash table */
void tst_reset_query();

/* processes a query file by 'process_query' on a new thread. 
The total number of running threads may not exceed MAX_N_THREADS.
Returns the number of the thread if successful (this number can be passed to 
'is_thread_alive'), -1 if the number of available threads has been exceeded. */
int thread_process_query(char* filename);

/* checks and returns the status of the thread 
started by 'thread_process_query'. */
int is_thread_alive(int thrNumber);

/* processes query files with names provided as separate lines
in a text file filenameBatch. Each query is processed on a separate thread 
by 'process_query'. The number of threads is set by set_n_threads.
Returns 0 if successful, 1 if filenameBatch was not found, 
2 if translation table has not been set,
3 if the reference has not been set. */
int batch_process_query(char* filenameBatch);

/* builds a query hash table from fasta nucleotide sequences contained in filename
by translating all possible 6 orfs. Every query kmer that is found
in the reference hash table is entered in the query table, and likelihoods 
for the codons encoding the kmer are updated. The bahavior of this function is affected by
set_translation_table and set_max_query_count.
Always returns 0 but reports errors to the terminal. 
This function is also called from the batch porcessing mode. */
void* process_query(void* filename);

/* calls _hquery_put_seq_uni_orf on the forward orf (no shift) */
int tst_process_query_orf1f(char* filename);

/* calls _hquery_put_seq_uni_orf on the reverse orf (no shift) */
int tst_process_query_orf1r(char* filename);

/* calls _hquery_put_seq_uni_orf on the forward orf (1nt shift) */
int tst_process_query_orf2f(char* filename);

/* calls _hquery_put_seq_uni_orf on the reverse orf (1nt shift) */
int tst_process_query_orf2r(char* filename);

/* calls _hquery_put_seq_uni_orf on the forward orf (2nt shift) */
int tst_process_query_orf3f(char* filename);

/* calls _hquery_put_seq_uni_orf on the reverse orf (2nt shift) */
int tst_process_query_orf3r(char* filename);

/* calculates decoding probabilities 
for the current global instances of likelihoods, totsum, and nConsensus */
void tst_calculate_decoding_probs();

/* writes the current global decoding probabilities to a json file */
int tst_write_decoding_probs(char* fnOut);

/* prints all reference kmers and the associated values to the terminal.
showExp_P flag indicates if the probability should be displayed as 
p (showExp_P=1) or -log(p) (showExp_P=0) */
void tst_print_ref_kmers(int showExp_P);

/* prints all query kmers and the associated values to the terminal.
showExp_P flag indicates if the probability should be displayed as 
p (showExp_P=1) or -log(p) (showExp_P=0) */
void tst_print_query_kmers();

/* copies the number of times each codon was entered into calculation
during query processing into 'ns'. */
void tst_get_n_consensus(uint64_t* ns);

/* copies decoding probabilities into 'probs'. */
void tst_get_decoding_probs(double* probs);

/* returns the number of kmers in the reference hash table */
unsigned long get_ref_size();

/* deletes all hash tables, releases all resources */
void free_resources();
