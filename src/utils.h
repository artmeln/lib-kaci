#ifndef _UTILS_H_
#define _UTILS_H_

#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "constants.h"
#include "operations128.h"


/* returns reverse complement of nucleotide n */
uint8_t inline _rc(uint8_t n) {
	if (n!=255) {
		n += 2;
		if (n>3) n-=4;
	}
	return n;
}

/* extracts familyName from description string.
description example: "proteinName;familyName;anything else" */
void _parse_family_name(char* description, char* familyName) {
	char* word;
	char* str = strdup(description);
	char* p = str;
	strsep(&str, ";"); // skip the first one
	word = strsep(&str, ";");
	strcpy(familyName, word);
	free(p);
}

/* extracts the file name (also eliminating .fna and .gz extensions)
from the full file path */
void _get_filename_base(char* input, char* output) {
	// get file name by splitting on '/'
	char* word;
	char* str = strdup(input);
	char* p = str;
	word = strsep(&str, "/");
	while (str!=NULL) {
		word = strsep(&str, "/");
	}
	// trim .gz and .fna extensions, if present
	int len;
	len = strlen(word);
	if (len>3) {
		if (word[len-3]=='.' && word[len-2]=='g' && word[len-1]=='z') word[len-3]='\0';
	}
	len = strlen(word);
	if (len>4) {
		if (word[len-4]=='.' && word[len-3]=='f' && word[len-2]=='n' && word[len-1]=='a') word[len-4]='\0';
	}
	// copy the output
	strcpy(output, word);
	free(p);
}

/* parses the description line of reference kmers stored as text.
The description consists of the family id followed by 20 probabilities
for the uncertain amino acid in the kmer.
An example description: "2;0.05,0.05,0.05,...,0.05" */
void _parse_ref_entry(char* description, int* famNumber, uint16_t* probs) {
	char* word;
	char* str = strdup(description);
	char* p = str;
	word = strsep(&str, ";"); // family name
	sscanf(word, "%d", famNumber);
	// read in probability array
	for (int jj=0; jj<AA_ALPHABET_SIZE; jj++) {
		word = strsep(&str, ",");
		sscanf(word, "%hd", probs);
		probs++;
	}
	free(p);
}

/* trims spaces, eol, etc from both ends of the string */
void _trim_string(char* input, char* output) {
	int end=strlen(input)-1;
	int ii;
	for (ii=0; ii<=end; ii++) {
		if (!isspace(input[ii])) break;
	}
	if (ii>end) {output[0]='\0'; return;}
	while (isspace(input[0])) input++;
	end=strlen(input)-1;
	while (isspace(input[end])) end--;
	memcpy(output, input, end+1);
	output[end+1] = '\0';
}

/* prints int64 encoded amino aicd kmer to the terminal as text */
void _print_aa_kmer(FILE* fOut, uint64_t ak, int kmer_length) {
	uint64_t mask = (1ULL<<kmer_length*BITS_PER_AA) - 1;
	for (int ii=0; ii<kmer_length; ii++) {
		fprintf(fOut, "%c", aa_table[(ak & mask) >> (kmer_length-1)*BITS_PER_AA]);
		ak = ak << BITS_PER_AA;
	}
}

/* prints int128 encoded codon kmer to the terminal as text */
void _print_cdn_kmer(FILE* fOut, uint128 ck, int kmer_length) {
	uint128 maskC;
	if (kmer_length>64/BITS_PER_CODON) {
		maskC.low = UINT64_MAX;
		maskC.high = (1ULL<<(kmer_length*BITS_PER_CODON-64)) - 1;
	} else {
		maskC.low = (1ULL<<(kmer_length*BITS_PER_CODON)) - 1;
		maskC.high = 0ULL;
	}
	for (int ii=0; ii<kmer_length; ii++) {
		fprintf(fOut, "%s", codon_table[_shift_right_128(_and_128(ck,maskC),(kmer_length-1)*BITS_PER_CODON).low]);
		ck = _shift_left_128(ck,BITS_PER_CODON);
	}
}

/* counts how many times each amino acid is present in the kmer. 
A kmer is divese if no single amino acid has a count that exceeds maxCountForOneAA */
int _is_aa_kmer_diverse(uint64_t ak, int maxCountForOneAA, int kmer_length) {
	// akEmpty is filled with '*'
	uint64_t akEmpty = 0ULL;
	for (int ss=0; ss<kmer_length; ss++) akEmpty = (akEmpty << BITS_PER_AA) + 20ULL;
	uint64_t maskLastAA = (1ULL<<BITS_PER_AA) - 1;
	uint64_t akTemp = ak;
	int aaCount = 0; // number of aa in this kmer
	for (int jj=0; jj<kmer_length; jj++) { // go over all aa in the kmer
		akTemp = ak >> BITS_PER_AA*jj;
		if ((akTemp & maskLastAA)==20) continue;
		uint64_t aaToCmpare = (akTemp & maskLastAA) << BITS_PER_AA*jj;
		int ii;
		for (ii=jj; ii<kmer_length; ii++) {
			uint64_t maskAii = maskLastAA << BITS_PER_AA*ii;
			if ((ak & maskAii)==aaToCmpare) {
				aaCount++;
				ak = (ak & ~maskAii) + (20ULL << BITS_PER_AA*ii);
			}
			aaToCmpare = aaToCmpare << BITS_PER_AA;
			if (ak==akEmpty) break;
		}
		if (ak==akEmpty) break;
		if (aaCount>maxCountForOneAA) return 0;
		aaCount = 0;
	}
	if (aaCount>maxCountForOneAA) return 0; else return 1;
}

// a different version that require a min number of aa in the kmer
/*int _is_aa_kmer_diverse(uint64_t ak, int minNumAA) {
	// akEmpty is filled with '*'
	uint64_t akEmpty = 0ULL;
	for (int ss=0; ss<klen; ss++) akEmpty = (akEmpty << BITS_PER_AA) + 20ULL;
	uint64_t maskLastAA = (1ULL<<BITS_PER_AA) - 1;
	uint64_t akTemp = ak;
	int aaCount = 0; // number of aa in this kmer
	for (int jj=0; jj<klen; jj++) { // go over all aa in the kmer
		akTemp = ak >> BITS_PER_AA*jj;
		if ((akTemp & maskLastAA)!=20) {
			aaCount++;
			if (aaCount>=minNumAA) return 1;
		} else {
			continue;
		}
		uint64_t aaToCmpare = (akTemp & maskLastAA) << BITS_PER_AA*jj;
		int ii;
		for (ii=jj; ii<klen; ii++) {
			uint64_t maskAii = maskLastAA << BITS_PER_AA*ii;
			if ((ak & maskAii)==aaToCmpare) ak = (ak & ~maskAii) + (20ULL << BITS_PER_AA*ii);
			aaToCmpare = aaToCmpare << BITS_PER_AA;
			if (ak==akEmpty) break;
		}
	}
	if (aaCount<minNumAA) return 0; else return 1;
}*/

/* mutates the amino acid at position 'pos' to 'aa'
('aa' is provided as a character from aa_table).
Both input and output kmers are int64 encoded. */
uint64_t _mutate_aa_kmer(uint64_t ak, int pos, char aa, int kmer_length) {
	uint64_t mask = 31ULL << (kmer_length-pos)*BITS_PER_AA;
	uint64_t ak_mut;
	ak_mut = (ak & ~mask) | ((uint64_t)aa_char2n[(uint8_t)aa] << (kmer_length-pos)*BITS_PER_AA);
	return ak_mut;
}

/* mutates the codon at position 'pos' to 'codon'.
('codon' is provided as its int128 encoding).
Both input and output kmers are int128 encoded. */
uint128 _mutate_cdn_kmer(uint128 ck, int pos, uint128 codon, int kmer_length) {
	uint128 mask = {127ULL,0}; // 2^BITS_PER_CODON - 1 = 127
	mask = _shift_left_128(mask, (kmer_length-pos)*BITS_PER_CODON);
	uint128 ck_mut;
	ck_mut = _and_128(ck,_not_128(mask));
	ck_mut = _or_128(ck_mut, _shift_left_128(codon,(kmer_length-pos)*BITS_PER_CODON));
	return ck_mut;
}

#endif
