import time
import os
import json
from subprocess import Popen, PIPE
from src.libkaci_wrapper import LibKaci

localPath = os.path.dirname(os.path.abspath(__file__)) + '/'
libPath = localPath+'src/libkaci.so'

tt = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
codons = ['TTT', 'TTC', 'TTA', 'TTG', 'TCT', 'TCC', 'TCA', 'TCG', 
          'TAT', 'TAC', 'TAA', 'TAG', 'TGT', 'TGC', 'TGA', 'TGG',
          'CTT', 'CTC', 'CTA', 'CTG', 'CCT', 'CCC', 'CCA', 'CCG', 
          'CAT', 'CAC', 'CAA', 'CAG', 'CGT', 'CGC', 'CGA', 'CGG', 
          'ATT', 'ATC', 'ATA', 'ATG', 'ACT', 'ACC', 'ACA', 'ACG', 
          'AAT', 'AAC', 'AAA', 'AAG', 'AGT', 'AGC', 'AGA', 'AGG', 
          'GTT', 'GTC', 'GTA', 'GTG', 'GCT', 'GCC', 'GCA', 'GCG', 
          'GAT', 'GAC', 'GAA', 'GAG', 'GGT', 'GGC', 'GGA', 'GGG']
aas = 'ACDEFGHIKLMNPQRSTVWY'

def compare_json_files(fnA, fnB):    
    fileA = open(fnA,'r')
    fileB = open(fnB,'r')
    dictA = json.load(fileA)
    dictB = json.load(fileB)
    for keyA in dictA:
        if keyA=="timestamp": continue
        if keyA=="version": continue
        valA = dictA[keyA]
        valB = dictB.get(keyA)
        if keyA=="referenceFile":
             valA = valA.split('/')[len(valA.split('/'))-1]
             valB = valB.split('/')[len(valB.split('/'))-1]
        if valA!=valB:
            print(f"At {keyA} the following two values are not equal:")
            print(valA)
            print(valB)
            return False
    for keyB in dictB:
        if keyB=="timestamp": continue
        if keyB=="version": continue
        valB = dictB[keyB]
        valA = dictA.get(keyB)
        if keyB=="referenceFile":
             valA = valA.split('/')[len(valA.split('/'))-1]
             valB = valB.split('/')[len(valB.split('/'))-1]
        if valB!=valA:
            print(f"At {keyB} the following two values are not equal:")
            print(valA)
            print(valB)
            return False
    fileA.close()
    fileB.close()
    return True
      

lib = LibKaci(libPath)

lib.set_verbose(True)

# test1: check the reference built from a simple example file.
# Includes tests for loading and comparing a stored reference file.

lib.set_kmer_size(4)

lib.build_ref(localPath + 'test/ref1.faa',False)
#lib.save_ref_txt('test/ref1_hashtable.faa.gz')
#os.system(f'gzip -d test/ref1_hashtable.faa.gz')
#lib.save_ref_bin('test/ref1_hashtable.faa.gz')

lib.tst_compare_ref_to_file_txt('test/ref1_hashtable.faa')

#lib.print_ref_kmers(True) 
# prints out 
#  - kmer
#  - family number
#  - probabilities

lib.load_ref_txt('test/ref1_hashtable.faa')
lib.tst_compare_ref_to_file_txt('test/ref1_hashtable.faa')

lib.set_bin_block_size(4)
lib.tst_compare_ref_to_file_bin('test/ref1_hashtable.bin')

lib.load_ref_bin('test/ref1_hashtable.bin')
lib.tst_compare_ref_to_file_txt('test/ref1_hashtable.faa')
lib.set_bin_block_size(4096)

print('Successfully completed test1.')


# test2: check that during reference construction
# - kmers with single amino acid count greater than aaPerKmer are omitted
# - kmers present in more than one family are given family number 0 and removed
# - only kmers over minFamrefKmerCount are kept
# - families with fewer than minProteinsPerFamily are skipped
aaPerKmer = 2
minFamrefKmerCount = 2
minProteinsPerFamily = 2
lib.set_kmer_size(4)
lib.set_max_familyref_count_for_one_aa_per_kmer(aaPerKmer)
lib.set_min_familyref_count_for_kmer(minFamrefKmerCount)
lib.set_min_proteins_per_family(minProteinsPerFamily)
lib.set_output_dir(localPath + 'test/')
lib.build_ref(localPath + 'test/ref2.faa',True)

#lib.save_ref_txt('test/ref2_hashtable.faa.gz')
#os.system(f'gzip -d test/ref2_hashtable.faa.gz')
#lib.print_ref_kmers(True) 
# prints out 
# kmer
# family number
# probabilities

#lib.save_parameters_ref("test/ref1_params.json")

lib.tst_compare_ref_to_file_txt('test/ref2_hashtable.faa')
print('Successfully completed test2.')


# test3 (may need a better test here)
# a test for equal query decoding probabilities for a simple collection of kmers
lib.set_kmer_size(7)
lib.set_verbose(False)
lib.set_report_query_details(True)
lib.set_max_familyref_count_for_one_aa_per_kmer(3)
lib.set_min_familyref_count_for_kmer(1) # count everything
lib.set_min_proteins_per_family(1)      # count everything
lib.set_max_query_count(1)
lib.build_ref(localPath + 'test/ref3.faa',False)
lib.set_translation_table(tt)
lib.tst_reset_query()
lib.tst_process_query_orf1f('test/query3.fna')
lib.tst_calculate_decoding_probs()

result = lib.tst_get_decoding_probs()
# TGC is present in seq1 only but based on the ref
# the decoding is uncertain, equally split between C, D, S, T
expectedProb_TGC_CDST = 0.073294  
measuredProb = result[codons.index('TGC')][aas.index('C')]
if abs(expectedProb_TGC_CDST-measuredProb)>1e-5:
    raise Exception(f'Test3 failed: expected probability {expectedProb_TGC_CDST} for TGC encoded by C, measured {measuredProb}')
measuredProb = result[codons.index('TGC')][aas.index('D')]
if abs(expectedProb_TGC_CDST-measuredProb)>1e-5:
    raise Exception(f'Test3 failed: expected probability {expectedProb_TGC_CDST} for TGC encoded by D, measured {measuredProb}')
measuredProb = result[codons.index('TGC')][aas.index('S')]
if abs(expectedProb_TGC_CDST-measuredProb)>1e-5:
    raise Exception(f'Test3 failed: expected probability {expectedProb_TGC_CDST} for TGC encoded by S, measured {measuredProb}')
measuredProb = result[codons.index('TGC')][aas.index('T')]
if abs(expectedProb_TGC_CDST-measuredProb)>1e-5:
    raise Exception(f'Test3 failed: expected probability {expectedProb_TGC_CDST} for TGC encoded by T, measured {measuredProb}')

# check that all probabilities are properly normalized
for cdn in range(0,len(codons)):
    s = sum(result[cdn])
    if abs(sum(result[cdn])-1.0)>1e-12:
         raise Exception(f'Test3 failed: decoding probabilities for {codons[cdn]} add up to {s}, not to 1.0,')
print('Successfully completed test3.')


# test4 - forward and reverse calculations should give identical results
lib.set_kmer_size(10)
lib.set_verbose(False)
lib.set_max_familyref_count_for_one_aa_per_kmer(3)
lib.set_min_familyref_count_for_kmer(1) # count everything
lib.set_min_proteins_per_family(1)      # count everything
lib.set_max_query_count(10)             # count everything
lib.build_ref(localPath + 'test/ref4.faa',False)

lib.set_translation_table(tt)

lib.tst_reset_query()
lib.tst_process_query_orf1f('test/query4.fna')
lib.tst_calculate_decoding_probs()
resultFwd = lib.tst_get_decoding_probs()

lib.tst_reset_query()
lib.tst_process_query_orf1r('test/query4_rc.fna')
lib.tst_calculate_decoding_probs()
resultRev = lib.tst_get_decoding_probs()

# compare results
for cc in range(len(resultFwd)):
     for aa in range(len(resultFwd[cc])):
          if abs(resultFwd[cc][aa]-resultRev[cc][aa])/resultFwd[cc][aa]>1e-6:
          #if resultFwd[cc][aa] != resultRev[cc][aa]:
               codon = codons[cc]
               amino_acid = aas[aa]
               raise Exception(f'Test4 failed: mismatch at codon={codon}, amino acid={amino_acid}; decoding probs fwd={resultFwd[cc][aa]}, rev={resultRev[cc][aa]}')
print('Successfully completed test4.')


# test5: check that
# - kmers with link number below minLinkNumber are omitted
lib.set_kmer_size(4)
lib.set_max_familyref_count_for_one_aa_per_kmer(4)
lib.set_min_familyref_count_for_kmer(1)
lib.set_min_proteins_per_family(1)

# everything should be counted
minLinkNumber = 2
lib.set_min_link_number(minLinkNumber)
lib.build_ref(localPath + 'test/ref5.faa',False)
res = lib.get_ref_size()
expected = 4
#lib.print_ref_kmers(True)
if res!=expected:
        raise Exception(f'Test5 failed: expected {expected} but found {res} kmers in the reference table for minLinkNumber = {minLinkNumber}.')
# the cutoff value at which kmers are not going to be counted
minLinkNumber = 3
lib.set_min_link_number(minLinkNumber)
lib.build_ref(localPath + 'test/ref5.faa',False)
res = lib.get_ref_size()
expected = 0
if res!=expected:
        raise Exception(f'Test5 failed: expected {expected} but found {res} kmers in the reference table for minLinkNumber = {minLinkNumber}.')
print('Successfully completed test5.')


# test6: check that process_query() is equivalent to sequentially
# calling tst_process_query_orf1f, ... , tst_process_query_orf3r
# and that the result agrees with the stored expectation
lib.set_kmer_size(11)
lib.set_max_familyref_count_for_one_aa_per_kmer(3)  
lib.set_min_familyref_count_for_kmer(2)
lib.set_min_proteins_per_family(1)
lib.set_min_link_number(10)
lib.set_max_query_count(2)
lib.set_output_dir(localPath + 'test/')
lib.build_ref(localPath + 'test/Sugar_tr.faa.gz',False)

# sequential calculation
lib.tst_process_query_orf1f(localPath + 'test/genomes/GCF_000146045.2_R64_genomic.fna.gz')
lib.tst_process_query_orf2f(localPath + 'test/genomes/GCF_000146045.2_R64_genomic.fna.gz')
lib.tst_process_query_orf3f(localPath + 'test/genomes/GCF_000146045.2_R64_genomic.fna.gz')
lib.tst_process_query_orf1r(localPath + 'test/genomes/GCF_000146045.2_R64_genomic.fna.gz')
lib.tst_process_query_orf2r(localPath + 'test/genomes/GCF_000146045.2_R64_genomic.fna.gz')
lib.tst_process_query_orf3r(localPath + 'test/genomes/GCF_000146045.2_R64_genomic.fna.gz')
lib.tst_calculate_decoding_probs()
lib.tst_write_decoding_probs('test/GCF_000146045.2_R64_genomic.json')
resultOfComparison = compare_json_files('test/GCF_000146045.2_R64_genomic.json','test/GCF_000146045.2_R64_genomic_expected.json')
if not resultOfComparison:
    raise Exception(f'Test6 failed: results of tst_process_query_orf1f/2f/3f/1r/2r/3r do not agree with the expectation.')

# single function call calculation
lib.process_query(localPath + 'test/genomes/GCF_000146045.2_R64_genomic.fna.gz')
resultOfComparison = compare_json_files('test/GCF_000146045.2_R64_genomic.json','test/GCF_000146045.2_R64_genomic_expected.json')
if not resultOfComparison:
    raise Exception(f'Test6 failed: results of process_query do not agree with the expectation.')

# batch calculation
lib.set_input_dir(localPath + 'test/genomes')
lib.set_n_threads(3)
lib.batch_process_query(localPath + 'test/genome_list.txt')
resultOfComparison = compare_json_files('test/GCF_000146045.2_R64_genomic.json','test/GCF_000146045.2_R64_genomic_expected.json')
if not resultOfComparison:
    raise Exception(f'Test6 failed: results of batch_process_query do not agree with the expectation for GCF_000146045.2_R64_genomic.')
os.system('rm test/GCF_000146045.2_R64_genomic.json')
resultOfComparison = compare_json_files('test/GCF_000182965.3_ASM18296v3_genomic.json','test/GCF_000182965.3_ASM18296v3_genomic_expected.json')
if not resultOfComparison:
    raise Exception(f'Test6 failed: results of batch_process_query do not agree with the expectation for GCF_000182965.3_ASM18296v3_genomic.')
os.system('rm test/GCF_000182965.3_ASM18296v3_genomic.json')
resultOfComparison = compare_json_files('test/GCA_001661245.1_Pacta1_2_genomic.json','test/GCA_001661245.1_Pacta1_2_genomic_expected.json')
if not resultOfComparison:
    raise Exception(f'Test6 failed: results of batch_process_query do not agree with the expectation for GCA_001661245.1_Pacta1_2_genomic.')
os.system('rm test/GCA_001661245.1_Pacta1_2_genomic.json')

# new thread calculation
if os.path.isfile('test/GCF_000146045.2_R64_genomic.json'):
    os.system(f'rm test/GCF_000146045.2_R64_genomic.json')
thrN = lib.thread_process_query(localPath + 'test/genomes/GCF_000146045.2_R64_genomic.fna.gz')
while lib.is_thread_alive(thrN):
    print("Waiting...")
    time.sleep(1)
resultOfComparison = compare_json_files('test/GCF_000146045.2_R64_genomic.json','test/GCF_000146045.2_R64_genomic_expected.json')
if not resultOfComparison:
    raise Exception(f'Test6 failed: results of process_query do not agree with the expectation.')
os.system('rm test/GCF_000146045.2_R64_genomic.json')
print('Successfully completed test6.')


# test7: check that details saved along with the kmers
# correctly identify kmer position in the sequence
genomeName = 'test/genomes/GCF_000146045.2_R64_genomic.fna.gz'
resultName = 'test/GCF_000146045.2_R64_genomic_kmers.csv'
fR64 = open(resultName)
line = fR64.readline().strip() # header
line = fR64.readline().strip()
seqDirection = int(line.split(',')[5])
# find a forward sequence
while seqDirection==0:
    line = fR64.readline().strip()
    seqDirection = int(line.split(',')[5])
seqNumber = int(line.split(',')[2])
seqPosition = int(line.split(',')[4])
ntKmer = line.split(',')[-1]
# extract the sequence using seqkit
cmd = f"seqkit head -n {seqNumber} {genomeName} | " \
    + f"seqkit range -r -1:-1 | " \
    + f"seqkit subseq -r {seqPosition}:{seqPosition+len(ntKmer)-1} | " \
    + "seqkit seq -s"
proc = Popen(cmd, shell=True,stdout=PIPE,stderr=PIPE)
stdout, stderr = proc.communicate()
seqkitKmer = stdout.decode('utf-8').strip()
if ntKmer!=seqkitKmer:
    print(cmd)
    raise Exception(f'Test7 failed: expected kmer {ntKmer} but seqkit extracted {seqkitKmer} with the above command.')

line = fR64.readline().strip()
seqDirection = int(line.split(',')[5])
# find a reverse sequence
while seqDirection!=0:
    line = fR64.readline().strip()
    seqDirection = int(line.split(',')[5])
seqNumber = int(line.split(',')[2])
seqPosition = int(line.split(',')[4])
ntKmer = line.split(',')[-1]
# extract the sequence using seqkit
cmd = f"seqkit head -n {seqNumber} {genomeName} | " \
    + f"seqkit range -r -1:-1 | seqkit subseq -r {seqPosition}:{seqPosition+len(ntKmer)-1} | " \
    + "seqkit seq -r -p -s "
proc = Popen(cmd, shell=True,stdout=PIPE,stderr=PIPE)
stdout, stderr = proc.communicate()
seqkitKmer = stdout.decode('utf-8').strip()
if ntKmer!=seqkitKmer:
    print(cmd)
    raise Exception(f'Test7 failed: expected kmer {ntKmer} but seqkit extracted {seqkitKmer} with the above command.')
print('Successfully completed test7.')


lib.free_resources()
