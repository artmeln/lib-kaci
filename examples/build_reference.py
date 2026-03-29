import os
from datetime import datetime
import sys
sys.path.append("..")
from  src.libkaci_wrapper import LibKaci

localPath = os.path.dirname(os.path.abspath(__file__)) + '/../'
libPath = localPath+'src/libkaci.so'
lib = LibKaci(libPath)

klen = 11
maxRefCount = round(klen*0.3)
minFamrefCountForKmer = 2
minProteinsPerFamily = 1
minLinkNumber = 10

dirName = f'{localPath}/examples/ref'
os.system(f'mkdir {dirName}') 

lib.set_kmer_size(klen)
lib.set_verbose(True)
lib.set_min_proteins_per_family(minProteinsPerFamily)
lib.set_max_familyref_count_for_one_aa_per_kmer(maxRefCount)
lib.set_min_familyref_count_for_kmer(minFamrefCountForKmer)
lib.set_min_link_number(minLinkNumber)
lib.set_output_dir(f'{localPath}/examples/ref')

# build the reference
lib.build_ref(f'{localPath}/test/Sugar_tr.faa.gz', True)

lib.save_parameters_ref(f'{dirName}/params.json')
lib.save_ref_bin(f'{dirName}/ref.bin')

lib.free_resources()