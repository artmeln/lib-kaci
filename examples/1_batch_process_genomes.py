import os
from datetime import datetime
import sys
sys.path.append("../../lib-kaci/")
from  src.libcodonprob_wrapper import LibCodonProb

#standardTT = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
standardTT = "FFLLSSSSYYAGCCSWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"

localPath = os.path.dirname(os.path.abspath(__file__)) + '/'
libPath = localPath+'../lib-kaci/src/libkaci.so'
lib = LibCodonProb(libPath)

nThreads = 7
klen = 11
maxQueryCount = 2
minLinkNumber = 20

dirGenomes = "genomes/"
dirResults = "results/"
referenceName = f"ref/ref_k11_link20/ref.bin"

# make genome list
f = open("genome_list.txt",'w')
fileList = os.listdir(dirGenomes)
for fileName in fileList:
    f.write(fileName+'\n')
f.close()


lib.set_input_dir(dirGenomes)
lib.set_n_threads(nThreads)
lib.set_translation_table(standardTT)

lib.set_kmer_size(klen)
lib.set_verbose(True)
lib.set_max_query_count(maxQueryCount)
lib.set_output_dir(dirResults)
lib.set_report_query_details(False)

# load the reference
lib.load_ref_bin(referenceName)

lib.batch_process_query(f"genome_list.txt")

lib.free_resources()
