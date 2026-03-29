import os
from datetime import datetime
import sys
sys.path.append("..")
from src.libkaci_wrapper import LibKaci

standardTT = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"

localPath = os.path.dirname(os.path.abspath(__file__)) + '/../'
libPath = localPath+'src/libkaci.so'
lib = LibKaci(libPath)

klen = 11
maxQueryCount = 2

lib.set_kmer_size(klen)
lib.set_translation_table(standardTT)
lib.set_verbose(True)
lib.set_report_query_details(True)
lib.set_input_dir(f'{localPath}/test/genomes')
lib.set_output_dir(f'{localPath}/examples/results')
lib.set_max_query_count(maxQueryCount)
lib.set_n_threads(6)
lib.load_ref_bin(f'{localPath}/examples/ref/ref.bin')

#lib.process_query(f'{localPath}/test/genomes/GCA_000146045.2_R64_genomic.fna.gz')
lib.batch_process_query(f'{localPath}/test/genomes/genome_list.txt')

lib.free_resources()