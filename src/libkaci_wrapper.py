import ctypes as ct
import numpy as np

class LibKaci():

    lib = None

    def __init__(self,libPath):
        self.lib = ct.CDLL(libPath)

    def set_kmer_size(self, klength: int):
        ret = self.lib.set_kmer_size(klength)
        if ret!=0:
            raise Exception(f'kmer size {klength} is not allowed.')

    def set_max_query_count(self, maxQueryCount: int):
        ret = self.lib.set_max_query_count(maxQueryCount)            
        if ret!=0:
            raise Exception(f'Max query count {maxQueryCount} is not allowed.')

    def set_report_query_details(self, value: bool):
        self.lib.set_report_query_details(value)            

    #def set_ref_min_aa_per_kmer(self, aaPerKmer: int):
    #    ret = self.lib.set_ref_min_aa_per_kmer(aaPerKmer)            
    #    if ret!=0:
    #        raise Exception(f'Min number of amino acids per reference kmer = {aaPerKmer} is not allowed.')

    def set_max_familyref_count_for_one_aa_per_kmer(self, aaPerKmer: int):
        ret = self.lib.set_max_familyref_count_for_one_aa_per_kmer(aaPerKmer)            
        if ret!=0:
            raise Exception(f'Max count {aaPerKmer} for one amino acids per reference kmer is not allowed.')

    def set_min_familyref_count_for_kmer(self, count: int):
        ret = self.lib.set_min_familyref_count_for_kmer(count)            
        if ret!=0:
            raise Exception(f'Minimum family reference kmer count {count} is not allowed.')

    def set_min_proteins_per_family(self, minProtsPerFamily: int):
        ret = self.lib.set_min_proteins_per_family(minProtsPerFamily)            
        if ret!=0:
            raise Exception(f'Min number of proteins {minProtsPerFamily} per family is not allowed.')

    def set_min_link_number(self, minLinkNumber: int):
        ret = self.lib.set_min_link_number(minLinkNumber)            
        if ret!=0:
            raise Exception(f'Min number of proteins {minLinkNumber} per family is not allowed.')

    def set_verbose(self, value: bool):
        self.lib.set_verbose(value)            

    def set_bin_block_size(self, binBlockSize: int):
        ret = self.lib.set_bin_block_size(binBlockSize)            
        if ret!=0:
            raise Exception(f'Binary block size  {binBlockSize} is not allowed.')

    def set_n_threads(self, nThreads: int):
        ret = self.lib.set_n_threads(nThreads)            
        if ret!=0:
            raise Exception(f'Number of threads {nThreads} is not allowed.')

    def set_input_dir(self, dirPath: str):
        self.lib.set_input_dir.argtype = ct.c_char_p
        self.lib.set_input_dir.restype = ct.c_int
        self.lib.set_input_dir(bytes(dirPath,'utf-8'))

    def set_output_dir(self, dirPath: str):
        self.lib.set_output_dir.argtype = ct.c_char_p
        self.lib.set_output_dir.restype = ct.c_int
        self.lib.set_output_dir(bytes(dirPath,'utf-8'))

    def set_translation_table(self, translation_table: str):
        if len(translation_table)!=64:
            raise Exception(f'Translation table must be 64 characters long, not {len(translation_table)}.')
        self.lib.set_translation_table.argtype = ct.c_char_p
        self.lib.set_translation_table.restype = ct.c_int
        self.lib.set_translation_table(bytes(translation_table,'utf-8'))

    def build_ref(self, filename: str, saveFamilyList: bool):
        self.lib.build_ref.argtypes = [ct.c_char_p, ct.c_int]
        self.lib.build_ref.restype = ct.c_int
        ret = self.lib.build_ref(bytes(filename,'utf-8'),saveFamilyList)
        match ret:
            case 0:
                pass
            case 1:
                raise Exception(f'Unable to open file {filename}.')
            case 2:
                raise Exception(f'Unable to write family list; check that the output directory exists.')
            case _:  # else
                raise Exception(f'Unknown error.')

    def save_parameters_ref(self, filename: str):
        self.lib.save_parameters_ref.argtype = ct.c_char_p
        self.lib.save_parameters_ref.restype = ct.c_int
        ret = self.lib.save_parameters_ref(bytes(filename,'utf-8'))
        match ret:
            case 0:
                pass
            case 1:
                raise Exception(f'Unable to open file {filename}.')
            case _:  # else
                raise Exception(f'Unknown error.')

    def save_ref_txt(self, filename: str):
        self.lib.save_ref_txt.argtype = ct.c_char_p
        self.lib.save_ref_txt.restype = ct.c_int
        ret = self.lib.save_ref_txt(bytes(filename,'utf-8'))
        match ret:
            case 0:
                pass
            case 1:
                raise Exception(f'Unable to open file {filename}.')
            case 2:
                raise Exception(f'Error with memory allocation in save_ref_txt.')
            case _:  # else
                raise Exception(f'Unknown error.')

    def save_ref_bin(self, filename: str):
        self.lib.save_ref_bin.argtype = ct.c_char_p
        self.lib.save_ref_bin.restype = ct.c_int
        ret = self.lib.save_ref_bin(bytes(filename,'utf-8'))
        match ret:
            case 0:
                pass
            case 1:
                raise Exception(f'Unable to open file {filename}.')
            case 2:
                raise Exception(f'Error with memory allocation in save_ref_txt.')
            case _:  # else
                raise Exception(f'Unknown error.')

    def load_ref_txt(self, filename: str):
        self.lib.load_ref_txt.argtype = ct.c_char_p
        self.lib.load_ref_txt.restype = ct.c_int
        ret = self.lib.load_ref_txt(bytes(filename,'utf-8'))
        match ret:
            case 0:
                pass
            case 1:
                raise Exception(f'Unable to open file {filename}.')
            case _:  # else
                raise Exception(f'Unknown error.')

    def load_ref_bin(self, filename: str):
        self.lib.load_ref_bin.argtype = ct.c_char_p
        self.lib.load_ref_bin.restype = ct.c_int
        ret = self.lib.load_ref_bin(bytes(filename,'utf-8'))
        match ret:
            case 0:
                pass
            case 1:
                raise Exception(f'Unable to open file {filename}.')
            case _:  # else
                raise Exception(f'Unknown error.')

    def tst_compare_ref_to_file_txt(self, filename: str):
        self.lib.tst_compare_ref_to_file_txt.argtype = ct.c_char_p
        self.lib.tst_compare_ref_to_file_txt.restype = ct.c_int
        ret = self.lib.tst_compare_ref_to_file_txt(bytes(filename,'utf-8'))
        match ret:
            case 0:
                pass
            case 1:
                raise Exception(f'Unable to open file {filename}.')
            case 4:
                raise Exception(f'Detected a mismatch while comparing reference to a file.')
            case 5:
                raise Exception(f'Preloaded reference contains kmers not present in {filename}.')
            case _:  # else
                raise Exception(f'Unknown error.')

    def tst_compare_ref_to_file_bin(self, filename: str):
        self.lib.tst_compare_ref_to_file_bin.argtype = ct.c_char_p
        self.lib.tst_compare_ref_to_file_bin.restype = ct.c_int
        ret = self.lib.tst_compare_ref_to_file_bin(bytes(filename,'utf-8'))
        match ret:
            case 0:
                pass
            case 1:
                raise Exception(f'Unable to open file {filename}.')
            case 2:
                raise Exception(f'Error parsing {filename}.')
            case 3:
                raise Exception(f'kmer size loaded from {filename} does not match the current settings.')
            case 4:
                raise Exception(f'Detected a mismatch while comparing reference to a file.')
            case 5:
                raise Exception(f'Preloaded reference contains kmers not present in {filename}.')
            case _:  # else
                raise Exception(f'Unknown error.')

    def tst_get_ref_probs_for_aa_kmer(self, kmer: str):
        buf = (ct.c_double*20)()
        self.lib.tst_get_ref_probs_for_aa_kmer.restype = ct.c_int
        ret = self.lib.tst_get_ref_probs_for_aa_kmer(bytes(kmer,'utf-8'), ct.byref(buf))
        match ret:
            case 0:
                pass
            case 1:
                raise Exception(f'Reference has not been set.')
            case 2:
                raise Exception(f'The size of {kmer} does not match the current kmer length setting.')
            case 3:
                raise Exception(f'Some of the characters in kmer {kmer} are not in the amino acid alphabet.')
            case 4:
                raise Exception(f'Kmer {kmer} is not in the reference.')
            case _:  # else
                raise Exception(f'Unknown error.')
        return np.array(buf)


    def tst_reset_query(self):
        self.lib.tst_reset_query()            

    def batch_process_query(self, filename: str):
        self.lib.batch_process_query.argtype = ct.c_char_p
        self.lib.batch_process_query.restype = ct.c_int
        ret = self.lib.batch_process_query(bytes(filename,'utf-8'))
        match ret:
            case 0:
                pass
            case 1:
                raise Exception(f'Unable to open file {filename}.')
            case 2:
                raise Exception(f'Translation table has not been set.')
            case 3:
                raise Exception(f'Reference has not been set.')
            case _:  # else
                raise Exception(f'Unknown error.')

    def thread_process_query(self, filename: str):
        self.lib.thread_process_query.argtype = ct.c_char_p
        self.lib.thread_process_query.restype = ct.c_int
        ret = self.lib.thread_process_query(bytes(filename,'utf-8'))
        match ret:
            case -1:
                raise Exception(f'All available threads are busy.')
            case _:  # else
                pass
        return ret

    def is_thread_alive(self, threadNumber: int):
        self.lib.is_thread_alive.argtype = ct.c_int
        self.lib.is_thread_alive.restype = ct.c_int
        ret = self.lib.is_thread_alive(threadNumber)            
        if ret==-1:
            raise Exception(f'The thread number {threadNumber} is invalid.')
        if ret==0:
            return False
        else:
            return True

    def process_query(self, filename: str):
        self.lib.process_query.argtype = ct.c_char_p
        self.lib.process_query(bytes(filename,'utf-8'))

    def tst_process_query_orf1f(self, filename: str):
        self.lib.tst_process_query_orf1f.argtype = ct.c_char_p
        self.lib.tst_process_query_orf1f.restype = ct.c_int
        ret = self.lib.tst_process_query_orf1f(bytes(filename,'utf-8'))
        if ret==1:
            raise Exception(f'Unable to open file {filename}.')
        elif ret==2:
            raise Exception(f'Translation table has not been set.')

    def tst_process_query_orf1r(self, filename: str):
        self.lib.tst_process_query_orf1r.argtype = ct.c_char_p
        self.lib.tst_process_query_orf1r.restype = ct.c_int
        ret = self.lib.tst_process_query_orf1r(bytes(filename,'utf-8'))
        if ret==1:
            raise Exception(f'Unable to open file {filename}.')
        elif ret==2:
            raise Exception(f'Translation table has not been set.')

    def tst_process_query_orf2f(self, filename: str):
        self.lib.tst_process_query_orf2f.argtype = ct.c_char_p
        self.lib.tst_process_query_orf2f.restype = ct.c_int
        ret = self.lib.tst_process_query_orf2f(bytes(filename,'utf-8'))
        if ret==1:
            raise Exception(f'Unable to open file {filename}.')
        elif ret==2:
            raise Exception(f'Translation table has not been set.')

    def tst_process_query_orf2r(self, filename: str):
        self.lib.tst_process_query_orf2r.argtype = ct.c_char_p
        self.lib.tst_process_query_orf2r.restype = ct.c_int
        ret = self.lib.tst_process_query_orf2r(bytes(filename,'utf-8'))
        if ret==1:
            raise Exception(f'Unable to open file {filename}.')
        elif ret==2:
            raise Exception(f'Translation table has not been set.')

    def tst_process_query_orf3f(self, filename: str):
        self.lib.tst_process_query_orf3f.argtype = ct.c_char_p
        self.lib.tst_process_query_orf3f.restype = ct.c_int
        ret = self.lib.tst_process_query_orf3f(bytes(filename,'utf-8'))
        if ret==1:
            raise Exception(f'Unable to open file {filename}.')
        elif ret==2:
            raise Exception(f'Translation table has not been set.')

    def tst_process_query_orf3r(self, filename: str):
        self.lib.tst_process_query_orf3r.argtype = ct.c_char_p
        self.lib.tst_process_query_orf3r.restype = ct.c_int
        ret = self.lib.tst_process_query_orf3r(bytes(filename,'utf-8'))
        if ret==1:
            raise Exception(f'Unable to open file {filename}.')
        elif ret==2:
            raise Exception(f'Translation table has not been set.')

    def tst_calculate_decoding_probs(self):
        self.lib.tst_calculate_decoding_probs()

    def tst_write_decoding_probs(self, filename: str):
        self.lib.tst_write_decoding_probs.argtype = ct.c_char_p
        self.lib.tst_write_decoding_probs.restype = ct.c_int
        ret = self.lib.tst_write_decoding_probs(bytes(filename,'utf-8'))
        if ret==1:
            raise Exception(f'Unable to write file {filename}.')
        
    def tst_print_ref_kmers(self, showExp_P: bool):
        self.lib.tst_print_ref_kmers(int(showExp_P))

    def tst_print_query_kmers(self):
        self.lib.tst_print_query_kmers()

    def tst_get_n_consensus(self):
        buf = (ct.c_uint64*64)()
        self.lib.tst_get_n_consensus.restype = None
        self.lib.tst_get_n_consensus(ct.byref(buf))
        return np.array(buf)

    def tst_get_decoding_probs(self):
        buf = (ct.c_double*21*64)()
        self.lib.tst_get_decoding_probs.restype = None
        self.lib.tst_get_decoding_probs(ct.byref(buf))
        return np.array(buf)

    def get_ref_size(self):
        self.lib.get_ref_size.restype = ct.c_ulong
        return self.lib.get_ref_size()

    def free_resources(self):
        self.lib.free_resources.argtype = None
        self.lib.free_resources.restype = None
        self.lib.free_resources()
