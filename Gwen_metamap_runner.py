#!/usr/bin/env python3

import sys
sys.path.insert(0, '/users/knarsinh/projects/clinical_trials/metamap/pymetamap') # for hypatia

# sys.path.insert(0, '/Volumes/TOSHIBA_EXT/ISB/clinical_trials/pymetamap-master') # for local

from pymetamap import MetaMap  # https://github.com/AnthonyMRios/pymetamap/blob/master/pymetamap/SubprocessBackend.py
import csv
import pandas as pd
import os
import subprocess
from time import sleep


metamap_base_dir = "/users/knarsinh/projects/clinical_trials/metamap/public_mm/"    # /users/knarsinh/projects/clinical_trials/metamap/public_mm 
metamap_bin_dir = 'bin/metamap20'

# metamap_base_dir = '/Volumes/TOSHIBA_EXT/ISB/clinical_trials/public_mm/' # for running on local
# metamap_bin_dir = 'bin/metamap18'

mm = MetaMap.get_instance(metamap_base_dir + metamap_bin_dir)

def start_metamap_servers(metamap_base_dir, metamap_bin_dir):
    global metamap_pos_server_dir
    global metamap_wsd_server_dir
    metamap_pos_server_dir = 'bin/skrmedpostctl' # Part of speech tagger
    metamap_wsd_server_dir = 'bin/wsdserverctl' # Word sense disambiguation 
    
    metamap_executable_path_pos = os.path.join(metamap_base_dir, metamap_pos_server_dir)
    metamap_executable_path_wsd = os.path.join(metamap_base_dir, metamap_wsd_server_dir)
    command_pos = [metamap_executable_path_pos, 'start']
    command_wsd = [metamap_executable_path_wsd, 'start']

    # Start servers, with open portion redirects output of metamap server printing output to NULL
    with open(os.devnull, "w") as fnull:
        result_post = subprocess.call(command_pos, stdout = fnull, stderr = fnull)
        result_wsd = subprocess.call(command_wsd, stdout = fnull, stderr = fnull)
    sleep(5)

def stop_metamap_servers(metamap_base_dir, metamap_bin_dir):
    metamap_executable_path_pos = os.path.join(metamap_base_dir, metamap_pos_server_dir)
    metamap_executable_path_wsd = os.path.join(metamap_base_dir, metamap_wsd_server_dir)
    command_pos = [metamap_executable_path_pos, 'stop']
    command_wsd = [metamap_executable_path_wsd, 'stop']
    
    # Stop servers, with open portion redirects output of metamap server printing output to NULL
    with open(os.devnull, "w") as fnull:
        result_post = subprocess.call(command_pos, stdout = fnull, stderr = fnull)
        result_wsd = subprocess.call(command_wsd, stdout = fnull, stderr = fnull)
    sleep(5)  




start_metamap_servers(metamap_base_dir, metamap_bin_dir)



with open("example_terms_list.txt", 'r') as termsfile:
    terms = termsfile.read().splitlines()
    terms = list(set(terms))
    for term in terms:
        sents = [term]
        concepts,error = mm.extract_concepts(sents)
        for concept in concepts:
            with open("mapped_terms_output.tsv", newline='', mode='a+', encoding="utf-8") as output:
                tsv_output = csv.writer(output, delimiter='\t')
                term_mm_list = []
                term_mm_list.extend([term, concept])
                tsv_output.writerow(term_mm_list)





# terms_to_map = pd.read_csv("example_terms_list.txt", sep ="\t", index_col=False, header=0, on_bad_lines = 'skip', encoding="utf-8")


# terms_to_map = terms_to_map.unique().tolist()
# for term in terms_to_map:
#     print(term)
#     sents = [term]
#     concepts,error = mm.extract_concepts(sents)
#     for concept in concepts:
#         print(concept)
#         with open("mapped_terms_output.tsv", mode='a+', encoding="utf-8") as output:
#             output.write(str(concept))
#             output.write("\n")

stop_metamap_servers(metamap_base_dir, metamap_bin_dir)




#             # open mapping cache to add mapped terms
#     mapping_filename = "mapping_cache.tsv"
#     if os.path.exists(mapping_filename):
#         output = open(mapping_filename, 'a', newline='', encoding="utf-8") 
#         output.close()
#     else:
#         output = open(mapping_filename, 'w+', newline='', encoding='utf-8')
#         col_names = ['mapping_tool', 'term_type', 'clintrial_term', 'input_term', 'mapping_tool_response', 'score']
#         csv_writer = csv.writer(output, delimiter='\t')
#         csv_writer.writerow(col_names)
#         output.close()




# start_metamap_servers(metamap_dirs)
# sents = ["haemodynamic effects of dexmedetomidine"]  # this has to be 1 sentence
# concepts,error = mm.extract_concepts(sents)
# for concept in concepts:
#     print(concept)
#     print("\n")
# stop_metamap_servers(metamap_dirs)










# # sents = ["haemodynamic effects of dexmedetomidine"]  # this has to be 1 sentence
# # concepts,error = mm.extract_concepts(sents)


# # for concept in concepts:
# #     print(concept)
# #     print("\n")
# # stop_metamap_servers(metamap_dirs)




# # def run_mappers(term_pair, params, term_type, mapping_filename):

# #     output = open("output_from_metamap.tsv", 'a', newline='', encoding="utf-8") 
# #     csv_writer = csv.writer(output, delimiter='\t')

# #     orig_term = term_pair[0]
# #     input_term = term_pair[1]
# #     from_mapper = []
# #     mm = MetaMap.get_instance(metamap_dirs["metamap_base_dir"] + metamap_dirs["metamap_bin_dir"])
