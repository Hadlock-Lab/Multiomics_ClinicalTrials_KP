#!/usr/bin/env python3

import os
import pandas as pd
# shell command: brew install postgresql 
import psycopg2
import pandas.io.sql as sqlio
# from config import API_KEY
from time import sleep
import re
import sys
import multiprocessing
import pickle
import json
from itertools import repeat


# adding Folder_2/subfolder to the system path
# sys.path.insert(0, '/Volumes/TOSHIBA_EXT/ISB/clinical_trials/pymetamap-master')  # for running on local
sys.path.insert(0, '/users/knarsinh/projects/clinical_trials/metamap')

from pymetamap import MetaMap

# uncomment when switching to Linux
# sys.path.insert(0, '/users/knarsinh/projects/clinical_trials/metamap/pymetamap/pymetamap')
# from pymetamap import MetaMap


# import requests
# from bs4 import BeautifulSoup


# Setup UMLS Server global vars
# metamap_base_dir = '/Volumes/TOSHIBA_EXT/ISB/clinical_trials/public_mm/' # for running on local
metamap_base_dir = "/users/knarsinh/projects/clinical_trials/metamap/public_mm"
metamap_bin_dir = 'bin/metamap18'
metamap_pos_server_dir = 'bin/skrmedpostctl'
metamap_wsd_server_dir = 'bin/wsdserverctl'


#####	----	****	----	----	****	----	----	****	----    #####

# ACCESSING DATA BY DOWNLOADING STATIC COPY OF DATABASE, AND INSTALLING POSTGRESQL SOFTWARE TO THEN POPULATE DATABASE ON MACHINE
# THEN CONNECTING TO DB

# connect to DB and get the column names of the table

# extract data from ClinicalTrials.gov
def get_trial_data():
    con = None
    con = psycopg2.connect(database="aact")
    con.rollback()
    cursor = con.cursor()

    con.autocommit = True # SQL statement is treated as a transaction and is automatically committed right after it is executed
    # grab the conditions
    sql = '''SELECT * FROM ctgov.conditions;'''
    cursor.execute(sql)
    column_names = [desc[0] for desc in cursor.description]
    tuples = cursor.fetchall()
    conditions_df = pd.DataFrame(tuples, columns=column_names)

    #grab the browse_conditions
    sql = '''SELECT * FROM ctgov.browse_conditions;'''
    cursor.execute(sql)
    column_names = [desc[0] for desc in cursor.description]
    tuples = cursor.fetchall()
    browse_conditions_df = pd.DataFrame(tuples, columns=column_names)

    #grab the interventions
    sql = '''SELECT * FROM ctgov.interventions;'''
    cursor.execute(sql)
    column_names = [desc[0] for desc in cursor.description]
    tuples = cursor.fetchall()
    interventions_df = pd.DataFrame(tuples, columns=column_names)

    #grab the browse_interventions
    sql = '''SELECT * FROM ctgov.browse_interventions;'''
    cursor.execute(sql)
    column_names = [desc[0] for desc in cursor.description]
    tuples = cursor.fetchall()
    browse_interventions_df = pd.DataFrame(tuples, columns=column_names)

    con.close()


    # rename and drop df relevant columns to prepare for merging
    interventions_df = interventions_df.rename(columns={'id': 'int_id',
                                                        'nct_id': 'int_nctid',
                                                        'intervention_type': 'int_type',
                                                        'name': 'int_name',
                                                        'description': 'int_description'})
    interventions_df = interventions_df.drop(columns=['int_id', 'int_description'])
    conditions_df = conditions_df.rename(columns={'id': 'con_id',
                                                  'nct_id': 'con_nctid',
                                                  'name': 'con_name',
                                                  'downcase_name': 'con_downcase_name'})
    conditions_df = conditions_df.drop(columns=['con_id', 'con_name'])
    browse_interventions_df = browse_interventions_df.rename(columns={'id': 'browseint_id',
                                                                      'nct_id': 'browseint_nctid',
                                                                      'mesh_term': 'browseint_meshterm',
                                                                      'downcase_mesh_term': 'browseint_meshterm_downcase',
                                                                      'mesh_type': 'browseint_meshtype'})

    browse_interventions_df = browse_interventions_df.drop(columns=['browseint_id', 'browseint_meshterm'])
    browse_conditions_df = browse_conditions_df.rename(columns={'id': 'browsecon_id',
                                                                'nct_id': 'browsecon_nctid',
                                                                'mesh_term': 'browsecon_meshterm',
                                                                'downcase_mesh_term': 'browsecon_meshterm_downcase',
                                                                'mesh_type': 'browsecon_meshtype'})
    browse_conditions_df = browse_conditions_df.drop(columns=['browsecon_id', 'browsecon_meshterm'])                                                                                                                          

    # merge conditions_df and interventions_df since they have relevant terms 
    df = pd.merge(conditions_df, interventions_df, left_on='con_nctid', right_on = 'int_nctid')
    df_dedup = df.drop_duplicates(subset = ['con_downcase_name', 'int_name'],
                                          keep = 'first').reset_index(drop = True)

    return df_dedup

    # with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
    #     display(df_dedup.head(100))


    # with pd.option_context('expand_frame_repr', False, 'display.max_rows', None):
    # 	print(df_dedup[['con_nctid', 'con_downcase_name', 'int_name']].head(100))


    # df_dedup[['con_downcase_name', 'int_name']].to_csv('ClinTrials_KG_nodes.csv', sep ='\t', index=False)

# the boss wants to see the properly formatted .tsv file output first, even if it's filled with dummy data...

# Massage data extracted from ClinicalTrials.gov
def preprocess_ct_data(ct_data):    
    # first get only relevant columns from DB
    ct_extract = pd.DataFrame(ct_data[['con_nctid', 'con_downcase_name', 'int_type', 'int_name']])
    ct_extract = ct_extract.rename(columns={'con_nctid': 'nctid'})
    # get CURIE column for nct_id column (https://bioregistry.io/registry/clinicaltrials)

    ct_extract['nctid_curie'] = ct_extract['nctid']
    ct_extract['nctid_curie'] = 'clinicaltrials:' + ct_extract['nctid'].astype(str)  # generate CURIEs for each clinical trials study

    ct_final = pd.DataFrame(columns=['subject','predicate','object', 'subject_name','object_name','category'])

    ct_final['subject'] = 'condition:' + ct_extract['con_downcase_name'].astype(str)
    ct_final['predicate'] = 'biolink:hypothesized_for'
    ct_final['object'] = 'intervention:' + ct_extract['int_name'].astype(str)   # this will not all be RxNorm CURIEs since some interventions are not drugs
    ct_final.subject_name = ct_extract.con_downcase_name
    ct_final.object_name = ct_extract.int_name
    ct_final.category = 'biolink:Association'
    ct_final['nctid'] = ct_extract['nctid']
    ct_final['nctid_curie'] = ct_extract['nctid_curie']

    return(ct_final)


def select_cols_to_preprocess(ct_preprocessed):
    # prep a list to hold dicts of de-asciied cols from the ClinTrials df
    processed_ct_cols = []
    metamap_compatible = {}
    

    # get unique lists of the concepts to MetaMap from the ClinTrials df
    conditions = list(set(list(ct_preprocessed['subject_name'])))
    interventions = list(filter(None, ct_preprocessed['object_name'].unique()))

    if not all(condition.isascii() for condition in conditions):
        print("Non-ASCII chars are detected in Conditions col from ClinTrial (not checking other cols), proceed with text processing")
        print("This step is unnecessary for MetaMap 2020+")
        # process lists to remove non-ascii chars (this is not required for MetaMap 2020!!!!)
        conditions_translated = dict((i, preprocess_cols_for_metamap(i)) for i in conditions)
        interventions_translated = dict((i, preprocess_cols_for_metamap(i)) for i in interventions)
        
        metamap_compatible["conditions"] = conditions_translated
        metamap_compatible["interventions"] = interventions_translated

    else:
        print("No non-ASCII chars detected or they are present, but we're not using MetaMap versions prior to 2020, no text processing required")
        metamap_compatible["conditions"] = conditions
        metamap_compatible["interventions"] = interventions

    return(metamap_compatible)
    

def preprocess_cols_for_metamap(text):
    non_ascii = "[^\x00-\x7F]"
    pattern = re.compile(r"[^\x00-\x7F]")
    non_ascii_text = re.sub(pattern, ' ', text)
    return non_ascii_text

# setup MetaMap
def start_metamap_servers():

    # uncomment when switching to Linux

    #  ***   -----   ******   -----   ******   -----   ******   -----   ******   -----   *** # 

    # # Start servers
    # os.system(metamap_base_dir + metamap_pos_server_dir + ' start') # Part of speech tagger
    # os.system(metamap_base_dir + metamap_wsd_server_dir + ' start') # Word sense disambiguation 
     
    # # Sleep a bit to give time for these servers to start up
    # sleep(40)
    #  ***   -----   ******   -----   ******   -----   ******   -----   ******   -----   *** # 


    # Start servers
    os.system(metamap_base_dir + metamap_pos_server_dir + ' start') # Part of speech tagger
    os.system(metamap_base_dir + metamap_wsd_server_dir + ' start') # Word sense disambiguation 
     
    # # Sleep a bit to give time for these servers to start up
    sleep(60)

def stop_metamap_servers():
    
    # Uncomment when switching to Linux
    #  ***   -----   ******   -----   ******   -----   ******   -----   ******   -----   *** # 

    # # Start servers
    # os.system(metamap_base_dir + metamap_pos_server_dir + ' stop') # Part of speech tagger
    # os.system(metamap_base_dir + metamap_wsd_server_dir + ' stop') # Word sense disambiguation 
    
    # # Sleep a bit to give time for these servers to start up
    #  ***   -----   ******   -----   ******   -----   ******   -----   ******   -----   *** # 

    # Stop servers
    os.system(metamap_base_dir + metamap_pos_server_dir + ' stop') # Part of speech tagger
    os.system(metamap_base_dir + metamap_wsd_server_dir + ' stop') # Word sense disambiguation 

def run_metamap(input_term, source_restriction):
    mm = MetaMap.get_instance(metamap_base_dir + metamap_bin_dir)
    concepts_dict = dict()
    if all(x is None for x in source_restriction):
        try:
            concepts,error = mm.extract_concepts([input_term],
                                                 word_sense_disambiguation = True,
                                                 prune = 10,
                                                 composite_phrase = 1)
            concepts_dict[input_term] = concepts
        except:
            concepts_dict[input_term] = None 
    else:
        try:
            concepts,error = mm.extract_concepts([input_term],
                                                 word_sense_disambiguation = True,
                                                 restrict_to_sources=source_restriction,
                                                 prune = 10,
                                                 composite_phrase = 1)
            concepts_dict[input_term] = concepts
        except:
            concepts_dict[input_term] = None 
    return(concepts_dict)


def write_output_out(output_to_store):
    # CAREFUL: Running WILL OVERWRITE METAMAPPED DISEASES AND INTERVENTIONS FILES, WHICH TAKE FOREVER TO GENERATE!

    # write mapped diseases and interventions to pickle file
    with open('metamapped_conditions.txt', 'wb') as file:
        pickle.dump(output_to_store.get("conditions"), file)
    with open('metamapped_interventions.txt', 'wb') as file:
        pickle.dump(output_to_store.get("interventions"), file)

    # write mapped diseases and interventions to human readable file
    with open('readable_metamapped_conditions.txt', 'w') as file:
        file.write(json.dumps(output_to_store.get("conditions"), indent=4))
    with open('readable_metamapped_interventions.txt', 'w') as file:
        file.write(json.dumps(output_to_store.get("interventions"), indent=4))


    # # write to pickle file
    # with open('metamapped_conditions.txt', 'wb') as file:
    #     pickle.dump(output_to_store, file)

    # # write to human readable file
    # with open('readable_metamapped_conditions.txt', 'w') as file:
    #     file.write(json.dumps(output_to_store, indent=4))


def get_mapped_ct_data(metamapped_concept, translation_dict):
    for key, val in metamapped_concept.items():
        for concept in val:
            concept = concept._asdict()
            row = [key]
            original_concept = (list(translation_dict.keys())[list(translation_dict.values()).index(key)])
            row.append(original_concept)
            row.extend([concept.get(k) for k in ['preferred_name', 'cui', 'score', 'semtypes']])
            return row # we only want the first mapped concept



def driver():
    # with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
    #     display(df_dedup.head(100))
    metamapped_dict = {}

    ct_data = get_trial_data()
    ct_preprocessed = preprocess_ct_data(ct_data)
    ct_processed = select_cols_to_preprocess(ct_preprocessed)
    # ct_processed.get("conditions")
    # print(ct_processed.get("interventions"))
    start_metamap_servers()
    
    import random
    conditions_proc = dict(random.sample(ct_processed.get("conditions").items(), 2000))
    interventions_proc = dict(random.sample(ct_processed.get("interventions").items(), 2000))

    # print(conditions_proc.values())

    # print(ct_processed.get("conditions").values())
    # print(type(ct_processed.get("conditions").values()))

    metamapped_conditions = multiprocessing.Pool(multiprocessing.cpu_count() - 1).starmap(run_metamap,
        zip(list(conditions_proc.values()),
            [['SNOMEDCT_US', 'SNOMEDCT_VET', 'ICD10CM', 'ICD10CM', 'ICD10PCS', 'ICD9CM', 'ICD9CM']]*len(list(conditions_proc.values()))))
    metamapped_interventions = multiprocessing.Pool(multiprocessing.cpu_count() - 1).starmap(run_metamap, 
                                                                             zip(list(interventions_proc.values()),
                                                                                 [[None]]*len(list(interventions_proc.values()))))

    # UNCOMMENT TO RUN ON FULL DATASET
    # metamapped_conditions = multiprocessing.Pool(multiprocessing.cpu_count() - 1).starmap(run_metamap,
    #     zip(list(ct_processed.get("conditions").values()),
    #         [['SNOMEDCT_US', 'SNOMEDCT_VET', 'ICD10CM', 'ICD10CM', 'ICD10PCS', 'ICD9CM', 'ICD9CM']]*len(list(ct_processed.get("conditions").values()))))

    # metamapped_interventions = multiprocessing.Pool(multiprocessing.cpu_count() - 1).starmap(run_metamap,
    #     zip(list(ct_processed.get("interventions").values()),
    #         [['SNOMEDCT_US', 'SNOMEDCT_VET', 'ICD10CM', 'ICD10CM', 'ICD10PCS', 'ICD9CM', 'ICD9CM']]*len(list(ct_processed.get("interventions").values()))))

    stop_metamap_servers()
    metamapped_dict["conditions"] = metamapped_conditions
    metamapped_dict["interventions"] = metamapped_interventions
    write_output_out(metamapped_dict)


    # with open('metamapped_conditions.txt', 'rb') as f:
    #     metamapped_conditions = pickle.load(f)
    
    # with open('metamapped_interventions.txt', 'rb') as f:
    #     metamapped_interventions = pickle.load(f)


    # print(metamapped_dict.get("conditions"))
    # print("\n\n\n\n")
    # print(ct_processed.get("conditions"))
    
    with multiprocessing.Pool(multiprocessing.cpu_count() - 1) as pool:
        final_metamapped_conditions = pool.starmap(get_mapped_ct_data,
            zip(metamapped_dict.get("conditions"),
                repeat(ct_processed.get("conditions"))))
        final_metamapped_interventions = pool.starmap(get_mapped_ct_data,
            zip(metamapped_dict.get("interventions"),
                repeat(ct_processed.get("interventions"))))

    master_condition_dict = {condition[0]:condition[1:] for condition in final_metamapped_conditions if condition}  
    master_intervention_dict = {intervention[0]:intervention[1:] for intervention in final_metamapped_interventions if intervention}    

    baseline_final_condition_dict = {}
    baseline_final_intervention_dict = {}

    for k, v in master_condition_dict.items():
        orig_con = k.lower()
        orig_con = orig_con.strip()
        mapped_con = str(v[1])
        mapped_con = mapped_con.lower()
        mapped_con = mapped_con.strip()
        if orig_con == mapped_con:
            baseline_final_condition_dict[k] = v[1:] 
            
    for k, v in master_intervention_dict.items():
        orig_int = k.lower()
        orig_int = orig_int.strip()
        mapped_int = str(v[1])
        mapped_int = mapped_int.lower()
        mapped_int = mapped_int.strip()
        if orig_int == mapped_int:
            baseline_final_intervention_dict[k] = v[1:] 

    # create CURIEs from values of baseline lists
    condition_curie_dict = {}
    intervention_curie_dict = {}

    for k, v in baseline_final_condition_dict.items():
        condition_curie_dict[k] = "UMLS:{}".format(v[1])    
    for k, v in baseline_final_intervention_dict.items():
        intervention_curie_dict[k] = "UMLS:{}".format(v[1])

    # print(intervention_curie_dict)

    ct_final_mapped = ct_preprocessed.copy()
    print(ct_final_mapped[:10])
    ct_final_mapped['subject'] = ct_preprocessed['subject_name'].map(condition_curie_dict)
    ct_final_mapped['object'] = ct_preprocessed['object_name'].map(intervention_curie_dict)

    ct_final_mapped = ct_final_mapped.dropna(subset=['subject','object'])
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
        print(ct_final_mapped.head(20))    
    
    # make dataframe for nodes file derived from edges dataframe (ct_final_mapped)
    ct_nodes_subjects = ct_final_mapped[['subject', 'subject_name']]
    ct_nodes_subjects.columns = ['id', 'name']
    ct_nodes_subjects = ct_nodes_subjects.assign(category='biolink:DiseaseOrPhenotypicFeature')

    ct_nodes_objects = ct_final_mapped[['object', 'object_name']]
    ct_nodes_objects.columns = ['id', 'name']
    ct_nodes_objects = ct_nodes_objects.assign(category='biolink:Treatment')

    ct_nodes = pd.concat([ct_nodes_subjects, ct_nodes_objects], ignore_index=True, axis=0)
    ct_nodes = ct_nodes.drop_duplicates('id')

    # write nodes and edges files

    ct_final_mapped.to_csv('ClinTrials_KG_edges_v01.csv', sep ='\t', index=False)
    ct_nodes.to_csv('ClinTrials_KG_nodes_v01.csv', sep ='\t', index=False)



    # print(baseline_final_intervention_dict)

    # print(master_intervention_dict) 
    # master_disease_dict = {disease[0]:disease[1:] for disease in final_metamapped_diseases if disease}  
    # master_intervention_dict = {intervention[0]:intervention[1:] for intervention in final_metamapped_interventions if intervention}    





    

    

    # print(metamapped_conditions)






    # with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
    #     print(ct_preprocessed.head(20))


driver()



#####	----	****	----	----	****	----	----	****	----    #####


