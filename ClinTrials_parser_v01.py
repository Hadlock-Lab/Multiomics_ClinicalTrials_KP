#!/usr/bin/env conda run -n ct_extract_env python

import pandas as pd
import os
import json

def parse_edges(data_folder):
    filename = "ClinTrials_KG_edges_v01_3.csv"
    filepath = os.path.join(data_folder, filename)

    edges = pd.read_csv(filepath, sep='\t')
    edges.rename(columns={'subject': 'disease', 'object': 'intervention', 'subject_name': 'disease_name', 'object_name': 'intervention_name'}, inplace=True)
    edges.rename(columns={'disease': 'object', 'intervention': 'subject', 'disease_name': 'object_name', 'intervention_name': 'subject_name'}, inplace=True)

    for index, row in edges.iterrows():
        id_dict = {}
        subject_dict = {}
        association_dict = {}
        object_dict = {}
        source_dict = {}

        id_dict["_id"] = "{}_{}_{}".format(row["nctid"].split("NCT")[1], row["subject"].split(':')[1], row["object"].split(':')[1])

        subject_dict["{}".format(row["subject"].split(':')[0])] = "{}".format(row["subject"].split(':')[1])
        subject_dict["name"] = row["subject_name"]
        # subject_dict["{}_semantic_types".format(row["subject"].split(':')[0])] = "TBD" # fix in next version
        subject_dict["type"] = "biolink:Treatment"

        association_dict["predicate"] = "{}".format(row["predicate"].split(':')[1])
        association_dict["edge_attributes"] = []
        association_dict["edge_attributes"].append(
            {"attribute_type_id":"clinicaltrials_id",
             "value":row["nctid"]
            }
        )
        association_dict["edge_attributes"].append(
            {"attribute_type_id":"biolink:aggregator_knowledge_source",
             "value":"infores:aact"}
            )
        association_dict["edge_attributes"].append(
            {"attribute_type_id": "biolink:primary_knowledge_source",
             "value": "infores:clinicaltrials"}
        )
        association_dict["edge_attributes"].append(
        {"attribute_type_id": "biolink:aggregator_knowledge_source",
         "value": "infores:biothings-multiomics-clinicaltrials"})
        
        association_dict["edge_attributes"].append(
            {
                "attribute_type_id": "biolink:knowledge_level",
                "value": "biolink:knowledge_assertion"
            }
        )
        association_dict["edge_attributes"].append(
            {
                "attribute_type_id": "biolink:agent_type",
                "value": "biolink:text_mining_agent"
            }
        )

        object_dict["{}".format(row["object"].split(':')[0])] = "{}".format(row["object"].split(':')[1])
        object_dict["name"] = row["object_name"]
        object_dict["type"] = "biolink:DiseaseorPhenotypicFeature"
        # object_dict["{}_semantic_types".format(row["object"].split(':')[0])] = "TBD" # fix in next version

        source_dict["edge_sources"] = []
        source_dict["edge_sources"].append(
            {
                "resource_id": "infores:biothings-multiomics-clinicaltrials",
                "resource_role": "aggregator_knowledge_source"
            }
        )

        source_dict["edge_sources"].append(
            {
                "resource_id": "infores:aact",
                "resource_role": "aggregator_knowledge_source"
            }
        )
        source_dict["edge_sources"].append(
            {
                "resource_id": "infores:clinicaltrials",
                "resource_role": "primary_knowledge_source"
            }
        )

        id_dict["subject"] = subject_dict
        id_dict["association"] = association_dict
        id_dict["object"] = object_dict 
        id_dict["source"] = source_dict

        # print(json.dumps(id_dict,sort_keys=True, indent=2))

        # yield the JSON one by one
        yield id_dict # comment for testing
        
def main():
    # data_folder = "../outputs/version_1" # uncomment for testing
    parse_edges(data_folder) # uncomment for testing

if __name__ == "__main__":
    main()


