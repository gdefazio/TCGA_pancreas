#!/usr/bin/env python
__name__ = 'TCGA querying'
__author__ = 'G. Defazio'
__version__ = 'draft'

import json
import sys
from datetime import datetime
import pandas as pd
import requests
import subprocess as sp
from shutil import rmtree
from shlex import split as shlit
import os
from os.path import join as pjoin
from functools import reduce
import argparse
import argcomplete
from gzip import open as gzopen
from concurrent.futures import ThreadPoolExecutor as TPE


def split_options():

    parser = argparse.ArgumentParser(
        description="""kMetaShot reference library builder""",
        prefix_chars="--")

    parser.add_argument("-p", "--primary_site", type=str,
                        nargs='+',
                        help="Primary site of cancer. Eg. breast",
                        action="store", required=False)
    parser.add_argument("-r", "--program_name", type=str,
                        nargs='+',
                        help="Program name Eg. TGCA",
                        action="store", required=False)
    parser.add_argument("-d", "--disease", type=str,
                        nargs='+',
                        help="Disease. Eg. breast cancer",
                        action="store", required=False)
    parser.add_argument('-o', '--output',
                        help="Output path. If not exists, it creates dir.",
                        action='store', required=True, type=str)
    parser.add_argument('-x', '--prefix',
                        help="Prefix to use in output dir and file names",
                        action='store', required=True, type=str)
    parser.add_argument('-s', '--store',
                        help="It retains sample files, if indicated.",
                        action='store_true')

    argcomplete.autocomplete(parser)
    return parser.parse_args()


def filter_gen(primary_site=None, program_name=None, disease=None):
    filters = {"op": "AND",
               "content": [{"op": "in",
                            "content": {
                                "field": "files.data_category",
                                "value": ["transcriptome profiling"]}}]}
    if primary_site is not None:
        if type(primary_site) is list:
            filters['content'].append({"op": "in",
                                       "content": {
                                           "field": "cases.primary_site",
                                           "value": primary_site}})
        elif type(primary_site) is str:
            primary_site = [primary_site]
            filters['content'].append({"op": "in",
                                       "content": {
                                           "field": "cases.primary_site",
                                           "value": primary_site}})
        else:
            raise NotImplementedError

    if program_name is not None:
        if type(program_name) is list:
            filters['content'].append({"op": "in",
                                       "content": {
                                           "field": "cases.project.program.name",
                                           "value": program_name}})
        elif type(program_name) is str:
            program_name = [program_name]
            filters['content'].append({"op": "in",
                                       "content": {
                                           "field": "cases.project.program.name",
                                           "value": program_name}})
        else:
            raise NotImplementedError

    if disease is not None:
        if type(disease) is list:
            filters['content'].append({"op": "in",
                                       "content": {
                                           "field": "cases.disease_type",
                                           "value": disease}})
        elif type(disease) is str:
            disease = [disease]
            filters['content'].append({"op": "in",
                                       "content": {
                                           "field": "cases.disease_type",
                                           "value": disease}})
        else:
            raise NotImplementedError

    # print(filters)
    return filters


def multi_merge(left: pd.DataFrame, right: pd.DataFrame) -> pd.DataFrame:
    mgd = left.merge(right, on=subcols[:-1])
    return mgd


def download(_id, _name):
    _scr = os.path.join(data_endpt, _id)
    _file = os.path.join(csv_tmp, "%s.gz" % _name)
    r = requests.get(_scr, allow_redirects=True)
    with gzopen(_file, 'wb') as csv_file:
        csv_file.write(r.content)
    # print(_file)


if __name__ == 'TCGA querying':
    print("""
    ###################################################################
    name: TCGA querying
    auth.: G.Defazio
    vers.: draft
    This program downloads open access (no token) rna_seq files 
    by selecting for primary site, program name and disease.
    It merges samples tables in a single gzipped CSV table.
    
    Exe date: %s
    ###################################################################
    """ % datetime.now())
    files_endpt = 'https://api.gdc.cancer.gov/files'
    data_endpt = 'https://api.gdc.cancer.gov/data'
    args = split_options()
    primary_site, project_name, disease, output, prefix, store = \
        args.primary_site, args.program_name, args.disease, args.output, args.prefix, args.store
    print("""
    Args.:
        -primary site: %s
        -program name: %s
        -disease: %s
        -output: %s
        -prefix: %s
        -store: %s
    """ % (primary_site, project_name, disease, output, prefix, store))
    try:
        os.mkdir(output)
    except FileExistsError:
        pass
    flds = ['data_format',
            'access',
            'file_name',
            'submitter_id',
            'data_category',
            'cases.case_id',
            'cases.project.name',
            'cases.primary_site',
            'cases.disease_type',
            'cases.demographic.race',
            'cases.samples.tissue_type',
            'cases.diagnoses.tissue_or_organ_of_origin',
            'cases.exposures.asbestos_exposure',
            'cases.samples.portions.portion_id']
    filt = filter_gen(primary_site=primary_site,program_name=project_name,disease=disease)
    params = {'fields': ",".join(flds),
              'size': 1,
              'filters': json.dumps(filt)
              }
    response = requests.get(files_endpt, params=params)
    # print(json.dumps(response.json(), indent=2))
    items = response.json()['data']['pagination']['total']
    params = {'fields': ",".join(flds),
              'size': items,
              'filters': json.dumps(filt)
              }
    response = requests.get(files_endpt, params=params)
    with open('/tmp/tmp.json', 'w') as tmp_json:
        aux = response.json()
        tmp_json.write(json.dumps(aux['data']['hits'], indent=2))
    if os.path.getsize('/tmp/tmp.json') <= 2:
        sys.exit('No files are available for these arguments')
    else:
        print(os.path.getsize('/tmp/tmp.json'))
        # exe = '/home/giuseppedefazio/Documenti/lavoro_unifi/Ramazzotti_TCGA_interrogation/scripts/'
        exe = os.path.dirname(__file__)
        # print(exe)
        metadata = pjoin(output, '%s_metadata.csv' % prefix)
        sp.run(shlit('python %s hits /tmp/tmp.json %s' % (pjoin(exe, 'json_to_csv.py'), metadata)))
        newtab = pd.read_csv(metadata)
        newtab['info_type'] = newtab['hits_file_name'].apply(lambda x: x.split('.')[1])
        newtab = newtab[newtab['info_type'] == 'rna_seq']
        newtab.to_csv(metadata, index=None)
        csv_tmp = pjoin(output, '%s_samples_csv' % prefix)
        try:
            os.mkdir(csv_tmp)
        except FileExistsError:
            rmtree(csv_tmp)
            os.mkdir(csv_tmp)

        to_download = list()
        gene_count = newtab[(newtab['info_type'] == 'rna_seq') & (newtab['hits_access'] == 'open')]
        print('Download for %s rna_seq files' % gene_count.shape[0])
        for cid in gene_count.hits_cases_0_case_id.unique():
            to_download.extend(gene_count[(gene_count['hits_cases_0_case_id'] == cid)
                               ][['hits_id', 'hits_file_name']].to_numpy().tolist())

        with TPE(10) as tpe:
            tpe.map(lambda x: download(x[0], x[1]), to_download)
        subcols = ['gene_id', 'gene_name', 'gene_type', 'unstranded']
        done = reduce(lambda x, y: multi_merge(x, y), [pd.read_csv(pjoin(csv_tmp, el),
                                                                   usecols=subcols,
                                                                   sep='\t',
                                                                   skiprows=[0,2,3,4,5])
                                                       for el in os.listdir(csv_tmp)])
        columns = list(done.columns)
        columns.__setitem__(slice(3, None),
                            [newtab[newtab['hits_file_name'] == el[:-3]]['hits_id'].iloc[0] for
                             el in os.listdir(csv_tmp)])
        # print(os.listdir(csv_tmp))
        done.columns = columns
        done.to_csv(pjoin(output, '%s_final_table.csv.gz' % prefix),
                    compression='gzip',
                    index=False)
        if store is False:
            rmtree(csv_tmp)

        print('All done %s' % datetime.now())


# filters = {"op": "AND",
#            "content":[
#                {"op": "in",
#                 "content": {
#                     "field": "cases.primary_site",
#                     "value": ["breast"]}},
#                {"op": "in",
#                 "content": {
#                     "field": "cases.project.program.name",
#                     "value": ["TCGA"]}},
#                {"op": "in",
#                 "content": {
#                     "field": "cases.disease_type",
#                     "value": ["Adenomas and Adenocarcinomas"]}},
#                {"op": "in",
#                 "content": {
#                     "field": "files.data_category",
#                     "value": ["transcriptome profiling"]}}
#            ]
#            }# #
# # ############################
#
#
#
#
# class QueryFields4File:
#     def __init__(self):
#         self.fields = [
#             'access',
#             'acl',
#             'created_datetime',
#             'data_category',
#             'data_format',
#             'data_type',
#             'error_type',
#             'experimental_strategy',
#             'file_id',
#             'file_name',
#             'file_size',
#             'file_state',
#             'md5sum',
#             'origin',
#             'platform',
#             'revision',
#             'state',
#             'state_comment',
#             'submitter_id',
#             'tags',
#             'type',
#             'updated_datetime',
#             'analysis.analysis_id',
#             'analysis.analysis_type',
#             'analysis.created_datetime',
#             'analysis.state',
#             'analysis.submitter_id',
#             'analysis.updated_datetime',
#             'analysis.workflow_end_datetime',
#             'analysis.workflow_link',
#             'analysis.workflow_start_datetime',
#             'analysis.workflow_type',
#             'analysis.workflow_version',
#             'analysis.input_files.access',
#             'analysis.input_files.created_datetime',
#             'analysis.input_files.data_category',
#             'analysis.input_files.data_format',
#             'analysis.input_files.data_type',
#             'analysis.input_files.error_type',
#             'analysis.input_files.experimental_strategy',
#             'analysis.input_files.file_id',
#             'analysis.input_files.file_name',
#             'analysis.input_files.file_size',
#             'analysis.input_files.file_state',
#             'analysis.input_files.md5sum',
#             'analysis.input_files.platform',
#             'analysis.input_files.revision',
#             'analysis.input_files.state',
#             'analysis.input_files.state_comment',
#             'analysis.input_files.submitter_id',
#             'analysis.input_files.updated_datetime',
#             'analysis.metadata.read_groups.adapter_name',
#             'analysis.metadata.read_groups.adapter_sequence',
#             'analysis.metadata.read_groups.base_caller_name',
#             'analysis.metadata.read_groups.base_caller_version',
#             'analysis.metadata.read_groups.created_datetime',
#             'analysis.metadata.read_groups.experiment_name',
#             'analysis.metadata.read_groups.flow_cell_barcode',
#             'analysis.metadata.read_groups.includes_spike_ins',
#             'analysis.metadata.read_groups.instrument_model',
#             'analysis.metadata.read_groups.is_paired_end',
#             'analysis.metadata.read_groups.library_name',
#             'analysis.metadata.read_groups.library_preparation_kit_catalog_number',
#             'analysis.metadata.read_groups.library_preparation_kit_name',
#             'analysis.metadata.read_groups.library_preparation_kit_vendor',
#             'analysis.metadata.read_groups.library_preparation_kit_version',
#             'analysis.metadata.read_groups.library_selection',
#             'analysis.metadata.read_groups.library_strand',
#             'analysis.metadata.read_groups.library_strategy',
#             'analysis.metadata.read_groups.platform',
#             'analysis.metadata.read_groups.read_group_id',
#             'analysis.metadata.read_groups.read_group_name',
#             'analysis.metadata.read_groups.read_length',
#             'analysis.metadata.read_groups.RIN',
#             'analysis.metadata.read_groups.sequencing_center',
#             'analysis.metadata.read_groups.sequencing_date',
#             'analysis.metadata.read_groups.size_selection_range',
#             'analysis.metadata.read_groups.spike_ins_concentration',
#             'analysis.metadata.read_groups.spike_ins_fasta',
#             'analysis.metadata.read_groups.state',
#             'analysis.metadata.read_groups.submitter_id',
#             'analysis.metadata.read_groups.target_capture_kit_catalog_number',
#             'analysis.metadata.read_groups.target_capture_kit_name',
#             'analysis.metadata.read_groups.target_capture_kit_target_region',
#             'analysis.metadata.read_groups.target_capture_kit_vendor',
#             'analysis.metadata.read_groups.target_capture_kit_version',
#             'analysis.metadata.read_groups.to_trim_adapter_sequence',
#             'analysis.metadata.read_groups.updated_datetime',
#             'analysis.metadata.read_groups.read_group_qcs.adapter_content',
#             'analysis.metadata.read_groups.read_group_qcs.basic_statistics',
#             'analysis.metadata.read_groups.read_group_qcs.created_datetime',
#             'analysis.metadata.read_groups.read_group_qcs.encoding',
#             'analysis.metadata.read_groups.read_group_qcs.fastq_name',
#             'analysis.metadata.read_groups.read_group_qcs.kmer_content',
#             'analysis.metadata.read_groups.read_group_qcs.overrepresented_sequences',
#             'analysis.metadata.read_groups.read_group_qcs.per_base_n_content',
#             'analysis.metadata.read_groups.read_group_qcs.per_base_sequence_content',
#             'analysis.metadata.read_groups.read_group_qcs.per_base_sequence_quality',
#             'analysis.metadata.read_groups.read_group_qcs.per_sequence_gc_content',
#             'analysis.metadata.read_groups.read_group_qcs.per_sequence_quality_score',
#             'analysis.metadata.read_groups.read_group_qcs.per_tile_sequence_quality',
#             'analysis.metadata.read_groups.read_group_qcs.percent_gc_content',
#             'analysis.metadata.read_groups.read_group_qcs.read_group_qc_id',
#             'analysis.metadata.read_groups.read_group_qcs.sequence_duplication_levels',
#             'analysis.metadata.read_groups.read_group_qcs.sequence_length_distribution',
#             'analysis.metadata.read_groups.read_group_qcs.state',
#             'analysis.metadata.read_groups.read_group_qcs.submitter_id',
#             'analysis.metadata.read_groups.read_group_qcs.total_sequences',
#             'analysis.metadata.read_groups.read_group_qcs.updated_datetime',
#             'analysis.metadata.read_groups.read_group_qcs.workflow_end_datetime',
#             'analysis.metadata.read_groups.read_group_qcs.workflow_link',
#             'analysis.metadata.read_groups.read_group_qcs.workflow_start_datetime',
#             'analysis.metadata.read_groups.read_group_qcs.workflow_type',
#             'analysis.metadata.read_groups.read_group_qcs.workflow_version',
#             'annotations.annotation_id',
#             'annotations.case_id',
#             'annotations.case_submitter_id',
#             'annotations.category',
#             'annotations.classification',
#             'annotations.created_datetime',
#             'annotations.creator',
#             'annotations.entity_id',
#             'annotations.entity_submitter_id',
#             'annotations.entity_type',
#             'annotations.legacy_created_datetime',
#             'annotations.legacy_updated_datetime',
#             'annotations.notes',
#             'annotations.state',
#             'annotations.status',
#             'annotations.submitter_id',
#             'annotations.updated_datetime',
#             'archive.archive_id',
#             'archive.created_datetime',
#             'archive.data_category',
#             'archive.data_format',
#             'archive.data_type',
#             'archive.error_type',
#             'archive.file_name',
#             'archive.file_size',
#             'archive.file_state',
#             'archive.md5sum',
#             'archive.revision',
#             'archive.state',
#             'archive.state_comment',
#             'archive.submitter_id',
#             'archive.updated_datetime',
#             'associated_entities.case_id',
#             'associated_entities.entity_id',
#             'associated_entities.entity_submitter_id',
#             'associated_entities.entity_type',
#             'cases.aliquot_ids',
#             'cases.analyte_ids',
#             'cases.case_id',
#             'cases.created_datetime',
#             'cases.days_to_index',
#             'cases.portion_ids',
#             'cases.sample_ids',
#             'cases.slide_ids',
#             'cases.state',
#             'cases.submitter_aliquot_ids',
#             'cases.submitter_analyte_ids',
#             'cases.submitter_id',
#             'cases.submitter_portion_ids',
#             'cases.submitter_sample_ids',
#             'cases.submitter_slide_ids',
#             'cases.updated_datetime',
#             'cases.annotations.annotation_id',
#             'cases.annotations.case_id',
#             'cases.annotations.case_submitter_id',
#             'cases.annotations.category',
#             'cases.annotations.classification',
#             'cases.annotations.created_datetime',
#             'cases.annotations.creator',
#             'cases.annotations.entity_id',
#             'cases.annotations.entity_submitter_id',
#             'cases.annotations.entity_type',
#             'cases.annotations.legacy_created_datetime',
#             'cases.annotations.legacy_updated_datetime',
#             'cases.annotations.notes',
#             'cases.annotations.state',
#             'cases.annotations.status',
#             'cases.annotations.submitter_id',
#             'cases.annotations.updated_datetime',
#             'cases.demographic.created_datetime',
#             'cases.demographic.demographic_id',
#             'cases.demographic.ethnicity',
#             'cases.demographic.gender',
#             'cases.demographic.race',
#             'cases.demographic.state'
#             'cases.demographic.submitter_id',
#             'cases.demographic.updated_datetime'
#             'cases.demographic.year_of_birth',
#             'cases.demographic.year_of_death',
#             'cases.diagnoses.age_at_diagnosis',
#             'cases.diagnoses.classification_of_tumor',
#             'cases.diagnoses.created_datetime',
#             'cases.diagnoses.days_to_birth',
#             'cases.diagnoses.days_to_death',
#             'cases.diagnoses.days_to_last_follow_up',
#             'cases.diagnoses.days_to_last_known_disease_status',
#             'cases.diagnoses.days_to_recurrence',
#             'cases.diagnoses.diagnosis_id',
#             'cases.diagnoses.last_known_disease_status',
#             'cases.diagnoses.morphology',
#             'cases.diagnoses.primary_diagnosis',
#             'cases.diagnoses.prior_malignancy',
#             'cases.diagnoses.progression_or_recurrence',
#             'cases.diagnoses.site_of_resection_or_biopsy',
#             'cases.diagnoses.state',
#             'cases.diagnoses.submitter_id',
#             'cases.diagnoses.tissue_or_organ_of_origin'
#             'cases.diagnoses.tumor_grade',
#             'cases.diagnoses.tumor_stage',
#             'cases.diagnoses.updated_datetime',
#             'cases.diagnoses.vital_status',
#             'cases.diagnoses.treatments.created_datetime',
#             'cases.diagnoses.treatments.days_to_treatment',
#             'cases.diagnoses.treatments.state',
#             'cases.diagnoses.treatments.submitter_id',
#             'cases.diagnoses.treatments.therapeutic_agents',
#             'cases.diagnoses.treatments.treatment_id',
#             'cases.diagnoses.treatments.treatment_intent_type',
#             'cases.diagnoses.treatments.treatment_or_therapy',
#             'cases.diagnoses.treatments.updated_datetime',
#             'cases.exposures.alcohol_history',
#             'cases.exposures.alcohol_intensity',
#             'cases.exposures.bmi',
#             'cases.exposures.cigarettes_per_day',
#             'cases.exposures.created_datetime',
#             'cases.exposures.exposure_id',
#             'cases.exposures.height',
#             'cases.exposures.state',
#             'cases.exposures.submitter_id',
#             'cases.exposures.updated_datetime',
#             'cases.exposures.weight',
#             'cases.exposures.years_smoked',
#             'cases.family_histories.created_datetime',
#             'cases.family_histories.family_history_id',
#             'cases.family_histories.relationship_age_at_diagnosis',
#             'cases.family_histories.relationship_gender',
#             'cases.family_histories.relationship_primary_diagnosis',
#             'cases.family_histories.relationship_type',
#             'cases.family_histories.relative_with_cancer_history',
#             'cases.family_histories.state',
#             'cases.family_histories.submitter_id',
#             'cases.family_histories.updated_datetime',
#             'cases.files.created_datetime',
#             'cases.files.error_type',
#             'cases.files.file_id',
#             'cases.files.file_name',
#             'cases.files.file_size',
#             'cases.files.file_state',
#             'cases.files.md5sum',
#             'cases.files.state',
#             'cases.files.state_comment',
#             'cases.files.submitter_id',
#             'cases.files.updated_datetime',
#             'cases.project.dbgap_accession_number',
#             'cases.project.disease_type',
#             'cases.project.name',
#             'cases.project.primary_site',
#             'cases.project.project_id',
#             'cases.project.released',
#             'cases.project.state',
#             'cases.project.program.dbgap_accession_number',
#             'cases.project.program.name',
#             'cases.project.program.program_id',
#             'cases.samples.composition',
#             'cases.samples.created_datetime',
#             'cases.samples.current_weight',
#             'cases.samples.days_to_collection',
#             'cases.samples.days_to_sample_procurement',
#             'cases.samples.freezing_method',
#             'cases.samples.initial_weight',
#             'cases.samples.intermediate_dimension',
#             'cases.samples.is_ffpe',
#             'cases.samples.longest_dimension',
#             'cases.samples.oct_embedded',
#             'cases.samples.pathology_report_uuid',
#             'cases.samples.preservation_method',
#             'cases.samples.sample_id',
#             'cases.samples.sample_type',
#             'cases.samples.sample_type_id',
#             'cases.samples.shortest_dimension',
#             'cases.samples.state',
#             'cases.samples.submitter_id',
#             'cases.samples.time_between_clamping_and_freezing',
#             'cases.samples.time_between_excision_and_freezing',
#             'cases.samples.tissue_type',
#             'cases.samples.tumor_code',
#             'cases.samples.tumor_code_id',
#             'cases.samples.tumor_descriptor',
#             'cases.samples.updated_datetime',
#             'cases.samples.annotations.annotation_id',
#             'cases.samples.annotations.case_id',
#             'cases.samples.annotations.case_submitter_id',
#             'cases.samples.annotations.category',
#             'cases.samples.annotations.classification',
#             'cases.samples.annotations.created_datetime',
#             'cases.samples.annotations.creator',
#             'cases.samples.annotations.entity_id',
#             'cases.samples.annotations.entity_submitter_id',
#             'cases.samples.annotations.entity_type',
#             'cases.samples.annotations.legacy_created_datetime',
#             'cases.samples.annotations.legacy_updated_datetime',
#             'cases.samples.annotations.notes',
#             'cases.samples.annotations.state',
#             'cases.samples.annotations.status',
#             'cases.samples.annotations.submitter_id',
#             'cases.samples.annotations.updated_datetime',
#             'cases.samples.portions.created_datetime',
#             'cases.samples.portions.creation_datetime',
#             'cases.samples.portions.is_ffpe',
#             'cases.samples.portions.portion_id',
#             'cases.samples.portions.portion_number',
#             'cases.samples.portions.state',
#             'cases.samples.portions.submitter_id',
#             'cases.samples.portions.updated_datetime',
#             'cases.samples.portions.weight',
#             'cases.samples.portions.analytes.a260_a280_ratio',
#             'cases.samples.portions.analytes.amount',
#             'cases.samples.portions.analytes.analyte_id',
#             'cases.samples.portions.analytes.analyte_type',
#             'cases.samples.portions.analytes.analyte_type_id',
#             'cases.samples.portions.analytes.concentration',
#             'cases.samples.portions.analytes.created_datetime',
#             'cases.samples.portions.analytes.spectrophotometer_method',
#             'cases.samples.portions.analytes.state',
#             'cases.samples.portions.analytes.submitter_id',
#             'cases.samples.portions.analytes.updated_datetime',
#             'cases.samples.portions.analytes.well_number',
#             'cases.samples.portions.analytes.aliquots.aliquot_id',
#             'cases.samples.portions.analytes.aliquots.amount',
#             'cases.samples.portions.analytes.aliquots.analyte_type',
#             'cases.samples.portions.analytes.aliquots.analyte_type_id',
#             'cases.samples.portions.analytes.aliquots.concentration',
#             'cases.samples.portions.analytes.aliquots.created_datetime',
#             'cases.samples.portions.analytes.aliquots.source_center',
#             'cases.samples.portions.analytes.aliquots.state',
#             'cases.samples.portions.analytes.aliquots.submitter_id',
#             'cases.samples.portions.analytes.aliquots.updated_datetime',
#             'cases.samples.portions.analytes.aliquots.annotations.annotation_id',
#             'cases.samples.portions.analytes.aliquots.annotations.case_id',
#             'cases.samples.portions.analytes.aliquots.annotations.case_submitter_id',
#             'cases.samples.portions.analytes.aliquots.annotations.category',
#             'cases.samples.portions.analytes.aliquots.annotations.classification',
#             'cases.samples.portions.analytes.aliquots.annotations.created_datetime',
#             'cases.samples.portions.analytes.aliquots.annotations.creator',
#             'cases.samples.portions.analytes.aliquots.annotations.entity_id',
#             'cases.samples.portions.analytes.aliquots.annotations.entity_submitter_id',
#             'cases.samples.portions.analytes.aliquots.annotations.entity_type',
#             'cases.samples.portions.analytes.aliquots.annotations.legacy_created_datetime',
#             'cases.samples.portions.analytes.aliquots.annotations.legacy_updated_datetime',
#             'cases.samples.portions.analytes.aliquots.annotations.notes',
#             'cases.samples.portions.analytes.aliquots.annotations.state',
#             'cases.samples.portions.analytes.aliquots.annotations.status',
#             'cases.samples.portions.analytes.aliquots.annotations.submitter_id',
#             'cases.samples.portions.analytes.aliquots.annotations.updated_datetime',
#             'cases.samples.portions.analytes.aliquots.center.center_id',
#             'cases.samples.portions.analytes.aliquots.center.center_type',
#             'cases.samples.portions.analytes.aliquots.center.code',
#             'cases.samples.portions.analytes.aliquots.center.name',
#             'cases.samples.portions.analytes.aliquots.center.namespace',
#             'cases.samples.portions.analytes.aliquots.center.short_name',
#             'cases.samples.portions.analytes.annotations.annotation_id',
#             'cases.samples.portions.analytes.annotations.case_id',
#             'cases.samples.portions.analytes.annotations.case_submitter_id',
#             'cases.samples.portions.analytes.annotations.category',
#             'cases.samples.portions.analytes.annotations.classification',
#             'cases.samples.portions.analytes.annotations.created_datetime',
#             'cases.samples.portions.analytes.annotations.creator',
#             'cases.samples.portions.analytes.annotations.entity_id',
#             'cases.samples.portions.analytes.annotations.entity_submitter_id',
#             'cases.samples.portions.analytes.annotations.entity_type',
#             'cases.samples.portions.analytes.annotations.legacy_created_datetime',
#             'cases.samples.portions.analytes.annotations.legacy_updated_datetime',
#             'cases.samples.portions.analytes.annotations.notes',
#             'cases.samples.portions.analytes.annotations.state',
#             'cases.samples.portions.analytes.annotations.status',
#             'cases.samples.portions.analytes.annotations.submitter_id',
#             'cases.samples.portions.analytes.annotations.updated_datetime',
#             'cases.samples.portions.annotations.annotation_id',
#             'cases.samples.portions.annotations.case_id',
#             'cases.samples.portions.annotations.case_submitter_id',
#             'cases.samples.portions.annotations.category',
#             'cases.samples.portions.annotations.classification',
#             'cases.samples.portions.annotations.created_datetime',
#             'cases.samples.portions.annotations.creator',
#             'cases.samples.portions.annotations.entity_id',
#             'cases.samples.portions.annotations.entity_submitter_id',
#             'cases.samples.portions.annotations.entity_type',
#             'cases.samples.portions.annotations.legacy_created_datetime',
#             'cases.samples.portions.annotations.legacy_updated_datetime',
#             'cases.samples.portions.annotations.notes',
#             'cases.samples.portions.annotations.state',
#             'cases.samples.portions.annotations.status',
#             'cases.samples.portions.annotations.submitter_id',
#             'cases.samples.portions.annotations.updated_datetime',
#             'cases.samples.portions.center.center_id',
#             'cases.samples.portions.center.center_type',
#             'cases.samples.portions.center.code',
#             'cases.samples.portions.center.name',
#             'cases.samples.portions.center.namespace',
#             'cases.samples.portions.center.short_name',
#             'cases.samples.portions.slides.created_datetime',
#             'cases.samples.portions.slides.number_proliferating_cells',
#             'cases.samples.portions.slides.percent_eosinophil_infiltration',
#             'cases.samples.portions.slides.percent_granulocyte_infiltration',
#             'cases.samples.portions.slides.percent_inflam_infiltration',
#             'cases.samples.portions.slides.percent_lymphocyte_infiltration',
#             'cases.samples.portions.slides.percent_monocyte_infiltration',
#             'cases.samples.portions.slides.percent_necrosis',
#             'cases.samples.portions.slides.percent_neutrophil_infiltration',
#             'cases.samples.portions.slides.percent_normal_cells',
#             'cases.samples.portions.slides.percent_stromal_cells',
#             'cases.samples.portions.slides.percent_tumor_cells',
#             'cases.samples.portions.slides.percent_tumor_nuclei',
#             'cases.samples.portions.slides.section_location',
#             'cases.samples.portions.slides.slide_id',
#             'cases.samples.portions.slides.state',
#             'cases.samples.portions.slides.submitter_id',
#             'cases.samples.portions.slides.updated_datetime',
#             'cases.samples.portions.slides.annotations.annotation_id',
#             'cases.samples.portions.slides.annotations.case_id',
#             'cases.samples.portions.slides.annotations.case_submitter_id',
#             'cases.samples.portions.slides.annotations.category',
#             'cases.samples.portions.slides.annotations.classification',
#             'cases.samples.portions.slides.annotations.created_datetime',
#             'cases.samples.portions.slides.annotations.creator',
#             'cases.samples.portions.slides.annotations.entity_id',
#             'cases.samples.portions.slides.annotations.entity_submitter_id',
#             'cases.samples.portions.slides.annotations.entity_type',
#             'cases.samples.portions.slides.annotations.legacy_created_datetime',
#             'cases.samples.portions.slides.annotations.legacy_updated_datetime',
#             'cases.samples.portions.slides.annotations.notes',
#             'cases.samples.portions.slides.annotations.state',
#             'cases.samples.portions.slides.annotations.status',
#             'cases.samples.portions.slides.annotations.submitter_id',
#             'cases.samples.portions.slides.annotations.updated_datetime',
#             'cases.summary.file_count',
#             'cases.summary.file_size',
#             'cases.summary.data_categories.data_category',
#             'cases.summary.data_categories.file_count',
#             'cases.summary.experimental_strategies.experimental_strategy',
#             'cases.summary.experimental_strategies.file_count',
#             'cases.tissue_source_site.bcr_id',
#             'cases.tissue_source_site.code',
#             'cases.tissue_source_site.name',
#             'cases.tissue_source_site.project',
#             'cases.tissue_source_site.tissue_source_site_id',
#             'center.center_id',
#             'center.center_type',
#             'center.code',
#             'center.name',
#             'center.namespace',
#             'center.short_name',
#             'downstream_analyses.analysis_id',
#             'downstream_analyses.analysis_type',
#             'downstream_analyses.created_datetime',
#             'downstream_analyses.state',
#             'downstream_analyses.submitter_id',
#             'downstream_analyses.updated_datetime',
#             'downstream_analyses.workflow_end_datetime',
#             'downstream_analyses.workflow_link',
#             'downstream_analyses.workflow_start_datetime',
#             'downstream_analyses.workflow_type',
#             'downstream_analyses.workflow_version',
#             'downstream_analyses.output_files.access',
#             'downstream_analyses.output_files.created_datetime',
#             'downstream_analyses.output_files.data_category',
#             'downstream_analyses.output_files.data_format',
#             'downstream_analyses.output_files.data_type',
#             'downstream_analyses.output_files.error_type',
#             'downstream_analyses.output_files.experimental_strategy',
#             'downstream_analyses.output_files.file_id',
#             'downstream_analyses.output_files.file_name',
#             'downstream_analyses.output_files.file_size',
#             'downstream_analyses.output_files.file_state',
#             'downstream_analyses.output_files.md5sum',
#             'downstream_analyses.output_files.platform',
#             'downstream_analyses.output_files.revision',
#             'downstream_analyses.output_files.state',
#             'downstream_analyses.output_files.state_comment',
#             'downstream_analyses.output_files.submitter_id',
#             'downstream_analyses.output_files.updated_datetime',
#             'index_files.access',
#             'index_files.created_datetime',
#             'index_files.data_category',
#             'index_files.data_format',
#             'index_files.data_type',
#             'index_files.error_type',
#             'index_files.experimental_strategy',
#             'index_files.file_id',
#             'index_files.file_name',
#             'index_files.file_size',
#             'index_files.file_state',
#             'index_files.md5sum',
#             'index_files.platform',
#             'index_files.revision',
#             'index_files.state',
#             'index_files.state_comment',
#             'index_files.submitter_id',
#             'index_files.updated_datetime',
#             'metadata_files.access',
#             'metadata_files.created_datetime',
#             'metadata_files.data_category',
#             'metadata_files.data_format',
#             'metadata_files.data_type',
#             'metadata_files.error_type',
#             'metadata_files.file_id',
#             'metadata_files.file_name',
#             'metadata_files.file_size',
#             'metadata_files.file_state',
#             'metadata_files.md5sum',
#             'metadata_files.state',
#             'metadata_files.state_comment',
#             'metadata_files.submitter_id',
#             'metadata_files.type',
#             'metadata_files.updated_datetime'
#         ]