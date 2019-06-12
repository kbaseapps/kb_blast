# -*- coding: utf-8 -*-
import gzip
import os
import re
import subprocess
import sys
import traceback
import uuid
from datetime import datetime
from pprint import pformat

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import requests
from requests_toolbelt import MultipartEncoder

# SDK Utils
from installed_clients.AbstractHandleClient import AbstractHandle
from installed_clients.KBaseDataObjectToFileUtilsClient import KBaseDataObjectToFileUtils
from installed_clients.DataFileUtilClient import DataFileUtil as DFUClient
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.WorkspaceClient import Workspace as workspaceService


class BlastUtil:

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "1.1.0"
    GIT_URL = "https://github.com/kbaseapps/kb_blast.git"
    GIT_COMMIT_HASH = "0722ff0b7d723e654ef9ebe470e2b515d13671bc"

    #BEGIN_CLASS_HEADER
    workspaceURL = None
    shockURL     = None
    handleURL    = None
    callbackURL  = None
    scratch      = None

    Make_BLAST_DB = '/kb/module/blast/bin/makeblastdb'
    BLASTn        = '/kb/module/blast/bin/blastn'
    BLASTp        = '/kb/module/blast/bin/blastp'
    BLASTx        = '/kb/module/blast/bin/blastx'
    tBLASTn       = '/kb/module/blast/bin/tblastn'
    tBLASTx       = '/kb/module/blast/bin/tblastx'
    psiBLAST      = '/kb/module/blast/bin/psiblast'


    # timestamp
    def now_ISO(self):
        now_timestamp = datetime.now()
        now_secs_from_epoch = (now_timestamp - datetime(1970,1,1)).total_seconds()
        now_timestamp_in_iso = datetime.fromtimestamp(int(now_secs_from_epoch)).strftime('%Y-%m-%d_%T')
        return now_timestamp_in_iso

    # message logging
    def log(self, target, message):
        message = '['+self.now_ISO()+'] '+message
        if target is not None:
            target.append(message)
        print(message)
        sys.stdout.flush()


    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config, ctx):
        #BEGIN_CONSTRUCTOR
        self.config = config
        self.ctx = ctx

        self.workspaceURL = config['workspace-url']
        self.shockURL = config['shock-url']
        self.handleURL = config['handle-service-url']
        self.serviceWizardURL = config['service-wizard-url']

#        self.callbackURL = os.environ['SDK_CALLBACK_URL'] if os.environ['SDK_CALLBACK_URL'] != None else 'https://kbase.us/services/njs_wrapper'
        self.callbackURL = os.environ.get('SDK_CALLBACK_URL')
        if self.callbackURL == None:
            raise ValueError ("SDK_CALLBACK_URL not set in environment")

        self.scratch = os.path.abspath(config['scratch'])
        if self.scratch == None:
            self.scratch = os.path.join('/kb','module','local_scratch')
        if not os.path.exists(self.scratch):
            os.makedirs(self.scratch)

        #END_CONSTRUCTOR
        pass


    #### Sequence Validation
    ##
    def validateSeq (self, seq_type, sequence_str, header_id):
        console = []
        PROT_pattern = re.compile("^[acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWYxX ]+$")
        NUC_pattern  = re.compile("^[acgtuACGTUnryNRY ]+$")   

        if header_id == None:  header_id = 'N/A'

        if seq_type.startswith('NUC'):
            if not NUC_pattern.match(sequence_str):
                self.log(console,"Not finding NUCLEOTIDE sequence for ID "+str(header_id)+" sequence: "+str(sequence_str))
                return False
        elif seq_type.startswith('PROT'): 
            if NUC_pattern.match(sequence_str):
                self.log(console,"Finding NUCLEOTIDE instead of PROTEIN sequence for ID "+str(header_id)+" sequence: "+str(sequence_str))
                return False
            elif not PROT_pattern.match(sequence_str):
                self.log(console,"Not finding PROTEIN sequence for ID "+str(header_id)+" sequence: "+str(sequence_str))
                return False

        return True


    #### Validate App input params
    ##
    def CheckBlastParams (self, params, app_name):

        # do some basic checks
        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
        if 'input_many_ref' not in params:
            raise ValueError('input_many_ref parameter is required')
        if 'output_filtered_name' not in params:
            raise ValueError('output_filtered_name parameter is required')

        # check query
        if app_name == 'psiBLAST':
            if 'input_msa_ref' not in params:
                raise ValueError('input_msa_ref parameter is required')

        else:  # need to store textarea query
            if ('output_one_name' not in params or params['output_one_name'] == None) \
                and ('input_one_sequence' in params and params['input_one_sequence'] != None):
                raise ValueError('output_one_name parameter required if input_one_sequence parameter is provided')
            if ('output_one_name' in params and params['output_one_name'] != None) \
                and ('input_one_sequence' not in params or params['input_one_sequence'] == None):
                raise ValueError('input_one_sequence parameter required if output_one_name parameter is provided')
            if ('input_one_ref' in params and params['input_one_ref'] != None) \
                and ('input_one_sequence' in params and params['input_one_sequence'] != None):
                raise ValueError('cannot have both input_one_sequence and input_one_ref parameter')
            if ('input_one_ref' in params and params['input_one_ref'] != None) \
                and ('output_one_name' in params and params['output_one_name'] != None):
                raise ValueError('cannot have both input_one_ref and output_one_name parameter')
            if ('input_one_ref' not in params or params['input_one_ref'] == None) \
                and ('input_one_sequence' not in params or params['input_one_sequence'] == None):
                raise ValueError('input_one_sequence or input_one_ref parameter is required')

        return True


    #### Store textarea input query sequence in object
    #
    def objectify_text_query (self, params, seq_type, search_tool_name):
        console = []
        invalid_msgs = []

        # Write the input_one_sequence to a SequenceSet object
        #
        if 'input_one_sequence' in params \
                and params['input_one_sequence'] != None \
                and params['input_one_sequence'] != "Optionally enter DNA sequence...":

            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=self.ctx['token'])
            ParseFastaStr_retVal = DOTFU.ParseFastaStr ({
                'fasta_str':    params['input_one_sequence'],
                'residue_type': seq_type,
                'case':         'UPPER',
                'console':      console,
                'invalid_msgs': invalid_msgs
                })
            header_id        = ParseFastaStr_retVal['id']
            header_desc      = ParseFastaStr_retVal['desc']
            sequence_str_buf = ParseFastaStr_retVal['seq']


            # load the method provenance from the context object
            #
            self.log(console,"SETTING PROVENANCE")  # DEBUG
            provenance = [{}]
            if 'provenance' in self.ctx:
                provenance = self.ctx['provenance']
            # add additional info to provenance here, in this case the input data object reference
            provenance[0]['input_ws_objects'] = []
            provenance[0]['service'] = 'kb_blast'
            provenance[0]['method'] = search_tool_name
                
            # Upload query
            #
            self.log(console,"UPLOADING OUTPUT QUERY OBJECT")  # DEBUG

            output_one_sequenceSet = { 'sequence_set_id': header_id,  
                                       'description': header_desc,
                                       'sequences': [ { 'sequence_id': header_id,
                                                        'description': header_desc,
                                                        'sequence': sequence_str_buf
                                                        }
                                                      ] 
                                       }

            try:
                ws = workspaceService(self.workspaceURL, token=self.ctx['token'])
                new_obj_info = ws.save_objects({
                            'workspace': params['workspace_name'],
                            'objects':[{
                                    'type': 'KBaseSequences.SequenceSet',
                                    'data': output_one_sequenceSet,
                                    'name': params['output_one_name'],
                                    'meta': {},
                                    'provenance': provenance
                                    }]
                            })[0]
                output_one_ref = str(new_obj_info[6])+'/'+str(new_obj_info[0])+'/'+str(new_obj_info[4])
            except Exception as e:
                raise ValueError('Unable to store output_one_name SequenceSet object from workspace: ' + str(e))
                #to get the full stack trace: traceback.format_exc()

            self.log(console, 'done')

        return output_one_ref


    # Handle overloading (input_one can be SequenceSet or FeatureSet)
    #
    def write_query_obj_to_file (self, params, input_one_ref, seq_type):
        console = []
        invalid_msgs = []
        appropriate_sequence_found_in_one_input = False

        # determine query object type
        #
        try:
            ws = workspaceService(self.workspaceURL, token=self.ctx['token'])
            #objects = ws.get_objects([{'ref': input_one_ref}])
            objects = ws.get_objects2({'objects':[{'ref': input_one_ref}]})['data']
            input_one_data = objects[0]['data']
            input_one_name = str(objects[0]['info'][1])
            info = objects[0]['info']
                                                             
            one_type_name = info[2].split('.')[1].split('-')[0]
        except Exception as e:
            raise ValueError('Unable to fetch input_one_ref object from workspace: ' + str(e))
        #to get the full stack trace: traceback.format_exc()


        # SequenceSet
        #
        if one_type_name == 'SequenceSet':
            if 'input_one_sequence' in params \
                and params['input_one_sequence'] != None \
                and params['input_one_sequence'] != "Optionally enter DNA sequence..." \
                and one_type_name != 'SequenceSet':

                self.log(invalid_msgs,"ERROR: Mismatched input type for Query Object: "+input_one_ref+" should be SequenceSet instead of: "+one_type_name)

            try:
                input_one_sequenceSet = input_one_data
            except Exception as e:
                print((traceback.format_exc()))
                raise ValueError('Unable to get sequenceSet object: ' + str(e))

            header_id = input_one_sequenceSet['sequences'][0]['sequence_id']
            sequence_str = input_one_data['sequences'][0]['sequence']

            if not self.validateSeq (seq_type, sequence_str, header_id):
                raise ValueError ("BAD record for sequence_id: "+header_id+"\n"+sequence_str+"\n")
            else:
                appropriate_sequence_found_in_one_input = True

            one_forward_reads_file_path = os.path.join(self.scratch, header_id+'.fasta')
            one_forward_reads_file_handle = open(one_forward_reads_file_path, 'w')
            self.log(console, 'writing reads file: '+str(one_forward_reads_file_path))
            one_forward_reads_file_handle.write('>'+header_id+"\n")
            one_forward_reads_file_handle.write(sequence_str+"\n")
            one_forward_reads_file_handle.close()

            self.log(console, 'done')

        # FeatureSet
        #
        elif one_type_name == 'FeatureSet':
            # retrieve sequences for features
            #input_one_featureSet = input_one_data
            one_forward_reads_file_dir = self.scratch
            one_forward_reads_file = input_one_name+".fasta"

            # DEBUG
            #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            FeatureSetToFASTA_params = {
                'featureSet_ref':      input_one_ref,
                'file':                one_forward_reads_file,
                'dir':                 one_forward_reads_file_dir,
                'console':             console,
                'invalid_msgs':        invalid_msgs,
                'residue_type':        seq_type,
                'feature_type':        'ALL',
                'record_id_pattern':   '%%genome_ref%%'+genome_id_feature_id_delim+'%%feature_id%%',
                'record_desc_pattern': '[%%genome_ref%%]',
                'case':                'upper',
                'linewrap':            50,
                'merge_fasta_files':   'TRUE'
                }

            #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=self.ctx['token'])
            FeatureSetToFASTA_retVal = DOTFU.FeatureSetToFASTA (FeatureSetToFASTA_params)
            one_forward_reads_file_path = FeatureSetToFASTA_retVal['fasta_file_path']
            if len(list(FeatureSetToFASTA_retVal['feature_ids_by_genome_ref'].keys())) > 0:
                appropriate_sequence_found_in_one_input = True


            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "FeatureSetToFasta() took "+str(end_time-beg_time)+" secs")


        return (one_forward_reads_file_path, appropriate_sequence_found_in_one_input)
