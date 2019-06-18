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

        try:
            self.wsClient = workspaceService(self.workspaceURL, token=self.ctx['token'])
        except:
            raise ValueError ("Failed to connect to workspace service")

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
                new_obj_info = self.wsClient.save_objects({
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
        feature_ids_by_genome_ref = None

        # determine query object type
        #
        try:
            ws = workspaceService(self.workspaceURL, token=self.ctx['token'])
            #objects = ws.get_objects([{'ref': input_one_ref}])
            objects = self.wsClient.get_objects2({'objects':[{'ref': input_one_ref}]})['data']
            input_one_data = objects[0]['data']
            input_one_name = str(objects[0]['info'][1])
            info = objects[0]['info']
                                                             
            query_type_name = info[2].split('.')[1].split('-')[0]
        except Exception as e:
            raise ValueError('Unable to fetch input_one_ref object from workspace: ' + str(e))
        #to get the full stack trace: traceback.format_exc()

        # DEBUG
        self.log (console,"QUERY_TYPE: '"+query_type_name+"'")


        # SequenceSet
        #
        if query_type_name == 'SequenceSet':
            if 'input_one_sequence' in params \
                and params['input_one_sequence'] != None \
                and params['input_one_sequence'] != "Optionally enter DNA sequence..." \
                and query_type_name != 'SequenceSet':

                self.log(invalid_msgs,"ERROR: Mismatched input type for Query Object: "+input_one_ref+" should be SequenceSet instead of: "+query_type_name)

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

            query_fasta_file_path = os.path.join(self.scratch, header_id+'.fasta')
            query_fasta_file_handle = open(query_fasta_file_path, 'w')
            self.log(console, 'writing reads file: '+str(query_fasta_file_path))
            query_fasta_file_handle.write('>'+header_id+"\n")
            query_fasta_file_handle.write(sequence_str+"\n")
            query_fasta_file_handle.close()

            self.log(console, 'done')

        # FeatureSet
        #
        elif query_type_name == 'FeatureSet':
            # retrieve sequences for features
            #input_one_featureSet = input_one_data
            query_fasta_file_dir = self.scratch
            query_fasta_file = input_one_name+".fasta"

            # DEBUG
            #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            FeatureSetToFASTA_params = {
                'featureSet_ref':      input_one_ref,
                'file':                query_fasta_file,
                'dir':                 query_fasta_file_dir,
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
            query_fasta_file_path = FeatureSetToFASTA_retVal['fasta_file_path']
            if len(list(FeatureSetToFASTA_retVal['feature_ids_by_genome_ref'].keys())) > 0:
                appropriate_sequence_found_in_one_input = True


            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "FeatureSetToFasta() took "+str(end_time-beg_time)+" secs")


        return ({ 'query_type_name': query_type_name,
                  'query_fasta_file_path': query_fasta_file_path,
                  'appropriate_sequence_found_in_one_input': appropriate_sequence_found_in_one_input,
                  'invalid_msgs': invalid_msgs
              })
                 

    #### Get the input_many object
    ##
    def write_target_obj_to_file (self, params, input_many_ref, seq_type):
        console = []
        invalid_msgs = []
        appropriate_sequence_found_in_many_input = False
        feature_ids = None
        feature_ids_by_genome_ref = None
        feature_ids_by_genome_id = None
        feature_id_to_function = None
        genome_ref_to_sci_name = None

        target_fasta_file_compression = None
        sequencing_tech = 'N/A'

        try:
            #objects = ws.get_objects([{'ref': input_many_ref}])
            objects = self.wsClient.get_objects2({'objects':[{'ref': input_many_ref}]})['data']
            input_many_data = objects[0]['data']
            info = objects[0]['info']
            input_many_name = str(info[1])
            target_type_name = info[2].split('.')[1].split('-')[0]

            if target_type_name == 'SingleEndLibrary':
                target_type_namespace = info[2].split('.')[0]
                if target_type_namespace == 'KBaseAssembly':
                    file_name = input_many_data['handle']['file_name']
                elif target_type_namespace == 'KBaseFile':
                    file_name = input_many_data['lib']['file']['file_name']
                else:
                    raise ValueError('bad data type namespace: '+target_type_namespace)
                #self.log(console, 'INPUT_MANY_FILENAME: '+file_name)  # DEBUG
                if file_name[-3:] == ".gz":
                    target_fasta_file_compression = 'gz'
                if 'sequencing_tech' in input_many_data:
                    sequencing_tech = input_many_data['sequencing_tech']

        except Exception as e:
            raise ValueError('Unable to fetch input_many_name object from workspace: ' + str(e))
            #to get the full stack trace: traceback.format_exc()


        # Handle overloading (input_many can be FeatureSet, Genome, or GenomeSet)
        #
        """
        if target_type_name == 'SequenceSet':
            try:
                input_many_sequenceSet = input_many_data
            except Exception as e:
                print((traceback.format_exc()))
                raise ValueError('Unable to get SequenceSet: ' + str(e))

            header_id = input_many_sequenceSet['sequences'][0]['sequence_id']
            target_fasta_file_path = os.path.join(self.scratch, header_id+'.fasta')
            target_fasta_file_handle = open(target_fasta_file_path, 'w')
            self.log(console, 'writing reads file: '+str(target_fasta_file_path))

            for seq_obj in input_many_sequenceSet['sequences']:
                header_id = seq_obj['sequence_id']
                sequence_str = seq_obj['sequence']

                if not self.validateSeq (seq_type, sequence_str, header_id):
                    raise ValueError ("BAD record for sequence_id: "+header_id+"\n"+sequence_str+"\n")
                else:
                    appropriate_sequence_found_in_many_input = True

                target_fasta_file_handle.write('>'+header_id+"\n")
                target_fasta_file_handle.write(sequence_str+"\n")
            target_fasta_file_handle.close();
            self.log(console, 'done')

        # SingleEndLibrary
        #
        elif target_type_name == 'SingleEndLibrary':

            # DEBUG
            #for k in data:
            #    self.log(console,"SingleEndLibrary ["+k+"]: "+str(data[k]))

            try:
                if 'lib' in input_many_data:
                    target_fasta = input_many_data['lib']['file']
                elif 'handle' in input_many_data:
                    target_fasta = input_many_data['handle']
                else:
                    self.log(console,"bad structure for 'target_fasta'")
                    raise ValueError("bad structure for 'target_fasta'")
                #if 'lib2' in data:
                #    reverse_reads = data['lib2']['file']
                #elif 'handle_2' in data:
                #    reverse_reads = data['handle_2']
                #else:
                #    reverse_reads={}

                ### NOTE: this section is what could be replaced by the transform services
                target_fasta_file_path = os.path.join(self.scratch,target_fasta['file_name'])
                target_fasta_file_handle = open(target_fasta_file_path, 'w')
                self.log(console, 'downloading reads file: '+str(target_fasta_file_path))
                headers = {'Authorization': 'OAuth '+self.ctx['token']}
                r = requests.get(target_fasta['url']+'/node/'+target_fasta['id']+'?download', stream=True, headers=headers)
                for chunk in r.iter_content(1024):
                    appropriate_sequence_found_in_many_input = True
                    target_fasta_file_handle.write(chunk)
                target_fasta_file_handle.close();
                self.log(console, 'done')
                ### END NOTE


                # remove carriage returns
                new_file_path = target_fasta_file_path+"-CRfree"
                new_file_handle = open(new_file_path, 'w')
                target_fasta_file_handle = open(target_fasta_file_path, 'r')
                for line in target_fasta_file_handle:
                    line = re.sub("\r","",line)
                    new_file_handle.write(line)
                target_fasta_file_handle.close();
                new_file_handle.close()
                target_fasta_file_path = new_file_path


                # convert FASTQ to FASTA (if necessary)
                new_file_path = target_fasta_file_path+".fna"
                new_file_handle = open(new_file_path, 'w')
                if target_fasta_file_compression == 'gz':
                    target_fasta_file_handle = gzip.open(target_fasta_file_path, 'r')
                else:
                    target_fasta_file_handle = open(target_fasta_file_path, 'r')
                header = None
                last_header = None
                last_seq_buf = None
                last_line_was_header = False
                was_fastq = False
                for line in target_fasta_file_handle:
                    if line.startswith('>'):
                        break
                    elif line.startswith('@'):
                        was_fastq = True
                        header = line[1:]
                        if last_header != None:
                            new_file_handle.write('>'+last_header)
                            new_file_handle.write(last_seq_buf)
                        last_seq_buf = None
                        last_header = header
                        last_line_was_header = True
                    elif last_line_was_header:
                        last_seq_buf = line
                        last_line_was_header = False
                    else:
                        continue
                if last_header != None:
                    new_file_handle.write('>'+last_header)
                    new_file_handle.write(last_seq_buf)

                new_file_handle.close()
                target_fasta_file_handle.close()
                if was_fastq:
                    target_fasta_file_path = new_file_path

            except Exception as e:
                print((traceback.format_exc()))
                raise ValueError('Unable to download single-end read library files: ' + str(e))
        """

        # FeatureSet
        #
        #elif target_type_name == 'FeatureSet':
        if target_type_name == 'FeatureSet':
            # retrieve sequences for features
            input_many_featureSet = input_many_data
            target_fasta_file_dir = self.scratch
            target_fasta_file = input_many_name+".fasta"

            # DEBUG
            #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            FeatureSetToFASTA_params = {
                'featureSet_ref':      input_many_ref,
                'file':                target_fasta_file,
                'dir':                 target_fasta_file_dir,
                'console':             console,
                'invalid_msgs':        invalid_msgs,
                'residue_type':        'nucleotide',
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
            target_fasta_file_path = FeatureSetToFASTA_retVal['fasta_file_path']
            feature_ids_by_genome_ref = FeatureSetToFASTA_retVal['feature_ids_by_genome_ref']
            if len(list(feature_ids_by_genome_ref.keys())) > 0:
                appropriate_sequence_found_in_many_input = True
            feature_id_to_function = FeatureSetToFASTA_retVal['feature_id_to_function']
            genome_ref_to_sci_name = FeatureSetToFASTA_retVal['genome_ref_to_sci_name']

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "FeatureSetToFasta() took "+str(end_time-beg_time)+" secs")


        # Genome
        #
        elif target_type_name == 'Genome':
            target_fasta_file_dir = self.scratch
            target_fasta_file = input_many_name+".fasta"

            # DEBUG
            #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            GenomeToFASTA_params = {
                'genome_ref':          input_many_ref,
                'file':                target_fasta_file,
                'dir':                 target_fasta_file_dir,
                'console':             console,
                'invalid_msgs':        invalid_msgs,
                'residue_type':        'nucleotide',
                'feature_type':        'ALL',
                'record_id_pattern':   '%%feature_id%%',
                'record_desc_pattern': '[%%genome_id%%]',
                'case':                'upper',
                'linewrap':            50
                }

            #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=self.ctx['token'])
            GenomeToFASTA_retVal = DOTFU.GenomeToFASTA (GenomeToFASTA_params)
            target_fasta_file_path = GenomeToFASTA_retVal['fasta_file_path']
            feature_ids = GenomeToFASTA_retVal['feature_ids']
            if len(feature_ids) > 0:
                appropriate_sequence_found_in_many_input = True
            feature_id_to_function = GenomeToFASTA_retVal['feature_id_to_function']
            genome_ref_to_sci_name = GenomeToFASTA_retVal['genome_ref_to_sci_name']

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "Genome2Fasta() took "+str(end_time-beg_time)+" secs")


        # GenomeSet
        #
        elif target_type_name == 'GenomeSet':
            input_many_genomeSet = input_many_data
            target_fasta_file_dir = self.scratch
            target_fasta_file = input_many_name+".fasta"

            # DEBUG
            #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            GenomeSetToFASTA_params = {
                'genomeSet_ref':       input_many_ref,
                'file':                target_fasta_file,
                'dir':                 target_fasta_file_dir,
                'console':             console,
                'invalid_msgs':        invalid_msgs,
                'residue_type':        'nucleotide',
                'feature_type':        'ALL',
                'record_id_pattern':   '%%genome_ref%%'+genome_id_feature_id_delim+'%%feature_id%%',
                'record_desc_pattern': '[%%genome_ref%%]',
                'case':                'upper',
                'linewrap':            50,
                'merge_fasta_files':   'TRUE'
                }

            #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=self.ctx['token'])
            GenomeSetToFASTA_retVal = DOTFU.GenomeSetToFASTA (GenomeSetToFASTA_params)
            target_fasta_file_path = GenomeSetToFASTA_retVal['fasta_file_path_list'][0]
            feature_ids_by_genome_id = GenomeSetToFASTA_retVal['feature_ids_by_genome_id']
            if len(list(feature_ids_by_genome_id.keys())) > 0:
                appropriate_sequence_found_in_many_input = True
            feature_id_to_function = GenomeToFASTA_retVal['feature_id_to_function']
            genome_ref_to_sci_name = GenomeToFASTA_retVal['genome_ref_to_sci_name']

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "FeatureSetToFasta() took "+str(end_time-beg_time)+" secs")


        # Missing proper input_target_type
        #
        else:
            raise ValueError('Cannot yet handle input_many type of: '+target_type_name)


        return ({ 'target_type_name': target_type_name,
                  'target_fasta_file_path': target_fasta_file_path,
                  'appropriate_sequence_found_in_many_input': appropriate_sequence_found_in_many_input,
                  'invalid_msgs': invalid_msgs,
                  'feature_ids': feature_ids,
                  'feature_ids_by_genome_ref': feature_ids_by_genome_ref,
                  'feature_ids_by_genome_id': feature_ids_by_genome_id,
                  'feature_id_to_function': feature_id_to_function,
                  'genome_ref_to_sci_name': genome_ref_to_sci_name
              })

