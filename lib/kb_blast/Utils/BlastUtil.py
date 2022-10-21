# -*- coding: utf-8 -*-
import os
import re
import subprocess
import sys
import traceback
import uuid
from datetime import datetime
from pprint import pformat

# SDK Utils
from installed_clients.KBaseDataObjectToFileUtilsClient import KBaseDataObjectToFileUtils
from installed_clients.kb_SetUtilitiesClient import kb_SetUtilities
from installed_clients.DataFileUtilClient import DataFileUtil as DFUClient
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.WorkspaceClient import Workspace as workspaceService


###############################################################################
# BlastUtil: methods to support Apps in kb_blast KBase module
###############################################################################

class BlastUtil:

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "1.8.0"
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

#        self.callbackURL = os.environ['SDK_CALLBACK_URL'] if os.environ['SDK_CALLBACK_URL'] is not None else 'https://kbase.us/services/njs_wrapper'
        self.callbackURL = os.environ.get('SDK_CALLBACK_URL')
        if self.callbackURL is None:
            raise ValueError ("SDK_CALLBACK_URL not set in environment")

        self.scratch = os.path.abspath(config['scratch'])
        if self.scratch is None:
            self.scratch = os.path.join('/kb','module','local_scratch')
        if not os.path.exists(self.scratch):
            os.makedirs(self.scratch)

        try:
            self.wsClient = workspaceService(self.workspaceURL, token=self.ctx['token'])
        except:
            raise ValueError ("Failed to connect to workspace service")
        try:
            REPORT_SERVICE_VER = 'release'
            self.reportClient = KBaseReport(self.callbackURL, token=self.ctx['token'], service_ver=REPORT_SERVICE_VER)
        except:
            raise ValueError ("Failed to instantiate KBaseReport client")
        try:
            SU_SERVICE_VER = 'release'
            self.set_util_Client = kb_SetUtilities(url=self.callbackURL, token=self.ctx['token'], service_ver=SU_SERVICE_VER)  # SDK Local
            #self.set_util_Client = kb_SetUtilities(url=self.serviceWizardURL, token=self.ctx['token'], service_ver=SU_SERVICE_VER)  # Service
        except:
            raise ValueError("unable to instantiate kb_SetUtilities client")
        try:
            DOTFU_SERVICE_VER = 'release'
            #DOTFU_SERVICE_VER = 'beta'  # DEBUG
            self.DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=self.ctx['token'], service_ver=DOTFU_SERVICE_VER)
        except:
            raise ValueError ("Failed to instantiate DataObjectToFileUtils client")

        self.genome_id_feature_id_delim = '.f:'


        #END_CONSTRUCTOR
        pass


    # _instantiate_provenance()
    #
    def _instantiate_provenance(self, 
                                method_name=None,
                                input_obj_refs=None):
        service = 'kb_blast'

        provenance = [{}]
        if 'provenance' in self.ctx:
            provenance = self.ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects'] = []
        if input_obj_refs:
            for input_ref in input_obj_refs:
                if '/' in input_ref:
                    provenance[0]['input_ws_objects'].append(input_ref)
        provenance[0]['service'] = service
        provenance[0]['method'] = method_name

        return provenance


    #### Sequence Validation
    ##
    def validateSeq (self, seq_type, sequence_str, header_id):
        console = []
        PROT_pattern = re.compile("^[acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWYxX ]+$")
        NUC_pattern  = re.compile("^[acgtuACGTUnryNRY ]+$")   

        if header_id is None:  header_id = 'N/A'

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
    def validate_BLAST_app_params (self, params, search_tool_name):

        # do some basic checks
        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
        if 'input_many_refs' not in params and params.get('gtdb_targets','none') == 'none':
            raise ValueError('either Search Targets or GTDB Targets are required')
        if 'output_filtered_name' not in params:
            raise ValueError('output_filtered_name parameter is required')
        if 'genome_disp_name_config' not in params:
            raise ValueError('genome_disp_name_config parameter is required')

        # check query
        if search_tool_name == 'psiBLAST':
            if 'input_msa_ref' not in params:
                raise ValueError('input_msa_ref parameter is required')

        else:  # need to store textarea query
            if params.get('output_one_name') is None and params.get('input_one_sequence') is not None:
                raise ValueError('output_one_name parameter required if input_one_sequence parameter is provided')
            if params.get('output_one_name') is not None and params.get('input_one_sequence') is None:
                raise ValueError('input_one_sequence parameter required if output_one_name parameter is provided')
            if params.get('input_one_ref') is not None and params.get('input_one_sequence') is not None:
                raise ValueError('cannot have both input_one_sequence and input_one_ref parameter')
            if params.get('input_one_ref') is not None and params.get('output_one_name') is not None:
                raise ValueError('cannot have both input_one_ref and output_one_name parameter')
            if params.get('input_one_ref') is None and params.get('input_one_sequence') is None:
                raise ValueError('input_one_sequence or input_one_ref parameter is required')

        return True


    #### Store textarea input query sequence in object
    #
    def objectify_text_query (self, params, seq_type, search_tool_name):
        console = []
        invalid_msgs = []
        method_name = search_tool_name+'_Search'

        # Write the input_one_sequence to a SequenceSet object
        #
        if params.get('input_one_sequence') is not None \
                and params['input_one_sequence'] != "Optionally enter DNA sequence...":

            ParseFastaStr_retVal = self.DOTFU.ParseFastaStr ({
                'fasta_str':    params['input_one_sequence'],
                'residue_type': seq_type,
                'case':         'UPPER',
                'console':      console,
                'invalid_msgs': invalid_msgs
                })
            header_id        = ParseFastaStr_retVal['id']
            header_desc      = ParseFastaStr_retVal['desc']
            sequence_str_buf = ParseFastaStr_retVal['seq']


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
                                    'provenance': self._instantiate_provenance (method_name=method_name)
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

        # defaults
        if not params.get('write_off_code_prot_seq'):
            params['write_off_code_prot_seq'] = 1
        params['write_off_code_prot_seq'] = int(params['write_off_code_prot_seq'])
        
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
            if params.get('input_one_sequence') is not None \
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
                'record_id_pattern':   '%%genome_ref%%'+self.genome_id_feature_id_delim+'%%feature_id%%',
                'record_desc_pattern': '[%%genome_ref%%]',
                'case':                'upper',
                'linewrap':            50,
                'id_len_limit':        49,
                'write_off_code_prot_seq': params['write_off_code_prot_seq'],
                'merge_fasta_files':   'TRUE'
                }

            FeatureSetToFASTA_retVal = self.DOTFU.FeatureSetToFASTA (FeatureSetToFASTA_params)
            query_fasta_file_path = FeatureSetToFASTA_retVal['fasta_file_path']
            #query_short_id_to_rec_id = FeatureSetToFASTA_retVal['short_id_to_rec_id']
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
        target_feature_info = { 'feature_ids': None,
                                'feature_ids_by_genome_ref': None,
                                'feature_ids_by_genome_id': None,
                                'feature_id_to_function': None,
                                'genome_ref_to_sci_name': None,
                                'genome_ref_to_obj_name': None,
                                'genome_id_to_genome_ref': None
        }
        target_fasta_file_compression = None
        sequencing_tech = 'N/A'

        # defaults
        if not params.get('write_off_code_prot_seq'):
            params['write_off_code_prot_seq'] = 1
        params['write_off_code_prot_seq'] = int(params['write_off_code_prot_seq'])
        
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


        # Handle overloading (input_many can be FeatureSet, Genome, GenomeSet, or SpeciesTree - not SequenceSet or SingleEndLibrary at this time)
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
                        if last_header is not None:
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
                if last_header is not None:
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
                'residue_type':        seq_type,
                'feature_type':        'ALL',
                'record_id_pattern':   '%%genome_ref%%'+self.genome_id_feature_id_delim+'%%feature_id%%',
                'record_desc_pattern': '[%%genome_ref%%]',
                'case':                'upper',
                'linewrap':            50,
                'id_len_limit':        49,
                'write_off_code_prot_seq': params['write_off_code_prot_seq'],
                'merge_fasta_files':   'TRUE'
                }

            FeatureSetToFASTA_retVal = self.DOTFU.FeatureSetToFASTA (FeatureSetToFASTA_params)
            target_fasta_file_path = FeatureSetToFASTA_retVal['fasta_file_path']
            target_feature_info['short_id_to_rec_id'] = FeatureSetToFASTA_retVal['short_id_to_rec_id']
            target_feature_info['feature_ids_by_genome_ref'] = FeatureSetToFASTA_retVal['feature_ids_by_genome_ref']
            if len(list(target_feature_info['feature_ids_by_genome_ref'].keys())) > 0:
                appropriate_sequence_found_in_many_input = True
            target_feature_info['feature_id_to_function'] = FeatureSetToFASTA_retVal['feature_id_to_function']
            target_feature_info['genome_ref_to_sci_name'] = FeatureSetToFASTA_retVal['genome_ref_to_sci_name']
            target_feature_info['genome_ref_to_obj_name'] = FeatureSetToFASTA_retVal['genome_ref_to_obj_name']

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
                'residue_type':        seq_type,
                'feature_type':        'ALL',
                'record_id_pattern':   '%%feature_id%%',
                'record_desc_pattern': '[%%genome_id%%]',
                'case':                'upper',
                'linewrap':            50,
                'id_len_limit':        49,
                'write_off_code_prot_seq': params['write_off_code_prot_seq']
                }

            GenomeToFASTA_retVal = self.DOTFU.GenomeToFASTA (GenomeToFASTA_params)
            target_fasta_file_path = GenomeToFASTA_retVal['fasta_file_path']
            target_feature_info['short_id_to_rec_id'] = GenomeToFASTA_retVal['short_id_to_rec_id']
            target_feature_info['feature_ids'] = GenomeToFASTA_retVal['feature_ids']
            if len(target_feature_info['feature_ids']) > 0:
                appropriate_sequence_found_in_many_input = True
            target_feature_info['feature_id_to_function'] = GenomeToFASTA_retVal['feature_id_to_function']
            target_feature_info['genome_ref_to_sci_name'] = GenomeToFASTA_retVal['genome_ref_to_sci_name']
            target_feature_info['genome_ref_to_obj_name'] = GenomeToFASTA_retVal['genome_ref_to_obj_name']
            

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
                'residue_type':        seq_type,
                'feature_type':        'ALL',
                'record_id_pattern':   '%%genome_ref%%'+self.genome_id_feature_id_delim+'%%feature_id%%',
                'record_desc_pattern': '[%%genome_ref%%]',
                'case':                'upper',
                'linewrap':            50,
                'id_len_limit':        49,
                'write_off_code_prot_seq': params['write_off_code_prot_seq'],
                'merge_fasta_files':   'TRUE'
                }

            GenomeSetToFASTA_retVal = self.DOTFU.GenomeSetToFASTA (GenomeSetToFASTA_params)
            target_fasta_file_path = GenomeSetToFASTA_retVal['fasta_file_path_list'][0]
            target_feature_info['short_id_to_rec_id'] = GenomeSetToFASTA_retVal['short_id_to_rec_id']
            target_feature_info['feature_ids_by_genome_id'] = GenomeSetToFASTA_retVal['feature_ids_by_genome_id']
            if len(list(target_feature_info['feature_ids_by_genome_id'].keys())) > 0:
                appropriate_sequence_found_in_many_input = True
            target_feature_info['feature_id_to_function'] = GenomeSetToFASTA_retVal['feature_id_to_function']
            target_feature_info['genome_ref_to_sci_name'] = GenomeSetToFASTA_retVal['genome_ref_to_sci_name']
            target_feature_info['genome_ref_to_obj_name'] = GenomeSetToFASTA_retVal['genome_ref_to_obj_name']

            target_feature_info['genome_id_to_genome_ref'] = dict()
            for genome_id in input_many_genomeSet['elements'].keys():
                genome_ref = input_many_genomeSet['elements'][genome_id]['ref']
                target_feature_info['genome_id_to_genome_ref'][genome_id] = genome_ref

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "GenomeSetToFasta() took "+str(end_time-beg_time)+" secs")


        # SpeciesTree
        #
        elif target_type_name == 'Tree':
            input_many_speciesTree = input_many_data
            target_fasta_file_dir = self.scratch
            target_fasta_file = input_many_name+".fasta"

            # DEBUG
            #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            SpeciesTreeToFASTA_params = {
                'tree_ref':            input_many_ref,
                'file':                target_fasta_file,
                'dir':                 target_fasta_file_dir,
                'console':             console,
                'invalid_msgs':        invalid_msgs,
                'residue_type':        seq_type,
                'feature_type':        'ALL',
                'record_id_pattern':   '%%genome_ref%%'+self.genome_id_feature_id_delim+'%%feature_id%%',
                'record_desc_pattern': '[%%genome_ref%%]',
                'case':                'upper',
                'linewrap':            50,
                'id_len_limit':        49,
                'write_off_code_prot_seq': params['write_off_code_prot_seq'],
                'merge_fasta_files':   'TRUE'
                }

            SpeciesTreeToFASTA_retVal = self.DOTFU.SpeciesTreeToFASTA (SpeciesTreeToFASTA_params)
            target_fasta_file_path = SpeciesTreeToFASTA_retVal['fasta_file_path_list'][0]
            target_feature_info['short_id_to_rec_id'] = SpeciesTreeToFASTA_retVal['short_id_to_rec_id']
            target_feature_info['feature_ids_by_genome_id'] = SpeciesTreeToFASTA_retVal['feature_ids_by_genome_id']
            if len(list(target_feature_info['feature_ids_by_genome_id'].keys())) > 0:
                appropriate_sequence_found_in_many_input = True
            target_feature_info['feature_id_to_function'] = SpeciesTreeToFASTA_retVal['feature_id_to_function']
            target_feature_info['genome_ref_to_sci_name'] = SpeciesTreeToFASTA_retVal['genome_ref_to_sci_name']
            target_feature_info['genome_ref_to_obj_name'] = SpeciesTreeToFASTA_retVal['genome_ref_to_obj_name']

            target_feature_info['genome_id_to_genome_ref'] = dict()
            for genome_id in input_many_speciesTree['ws_refs'].keys():
                genome_ref = input_many_speciesTree['ws_refs'][genome_id]['g'][0]
                target_feature_info['genome_id_to_genome_ref'][genome_id] = genome_ref

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "SpeciesTreeToFasta() took "+str(end_time-beg_time)+" secs")


        # AnnotatedMetagenomeAssembly
        #
        elif target_type_name == 'AnnotatedMetagenomeAssembly':
            target_fasta_file_dir = self.scratch
            target_fasta_file = input_many_name+".fasta"

            # DEBUG
            #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            AnnotatedMetagenomeAssemblyToFASTA_params = {
                'ama_ref':             input_many_ref,
                'file':                target_fasta_file,
                'dir':                 target_fasta_file_dir,
                'console':             console,
                'invalid_msgs':        invalid_msgs,
                'residue_type':        seq_type,
                'feature_type':        'ALL',
                'record_id_pattern':   '%%feature_id%%',
                'record_desc_pattern': '[%%genome_id%%]',
                'case':                'upper',
                'linewrap':            50,
                'write_off_code_prot_seq': params['write_off_code_prot_seq'],
                'id_len_limit':        49
                }

            AnnotatedMetagenomeAssemblyToFASTA_retVal = self.DOTFU.AnnotatedMetagenomeAssemblyToFASTA (AnnotatedMetagenomeAssemblyToFASTA_params)
            target_fasta_file_path = AnnotatedMetagenomeAssemblyToFASTA_retVal['fasta_file_path']
            target_feature_info['short_id_to_rec_id'] = AnnotatedMetagenomeAssemblyToFASTA_retVal['short_id_to_rec_id']
            target_feature_info['feature_ids'] = AnnotatedMetagenomeAssemblyToFASTA_retVal['feature_ids']
            if len(target_feature_info['feature_ids']) > 0:
                appropriate_sequence_found_in_many_input = True
            target_feature_info['feature_id_to_function'] = AnnotatedMetagenomeAssemblyToFASTA_retVal['feature_id_to_function']
            target_feature_info['ama_ref_to_obj_name'] = AnnotatedMetagenomeAssemblyToFASTA_retVal['ama_ref_to_obj_name']

            # DEBUG
            #with open (target_fasta_file_path, 'r') as fasta_handle:
            #    for fasta_line in fasta_handle.readlines():
            #        print ("FASTA_LINE: '"+fasta_line)
            

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "Genome2Fasta() took "+str(end_time-beg_time)+" secs")


        # Missing proper input_target_type
        #
        else:
            raise ValueError('Cannot yet handle input_many type of: '+target_type_name)


        return ({ 'target_name': input_many_name,
                  'target_type_name': target_type_name,
                  'target_fasta_file_path': target_fasta_file_path,
                  'appropriate_sequence_found_in_many_input': appropriate_sequence_found_in_many_input,
                  'invalid_msgs': invalid_msgs,
                  'target_feature_info': target_feature_info
              })


    # input data failed validation.  Need to return
    #
    def save_error_report_with_invalid_msgs (self, invalid_msgs, input_one_ref, input_many_refs, method_name):
        console = []

        # build output report object
        #
        self.log(console,"BUILDING REPORT")  # DEBUG
        report += "FAILURE:\n\n"+"\n".join(invalid_msgs)+"\n"
        reportObj = {
            'objects_created':[],
            'text_message':report
        }

        reportName = 'blast_report_'+str(uuid.uuid4())
        report_obj_info = self.wsClient.save_objects({
            #'id':info[6],
            'workspace':params['workspace_name'],
            'objects':[
                {
                    'type':'KBaseReport.Report',
                    'data':reportObj,
                    'name':reportName,
                    'meta':{},
                    'hidden':1,
                    'provenance': self._instantiate_provenance(method_name = method_name,
                                                               input_obj_refs = [input_one_ref]+input_many_refs)
                }
            ]
        })[0]

        error_report_info = { 'name': reportName,
                              'ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4])
                          }
        return error_report_info


    #### format_BLAST_db()
    ##
    def format_BLAST_db (self, search_tool_name, target_fasta_file_path):
        console = []
        BLAST_FORMAT_successful = True

        # set seq type
        if search_tool_name == 'BLASTn' or search_tool_name == 'tBLASTn':
            seq_type = 'nucl'
        else:
            seq_type = 'prot'
            

        # FORMAT DB
        #
        # OLD SYNTAX: formatdb -i $database -o T -p F -> $database.nsq or $database.00.nsq
        # NEW SYNTAX: makeblastdb -in $database -parse_seqids -dbtype prot/nucl -out <basename>
        makeblastdb_cmd = [self.Make_BLAST_DB]

        # check for necessary files
        if not os.path.isfile(self.Make_BLAST_DB):
            self.log(console,"no such file '"+self.Make_BLAST_DB+"'")
            BLAST_DB_FORMAT_successful = False
        if not os.path.isfile(target_fasta_file_path):
            self.log(console,"no such file '"+target_fasta_file_path+"'")
            BLAST_DB_FORMAT_successful = False
        elif not os.path.getsize(target_fasta_file_path) > 0:
            self.log(console,"empty file '"+target_fasta_file_path+"'")
            BLAST_DB_FORMAT_successful = False

        makeblastdb_cmd.append('-in')
        makeblastdb_cmd.append(target_fasta_file_path)
        makeblastdb_cmd.append('-parse_seqids')
        makeblastdb_cmd.append('-dbtype')
        makeblastdb_cmd.append(seq_type)
        makeblastdb_cmd.append('-out')
        makeblastdb_cmd.append(target_fasta_file_path)

        # Run Make_BLAST_DB, capture output as it happens
        #
        self.log(console, 'RUNNING Make_BLAST_DB:')
        self.log(console, '    '+' '.join(makeblastdb_cmd))
#        report += "\n"+'running Make_BLAST_DB:'+"\n"
#        report += '    '+' '.join(makeblastdb_cmd)+"\n"

        p = subprocess.Popen(makeblastdb_cmd, \
                             cwd = self.scratch, \
                             stdout = subprocess.PIPE, \
                             stderr = subprocess.STDOUT, \
                             shell = False)

        while True:
            line = p.stdout.readline().decode()
            if not line: break
            self.log(console, line.replace('\n', ''))

        p.stdout.close()
        p.wait()
        self.log(console, 'return code: ' + str(p.returncode))
        if p.returncode != 0:
            self.log(console,'Error running makeblastdb, return code: '+str(p.returncode) + '\n\n')
            #'\n\n'+ '\n'.join(console))
            BLAST_DB_FORMAT_successful = False

        # Check for db output
        if seq_type.lower().startswith('n'):
            db_ext = 'nsq'
        else:
            db_ext = 'psq'
        if not os.path.isfile(target_fasta_file_path+"."+db_ext) and not os.path.isfile(target_fasta_file_path+".00."+db_ext):
            self.log(console,"makeblastdb failed to create DB file '"+target_fasta_file_path+"."+db_ext+"'")
            BLAST_DB_FORMAT_successful = False
        elif not os.path.getsize(target_fasta_file_path+"."+db_ext) > 0 and not os.path.getsize(target_fasta_file_path+".00."+db_ext) > 0:
            self.log(console,"makeblastdb created empty DB file '"+target_fasta_file_path+"."+db_ext+"'")
            BLAST_DB_FORMAT_successful = False

        return BLAST_FORMAT_successful


    # _check_BLAST_input_ready()
    #
    def _check_BLAST_input_ready (self, blast_bin, query_fasta_file_path, target_fasta_file_path):
        console = []
        BLAST_ready = True

        #db_file = re.sub('.fasta$','',target_fasta_file_path)
        db_file = re.sub('.faa$','',target_fasta_file_path)
        db_file_prot_marker = db_file+'.pdb'
        db_file_nuc_marker = db_file+'.ndb'
        db_file_prot_00_marker = db_file+'.00.pdb'
        db_file_nuc_00_marker = db_file+'.00.ndb'
        
        # check for necessary files
        if not os.path.isfile(blast_bin):
            self.log(console, "no such file '"+blast_bin+"'")
            BLAST_ready = False
        if not os.path.isfile(query_fasta_file_path):
            self.log(console, "no such file '"+query_fasta_file_path+"'")
            BLAST_ready = False
        elif not os.path.getsize(query_fasta_file_path) > 0:
            self.log(console, "empty file '"+query_fasta_file_path+"'")
            BLAST_ready = False
        if not os.path.isfile(db_file_prot_marker) and \
           not os.path.isfile(db_file_prot_00_marker) and \
           not os.path.isfile(db_file_nuc_marker) and \
           not os.path.isfile(db_file_nuc_00_marker):
            self.log(console, "no such db '"+db_file+"'")
            BLAST_ready = False
        elif os.path.isfile(db_file_prot_marker) and \
             not os.path.getsize(db_file_prot_marker) > 0:
            self.log(console, "empty db '"+db_file+"'")
            BLAST_ready = False
        elif os.path.isfile(db_file_prot_00_marker) and \
             not os.path.getsize(db_file_prot_00_marker) > 0:
            self.log(console, "empty db '"+db_file+"'")
            BLAST_ready = False
        elif os.path.isfile(db_file_nuc_marker) and \
             not os.path.getsize(db_file_nuc_marker) > 0:
            self.log(console, "empty db '"+db_file+"'")
            BLAST_ready = False
        elif os.path.isfile(db_file_nuc_00_marker) and \
             not os.path.getsize(db_file_nuc_00_marker) > 0:
            self.log(console, "empty db '"+db_file+"'")
            BLAST_ready = False

        return BLAST_ready


    # _set_BLAST_bin()
    #
    def _set_BLAST_bin (self, search_tool_name):
        blast_bin = { 'BLASTn':   self.BLASTn,
                      'BLASTp':   self.BLASTp,
                      'BLASTx':   self.BLASTx,
                      'tBLASTn':  self.tBLASTn,
                      'tBLASTx':  self.tBLASTx,
                      'psiBLAST': self.psiBLAST
                  }

        return blast_bin[search_tool_name]

    # _set_BLAST_seq_types()
    #
    def _set_BLAST_seq_types (self, search_tool_name):
        query_seq_type = { 'BLASTn': 'NUC',
                           'BLASTp': 'PRO',
                           'BLASTx': 'NUC',
                           'tBLASTn': 'PRO',
                           'tBLASTx': 'NUC',
                           'psiBLAST': 'PRO'
                       }
        target_seq_type = { 'BLASTn': 'NUC',
                            'BLASTp': 'PRO',
                            'BLASTx': 'PRO',
                            'tBLASTn': 'NUC',
                            'tBLASTx': 'NUC',
                            'psiBLAST': 'PRO'
                        }
        return (query_seq_type[search_tool_name], target_seq_type[search_tool_name])
        
    # _set_BLAST_output_path()
    #
    def _set_BLAST_output_path (self, BLAST_output_format_str):
        timestamp = int((datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()*1000)
        output_dir = os.path.join(self.scratch,'output.'+str(timestamp))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        return os.path.join(output_dir, 'alnout_m='+BLAST_output_format_str+'.txt');
        #output_filtered_fasta_file_path = os.path.join(output_dir, 'output_filtered.fna');  # only for SingleEndLibrary


    # _set_HTML_outdir()
    #
    def _set_HTML_outdir (self):
        timestamp = int((datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()*1000)
        output_dir = os.path.join(self.scratch,'output.'+str(timestamp))
        html_output_dir = os.path.join(output_dir,'html')
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        if not os.path.exists(html_output_dir):
            os.makedirs(html_output_dir)

        return (html_output_dir)


    # _set_HTML_file_path()
    #
    def _set_HTML_file_path (self, html_output_dir, html_file):
        html_file_path = os.path.join (html_output_dir, html_file)
        return (html_file_path)


    # _build_BLAST_cmd()
    #
    def _build_BLAST_cmd (self, 
                          search_tool_name=None, 
                          query_fasta_file_path=None,
                          target_fasta_file_path=None,
                          output_aln_file_path=None,
                          BLAST_output_format_str=None,
                          e_value=None,
                          maxaccepts=None):

        # set BLAST bin
        BLAST_bin = self._set_BLAST_bin (search_tool_name)

        # check if ready
        if not self._check_BLAST_input_ready (BLAST_bin, query_fasta_file_path, target_fasta_file_path):
            raise ValueError ("Not ready to run BLAST")


        # OLD SYNTAX: $blast -q $q -G $G -E $E -m $m -e $e_value -v $limit -b $limit -K $limit -p blastn -i $fasta_file -d $database -o $out_file
        # NEW SYNTAX: blastn -query <queryfile> -db <basename> -out <out_aln_file> -outfmt 0/7 (8 became 7) -evalue <e_value> -dust no (DNA) -seg no (AA) -num_threads <num_cores>
        blast_cmd = [BLAST_bin]
        blast_cmd.append('-query')
        blast_cmd.append(query_fasta_file_path)
        blast_cmd.append('-db')
        blast_cmd.append(target_fasta_file_path)
        blast_cmd.append('-out')
        blast_cmd.append(output_aln_file_path)
        #blast_cmd.append('-html')  # HTML is a flag so doesn't get an arg val
        blast_cmd.append('-outfmt')
        blast_cmd.append(BLAST_output_format_str)
        blast_cmd.append('-evalue')
        blast_cmd.append(str(e_value))
        # options (not allowed for format 0)
        if BLAST_output_format_str != '0' and maxaccepts is not None:
            blast_cmd.append('-max_target_seqs')
            blast_cmd.append(str(maxaccepts))

        return blast_cmd


    # _exec_BLAST()
    #
    def _exec_BLAST (self, BLAST_cmd):
        console = []

        # Run BLAST, capture output as it happens
        self.log(console, 'RUNNING BLAST:')
        self.log(console, ' '.join(BLAST_cmd))
        self.log(console, '--------------------------------------')

        p = subprocess.Popen(BLAST_cmd, \
                             cwd = self.scratch, \
                             stdout = subprocess.PIPE, \
                             stderr = subprocess.STDOUT, \
                             shell = False)

        while True:
            line = p.stdout.readline().decode()
            if not line: break
            self.log(console, line.replace('\n', ''))

        p.stdout.close()
        p.wait()
        self.log(console, 'return code: ' + str(p.returncode))
        if p.returncode != 0:
            return 'Error running BLAST, return code: '+str(p.returncode) + '\n\n'+ '\n'.join(console)

        return 'Success'


    #### run_BLAST(): actual invocation
    ##
    def run_BLAST (self, 
                   search_tool_name = None,
                   query_fasta_file_path = None,
                   target_fasta_file_path = None,
                   e_value = None, 
                   maxaccepts = None,
                   BLAST_output_format_str = None):
        console = []
        
        # set the output path
        output_aln_file_path = self._set_BLAST_output_path (BLAST_output_format_str)


        # construct the BLAST command
        BLAST_cmd = self._build_BLAST_cmd (search_tool_name = search_tool_name,
                                           query_fasta_file_path = query_fasta_file_path,
                                           target_fasta_file_path = target_fasta_file_path,
                                           output_aln_file_path = output_aln_file_path,
                                           BLAST_output_format_str = BLAST_output_format_str,
                                           e_value = e_value,
                                           maxaccepts = maxaccepts)


        # execute BLAST
        BLAST_exec_return_msg = self._exec_BLAST (BLAST_cmd)
        if BLAST_exec_return_msg != 'Success':
            self.log(console, BLAST_exec_return_msg)
            raise ValueError ("FAILURE executing BLAST with command: \n\n"+"\n".join(BLAST_cmd))


        # upload BLAST output
        dfu = DFUClient(self.callbackURL)
        try:
            bulk_save_info = dfu.file_to_shock({'file_path': output_aln_file_path,
                                                 # DEBUG
                                                 # 'make_handle': 0,
                                                 # 'pack': 'zip'})
                                                 'make_handle': 0})
        except:
            raise ValueError ('error uploading '+output_aln_file_path+' file')


        # return info
        return {
            'output_aln_file_path': output_aln_file_path,
            'bulk_save_info': bulk_save_info
        }


    #### get_query_len()
    ##
    def get_query_len (self, query_fasta_file_path):
        query_len = 0
        with open(query_fasta_file_path, 'r') as query_file_handle:
            for line in query_file_handle:
                if line.startswith('>'):
                    continue
                query_len += len(re.sub(" ","", line.rstrip())) 

        return query_len


    #### parse_BLAST_tab_output()
    ##
    def parse_BLAST_tab_output (self, 
                                output_aln_file_path = None, 
                                search_tool_name = None,
                                params = None, 
                                query_len = None,
                                num_targets = 1,
                                target_ref = None,
                                target_name = None,
                                target_type_name = None,
                                target_feature_info = None):
        console = []
        invalid_msgs = []
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = list(range(11))  # object_info tuple
        method_name = search_tool_name+'_Search'
        (q_seq_type, t_seq_type) = self._set_BLAST_seq_types (search_tool_name)
        output_featureSet_ref = None
        
        # Parse the BLAST tabular output and store ids to filter many set to make filtered object to save back to KBase
        #
        self.log(console, 'PARSING BLAST ALIGNMENT OUTPUT')
        if not os.path.isfile(output_aln_file_path):
            raise ValueError("failed to create BLAST output: "+output_aln_file_path)
        elif not os.path.getsize(output_aln_file_path) > 0:
            raise ValueError("created empty file for BLAST output: "+output_aln_file_path)
        hit_seq_ids = dict()
        accept_fids = dict()
        output_aln_file_handle = open (output_aln_file_path, 'r')
        output_aln_buf = output_aln_file_handle.readlines()
        output_aln_file_handle.close()

        # DEBUG
        self.log(console, "BLAST_OUTPUT:")
        self.log(console, "".join(output_aln_buf))

        hit_total = 0
        high_bitscore_line = dict()
        high_bitscore_score = dict()
        high_bitscore_ident = dict()
        high_bitscore_alnlen = dict()
        hit_order = []
        hit_buf = []
        header_done = False
        for line in output_aln_buf:

            if line.startswith('#'):
                if not header_done:
                    hit_buf.append(line)
                continue
            header_done = True
            hit_info = line.split("\t")
            sh_hit_seq_id  = hit_info[1]
            if target_feature_info is not None:
                if sh_hit_seq_id in target_feature_info['short_id_to_rec_id']:
                    hit_seq_id = target_feature_info['short_id_to_rec_id'][sh_hit_seq_id]
                    hit_info[1] = hit_seq_id
                    line = "\t".join(hit_info)
                else:
                    hit_seq_id = sh_hit_seq_id
            else:
                hit_seq_id = sh_hit_seq_id
            hit_ident      = float(hit_info[2]) / 100.0
            hit_aln_len    = hit_info[3]
            hit_mismatches = hit_info[4]
            hit_gaps       = hit_info[5]
            hit_q_beg      = hit_info[6]
            hit_q_end      = hit_info[7]
            hit_t_beg      = hit_info[8]
            hit_t_end      = hit_info[9]
            hit_e_value    = hit_info[10]
            hit_bitscore   = hit_info[11]

            if q_seq_type != t_seq_type:
                hit_aln_len *= 3
            
            # BLAST SOMETIMES ADDS THIS TO IDs.  NO IDEA WHY, BUT GET RID OF IT!
            if hit_seq_id.startswith('gnl\|'):
                hit_seq_id = hit_seq_id[4:]

            try:
                if float(hit_bitscore) > float(high_bitscore_score[hit_seq_id]):
                    high_bitscore_score[hit_seq_id] = hit_bitscore
                    high_bitscore_ident[hit_seq_id] = hit_ident
                    high_bitscore_alnlen[hit_seq_id] = hit_aln_len
                    high_bitscore_line[hit_seq_id] = line
            except:
                hit_order.append(hit_seq_id)
                high_bitscore_score[hit_seq_id] = hit_bitscore
                high_bitscore_ident[hit_seq_id] = hit_ident
                high_bitscore_alnlen[hit_seq_id] = hit_aln_len
                high_bitscore_line[hit_seq_id] = line

        filtering_fields = dict()
        for hit_seq_id in hit_order:
            hit_buf.append(high_bitscore_line[hit_seq_id])
            filtering_fields[hit_seq_id] = dict()

            filter = False
            if 'ident_thresh' in params and float(params['ident_thresh']) > 100*float(high_bitscore_ident[hit_seq_id]):
                filter = True
                filtering_fields[hit_seq_id]['ident_thresh'] = True
                self.log (console, "FILTERING "+hit_seq_id+" on IDENT")  # DEBUG
            if 'bitscore' in params and float(params['bitscore']) > float(high_bitscore_score[hit_seq_id]):
                filter = True
                filtering_fields[hit_seq_id]['bitscore'] = True
                self.log (console, "FILTERING "+hit_seq_id+" on BITSCORE")  # DEBUG
            if 'overlap_fraction' in params and float(params['overlap_fraction']) > 100*float(high_bitscore_alnlen[hit_seq_id])/float(query_len):
                filter = True
                filtering_fields[hit_seq_id]['overlap_fraction'] = True
                self.log (console, "FILTERING "+hit_seq_id+" on OVERLAP")  # DEBUG

            if filter:
                continue
            
            hit_total += 1
            hit_seq_ids[hit_seq_id] = True
            self.log(console, "HIT: '"+hit_seq_id+"'")  # DEBUG
        

        self.log(console, 'EXTRACTING HITS FROM INPUT')
        #self.log(console, 'MANY_TYPE_NAME: '+many_type_name)  # DEBUG


        # Not handling SequenceSet nor SingleEndLibrary at this time
        #
        """
        # SequenceSet input -> SequenceSet output
        #
        if target_type_name == 'SequenceSet':
            seq_total = len(input_many_sequenceSet['sequences'])

            output_sequenceSet = dict()

            if 'sequence_set_id' in input_many_sequenceSet and input_many_sequenceSet['sequence_set_id'] is not None:
                output_sequenceSet['sequence_set_id'] = input_many_sequenceSet['sequence_set_id'] + "."+search_tool_name+"_Search_filtered"
            else:
                output_sequenceSet['sequence_set_id'] = search_tool_name+"_Search_filtered"
            if 'description' in input_many_sequenceSet and input_many_sequenceSet['description'] is not None:
                output_sequenceSet['description'] = input_many_sequenceSet['description'] + " - "+search_tool_name+"_Search filtered"
            else:
                output_sequenceSet['description'] = search_tool_name+"_Search filtered"

            self.log(console,"ADDING SEQUENCES TO SEQUENCESET")
            output_sequenceSet['sequences'] = []

            for seq_obj in input_many_sequenceSet['sequences']:
                header_id = seq_obj['sequence_id']
                #header_desc = seq_obj['description']
                #sequence_str = seq_obj['sequence']

                id_untrans = header_id
                id_trans = re.sub ('\|',':',id_untrans)  # BLAST seems to make this translation now when id format has simple 'kb|blah' format
                if id_trans in hit_seq_ids or id_untrans in hit_seq_ids:
                    #self.log(console, 'FOUND HIT '+header_id)  # DEBUG
                    accept_fids[id_untrans] = True
                    output_sequenceSet['sequences'].append(seq_obj)

        # SingleEndLibrary input -> SingleEndLibrary output
        #
        elif target_type_name == 'SingleEndLibrary':

            #  Note: don't use SeqIO.parse because loads everything into memory
            #
#            with open(target_fasta_file_path, 'r', -1) as target_fasta_file_handle, open(output_filtered_fasta_file_path, 'w', -1) as output_filtered_fasta_file_handle:
            output_filtered_fasta_file_handle = open(output_filtered_fasta_file_path, 'w', -1)
            if target_fasta_file_compression == 'gz':
                target_fasta_file_handle = gzip.open(target_fasta_file_path, 'r', -1)
            else:
                target_fasta_file_handle = open(target_fasta_file_path, 'r', -1)

            seq_total = 0;
            filtered_seq_total = 0
            last_seq_buf = []
            last_seq_id = None
            last_header = None
            pattern = re.compile('^\S*')
            for line in target_fasta_file_handle:
                if line.startswith('>'):
                    #self.log(console, 'LINE: '+line)  # DEBUG
                    seq_total += 1
                    seq_id = line[1:]  # removes '>'
                    seq_id = pattern.findall(seq_id)[0]

                    if last_seq_id is not None:
                        #self.log(console, 'ID: '+last_seq_id)  # DEBUG
                        id_untrans = last_seq_id
                        id_trans = re.sub ('\|',':',id_untrans)  # BLAST seems to make this translation now when id format has simple 'kb|blah' format
                        if id_trans in hit_seq_ids or id_untrans in hit_seq_ids:
                            #self.log(console, 'FOUND HIT '+header_id)  # DEBUG
                            accept_fids[id_untrans] = True
                            filtered_seq_total += 1
                            output_filtered_fasta_file_handle.write(last_header)
                            output_filtered_fasta_file_handle.writelines(last_seq_buf)
                    last_seq_buf = []
                    last_seq_id = seq_id
                    last_header = line
                else:
                    last_seq_buf.append(line)

            if last_seq_id is not None:
                #self.log(console, 'ID: '+last_seq_id)  # DEBUG
                id_untrans = last_seq_id
                id_trans = re.sub ('\|',':',id_untrans)  # BLAST seems to make this translation now when id format has simple 'kb|blah' format
                if id_trans in hit_seq_ids or id_untrans in hit_seq_ids:
                    #self.log(console, 'FOUND HIT '+header_id)  # DEBUG
                    accept_fids[id_untrans] = True
                    filtered_seq_total += 1
                    output_filtered_fasta_file_handle.write(last_header)
                    output_filtered_fasta_file_handle.writelines(last_seq_buf)
            last_seq_buf = []
            last_seq_id = None
            last_header = None

            target_fasta_file_handle.close()
            output_filtered_fasta_file_handle.close()

            if filtered_seq_total != hit_total:
                self.log(console,'hits in BLAST alignment output '+str(hit_total)+' != '+str(filtered_seq_total)+' matched sequences in input file')
                raise ValueError('hits in BLAST alignment output '+str(hit_total)+' != '+str(filtered_seq_total)+' matched sequences in input file')
        """


        # FeatureSet input -> FeatureSet output
        #
        #elif target_type_name == 'FeatureSet':
        if target_type_name == 'FeatureSet':
            #seq_total = len(list(input_many_featureSet['elements'].keys()))
            seq_total = 0
            for genome_ref in target_feature_info['feature_ids_by_genome_ref'].keys():
                for fId in target_feature_info['feature_ids_by_genome_ref'][genome_ref]:
                    seq_total += 1

            output_featureSet = dict()
            #if 'description' in input_many_featureSet and input_many_featureSet['description'] is not None:
            #    output_featureSet['description'] = input_many_featureSet['description'] + " - "+search_tool_name+"_Search filtered"
            #else:
            #    output_featureSet['description'] = search_tool_name+"_Search filtered"
            output_featureSet['description'] = search_tool_name+"_Search filtered"
            output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()

            #fId_list = list(input_many_featureSet['elements'].keys())
            self.log(console,"ADDING FEATURES TO FEATURESET")
            for genome_ref in target_feature_info['feature_ids_by_genome_ref'].keys():
                for fId in target_feature_info['feature_ids_by_genome_ref'][genome_ref]:
                    id_untrans = genome_ref+self.genome_id_feature_id_delim+fId
                    id_trans = re.sub ('\|',':',id_untrans)  # BLAST seems to make this translation now when id format has simple 'kb|blah' format
                    if id_trans in hit_seq_ids or id_untrans in hit_seq_ids:
                        #self.log(console, 'FOUND HIT '+fId)  # DEBUG
                        accept_fids[id_untrans] = True
                        #fId = id_untrans  # don't change fId for output FeatureSet
                        try:
                            this_genome_ref_list = output_featureSet['elements'][fId]
                        except:
                            output_featureSet['elements'][fId] = []
                            output_featureSet['element_ordering'].append(fId)
                        output_featureSet['elements'][fId].append(genome_ref)

        # Parse Genome hits into FeatureSet
        #
        elif target_type_name == 'Genome':
            seq_total = 0
            output_featureSet = dict()
#            if 'scientific_name' in input_many_genome and input_many_genome['scientific_name'] is not None:
#                output_featureSet['description'] = input_many_genome['scientific_name'] + " - "+search_tool_name+"_Search filtered"
#            else:
#                output_featureSet['description'] = search_tool_name+"_Search filtered"
            output_featureSet['description'] = search_tool_name+"_Search filtered"
            output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()
            for fid in target_feature_info['feature_ids']:
                #if fid == 'AWN69_RS07145' or fid == 'AWN69_RS13375':
                #    self.log(console, 'CHECKING FID '+fid)  # DEBUG
                seq_total += 1
                id_untrans = fid
                id_trans = re.sub ('\|',':',id_untrans)  # BLAST seems to make this translation now when id format has simple 'kb|blah' format
                if id_trans in hit_seq_ids or id_untrans in hit_seq_ids:
                    self.log(console, 'FOUND HIT '+fid)  # DEBUG
                    #output_featureSet['element_ordering'].append(fid)
                    accept_fids[id_untrans] = True
                    #fid = input_many_ref+self.genome_id_feature_id_delim+id_untrans  # don't change fId for output FeatureSet
                    genome_ref = target_ref
                    output_featureSet['element_ordering'].append(fid)
                    output_featureSet['elements'][fid] = [genome_ref]


        # Parse GenomeSet or SpeciesTree hits into FeatureSet
        #
        elif target_type_name == 'GenomeSet' or target_type_name == 'Tree':
            seq_total = 0

            output_featureSet = dict()
            #if 'description' in input_many_genomeSet and input_many_genomeSet['description'] is not None:
            #    output_featureSet['description'] = input_many_genomeSet['description'] + " - "+search_tool_name+"_Search filtered"
            #else:
            #    output_featureSet['description'] = search_tool_name+"_Search filtered"
            output_featureSet['description'] = search_tool_name+"_Search filtered"
            output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()

            self.log(console,"READING HITS FOR GENOMES")  # DEBUG
            for genome_id in list(target_feature_info['feature_ids_by_genome_id'].keys()):
                self.log(console,"READING HITS FOR GENOME "+genome_id)  # DEBUG
                #genome_ref = input_many_genomeSet['elements'][genome_id]['ref']
                genome_ref = target_feature_info['genome_id_to_genome_ref'][genome_id]
                for feature_id in target_feature_info['feature_ids_by_genome_id'][genome_id]:
                    seq_total += 1
                    id_untrans = genome_ref+self.genome_id_feature_id_delim+feature_id
                    id_trans = re.sub ('\|',':',id_untrans)  # BLAST seems to make this translation now when id format has simple 'kb|blah' format
                    if id_trans in hit_seq_ids or id_untrans in hit_seq_ids:
                        #self.log(console, 'FOUND HIT '+fId)  # DEBUG
                        accept_fids[id_untrans] = True
                        #feature_id = id_untrans  # don't change fId for output FeatureSet
                        try:
                            this_genome_ref_list = output_featureSet['elements'][feature_id]
                        except:
                            output_featureSet['elements'][feature_id] = []
                            output_featureSet['element_ordering'].append(feature_id)
                        output_featureSet['elements'][feature_id].append(genome_ref)


        # Parse AnnotatedMetagenomeAssembly hits into FeatureSet
        #
        elif target_type_name == 'AnnotatedMetagenomeAssembly':
            seq_total = 0
            output_featureSet = dict()
#            if 'scientific_name' in input_many_genome and input_many_genome['scientific_name'] is not None:
#                output_featureSet['description'] = input_many_genome['scientific_name'] + " - "+search_tool_name+"_Search filtered"
#            else:
#                output_featureSet['description'] = search_tool_name+"_Search filtered"
            output_featureSet['description'] = search_tool_name+"_Search filtered"
            output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()
            for fid in target_feature_info['feature_ids']:
                #if fid == 'AWN69_RS07145' or fid == 'AWN69_RS13375':
                #    self.log(console, 'CHECKING FID '+fid)  # DEBUG
                seq_total += 1
                id_untrans = fid
                id_trans = re.sub ('\|',':',id_untrans)  # BLAST seems to make this translation now when id format has simple 'kb|blah' format
                #print ("TESTING FEATURES: ID_UNTRANS: '"+id_untrans+"'")  # DEBUG
                #print ("TESTING FEATURES: ID_TRANS: '"+id_trans+"'")  # DEBUG
                if id_trans in hit_seq_ids or id_untrans in hit_seq_ids:
                    self.log(console, 'FOUND HIT '+fid)  # DEBUG
                    #output_featureSet['element_ordering'].append(fid)
                    accept_fids[id_untrans] = True
                    #fid = input_many_ref+self.genome_id_feature_id_delim+id_untrans  # don't change fId for output FeatureSet
                    ama_ref = target_ref
                    output_featureSet['element_ordering'].append(fid)
                    output_featureSet['elements'][fid] = [ama_ref]

        # Parse RefData hits into FeatureSet
        #
        elif target_type_name == 'RefData':
            seq_total = 0

            output_featureSet = dict()
            output_featureSet['description'] = search_tool_name+"_Search filtered"
            output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()

            self.log(console,"READING HITS FOR REFDATA")  # DEBUG
            for id_untrans in list(hit_seq_ids.keys()):
                accept_fids[id_untrans] = True
                [genome_ref, feature_id] = id_untrans.split(self.genome_id_feature_id_delim)
                try:
                    this_genome_ref_list = output_featureSet['elements'][feature_id]
                except:
                    output_featureSet['elements'][feature_id] = []
                    output_featureSet['element_ordering'].append(feature_id)
                output_featureSet['elements'][feature_id].append(genome_ref)


        # Upload results
        #
        output_featureSet_ref = None
        if len(invalid_msgs) == 0 and len(list(hit_seq_ids.keys())) > 0:
            self.log(console,"UPLOADING RESULTS")  # DEBUG

            # we are now making FeatureSets with AMA feature
            if target_type_name != 'SingleEndLibrary' and target_type_name != 'SequenceSet':  
                if num_targets == 1:
                    output_featureSet_name = params['output_filtered_name']
                else:
                    output_featureSet_name = params['output_filtered_name']+'-'+target_name
                new_obj_info = self.wsClient.save_objects({
                            'workspace': params['workspace_name'],
                            'objects':[{
                                    'type': 'KBaseCollections.FeatureSet',
                                    'data': output_featureSet,
                                    'name': output_featureSet_name,
                                    'meta': {},
                                    'hidden': 1,
                                    'provenance': self._instantiate_provenance (method_name=method_name,
                                                                                input_obj_refs=[params['input_one_ref'],target_ref])
                                }]
                        })[0]
                output_featureSet_ref = '/'.join([str(new_obj_info[WSID_I]),
                                                  str(new_obj_info[OBJID_I]),
                                                  str(new_obj_info[VERSION_I])])
            else:
                raise ValueError ("Not currently supporting SingleEndLibrary nor SequenceSet as target type")
                
        return {
            'accept_fids': accept_fids,
            'filtering_fields': filtering_fields,
            'seq_total': seq_total,
            'hit_order': hit_order,
            'hit_total': hit_total,
            'hit_buf': hit_buf,
            'output_featureSet_ref': output_featureSet_ref
        }

    # _get_html_file_name
    #
    def _get_html_file_name(self, this_target_name, search_tool_name):
        return this_target_name+'-'+search_tool_name+'_Search.html'
    
    # _add_html_tabs()
    #
    def _add_html_tabs (self, search_tool_name, targets_name, input_many_refs, input_many_ref):
        # TODO: Add RefData tabs

        html_tabs = []
        for this_input_ref in input_many_refs:
            this_target_name = targets_name[this_input_ref]
            if this_input_ref == input_many_ref:
                html_tabs += [' <b>'+this_target_name+'</b> ' ]
            else:
                html_file_name = self._get_html_file_name(this_target_name, search_tool_name)
                html_tabs += [' <a href="'+html_file_name+'">'+this_target_name+'</a> ']
        
        return [" | ".join(html_tabs)]
    

    # _write_HTML_report()
    #
    def _write_HTML_report(self,
                           search_tool_name = None,
                           input_many_names = None,
                           input_many_refs = None,
                           targets_name = None,
                           targets_type_name = None,
                           targets_feature_info = None,
                           genome_disp_name_config = None,
                           query_len = None,
                           all_parsed_BLAST_results = None):
        html_dir = None
        html_file_path = None
        html_files = []
        console = []
        (q_seq_type, t_seq_type) = self._set_BLAST_seq_types (search_tool_name)

        # config
        head_color = "#eeeeff"
        border_head_color = "#ffccff"
        accept_row_color = 'white'
        #reject_row_color = '#ffeeee'
        reject_row_color = '#eeeeee'
        reject_cell_color = '#ffcccc'
        text_fontsize = "2"
        text_color = '#606060'
        border_body_color = "#cccccc"
        bar_width = 100
        bar_height = 15
        bar_color = "lightblue"
        bar_line_color = "#cccccc"
        bar_fontsize = "1"
        bar_char = "."
        cellpadding = "3"
        cellspacing = "2"
        border = "0"

        # build html page for each target 
        for input_many_ref in input_many_refs:
            
            # set characteristics for this target
            target_name = targets_name[input_many_ref]
            target_type_name = targets_type_name[input_many_ref]
            target_feature_info = targets_feature_info[input_many_ref]
            
            # get blast results for this target
            accept_fids = all_parsed_BLAST_results[input_many_ref]['accept_fids']
            filtering_fields = all_parsed_BLAST_results[input_many_ref]['filtering_fields']
            seq_total = all_parsed_BLAST_results[input_many_ref]['seq_total']
            hit_order = all_parsed_BLAST_results[input_many_ref]['hit_order']
            hit_total = all_parsed_BLAST_results[input_many_ref]['hit_total']
            hit_buf = all_parsed_BLAST_results[input_many_ref]['hit_buf']
            
            # begin buffer and table header
            html_report_lines = []
            html_report_lines += ['<html>']
            html_report_lines += ['<body bgcolor="white">']
            html_report_lines += self._add_html_tabs (search_tool_name, targets_name, input_many_refs, input_many_ref)
            html_report_lines += ['<p>']

            html_report_lines += ['<table cellpadding='+cellpadding+' cellspacing = '+cellspacing+' border='+border+'>']
            html_report_lines += ['<tr bgcolor="'+head_color+'">']
            html_report_lines += ['<td style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'ALIGNMENT COVERAGE'+'</font></td>']
            html_report_lines += ['<td style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'GENE ID'+'</font></td>']
            html_report_lines += ['<td style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'FUNCTION'+'</font></td>']
            html_report_lines += ['<td style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'GENOME'+'</font></td>']
            html_report_lines += ['<td align=center style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'IDENT'+'%</font></td>']
            html_report_lines += ['<td align=center  style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'ALN_LEN'+'</font></td>']
            html_report_lines += ['<td align=center  style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'E-VALUE'+'</font></td>']
            html_report_lines += ['<td align=center  style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'BIT SCORE'+'</font></td>']
            html_report_lines += ['<td align=center  style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'<nobr>Q_BEG-Q_END:</nobr> <nobr>H_BEG-H_END</nobr>'+'</font></td>']
            html_report_lines += ['<td align=center  style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'MIS MATCH'+'</font></td>']
            html_report_lines += ['<td align=center  style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'GAP OPEN'+'</font></td>']
            html_report_lines += ['</tr>']

            # add in hits
            for line in hit_buf:
                line = line.strip()
                if line == '' or line.startswith('#'):
                    continue

                [query_id, hit_id, identity, aln_len, mismatches, gap_openings, q_beg, q_end, h_beg, h_end, e_value, bit_score] = line.split("\t")[0:12]
                aln_len_calc = aln_len
                if q_seq_type != t_seq_type:
                    aln_len_calc *= 3
                aln_len_perc = round (100.0*float(aln_len_calc)/float(query_len), 1)
                identity = str(round(float(identity), 1))
                if identity == '100.0':  identity = '100'

                #if target_type_name == 'SingleEndLibrary':
                #    pass
                #elif target_type_name == 'SequenceSet':
                if target_type_name == 'SequenceSet':
                    pass
                elif target_type_name == 'Genome' or \
                     target_type_name == 'GenomeSet' or \
                     target_type_name == 'Tree' or \
                     target_type_name == 'FeatureSet' or \
                     target_type_name == 'AnnotatedMetagenomeAssembly':

                    if target_type_name != 'Genome' and target_type_name != 'AnnotatedMetagenomeAssembly':
                        [genome_ref, hit_fid] = hit_id.split(self.genome_id_feature_id_delim)
                    else:
                        genome_ref = input_many_ref
                        hit_fid = hit_id

                    # can't just use hit_fid because may have pipes translated and can't translate back
                    fid_lookup = None
                    for fid in list(target_feature_info['feature_id_to_function'][genome_ref].keys()):
                        id_untrans = fid
                        id_trans = re.sub ('\|',':',id_untrans)  # BLAST seems to make this translation now when id format has simple 'kb|blah' format

                        #self.log (console, "SCANNING FIDS.  HIT_FID: '"+str(hit_fid)+"' FID: '"+str(fid)+"' TRANS: '"+str(id_trans)+"'")  # DEBUG

                        if id_untrans == hit_fid or id_trans == hit_fid:
                            #self.log (console, "GOT ONE!")  # DEBUG
                            if target_type_name == 'Genome' or target_type_name == 'AnnotatedMetagenomeAssembly':
                                accept_id = fid
                            elif target_type_name == 'GenomeSet' or target_type_name == 'Tree' or target_type_name == 'FeatureSet':
                                accept_id = genome_ref+self.genome_id_feature_id_delim+fid
                            if accept_id in accept_fids:
                                row_color = accept_row_color
                            else:
                                row_color = reject_row_color
                            fid_lookup = fid
                            break
                    #self.log (console, "HIT_FID: '"+str(hit_fid)+"' FID_LOOKUP: '"+str(fid_lookup)+"'")  # DEBUG
                    if fid_lookup is None:
                        raise ValueError ("unable to find fid for hit_fid: '"+str(hit_fid))
                    elif fid_lookup not in target_feature_info['feature_id_to_function'][genome_ref]:
                        raise ValueError ("unable to find function for fid: '"+str(fid_lookup))
                    fid_disp = re.sub (r"^.*\.([^\.]+)\.([^\.]+)$", r"\1.\2", fid_lookup)

                    func_disp = target_feature_info['feature_id_to_function'][genome_ref][fid_lookup]

                    # set genome_disp_name
                    if target_type_name == 'AnnotatedMetagenomeAssembly':
                        genome_disp_name = target_feature_info['ama_ref_to_obj_name'][genome_ref]
                    else:
                        genome_obj_name = target_feature_info['genome_ref_to_obj_name'][genome_ref]
                        genome_sci_name = target_feature_info['genome_ref_to_sci_name'][genome_ref]
                        [ws_id, obj_id, genome_obj_version] = genome_ref.split('/')
                        genome_disp_name = ''
                        if 'obj_name' in genome_disp_name_config:
                            genome_disp_name = genome_obj_name
                        if 'ver' in genome_disp_name_config:
                            genome_disp_name += '.v'+str(genome_obj_version)
                        if 'sci_name' in genome_disp_name_config:
                            if genome_disp_name != '':
                                genome_disp_name += ': '+genome_sci_name
                            else:
                                genome_disp_name = genome_sci_name
                        if genome_disp_name == '':
                            genome_disp_name = genome_obj_name

                    #if 'overlap_fraction' in params and float(params['overlap_fraction']) > float(high_bitscore_alnlen[hit_seq_id])/float(query_len):

                    html_report_lines += ['<tr bgcolor="'+row_color+'">']
                    #html_report_lines += ['<tr bgcolor="'+'white'+'">']  # DEBUG
                    # add overlap bar

                    # coverage graphic
                    html_report_lines += ['<td valign=middle align=center style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'">']
                    html_report_lines += ['<table style="height:'+str(bar_height)+'px; width:'+str(bar_width)+'px" border=0 cellpadding=0 cellspacing=0>']
                    full_len_pos = bar_width
                    aln_beg_pos = int (float(bar_width) * float(int(q_beg)-1)/float(int(query_len)-1))
                    aln_end_pos = int (float(bar_width) * float(int(q_end)-1)/float(int(query_len)-1))
                    cell_pix_height = str(int(round(float(bar_height)/3.0, 0)))

                    cell_color = ['','','']
                    cell_width = []
                    cell_width.append(aln_beg_pos)
                    cell_width.append(aln_end_pos-aln_beg_pos)
                    cell_width.append(bar_width-aln_end_pos)

                    for row_i in range(3):
                        html_report_lines += ['<tr style="height:'+cell_pix_height+'px">']
                        unalign_color = row_color
                        if row_i == 1:
                            unalign_color = bar_line_color
                        cell_color[0] = unalign_color
                        cell_color[1] = bar_color
                        cell_color[2] = unalign_color

                        for col_i in range(3):
                            cell_pix_width = str(cell_width[col_i])
                            cell_pix_color = cell_color[col_i]
                            html_report_lines += ['<td style="height:'+cell_pix_height+'px; width:'+cell_pix_width+'px" bgcolor="'+cell_pix_color+'"></td>']
                        html_report_lines += ['</tr>']
                    html_report_lines += ['</table>']
                    html_report_lines += ['</td>']

                    # add other cells
                    # fid
                    html_report_lines += ['<td style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+str(fid_disp)+'</font></td>']
                    # func
                    html_report_lines += ['<td style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+func_disp+'</font></td>']
                    # sci name
                    html_report_lines += ['<td style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+genome_disp_name+'</font></td>']
                    # ident
                    if 'ident_thresh' in filtering_fields[hit_id]:
                        this_cell_color = reject_cell_color
                    else:
                        this_cell_color = row_color
                    html_report_lines += ['<td align=center bgcolor="'+this_cell_color+'" style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+str(identity)+'%</font></td>']
                    # aln len
                    if 'overlap_fraction' in filtering_fields[hit_id]:
                        this_cell_color = reject_cell_color
                    else:
                        this_cell_color = row_color
                    html_report_lines += ['<td align=center bgcolor="'+this_cell_color+'" style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+str(aln_len)+' ('+str(aln_len_perc)+'%)</font></td>']
                    # evalue
                    html_report_lines += ['<td align=center style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'><nobr>'+str(e_value)+'</nobr></font></td>']
                    # bitscore
                    if 'bitscore' in filtering_fields[hit_id]:
                        this_cell_color = reject_cell_color
                    else:
                        this_cell_color = row_color
                    html_report_lines += ['<td align=center bgcolor="'+this_cell_color+'" style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+str(bit_score)+'</font></td>']
                    # aln coords
                    html_report_lines += ['<td align=center style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'><nobr>'+str(q_beg)+'-'+str(q_end)+':</nobr> <nobr>'+str(h_beg)+'-'+str(h_end)+'</nobr></font></td>']
                    # mismatches
                    html_report_lines += ['<td align=center style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+str(mismatches)+'</font></td>']
                    # gaps
                    html_report_lines += ['<td align=center style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+str(gap_openings)+'</font></td>']
                    html_report_lines += ['</tr>']

            html_report_lines += ['</table>']
            html_report_lines += ['</body>']
            html_report_lines += ['</html>']

            # write html to file and upload
            html_report_str = "\n".join(html_report_lines)
            html_file = self._get_html_file_name(target_name, search_tool_name)
            if not html_dir:
                html_dir = self._set_HTML_outdir()
            html_path = self._set_HTML_file_path (html_dir, html_file)
            with open (html_path, 'w') as html_handle:
                html_handle.write(html_report_str)
            html_files.append(html_file)
            
        return (html_dir, html_files)


    #### build output report
    ##
    def build_BLAST_report (self, 
                            search_tool_name = None,
                            params = None,
                            targets_name = None,
                            targets_type_name = None,
                            targets_feature_info = None,
                            base_bulk_save_infos = None,
                            extra_bulk_save_infos = None,
                            query_len = None,
                            all_parsed_BLAST_results = None,
                            objects_created = None):

        # init
        method_name = search_tool_name+'_Search'
        console = []
        invalid_msgs = []
        report = ''
        self.log(console,"BUILDING REPORT")  # DEBUG

        input_many_refs = params['input_many_refs']

        # check output before proceeding
        all_hit_total = 0
        all_hit_order = []
        for input_many_ref in input_many_refs:

            target_name = targets_name[input_many_ref]
            
            # pass args
            accept_fids = all_parsed_BLAST_results[input_many_ref]['accept_fids']
            filtering_fields = all_parsed_BLAST_results[input_many_ref]['filtering_fields']
            seq_total = all_parsed_BLAST_results[input_many_ref]['seq_total']
            hit_order = all_parsed_BLAST_results[input_many_ref]['hit_order']
            hit_total = all_parsed_BLAST_results[input_many_ref]['hit_total']
            hit_buf   = all_parsed_BLAST_results[input_many_ref]['hit_buf']

            all_hit_total += hit_total
            all_hit_order.extend(hit_order)
            
            # text report
            report += 'TARGET: '+target_name+"\n"
            report += 'sequences in search db: '+str(seq_total)+"\n"
            report += 'sequences in hit set: '+str(len(hit_order))+"\n"
            report += 'sequences in accepted hit set: '+str(hit_total)+"\n"
            report += "\n"
            for line in hit_buf:
                report += line
            self.log (console, report)

        # don't waste time if no hits
        if all_hit_total == 0 and len(all_hit_order) == 0:
            report += "No hits were found\n"
            reportObj = {
                'objects_created':[],
                'text_message':report
                }

            reportName = 'blast_report_'+str(uuid.uuid4())
            report_obj_info = self.wsClient.save_objects({
                    #                'id':info[6],
                    'workspace':params['workspace_name'],
                    'objects':[
                        {
                            'type':'KBaseReport.Report',
                            'data':reportObj,
                            'name':reportName,
                            'meta':{},
                            'hidden':1,
                            'provenance':self._instantiate_provenance(method_name=method_name,
                                                                      input_obj_refs=[params['input_one_ref']]+params['input_many_refs'])
                            }
                        ]
                    })[0]
            report_info = dict()
            report_info['name'] = report_obj_info[1]
            report_info['ref'] = str(report_obj_info[6])+'/'+str(report_obj_info[0])+'/'+str(report_obj_info[4])
            
            return report_info

        # build html report
        (html_dir, html_file_names) = self._write_HTML_report (search_tool_name = search_tool_name,
                                                input_many_refs = params['input_many_refs'],
                                                targets_name = targets_name,
                                                targets_type_name = targets_type_name,
                                                targets_feature_info = targets_feature_info,
                                                genome_disp_name_config = params['genome_disp_name_config'],
                                                query_len = query_len,
                                                all_parsed_BLAST_results = all_parsed_BLAST_results)

        # upload html report
        dfu = DFUClient(self.callbackURL)
        try:
            html_upload_ret = dfu.file_to_shock({'file_path': html_dir,
                                                 'make_handle': 0,
                                                 'pack': 'zip'})
        except:
            raise ValueError ('Logging exception loading html_report to shock')


        # create report object
        reportName = 'blast_report_'+str(uuid.uuid4())
        reportObj = {'objects_created': [],
                     #'text_message': '',  # or is it 'message'?
                     'message': '',  # or is it 'text_message'?
                     'direct_html': '',
                     'direct_html_link_index': None,
                     'file_links': [],
                     'html_links': [],
                     'workspace_name': params['workspace_name'],
                     'report_object_name': reportName
        }
        #html_buf_lim = 16000  # really 16KB, but whatever
        #if len(html_report_str) <= html_buf_lim:
        #    reportObj['direct_html'] = html_report_str
        #else:
        #    reportObj['direct_html_link_index'] = 0

        reportObj['direct_html_link_index'] = 0
        reportObj['html_links'] = [{'shock_id': html_upload_ret['shock_id'],
                                    #'name': search_tool_name+'_results.html',
                                    'name': html_file_names[0],
                                    'label': search_tool_name+' Results'}
        ]
        reportObj['file_links'] = []
        for input_many_ref in input_many_refs:
            target_name = targets_name[input_many_ref]
            base_bulk_save_info = base_bulk_save_infos[input_many_ref]
            if input_many_ref in extra_bulk_save_infos:
                extra_bulk_save_info = extra_bulk_save_infos[input_many_ref]
            else:
                extra_bulk_save_info = None
                
            reportObj['file_links'].append({'shock_id': base_bulk_save_info['shock_id'],
                                            'name': target_name+'-'+search_tool_name+'_Search-m'+'7'+'.txt',
                                            'label': target_name+'-'+search_tool_name+' Results: m'+'7'})

            if extra_bulk_save_info is not None:
                extension = 'txt'
                if params['output_extra_format'] == '5':
                    extension = 'xml'
                elif params['output_extra_format'] == '8':
                    extension = 'asn1txt'
                elif params['output_extra_format'] == '9':
                    extension = 'asn1bin'
                elif params['output_extra_format'] == '10':
                    extension = 'csv'
                elif params['output_extra_format'] == '11':
                    extension = 'asn1arc'
                reportObj['file_links'].append({'shock_id': extra_bulk_save_info['shock_id'],
                                                'name': target_name+'-'+search_tool_name+'_Search-m'+str(params['output_extra_format'])+'.'+extension,
                                                'label': target_name+'-'+search_tool_name+' Results: m'+str(params['output_extra_format'])})
                            
                            
        # complete report
        objects_created.reverse()  # want merged featureset at position 0
        reportObj['objects_created'] = objects_created
        
        ##reportObj['message'] = report

        # save report object
        report_info = self.reportClient.create_extended_report(reportObj)

        return report_info


    #### run_BLAST_App(): top-level method
    ##
    def run_BLAST_App (self, search_tool_name, params):
        console = []
        invalid_msgs = []
        method_name = search_tool_name+'_Search()'
        self.log(console,'Running '+search_tool_name+'_Search with params=')
        self.log(console, "\n"+pformat(params))
        report = ''
#        report = 'Running '+search_tool_name+'_Search with params='
#        report += "\n"+pformat(params)
        appropriate_sequence_found_in_one_input = False
        appropriate_sequence_found_in_many_input = False
        base_BLAST_output_format = '7'
        (q_seq_type, t_seq_type) = self._set_BLAST_seq_types (search_tool_name)
        
        #### Validate App input params
        #
        if not params.get('output_one_name'):
            if not params.get('input_one_ref'):
                params['output_one_name'] = 'query-'+params['output_filtered_name']+'.Seq'
                
        if not self.validate_BLAST_app_params (params, method_name):
            raise ValueError('App input validation failed in validate_BLAST_app_params() for App ' + method_name)

        # Get input obj refs
        #
        input_many_refs = params.get('input_many_refs',[])

        # Get RefData targets
        #
        refdata_targets = []
        refdata_names = {
            'Archaea-RS': 'GTDB_Species_Reps-Archaea-RefSeq',
            'Archaea-GB': 'GTDB_Species_Reps-Archaea-GenBank',
            'Bacteria-RS': 'GTDB_Species_Reps-Bacteria-RefSeq',
            'Bacteria-GB': 'GTDB_Species_Reps-Bacteria-GenBank'
            }
        if params.get('gtdb_targets','none') != 'none':
            if params['gtdb_targets'] == 'all':
                refdata_targets.extend(['Archaea-RS', 'Archaea-GB', 'Bacteria-RS', 'Bacteria-GB'])
            elif params['gtdb_targets'] == 'archaea-RS+GB':
                refdata_targets.extend(['Archaea-RS', 'Archaea-GB'])
            elif params['gtdb_targets'] == 'archaea-RS':
                refdata_targets.extend(['Archaea-RS'])
            elif params['gtdb_targets'] == 'archaea-GB':
                refdata_targets.extend(['Archaea-GB'])
            elif params['gtdb_targets'] == 'bacteria-RS+GB':
                refdata_targets.extend(['Bacteria-RS', 'Bacteria-GB'])
            elif params['gtdb_targets'] == 'bacteria-RS':
                refdata_targets.extend(['Bacteria-RS'])
            elif params['gtdb_targets'] == 'bacteria-GB':
                refdata_targets.extend(['Bacteria-GB'])
            elif params['gtdb_targets'] == 'archaea-bacteria-RS':
                refdata_targets.extend(['Archaea-RS', 'Bacteria-RS'])
            elif params['gtdb_targets'] == 'archaea-bacteria-GB':
                refdata_targets.extend(['Archaea-GB', 'Bacteria-GB'])
                
        # Get query
        #
        if params.get('input_one_sequence') is not None \
                and not params['input_one_sequence'].startswith("Optionally enter"):
            input_one_ref = self.objectify_text_query (params, q_seq_type, method_name)
            params['input_one_ref'] = input_one_ref
        else:
            input_one_ref = params['input_one_ref']


        # Write query obj to fasta file (can be Feature, SequenceSet, or FeatureSet)
        #
        write_query_obj_to_file_result = self.write_query_obj_to_file (params, input_one_ref, q_seq_type)
        query_type_name = write_query_obj_to_file_result['query_type_name']
        query_fasta_file_path = write_query_obj_to_file_result['query_fasta_file_path']
        appropriate_sequence_found_in_one_input = write_query_obj_to_file_result['appropriate_sequence_found_in_one_input']
        invalid_msgs.extend(write_query_obj_to_file_result['invalid_msgs'])


        # Write target obj to fasta file (can be Feature, SequenceSet, or FeatureSet)
        #
        targets_name = dict()
        targets_type_name = dict()
        targets_fasta_file_path = dict()
        appropriate_sequence_found_in_many_inputs = dict()
        targets_feature_info = dict()
        for input_many_ref in input_many_refs:
            write_target_obj_to_file_result = self.write_target_obj_to_file (params, input_many_ref, t_seq_type)
            targets_name[input_many_ref] = write_target_obj_to_file_result['target_name']
            targets_type_name[input_many_ref] = write_target_obj_to_file_result['target_type_name']
            targets_fasta_file_path[input_many_ref] = write_target_obj_to_file_result['target_fasta_file_path']
            appropriate_sequence_found_in_many_inputs[input_many_ref] = write_target_obj_to_file_result['appropriate_sequence_found_in_many_input']
            invalid_msgs.extend(write_target_obj_to_file_result['invalid_msgs'])

            targets_feature_info[input_many_ref] = write_target_obj_to_file_result['target_feature_info']


        # Add refdata info to targets
        #
        for refdata_target in refdata_targets:
            targets_name[refdata_target] = refdata_names[refdata_target]
            targets_type_name[refdata_target] = 'RefData'
            targets_fasta_file_path[refdata_target] = os.path.join(os.sep,'data','blast_dbs',refdata_target)
            appropriate_sequence_found_in_many_inputs[refdata_target] = True
            targets_feature_info[refdata_target] = None
            

        # check for failed input file creation
        #
        if not appropriate_sequence_found_in_one_input:
            self.log(invalid_msgs,"no "+q_seq_type+" sequence found in '"+input_one_name+"'")
        for input_many_ref in input_many_refs + refdata_targets:
            if not appropriate_sequence_found_in_many_inputs[input_many_ref]:
                self.log(invalid_msgs,"no "+t_seq_type+" sequences found in '"+input_many_name+"'")

        if len(invalid_msgs) > 0:
            error_report_info = self.save_error_report_with_invalid_msgs (invalid_msgs, input_one_ref, input_many_refs, method_name)
            returnVal = { 'report_name': report_info['name'],
                          'report_ref': report_info['ref']
                      }        
            return [returnVal]


        #### FORMAT DB
        ##
        for input_many_ref in input_many_refs:  # Note: RefData dbs already formatted
            if not self.format_BLAST_db (search_tool_name, targets_fasta_file_path[input_many_ref]):
                raise ValueError ("failed to format BLAST db for "+input_many_ref)
            

        #### Run BLAST for base format
        ##
        output_aln_file_paths = dict()
        base_bulk_save_infos = dict()
        for input_many_ref in input_many_refs + refdata_targets:
            if input_many_ref in refdata_targets:
                max_accepts = str(params['refdata_maxaccepts'])
            else:
                max_accepts = str(params['maxaccepts'])
            BLAST_output_results = self.run_BLAST (search_tool_name = search_tool_name, 
                                                   query_fasta_file_path = query_fasta_file_path, 
                                                   target_fasta_file_path = targets_fasta_file_path[input_many_ref],
                                                   e_value = str(params['e_value']),
                                                   maxaccepts = max_accepts,
                                                   BLAST_output_format_str = str(base_BLAST_output_format)
            )
            output_aln_file_paths[input_many_ref] = BLAST_output_results['output_aln_file_path']
            base_bulk_save_infos[input_many_ref] = BLAST_output_results['bulk_save_info']


        #### Run BLAST for extra format
        ##
        output_extra_aln_file_paths = dict()
        extra_bulk_save_infos = dict()
        for input_many_ref in input_many_refs + refdata_targets:
            if input_many_ref in refdata_targets:
                max_accepts = str(params['refdata_maxaccepts'])
            else:
                max_accepts = str(params['maxaccepts'])
            if str(params.get('output_extra_format')) and str(params.get('output_extra_format')) != 'none':
                BLAST_extra_output_results = self.run_BLAST (search_tool_name = search_tool_name,
                                                             query_fasta_file_path = query_fasta_file_path, 
                                                             target_fasta_file_path = targets_fasta_file_path[input_many_ref],
                                                             e_value = str(params['e_value']),
                                                             maxaccepts = max_accepts,
                                                             BLAST_output_format_str = str(params['output_extra_format'])
                )

                output_extra_aln_file_paths[input_many_ref] = BLAST_extra_output_results['output_aln_file_path']
                extra_bulk_save_infos[input_many_ref] = BLAST_extra_output_results['bulk_save_info']


        # get query_len for filtering and reporting later
        #
        query_len = self.get_query_len (query_fasta_file_path)


        # Parse the BLAST tabular output and store ids to filter many set to make filtered object to save back to KBase
        #
        all_parsed_BLAST_results = dict()
        objects_created = []
        output_featureSet_refs = []
        num_targets = len(input_many_refs) + len(refdata_targets)
        for input_many_ref in input_many_refs + refdata_targets:
            this_parsed_BLAST_results = \
                self.parse_BLAST_tab_output (output_aln_file_path = output_aln_file_paths[input_many_ref],
                                             search_tool_name = search_tool_name,
                                             params = params,
                                             query_len = query_len,
                                             num_targets = num_targets,
                                             target_ref = input_many_ref,
                                             target_name = targets_name[input_many_ref],
                                             target_type_name = targets_type_name[input_many_ref],
                                             target_feature_info = targets_feature_info[input_many_ref])

            all_parsed_BLAST_results[input_many_ref] = this_parsed_BLAST_results
            if this_parsed_BLAST_results.get('output_featureSet_ref'):
                objects_created.append({'ref':this_parsed_BLAST_results['output_featureSet_ref'],'description':targets_name[input_many_ref]+" "+search_tool_name+' hits'})
                output_featureSet_refs.append(this_parsed_BLAST_results['output_featureSet_ref'])

        # Merge FeatureSets into one output
        if len(output_featureSet_refs) > 1:
            self.log(console, "CREATING MERGED OUTPUT FEATURESET")
            merge_featureSet_params = {'workspace_name': params['workspace_name'],
                                       'input_refs': output_featureSet_refs,
                                       'desc': 'Merged FeatureSets from '+search_tool_name+' Search',
                                       'output_name': params['output_filtered_name']
                                       }
            merge_retVal = self.set_util_Client.KButil_Merge_FeatureSet_Collection(merge_featureSet_params)
            merge_report_obj = self.wsClient.get_objects2({'objects':[{'ref':merge_retVal['report_ref']}]})['data'][0]['data']
            merged_featureSet_ref = merge_report_obj['objects_created'][0]['ref']
            objects_created.append({'ref': merged_featureSet_ref,'description':'ALL '+search_tool_name+' hits'})

            
        # build output report object
        #
        if len(invalid_msgs) == 0:
            report_info = self.build_BLAST_report (search_tool_name = search_tool_name,
                                                   params = params,
                                                   targets_name = targets_name,
                                                   targets_type_name = targets_type_name,
                                                   targets_feature_info = targets_feature_info,
                                                   base_bulk_save_infos = base_bulk_save_infos,
                                                   extra_bulk_save_infos = extra_bulk_save_infos,
                                                   query_len = query_len, 
                                                   all_parsed_BLAST_results = all_parsed_BLAST_results,
                                                   objects_created = objects_created)
        else:
            report_info = error_report_info


        # return
        #
        self.log(console,search_tool_name+"_Search DONE")
        returnVal = { 'report_name': report_info['name'],
                      'report_ref': report_info['ref']
                      }
        return returnVal
