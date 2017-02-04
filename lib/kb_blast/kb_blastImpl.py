# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os
import sys
import shutil
import hashlib
import subprocess
import requests
import re
import traceback
import uuid
from datetime import datetime
from pprint import pprint, pformat
import numpy as np
import gzip

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from biokbase.workspace.client import Workspace as workspaceService
from requests_toolbelt import MultipartEncoder
from biokbase.AbstractHandle.Client import AbstractHandle as HandleService

# SDK Utils
from KBaseDataObjectToFileUtils.KBaseDataObjectToFileUtilsClient import KBaseDataObjectToFileUtils
from KBaseReport.KBaseReportClient import KBaseReport

# silence whining
import requests
requests.packages.urllib3.disable_warnings()

#END_HEADER


class kb_blast:
    '''
    Module Name:
    kb_blast

    Module Description:
    ** A KBase module: kb_blast
**
** This module contains 6 methods from BLAST+: BLASTn, BLASTp, BLASTx, tBLASTx, tBLASTn, and PSI-BLAST
**
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.2"
    GIT_URL = "https://github.com/kbaseapps/kb_blast.git"
    GIT_COMMIT_HASH = "dc738060c636bacef584c58816e67aa4ad2803cf"

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

    # target is a list for collecting log messages
    def log(self, target, message):
        # we should do something better here...
        if target is not None:
            target.append(message)
        print(message)
        sys.stdout.flush()


    # Helper script borrowed from the transform service, logger removed
    #
    def upload_file_to_shock(self,
                             console,  # DEBUG
                             shock_service_url = None,
                             filePath = None,
                             ssl_verify = True,
                             token = None):
        """
        Use HTTP multi-part POST to save a file to a SHOCK instance.
        """
        self.log(console,"UPLOADING FILE "+filePath+" TO SHOCK")

        if token is None:
            raise Exception("Authentication token required!")

        #build the header
        header = dict()
        header["Authorization"] = "Oauth {0}".format(token)
        if filePath is None:
            raise Exception("No file given for upload to SHOCK!")

        dataFile = open(os.path.abspath(filePath), 'rb')
        m = MultipartEncoder(fields={'upload': (os.path.split(filePath)[-1], dataFile)})
        header['Content-Type'] = m.content_type

        #logger.info("Sending {0} to {1}".format(filePath,shock_service_url))
        try:
            response = requests.post(shock_service_url + "/node", headers=header, data=m, allow_redirects=True, verify=ssl_verify)
            dataFile.close()
        except:
            dataFile.close()
            raise
        if not response.ok:
            response.raise_for_status()
        result = response.json()
        if result['error']:
            raise Exception(result['error'][0])
        else:
            return result["data"]


    def upload_SingleEndLibrary_to_shock_and_ws (self,
                                                 ctx,
                                                 console,  # DEBUG
                                                 workspace_name,
                                                 obj_name,
                                                 file_path,
                                                 provenance,
                                                 sequencing_tech):

        self.log(console,'UPLOADING FILE '+file_path+' TO '+workspace_name+'/'+obj_name)

        # 1) upload files to shock
        token = ctx['token']
        forward_shock_file = self.upload_file_to_shock(
            console,  # DEBUG
            shock_service_url = self.shockURL,
            filePath = file_path,
            token = token
            )
        #pprint(forward_shock_file)
        self.log(console,'SHOCK UPLOAD DONE')

        # 2) create handle
        self.log(console,'GETTING HANDLE')
        hs = HandleService(url=self.handleURL, token=token)
        forward_handle = hs.persist_handle({
                                        'id' : forward_shock_file['id'], 
                                        'type' : 'shock',
                                        'url' : self.shockURL,
                                        'file_name': forward_shock_file['file']['name'],
                                        'remote_md5': forward_shock_file['file']['checksum']['md5']})

        
        # 3) save to WS
        self.log(console,'SAVING TO WORKSPACE')
        single_end_library = {
            'lib': {
                'file': {
                    'hid':forward_handle,
                    'file_name': forward_shock_file['file']['name'],
                    'id': forward_shock_file['id'],
                    'url': self.shockURL,
                    'type':'shock',
                    'remote_md5':forward_shock_file['file']['checksum']['md5']
                },
                'encoding':'UTF8',
                'type':'fasta',
                'size':forward_shock_file['file']['size']
            },
            'sequencing_tech':sequencing_tech
        }
        self.log(console,'GETTING WORKSPACE SERVICE OBJECT')
        ws = workspaceService(self.workspaceURL, token=ctx['token'])
        self.log(console,'SAVE OPERATION...')
        new_obj_info = ws.save_objects({
                        'workspace':workspace_name,
                        'objects':[
                            {
                                'type':'KBaseFile.SingleEndLibrary',
                                'data':single_end_library,
                                'name':obj_name,
                                'meta':{},
                                'provenance':provenance
                            }]
                        })[0]
        self.log(console,'SAVED TO WORKSPACE')

        return new_obj_info[0]

    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
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


    def BLASTn_Search(self, ctx, params):
        """
        Methods for BLAST of various flavors of one sequence against many sequences 
        **
        **    overloading as follows:
        **        input_one_type: SequenceSet, Feature, FeatureSet
        **        input_many_type: SequenceSet, SingleEndLibrary, FeatureSet, Genome, GenomeSet
        **        output_type: SequenceSet (if input_many is SS), SingleEndLibrary (if input_many is SELib), (else) FeatureSet
        :param params: instance of type "BLAST_Params" (BLAST Input Params)
           -> structure: parameter "workspace_name" of type "workspace_name"
           (** The workspace object refs are of form: ** **    objects =
           ws.get_objects([{'ref':
           params['workspace_id']+'/'+params['obj_name']}]) ** ** "ref" means
           the entire name combining the workspace id and the object name **
           "id" is a numerical identifier of the workspace or object, and
           should just be used for workspace ** "name" is a string identifier
           of a workspace or object.  This is received from Narrative.),
           parameter "input_one_sequence" of type "sequence", parameter
           "input_one_ref" of type "data_obj_ref", parameter "input_many_ref"
           of type "data_obj_ref", parameter "input_msa_ref" of type
           "data_obj_ref", parameter "output_one_name" of type
           "data_obj_name", parameter "output_filtered_name" of type
           "data_obj_name", parameter "ident_thresh" of Double, parameter
           "e_value" of Double, parameter "bitscore" of Double, parameter
           "overlap_fraction" of Double, parameter "maxaccepts" of Double,
           parameter "rounds" of Double
        :returns: instance of type "BLAST_Output" (BLAST Output) ->
           structure: parameter "report_name" of type "data_obj_name",
           parameter "report_ref" of type "data_obj_ref"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN BLASTn_Search
        console = []
        invalid_msgs = []
        search_tool_name = 'BLASTn'
        self.log(console,'Running '+search_tool_name+'_Search with params=')
        self.log(console, "\n"+pformat(params))
        report = ''
#        report = 'Running '+search_tool_name+'_Search with params='
#        report += "\n"+pformat(params)
        appropriate_sequence_found_in_one_input = False
        appropriate_sequence_found_in_many_input = False


        #### do some basic checks
        #
        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
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
#        if 'input_one_sequence' not in params:
#            raise ValueError('input_one_sequence parameter is required')
#        if 'input_one_ref' not in params:
#            raise ValueError('input_one_ref parameter is required')
        if 'input_many_ref' not in params:
            raise ValueError('input_many_ref parameter is required')
        if 'output_filtered_name' not in params:
            raise ValueError('output_filtered_name parameter is required')


        # set local names
        input_one_ref  = None
        output_one_ref = None
        input_many_ref = params['input_many_ref']
        

        # Write the input_one_sequence to a SequenceSet object
        #
        if 'input_one_sequence' in params \
                and params['input_one_sequence'] != None \
                and params['input_one_sequence'] != "Optionally enter DNA sequence...":

            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'])
            ParseFastaStr_retVal = DOTFU.ParseFastaStr ({
                'fasta_str':     params['input_one_sequence'],
                'residue_type': 'NUC',
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
            if 'provenance' in ctx:
                provenance = ctx['provenance']
            # add additional info to provenance here, in this case the input data object reference
            provenance[0]['input_ws_objects'] = []
            provenance[0]['service'] = 'kb_blast'
            provenance[0]['method'] = search_tool_name+'_Search'
                
            # Upload results
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
                ws = workspaceService(self.workspaceURL, token=ctx['token'])
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


        #### Get the input_one object
        ##
        if 'output_one_name' in params and output_one_ref:
            input_one_ref = output_one_ref
        else:
            input_one_ref = params['input_one_ref']

        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            objects = ws.get_objects([{'ref': input_one_ref}])
            input_one_data = objects[0]['data']
            input_one_name = str(objects[0]['info'][1])
            info = objects[0]['info']
                                                             
            one_type_name = info[2].split('.')[1].split('-')[0]
        except Exception as e:
            raise ValueError('Unable to fetch input_one_ref object from workspace: ' + str(e))
        #to get the full stack trace: traceback.format_exc()

        if 'input_one_sequence' in params \
                and params['input_one_sequence'] != None \
                and params['input_one_sequence'] != "Optionally enter DNA sequence..." \
                and one_type_name != 'SequenceSet':

            self.log(invalid_msgs,"ERROR: Mismatched input type for Query Object: "+input_one_ref+" should be SequenceSet instead of: "+one_type_name)


        # Handle overloading (input_one can be Feature, SequenceSet, or FeatureSet)
        #

        # SequenceSet
        #
        if one_type_name == 'SequenceSet':
            try:
                input_one_sequenceSet = input_one_data
            except Exception as e:
                print(traceback.format_exc())
                raise ValueError('Unable to get sequenceSet object: ' + str(e))

            header_id = input_one_sequenceSet['sequences'][0]['sequence_id']
            sequence_str = input_one_data['sequences'][0]['sequence']

            #PROT_pattern = re.compile("^[acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWYxX ]+$")
            DNA_pattern  = re.compile("^[acgtuACGTUnryNRY ]+$")   
            if not DNA_pattern.match(sequence_str):
                self.log(invalid_msgs,"BAD record for sequence_id: "+header_id+"\n"+sequence_str+"\n")
            else:
                appropriate_sequence_found_in_one_input = True

            one_forward_reads_file_path = os.path.join(self.scratch, header_id+'.fasta')
            one_forward_reads_file_handle = open(one_forward_reads_file_path, 'w', 0)
            self.log(console, 'writing reads file: '+str(one_forward_reads_file_path))
            one_forward_reads_file_handle.write('>'+header_id+"\n")
            one_forward_reads_file_handle.write(sequence_str+"\n")
            one_forward_reads_file_handle.close();

            self.log(console, 'done')

        # FeatureSet
        #
        elif one_type_name == 'FeatureSet':
            # retrieve sequences for features
            #input_one_featureSet = input_one_data
            genome_id_feature_id_delim = '.f:'
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
                'residue_type':        'nucleotide',
                'feature_type':        'ALL',
                'record_id_pattern':   '%%genome_ref%%'+genome_id_feature_id_delim+'%%feature_id%%',
                'record_desc_pattern': '[%%genome_ref%%]',
                'case':                'upper',
                'linewrap':            50,
                'merge_fasta_files':   'TRUE'
                }

            #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'])
            FeatureSetToFASTA_retVal = DOTFU.FeatureSetToFASTA (FeatureSetToFASTA_params)
            one_forward_reads_file_path = FeatureSetToFASTA_retVal['fasta_file_path']
            if len(FeatureSetToFASTA_retVal['feature_ids_by_genome_ref'].keys()) > 0:
                appropriate_sequence_found_in_one_input = True


            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "FeatureSetToFasta() took "+str(end_time-beg_time)+" secs")

        # Feature
        #
        elif one_type_name == 'Feature':
            # export feature to FASTA file
            feature = input_one_data
            one_forward_reads_file_path = os.path.join(self.scratch, input_one_name+".fasta")
            self.log(console, 'writing fasta file: '+one_forward_reads_file_path)
            # BLASTn is nuc-nuc
            if feature['dna_sequence'] != None:
                record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description='['+feature['genome_id']+']'+' '+feature['function'])
            #record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description='['+feature['genome_id']+']'+' '+feature['function'])
                SeqIO.write([record], one_forward_reads_file_path, "fasta")
                appropriate_sequence_found_in_one_input = True
        else:
            raise ValueError('Cannot yet handle input_one type of: '+type_name)


        #### Get the input_many object
        ##
        many_forward_reads_file_compression = None
        sequencing_tech = 'N/A'
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            objects = ws.get_objects([{'ref': input_many_ref}])
            input_many_data = objects[0]['data']
            info = objects[0]['info']
            input_many_name = str(info[1])
            many_type_name = info[2].split('.')[1].split('-')[0]

            if many_type_name == 'SingleEndLibrary':
                many_type_namespace = info[2].split('.')[0]
                if many_type_namespace == 'KBaseAssembly':
                    file_name = input_many_data['handle']['file_name']
                elif many_type_namespace == 'KBaseFile':
                    file_name = input_many_data['lib']['file']['file_name']
                else:
                    raise ValueError('bad data type namespace: '+many_type_namespace)
                #self.log(console, 'INPUT_MANY_FILENAME: '+file_name)  # DEBUG
                if file_name[-3:] == ".gz":
                    many_forward_reads_file_compression = 'gz'
                if 'sequencing_tech' in input_many_data:
                    sequencing_tech = input_many_data['sequencing_tech']

        except Exception as e:
            raise ValueError('Unable to fetch input_many_name object from workspace: ' + str(e))
            #to get the full stack trace: traceback.format_exc()


        # Handle overloading (input_many can be SequenceSet, SingleEndLibrary, FeatureSet, Genome, or GenomeSet)
        #
        if many_type_name == 'SequenceSet':
            try:
                input_many_sequenceSet = input_many_data
            except Exception as e:
                print(traceback.format_exc())
                raise ValueError('Unable to get SequenceSet: ' + str(e))

            header_id = input_many_sequenceSet['sequences'][0]['sequence_id']
            many_forward_reads_file_path = os.path.join(self.scratch, header_id+'.fasta')
            many_forward_reads_file_handle = open(many_forward_reads_file_path, 'w', 0)
            self.log(console, 'writing reads file: '+str(many_forward_reads_file_path))

            for seq_obj in input_many_sequenceSet['sequences']:
                header_id = seq_obj['sequence_id']
                sequence_str = seq_obj['sequence']

                #PROT_pattern = re.compile("^[acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWYxX ]+$")
                DNA_pattern = re.compile("^[acgtuACGTUnryNRY ]+$")   
                if not DNA_pattern.match(sequence_str):
                    self.log(invalid_msgs,"BAD record for sequence_id: "+header_id+"\n"+sequence_str+"\n")
                    continue
                appropriate_sequence_found_in_many_input = True
                many_forward_reads_file_handle.write('>'+header_id+"\n")
                many_forward_reads_file_handle.write(sequence_str+"\n")
            many_forward_reads_file_handle.close();
            self.log(console, 'done')

        # SingleEndLibrary
        #
        elif many_type_name == 'SingleEndLibrary':

            # DEBUG
            #for k in data:
            #    self.log(console,"SingleEndLibrary ["+k+"]: "+str(data[k]))

            try:
                if 'lib' in input_many_data:
                    many_forward_reads = input_many_data['lib']['file']
                elif 'handle' in input_many_data:
                    many_forward_reads = input_many_data['handle']
                else:
                    self.log(console,"bad structure for 'many_forward_reads'")
                    raise ValueError("bad structure for 'many_forward_reads'")
                #if 'lib2' in data:
                #    reverse_reads = data['lib2']['file']
                #elif 'handle_2' in data:
                #    reverse_reads = data['handle_2']
                #else:
                #    reverse_reads={}

                ### NOTE: this section is what could be replaced by the transform services
                many_forward_reads_file_path = os.path.join(self.scratch,many_forward_reads['file_name'])
                many_forward_reads_file_handle = open(many_forward_reads_file_path, 'w', 0)
                self.log(console, 'downloading reads file: '+str(many_forward_reads_file_path))
                headers = {'Authorization': 'OAuth '+ctx['token']}
                r = requests.get(many_forward_reads['url']+'/node/'+many_forward_reads['id']+'?download', stream=True, headers=headers)
                for chunk in r.iter_content(1024):
                    appropriate_sequence_found_in_many_input = True
                    many_forward_reads_file_handle.write(chunk)
                many_forward_reads_file_handle.close();
                self.log(console, 'done')
                ### END NOTE


                # remove carriage returns
                new_file_path = many_forward_reads_file_path+"-CRfree"
                new_file_handle = open(new_file_path, 'w', 0)
                many_forward_reads_file_handle = open(many_forward_reads_file_path, 'r', 0)
                for line in many_forward_reads_file_handle:
                    line = re.sub("\r","",line)
                    new_file_handle.write(line)
                many_forward_reads_file_handle.close();
                new_file_handle.close()
                many_forward_reads_file_path = new_file_path


                # convert FASTQ to FASTA (if necessary)
                new_file_path = many_forward_reads_file_path+".fna"
                new_file_handle = open(new_file_path, 'w', 0)
                if many_forward_reads_file_compression == 'gz':
                    many_forward_reads_file_handle = gzip.open(many_forward_reads_file_path, 'r', 0)
                else:
                    many_forward_reads_file_handle = open(many_forward_reads_file_path, 'r', 0)
                header = None
                last_header = None
                last_seq_buf = None
                last_line_was_header = False
                was_fastq = False
                for line in many_forward_reads_file_handle:
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
                many_forward_reads_file_handle.close()
                if was_fastq:
                    many_forward_reads_file_path = new_file_path

            except Exception as e:
                print(traceback.format_exc())
                raise ValueError('Unable to download single-end read library files: ' + str(e))

        # FeatureSet
        #
        elif many_type_name == 'FeatureSet':
            # retrieve sequences for features
            input_many_featureSet = input_many_data
            genome_id_feature_id_delim = '.f:'
            many_forward_reads_file_dir = self.scratch
            many_forward_reads_file = input_many_name+".fasta"

            # DEBUG
            #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            FeatureSetToFASTA_params = {
                'featureSet_ref':      input_many_ref,
                'file':                many_forward_reads_file,
                'dir':                 many_forward_reads_file_dir,
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
            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'])
            FeatureSetToFASTA_retVal = DOTFU.FeatureSetToFASTA (FeatureSetToFASTA_params)
            many_forward_reads_file_path = FeatureSetToFASTA_retVal['fasta_file_path']
            feature_ids_by_genome_ref = FeatureSetToFASTA_retVal['feature_ids_by_genome_ref']
            if len(feature_ids_by_genome_ref.keys()) > 0:
                appropriate_sequence_found_in_many_input = True

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "FeatureSetToFasta() took "+str(end_time-beg_time)+" secs")


        # Genome
        #
        elif many_type_name == 'Genome':
            many_forward_reads_file_dir = self.scratch
            many_forward_reads_file = input_many_name+".fasta"

            # DEBUG
            #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            GenomeToFASTA_params = {
                'genome_ref':          input_many_ref,
                'file':                many_forward_reads_file,
                'dir':                 many_forward_reads_file_dir,
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
            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'])
            GenomeToFASTA_retVal = DOTFU.GenomeToFASTA (GenomeToFASTA_params)
            many_forward_reads_file_path = GenomeToFASTA_retVal['fasta_file_path']
            feature_ids = GenomeToFASTA_retVal['feature_ids']
            if len(feature_ids) > 0:
                appropriate_sequence_found_in_many_input = True

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "Genome2Fasta() took "+str(end_time-beg_time)+" secs")


        # GenomeSet
        #
        elif many_type_name == 'GenomeSet':
            input_many_genomeSet = input_many_data
            genome_id_feature_id_delim = '.f:'
            many_forward_reads_file_dir = self.scratch
            many_forward_reads_file = input_many_name+".fasta"

            # DEBUG
            #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            GenomeSetToFASTA_params = {
                'genomeSet_ref':       input_many_ref,
                'file':                many_forward_reads_file,
                'dir':                 many_forward_reads_file_dir,
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
            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'])
            GenomeSetToFASTA_retVal = DOTFU.GenomeSetToFASTA (GenomeSetToFASTA_params)
            many_forward_reads_file_path = GenomeSetToFASTA_retVal['fasta_file_path_list'][0]
            feature_ids_by_genome_id = GenomeSetToFASTA_retVal['feature_ids_by_genome_id']
            if len(feature_ids_by_genome_id.keys()) > 0:
                appropriate_sequence_found_in_many_input = True

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "FeatureSetToFasta() took "+str(end_time-beg_time)+" secs")


        # Missing proper input_many_type
        #
        else:
            raise ValueError('Cannot yet handle input_many type of: '+type_name)


        # check for failed input file creation
        #
        if not appropriate_sequence_found_in_one_input:
            self.log(invalid_msgs,"no dna sequences found in '"+input_one_name+"'")
        if not appropriate_sequence_found_in_many_input:
            self.log(invalid_msgs,"no dna sequences found in '"+input_many_name+"'")


        # input data failed validation.  Need to return
        #
        if len(invalid_msgs) > 0:

            # load the method provenance from the context object
            #
            self.log(console,"SETTING PROVENANCE")  # DEBUG
            provenance = [{}]
            if 'provenance' in ctx:
                provenance = ctx['provenance']
            # add additional info to provenance here, in this case the input data object reference
            provenance[0]['input_ws_objects'] = []
            provenance[0]['input_ws_objects'].append(input_one_ref)
            provenance[0]['input_ws_objects'].append(input_many_ref)
            provenance[0]['service'] = 'kb_blast'
            provenance[0]['method'] = search_tool_name+'_Search'


            # build output report object
            #
            self.log(console,"BUILDING REPORT")  # DEBUG
            report += "FAILURE:\n\n"+"\n".join(invalid_msgs)+"\n"
            reportObj = {
                'objects_created':[],
                'text_message':report
                }

            reportName = 'blast_report_'+str(uuid.uuid4())
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            report_obj_info = ws.save_objects({
                    #'id':info[6],
                    'workspace':params['workspace_name'],
                    'objects':[
                        {
                        'type':'KBaseReport.Report',
                        'data':reportObj,
                        'name':reportName,
                        'meta':{},
                        'hidden':1,
                        'provenance':provenance  # DEBUG
                        }
                        ]
                    })[0]

            self.log(console,"BUILDING RETURN OBJECT")
            returnVal = { 'report_name': reportName,
                      'report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
                      }
            self.log(console,search_tool_name+"_Search DONE")
            return [returnVal]


        # FORMAT DB
        #
        # OLD SYNTAX: formatdb -i $database -o T -p F -> $database.nsq or $database.00.nsq
        # NEW SYNTAX: makeblastdb -in $database -parse_seqids -dbtype prot/nucl -out <basename>
        makeblastdb_cmd = [self.Make_BLAST_DB]

        # check for necessary files
        if not os.path.isfile(self.Make_BLAST_DB):
            raise ValueError("no such file '"+self.Make_BLAST_DB+"'")
        if not os.path.isfile(many_forward_reads_file_path):
            raise ValueError("no such file '"+many_forward_reads_file_path+"'")
        elif not os.path.getsize(many_forward_reads_file_path) > 0:
            raise ValueError("empty file '"+many_forward_reads_file_path+"'")

        makeblastdb_cmd.append('-in')
        makeblastdb_cmd.append(many_forward_reads_file_path)
        makeblastdb_cmd.append('-parse_seqids')
        makeblastdb_cmd.append('-dbtype')
        makeblastdb_cmd.append('nucl')
        makeblastdb_cmd.append('-out')
        makeblastdb_cmd.append(many_forward_reads_file_path)

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
            line = p.stdout.readline()
            if not line: break
            self.log(console, line.replace('\n', ''))

        p.stdout.close()
        p.wait()
        self.log(console, 'return code: ' + str(p.returncode))
        if p.returncode != 0:
            raise ValueError('Error running makeblastdb, return code: '+str(p.returncode) + 
                '\n\n'+ '\n'.join(console))

        # Check for db output
        if not os.path.isfile(many_forward_reads_file_path+".nsq") and not os.path.isfile(many_forward_reads_file_path+".00.nsq"):
            raise ValueError("makeblastdb failed to create DB file '"+many_forward_reads_file_path+".nsq'")
        elif not os.path.getsize(many_forward_reads_file_path+".nsq") > 0 and not os.path.getsize(many_forward_reads_file_path+".00.nsq") > 0:
            raise ValueError("makeblastdb created empty DB file '"+many_forward_reads_file_path+".nsq'")


        ### Construct the BLAST command
        #
        # OLD SYNTAX: $blast -q $q -G $G -E $E -m $m -e $e_value -v $limit -b $limit -K $limit -p blastn -i $fasta_file -d $database -o $out_file
        # NEW SYNTAX: blastn -query <queryfile> -db <basename> -out <out_aln_file> -outfmt 0/7 (8 became 7) -evalue <e_value> -dust no (DNA) -seg no (AA) -num_threads <num_cores>
        #
        blast_bin = self.BLASTn
        blast_cmd = [blast_bin]

        # check for necessary files
        if not os.path.isfile(blast_bin):
            raise ValueError("no such file '"+blast_bin+"'")
        if not os.path.isfile(one_forward_reads_file_path):
            raise ValueError("no such file '"+one_forward_reads_file_path+"'")
        elif not os.path.getsize(one_forward_reads_file_path) > 0:
            raise ValueError("empty file '"+one_forward_reads_file_path+"'")
        if not os.path.isfile(many_forward_reads_file_path):
            raise ValueError("no such file '"+many_forward_reads_file_path+"'")
        elif not os.path.getsize(many_forward_reads_file_path):
            raise ValueError("empty file '"+many_forward_reads_file_path+"'")

        # set the output path
        timestamp = int((datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()*1000)
        output_dir = os.path.join(self.scratch,'output.'+str(timestamp))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        output_aln_file_path = os.path.join(output_dir, 'alnout.txt');
        output_filtered_fasta_file_path = os.path.join(output_dir, 'output_filtered.fna');

        # this is command for basic search mode
        blast_cmd.append('-query')
        blast_cmd.append(one_forward_reads_file_path)
        blast_cmd.append('-db')
        blast_cmd.append(many_forward_reads_file_path)
        blast_cmd.append('-out')
        blast_cmd.append(output_aln_file_path)
        blast_cmd.append('-outfmt')
        blast_cmd.append('7')
        #blast_cmd.append('0')  # DEBUG
        blast_cmd.append('-evalue')
        blast_cmd.append(str(params['e_value']))

        # options
        if 'maxaccepts' in params:
            if params['maxaccepts']:
                blast_cmd.append('-max_target_seqs')
                blast_cmd.append(str(params['maxaccepts']))

        # Run BLAST, capture output as it happens
        #
        self.log(console, 'RUNNING BLAST:')
        self.log(console, '    '+' '.join(blast_cmd))
#        report += "\n"+'running BLAST:'+"\n"
#        report += '    '+' '.join(blast_cmd)+"\n"

        p = subprocess.Popen(blast_cmd, \
                             cwd = self.scratch, \
                             stdout = subprocess.PIPE, \
                             stderr = subprocess.STDOUT, \
                             shell = False)

        while True:
            line = p.stdout.readline()
            if not line: break
            self.log(console, line.replace('\n', ''))

        p.stdout.close()
        p.wait()
        self.log(console, 'return code: ' + str(p.returncode))
        if p.returncode != 0:
            raise ValueError('Error running BLAST, return code: '+str(p.returncode) + 
                '\n\n'+ '\n'.join(console))


        # get query_len for filtering later
        #
        query_len = 0
        with open(one_forward_reads_file_path, 'r', 0) as query_file_handle:
            for line in query_file_handle:
                if line.startswith('>'):
                    continue
                query_len += len(re.sub(" ","", line.rstrip())) 
        

        # DEBUG
        #one_forward_reads_file_handle = open(one_forward_reads_file_path, 'r', 0)
        #self.log(console, 'reading QUERY reads file: '+str(one_forward_reads_file_path))
        #self.log(console, "".join(one_forward_reads_file_handle.readlines()))
        #one_forward_reads_file_handle.close();

        #many_forward_reads_file_handle = open(many_forward_reads_file_path, 'r', 0)
        #self.log(console, 'reading TARGET reads file: '+str(many_forward_reads_file_path))
        #self.log(console, "".join(many_forward_reads_file_handle.readlines()))

        #in_target = False
        #for line in many_forward_reads_file_handle.readlines():
        #    if line.startswith('>'+'WP_053463618.1'):
        #        in_target = True
        #        self.log(console, line)
        #        continue
        #    elif line.startswith('>'):
        #        if in_target:
        #            in_target = False
        #        continue
        #    elif in_target:
        #        self.log(console, line)
        #many_forward_reads_file_handle.close();


        # Parse the BLAST tabular output and store ids to filter many set to make filtered object to save back to KBase
        #
        self.log(console, 'PARSING BLAST ALIGNMENT OUTPUT')
        if not os.path.isfile(output_aln_file_path):
            raise ValueError("failed to create BLAST output: "+output_aln_file_path)
        elif not os.path.getsize(output_aln_file_path) > 0:
            raise ValueError("created empty file for BLAST output: "+output_aln_file_path)
        hit_seq_ids = dict()
        output_aln_file_handle = open (output_aln_file_path, "r", 0)
        output_aln_buf = output_aln_file_handle.readlines()
        output_aln_file_handle.close()
        hit_total = 0
        high_bitscore_line = dict()
        high_bitscore_score = dict()
        high_bitscore_ident = dict()
        high_bitscore_alnlen = dict()
        hit_order = []
        hit_buf = []
        header_done = False
        for line in output_aln_buf:
            #self.log(console, "HIT_LINE: '"+line+"'")  # DEBUG

            if line.startswith('#'):
                if not header_done:
                    hit_buf.append(line)
                continue
            header_done = True
            hit_info = line.split("\t")
            hit_seq_id     = hit_info[1]
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

            # BLAST SOMETIMES ADDS THIS TO IDs.  NO IDEA WHY, BUT GET RID OF IT!
            if hit_seq_id.startswith('gnl|'):
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

        for hit_seq_id in hit_order:
            hit_buf.append(high_bitscore_line[hit_seq_id])

            if 'ident_thresh' in params and float(params['ident_thresh']) > float(high_bitscore_ident[hit_seq_id]):
                continue
            if 'bitscore' in params and float(params['bitscore']) > float(high_bitscore_score[hit_seq_id]):
                continue
            if 'overlap_fraction' in params and float(params['overlap_fraction']) > float(high_bitscore_alnlen[hit_seq_id])/float(query_len):
                continue
            
            hit_total += 1
            hit_seq_ids[hit_seq_id] = True
            self.log(console, "HIT: '"+hit_seq_id+"'")  # DEBUG
        

        self.log(console, 'EXTRACTING HITS FROM INPUT')
        #self.log(console, 'MANY_TYPE_NAME: '+many_type_name)  # DEBUG


        # SequenceSet input -> SequenceSet output
        #
        if many_type_name == 'SequenceSet':
            seq_total = len(input_many_sequenceSet['sequences'])

            output_sequenceSet = dict()

            if 'sequence_set_id' in input_many_sequenceSet and input_many_sequenceSet['sequence_set_id'] != None:
                output_sequenceSet['sequence_set_id'] = input_many_sequenceSet['sequence_set_id'] + "."+search_tool_name+"_Search_filtered"
            else:
                output_sequenceSet['sequence_set_id'] = search_tool_name+"_Search_filtered"
            if 'description' in input_many_sequenceSet and input_many_sequenceSet['description'] != None:
                output_sequenceSet['description'] = input_many_sequenceSet['description'] + " - "+search_tool_name+"_Search filtered"
            else:
                output_sequenceSet['description'] = search_tool_name+"_Search filtered"

            self.log(console,"ADDING SEQUENCES TO SEQUENCESET")
            output_sequenceSet['sequences'] = []

            for seq_obj in input_many_sequenceSet['sequences']:
                header_id = seq_obj['sequence_id']
                #header_desc = seq_obj['description']
                #sequence_str = seq_obj['sequence']

                try:
                    #self.log(console,"checking '"+header_id+"'")
                    id_trans = re.sub ('|',':',header_id)  # BLAST seems to make this translation now
                    in_filtered_set = hit_seq_ids[id_trans]
                    #self.log(console, 'FOUND HIT '+header_id)  # DEBUG
                    output_sequenceSet['sequences'].append(seq_obj)
                except:
                    pass

        # SingleEndLibrary input -> SingleEndLibrary output
        #
        elif many_type_name == 'SingleEndLibrary':

            #  Note: don't use SeqIO.parse because loads everything into memory
            #
#            with open(many_forward_reads_file_path, 'r', -1) as many_forward_reads_file_handle, open(output_filtered_fasta_file_path, 'w', -1) as output_filtered_fasta_file_handle:
            output_filtered_fasta_file_handle = open(output_filtered_fasta_file_path, 'w', -1)
            if many_forward_reads_file_compression == 'gz':
                many_forward_reads_file_handle = gzip.open(many_forward_reads_file_path, 'r', -1)
            else:
                many_forward_reads_file_handle = open(many_forward_reads_file_path, 'r', -1)

            seq_total = 0;
            filtered_seq_total = 0
            last_seq_buf = []
            last_seq_id = None
            last_header = None
            pattern = re.compile('^\S*')
            for line in many_forward_reads_file_handle:
                if line.startswith('>'):
                    #self.log(console, 'LINE: '+line)  # DEBUG
                    seq_total += 1
                    seq_id = line[1:]  # removes '>'
                    seq_id = pattern.findall(seq_id)[0]

                    if last_seq_id != None:
                        #self.log(console, 'ID: '+last_seq_id)  # DEBUG
                        try:
                            id_trans = re.sub('|',':',last_seq_id)  # BLAST seems to make this translation now
                            in_filtered_set = hit_seq_ids[id_trans]
                            #self.log(console, 'FOUND HIT '+last_seq_id)  # DEBUG
                            filtered_seq_total += 1
                            output_filtered_fasta_file_handle.write(last_header)
                            output_filtered_fasta_file_handle.writelines(last_seq_buf)
                        except:
                            pass
                        
                    last_seq_buf = []
                    last_seq_id = seq_id
                    last_header = line
                else:
                    last_seq_buf.append(line)

            if last_seq_id != None:
                #self.log(console, 'ID: '+last_seq_id)  # DEBUG
                try:
                    id_trans = re.sub ('|',':',last_seq_id)  # BLAST seems to make this translation now
                    in_filtered_set = hit_seq_ids[id_trans]
                    #self.log(console, 'FOUND HIT: '+last_seq_id)  # DEBUG
                    filtered_seq_total += 1
                    output_filtered_fasta_file_handle.write(last_header)
                    output_filtered_fasta_file_handle.writelines(last_seq_buf)
                except:
                    pass
                
            last_seq_buf = []
            last_seq_id = None
            last_header = None

            many_forward_reads_file_handle.close()
            output_filtered_fasta_file_handle.close()

            if filtered_seq_total != hit_total:
                self.log(console,'hits in BLAST alignment output '+str(hit_total)+' != '+str(filtered_seq_total)+' matched sequences in input file')
                raise ValueError('hits in BLAST alignment output '+str(hit_total)+' != '+str(filtered_seq_total)+' matched sequences in input file')


        # FeatureSet input -> FeatureSet output
        #
        elif many_type_name == 'FeatureSet':
            seq_total = len(input_many_featureSet['elements'].keys())

            output_featureSet = dict()
            if 'description' in input_many_featureSet and input_many_featureSet['description'] != None:
                output_featureSet['description'] = input_many_featureSet['description'] + " - "+search_tool_name+"_Search filtered"
            else:
                output_featureSet['description'] = search_tool_name+"_Search filtered"
            #output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()

            fId_list = input_many_featureSet['elements'].keys()
            self.log(console,"ADDING FEATURES TO FEATURESET")
            for fId in sorted(fId_list):
                for genome_ref in input_many_featureSet['elements'][fId]:
                    try:
                        #self.log(console,"checking '"+fId+"'")
                        #in_filtered_set = hit_seq_ids[fId]
                        id_trans = re.sub ('|',':',genome_ref+genome_id_feature_id_delim+fId)  # BLAST seems to make this translation now
                        in_filtered_set = hit_seq_ids[id_trans]
                        #self.log(console, 'FOUND HIT '+fId)  # DEBUG
                        try:
                            this_genome_ref_list = output_featureSet['elements'][fId]
                        except:
                            output_featureSet['elements'][fId] = []
                        output_featureSet['elements'][fId].append(genome_ref)
                    except:
                        pass

        # Parse Genome hits into FeatureSet
        #
        elif many_type_name == 'Genome':
            seq_total = 0
            output_featureSet = dict()
#            if 'scientific_name' in input_many_genome and input_many_genome['scientific_name'] != None:
#                output_featureSet['description'] = input_many_genome['scientific_name'] + " - "+search_tool_name+"_Search filtered"
#            else:
#                output_featureSet['description'] = search_tool_name+"_Search filtered"
            output_featureSet['description'] = search_tool_name+"_Search filtered"
            #output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()
            for fid in feature_ids:
                seq_total += 1
                try:
                    id_trans = re.sub ('|',':',fid)  # BLAST seems to make this translation now
                    in_filtered_set = hit_seq_ids[id_trans]
                    #output_featureSet['element_ordering'].append(fid)
                    output_featureSet['elements'][fid] = [input_many_ref]
                except:
                    pass

        # Parse GenomeSet hits into FeatureSet
        #
        elif many_type_name == 'GenomeSet':
            seq_total = 0

            output_featureSet = dict()
            if 'description' in input_many_genomeSet and input_many_genomeSet['description'] != None:
                output_featureSet['description'] = input_many_genomeSet['description'] + " - "+search_tool_name+"_Search filtered"
            else:
                output_featureSet['description'] = search_tool_name+"_Search filtered"
            #output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()

            self.log(console,"READING HITS FOR GENOMES")  # DEBUG
            for genome_id in feature_ids_by_genome_id.keys():
                self.log(console,"READING HITS FOR GENOME "+genome_id)  # DEBUG
                genome_ref = input_many_genomeSet['elements'][genome_id]['ref']
                for feature_id in feature_ids_by_genome_id[genome_id]:
                    seq_total += 1
                    try:
                        #in_filtered_set = hit_seq_ids[feature['id']]
                        id_trans = re.sub('|',':',genome_ref+genome_id_feature_id_delim+feature_id)
                        in_filtered_set = hit_seq_ids[id_trans]
                        #in_filtered_set = hit_seq_ids[feature_id]
                        #self.log(console, 'FOUND HIT: '+feature['id'])  # DEBUG
                        #output_featureSet['element_ordering'].append(feature['id'])
                        try:
                            this_genome_ref_list = output_featureSet['elements'][feature_id]
                        except:
                            output_featureSet['elements'][feature_id] = []
                        output_featureSet['elements'][feature_id].append(genome_ref)
                    except:
                        pass


        # load the method provenance from the context object
        #
        self.log(console,"SETTING PROVENANCE")  # DEBUG
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects'] = []
        provenance[0]['input_ws_objects'].append(input_one_ref)
        provenance[0]['input_ws_objects'].append(input_many_ref)
        provenance[0]['service'] = 'kb_blast'
        provenance[0]['method'] = search_tool_name+'_Search'


        # Upload results
        #
        if len(invalid_msgs) == 0:
            self.log(console,"UPLOADING RESULTS")  # DEBUG

            # input many SingleEndLibrary -> upload SingleEndLibrary
            #
            if many_type_name == 'SingleEndLibrary':
            
                self.upload_SingleEndLibrary_to_shock_and_ws (ctx,
                                                          console,  # DEBUG
                                                          params['workspace_name'],
                                                          params['output_filtered_name'],
                                                          output_filtered_fasta_file_path,
                                                          provenance,
                                                          sequencing_tech
                                                         )

            # input many SequenceSet -> save SequenceSet
            #
            elif many_type_name == 'SequenceSet':
                new_obj_info = ws.save_objects({
                            'workspace': params['workspace_name'],
                            'objects':[{
                                    'type': 'KBaseSequences.SequenceSet',
                                    'data': output_sequenceSet,
                                    'name': params['output_filtered_name'],
                                    'meta': {},
                                    'provenance': provenance
                                }]
                        })[0]

            else:  # input many FeatureSet, Genome, and GenomeSet -> upload FeatureSet output
                new_obj_info = ws.save_objects({
                            'workspace': params['workspace_name'],
                            'objects':[{
                                    'type': 'KBaseCollections.FeatureSet',
                                    'data': output_featureSet,
                                    'name': params['output_filtered_name'],
                                    'meta': {},
                                    'provenance': provenance
                                }]
                        })[0]

        # build output report object
        #
        self.log(console,"BUILDING REPORT")  # DEBUG
        if len(invalid_msgs) == 0:
            report += 'sequences in many set: '+str(seq_total)+"\n"
            report += 'sequences in hit set:  '+str(hit_total)+"\n"
            report += "\n"
            for line in hit_buf:
                report += line
            reportObj = {
                'objects_created':[{'ref':params['workspace_name']+'/'+params['output_filtered_name'], 'description':search_tool_name+'_Search hits'}],
                'text_message':report
                }
        else:
            report += "FAILURE\n\n"+"\n".join(invalid_msgs)+"\n"
            reportObj = {
                'objects_created':[],
                'text_message':report
                }

        reportName = 'blast_report_'+str(uuid.uuid4())
        report_obj_info = ws.save_objects({
#                'id':info[6],
                'workspace':params['workspace_name'],
                'objects':[
                    {
                        'type':'KBaseReport.Report',
                        'data':reportObj,
                        'name':reportName,
                        'meta':{},
                        'hidden':1,
                        'provenance':provenance
                    }
                ]
            })[0]

        self.log(console,"BUILDING RETURN OBJECT")
#        returnVal = { 'output_report_name': reportName,
#                      'output_report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
#                      'output_filtered_ref': params['workspace_name']+'/'+params['output_filtered_name']
#                      }
        returnVal = { 'report_name': reportName,
                      'report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
                      }
        self.log(console,search_tool_name+"_Search DONE")
        #END BLASTn_Search

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method BLASTn_Search return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def BLASTp_Search(self, ctx, params):
        """
        :param params: instance of type "BLAST_Params" (BLAST Input Params)
           -> structure: parameter "workspace_name" of type "workspace_name"
           (** The workspace object refs are of form: ** **    objects =
           ws.get_objects([{'ref':
           params['workspace_id']+'/'+params['obj_name']}]) ** ** "ref" means
           the entire name combining the workspace id and the object name **
           "id" is a numerical identifier of the workspace or object, and
           should just be used for workspace ** "name" is a string identifier
           of a workspace or object.  This is received from Narrative.),
           parameter "input_one_sequence" of type "sequence", parameter
           "input_one_ref" of type "data_obj_ref", parameter "input_many_ref"
           of type "data_obj_ref", parameter "input_msa_ref" of type
           "data_obj_ref", parameter "output_one_name" of type
           "data_obj_name", parameter "output_filtered_name" of type
           "data_obj_name", parameter "ident_thresh" of Double, parameter
           "e_value" of Double, parameter "bitscore" of Double, parameter
           "overlap_fraction" of Double, parameter "maxaccepts" of Double,
           parameter "rounds" of Double
        :returns: instance of type "BLAST_Output" (BLAST Output) ->
           structure: parameter "report_name" of type "data_obj_name",
           parameter "report_ref" of type "data_obj_ref"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN BLASTp_Search
        console = []
        invalid_msgs = []
        search_tool_name = 'BLASTp'
        self.log(console,'Running '+search_tool_name+'_Search with params=')
        self.log(console, "\n"+pformat(params))
        report = ''
#        report = 'Running '+search_tool_name+'_Search with params='
#        report += "\n"+pformat(params)
        appropriate_sequence_found_in_one_input = False
        appropriate_sequence_found_in_many_input = False


        #### do some basic checks
        #
        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
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
#        if 'input_one_sequence' not in params:
#            raise ValueError('input_one_sequence parameter is required')
#        if 'input_one_ref' not in params:
#            raise ValueError('input_one_ref parameter is required')
        if 'input_many_ref' not in params:
            raise ValueError('input_many_ref parameter is required')
        if 'output_filtered_name' not in params:
            raise ValueError('output_filtered_name parameter is required')


        # set local names
        input_one_ref  = None
        output_one_ref = None
        input_many_ref = params['input_many_ref']
        

        # Write the input_one_sequence to file
        #
        if 'input_one_sequence' in params \
                and params['input_one_sequence'] != None \
                and params['input_one_sequence'] != "Optionally enter PROTEIN sequence...":

            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'])
            ParseFastaStr_retVal = DOTFU.ParseFastaStr ({
                'fasta_str':     params['input_one_sequence'],
                'residue_type': 'PROT',
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
            if 'provenance' in ctx:
                provenance = ctx['provenance']
            # add additional info to provenance here, in this case the input data object reference
            provenance[0]['input_ws_objects'] = []
            provenance[0]['service'] = 'kb_blast'
            provenance[0]['method'] = search_tool_name+'_Search'
                
            # Upload results
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
                ws = workspaceService(self.workspaceURL, token=ctx['token'])
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


        #### Get the input_one object
        ##
        if 'output_one_name' in params and output_one_ref:
            input_one_ref = output_one_ref
        else:
            input_one_ref = params['input_one_ref']

        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            objects = ws.get_objects([{'ref': input_one_ref}])
            input_one_data = objects[0]['data']
            input_one_name = str(objects[0]['info'][1])
            info = objects[0]['info']
            
            one_type_name = info[2].split('.')[1].split('-')[0]
        except Exception as e:
            raise ValueError('Unable to fetch input_one_ref object from workspace: ' + str(e))
            #to get the full stack trace: traceback.format_exc()

        if 'input_one_sequence' in params \
                and params['input_one_sequence'] != None \
                and params['input_one_sequence'] != "Optionally enter PROTEIN sequence..." \
                and one_type_name != 'SequenceSet':

            self.log(invalid_msgs,"ERROR: Mismatched input type for Query Object: "+input_one_ref+" should be SequenceSet instead of: "+one_type_name)


        # Handle overloading (input_one can be SequenceSet, Feature, or FeatureSet)
        #
        if one_type_name == 'SequenceSet':
            try:
                input_one_sequenceSet = input_one_data
            except Exception as e:
                print(traceback.format_exc())
                raise ValueError('Unable to get sequenceSet object: ' + str(e))

            header_id = input_one_sequenceSet['sequences'][0]['sequence_id']
            sequence_str = input_one_data['sequences'][0]['sequence']

            PROT_pattern = re.compile("^[acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWYxX ]+$")
            #DNA_pattern = re.compile("^[acgtuACGTUnryNRY ]+$")
            if not PROT_pattern.match(sequence_str):
                self.log(invalid_msgs,"BAD record for sequence_id: "+header_id+"\n"+sequence_str+"\n")
            else:
                appropriate_sequence_found_in_one_input = True

            one_forward_reads_file_path = os.path.join(self.scratch, header_id+'.fasta')
            one_forward_reads_file_handle = open(one_forward_reads_file_path, 'w', 0)
            self.log(console, 'writing reads file: '+str(one_forward_reads_file_path))
            one_forward_reads_file_handle.write('>'+header_id+"\n")
            one_forward_reads_file_handle.write(sequence_str+"\n")
            one_forward_reads_file_handle.close();
            self.log(console, 'done')

        # FeatureSet
        #
        elif one_type_name == 'FeatureSet':
            # retrieve sequences for features
            #input_one_featureSet = input_one_data
            genome_id_feature_id_delim = '.f:'
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
                    'residue_type':        'protein',
                    'feature_type':        'CDS',
                    'record_id_pattern':   '%%genome_ref%%'+genome_id_feature_id_delim+'%%feature_id%%',
                    'record_desc_pattern': '[%%genome_ref%%]',
                    'case':                'upper',
                    'linewrap':            50,
                    'merge_fasta_files':   'TRUE'
                    }

            #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'])
            FeatureSetToFASTA_retVal = DOTFU.FeatureSetToFASTA (FeatureSetToFASTA_params)
            one_forward_reads_file_path = FeatureSetToFASTA_retVal['fasta_file_path']
            if len(FeatureSetToFASTA_retVal['feature_ids_by_genome_ref'].keys()) > 0:
                appropriate_sequence_found_in_one_input = True

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "FeatureSetToFasta() took "+str(end_time-beg_time)+" secs")


        # Feature
        #
        elif one_type_name == 'Feature':
            # export feature to FASTA file
            feature = input_one_data
            one_forward_reads_file_path = os.path.join(self.scratch, input_one_name+".fasta")
            self.log(console, 'writing fasta file: '+one_forward_reads_file_path)
            # BLASTp is prot-prot
            #record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description='['+feature['genome_id']+']'+' '+feature['function'])
            #if feature['type'] != 'CDS':
            #    self.log(console,params['input_one_ref']+" feature type must be CDS")
            #    self.log(invalid_msgs,params['input_one_ref']+" feature type must be CDS")
            if 'protein_translation' not in feature or feature['protein_translation'] == None:
                #self.log(console,"bad CDS Feature "+params['input_one_ref']+": no protein_translation found")
                #raise ValueError("bad CDS Feature "+params['input_one_ref']+": no protein_translation found")
                self.log(console,params['input_one_ref']+" feature type must be CDS")
                self.log(invalid_msgs,params['input_one_ref']+" feature type must be CDS")
            else:
                appropriate_sequence_found_in_one_input = True
                record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description='['+feature['genome_id']+']'+' '+feature['function'])
                SeqIO.write([record], one_forward_reads_file_path, "fasta")
                appropriate_sequence_found_in_one_input = True
        else:
            raise ValueError('Cannot yet handle input_one type of: '+type_name)            

        #### Get the input_many object
        ##
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            objects = ws.get_objects([{'ref': input_many_ref}])
            input_many_data = objects[0]['data']
            info = objects[0]['info']
            input_many_name = str(info[1])
            many_type_name = info[2].split('.')[1].split('-')[0]

        except Exception as e:
            raise ValueError('Unable to fetch input_many_name object from workspace: ' + str(e))
            #to get the full stack trace: traceback.format_exc()


        # Handle overloading (input_many can be SequenceSet, FeatureSet, Genome, or GenomeSet)
        #
        if many_type_name == 'SequenceSet':
            try:
                input_many_sequenceSet = input_many_data
            except Exception as e:
                print(traceback.format_exc())
                raise ValueError('Unable to get SequenceSet: ' + str(e))

            header_id = input_many_sequenceSet['sequences'][0]['sequence_id']
            many_forward_reads_file_path = os.path.join(self.scratch, header_id+'.fasta')
            many_forward_reads_file_handle = open(many_forward_reads_file_path, 'w', 0)
            self.log(console, 'writing reads file: '+str(many_forward_reads_file_path))

            for seq_obj in input_many_sequenceSet['sequences']:
                header_id = seq_obj['sequence_id']
                sequence_str = seq_obj['sequence']

                PROT_pattern = re.compile("^[acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWYxX ]+$")
                #DNA_pattern = re.compile("^[acgtuACGTUnryNRY ]+$")
                if not PROT_pattern.match(sequence_str):
                    self.log(invalid_msgs,"BAD record for sequence_id: "+header_id+"\n"+sequence_str+"\n")
                    continue
                appropriate_sequence_found_in_many_input = True
                many_forward_reads_file_handle.write('>'+header_id+"\n")
                many_forward_reads_file_handle.write(sequence_str+"\n")
            many_forward_reads_file_handle.close();
            self.log(console, 'done')


        # FeatureSet
        #
        elif many_type_name == 'FeatureSet':
            # retrieve sequences for features
            input_many_featureSet = input_many_data
            genome_id_feature_id_delim = '.f:'
            many_forward_reads_file_dir = self.scratch
            many_forward_reads_file = input_many_name+".fasta"

            # DEBUG
            #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            FeatureSetToFASTA_params = {
                'featureSet_ref':      input_many_ref,
                'file':                many_forward_reads_file,
                'dir':                 many_forward_reads_file_dir,
                'console':             console,
                'invalid_msgs':        invalid_msgs,
                'residue_type':        'protein',
                'feature_type':        'CDS',
                'record_id_pattern':   '%%genome_ref%%'+genome_id_feature_id_delim+'%%feature_id%%',
                'record_desc_pattern': '[%%genome_ref%%]',
                'case':                'upper',
                'linewrap':            50,
                'merge_fasta_files':   'TRUE'
                }

            #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'])
            FeatureSetToFASTA_retVal = DOTFU.FeatureSetToFASTA (FeatureSetToFASTA_params)
            many_forward_reads_file_path = FeatureSetToFASTA_retVal['fasta_file_path']
            feature_ids_by_genome_ref = FeatureSetToFASTA_retVal['feature_ids_by_genome_ref']
            if len(feature_ids_by_genome_ref.keys()) > 0:
                appropriate_sequence_found_in_many_input = True

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "FeatureSetToFasta() took "+str(end_time-beg_time)+" secs")


        # Genome
        #
        elif many_type_name == 'Genome':
            many_forward_reads_file_dir = self.scratch
            many_forward_reads_file = input_many_name+".fasta"

            # DEBUG
            #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            GenomeToFASTA_params = {
                'genome_ref':          input_many_ref,
                'file':                many_forward_reads_file,
                'dir':                 many_forward_reads_file_dir,
                'console':             console,
                'invalid_msgs':        invalid_msgs,
                'residue_type':        'protein',
                'feature_type':        'CDS',
                'record_id_pattern':   '%%feature_id%%',
                'record_desc_pattern': '[%%genome_ref%%]',
                'case':                'upper',
                'linewrap':            50
                }

            #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'])
            GenomeToFASTA_retVal = DOTFU.GenomeToFASTA (GenomeToFASTA_params)
            many_forward_reads_file_path = GenomeToFASTA_retVal['fasta_file_path']
            feature_ids = GenomeToFASTA_retVal['feature_ids']
            if len(feature_ids) > 0:
                appropriate_sequence_found_in_many_input = True

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "Genome2Fasta() took "+str(end_time-beg_time)+" secs")


        # GenomeSet
        #
        elif many_type_name == 'GenomeSet':
            input_many_genomeSet = input_many_data
            genome_id_feature_id_delim = '.f:'
            many_forward_reads_file_dir = self.scratch
            many_forward_reads_file = input_many_name+".fasta"

            # DEBUG
            #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            GenomeSetToFASTA_params = {
                'genomeSet_ref':       input_many_ref,
                'file':                many_forward_reads_file,
                'dir':                 many_forward_reads_file_dir,
                'console':             console,
                'invalid_msgs':        invalid_msgs,
                'residue_type':        'protein',
                'feature_type':        'CDS',
                'record_id_pattern':   '%%genome_ref%%'+genome_id_feature_id_delim+'%%feature_id%%',
                'record_desc_pattern': '[%%genome_ref%%]',
                'case':                'upper',
                'linewrap':            50,
                'merge_fasta_files':   'TRUE'
                }

            #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'])
            GenomeSetToFASTA_retVal = DOTFU.GenomeSetToFASTA (GenomeSetToFASTA_params)
            many_forward_reads_file_path = GenomeSetToFASTA_retVal['fasta_file_path_list'][0]
            feature_ids_by_genome_id = GenomeSetToFASTA_retVal['feature_ids_by_genome_id']
            if len(feature_ids_by_genome_id.keys()) > 0:
                appropriate_sequence_found_in_many_input = True

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "FeatureSetToFasta() took "+str(end_time-beg_time)+" secs")


        # Missing proper input_many_type
        #
        else:
            raise ValueError('Cannot yet handle input_many type of: '+type_name)            


        # check for failed input file creation
        #
        if not appropriate_sequence_found_in_one_input:
            self.log(invalid_msgs,"no protein sequences found in '"+input_one_name+"'")
        if not appropriate_sequence_found_in_many_input:
            self.log(invalid_msgs,"no protein sequences found in '"+input_many_name+"'")


        # input data failed validation.  Need to return
        #
        if len(invalid_msgs) > 0:

            # load the method provenance from the context object
            #
            self.log(console,"SETTING PROVENANCE")  # DEBUG
            provenance = [{}]
            if 'provenance' in ctx:
                provenance = ctx['provenance']
            # add additional info to provenance here, in this case the input data object reference
            provenance[0]['input_ws_objects'] = []
            provenance[0]['input_ws_objects'].append(input_one_ref)
            provenance[0]['input_ws_objects'].append(input_many_ref)
            provenance[0]['service'] = 'kb_blast'
            provenance[0]['method'] = search_tool_name+'_Search'


            # build output report object
            #
            self.log(console,"BUILDING REPORT")  # DEBUG
            report += "FAILURE:\n\n"+"\n".join(invalid_msgs)+"\n"
            reportObj = {
                'objects_created':[],
                'text_message':report
                }

            reportName = 'blast_report_'+str(uuid.uuid4())
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            report_obj_info = ws.save_objects({
                    #'id':info[6],
                    'workspace':params['workspace_name'],
                    'objects':[
                        {
                        'type':'KBaseReport.Report',
                        'data':reportObj,
                        'name':reportName,
                        'meta':{},
                        'hidden':1,
                        'provenance':provenance  # DEBUG
                        }
                        ]
                    })[0]

            self.log(console,"BUILDING RETURN OBJECT")
            returnVal = { 'report_name': reportName,
                      'report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
                      }
            self.log(console,search_tool_name+"_Search DONE")
            return [returnVal]


        # FORMAT DB
        #
        # OLD SYNTAX: formatdb -i $database -o T -p F -> $database.nsq or $database.00.nsq
        # NEW SYNTAX: makeblastdb -in $database -parse_seqids -dbtype prot/nucl -out <basename>
        makeblastdb_cmd = [self.Make_BLAST_DB]

        # check for necessary files
        if not os.path.isfile(self.Make_BLAST_DB):
            raise ValueError("no such file '"+self.Make_BLAST_DB+"'")
        if not os.path.isfile(many_forward_reads_file_path):
            raise ValueError("no such file '"+many_forward_reads_file_path+"'")
        elif not os.path.getsize(many_forward_reads_file_path) > 0:
            raise ValueError("empty file '"+many_forward_reads_file_path+"'")

        makeblastdb_cmd.append('-in')
        makeblastdb_cmd.append(many_forward_reads_file_path)
        makeblastdb_cmd.append('-parse_seqids')
        makeblastdb_cmd.append('-dbtype')
        makeblastdb_cmd.append('prot')
        makeblastdb_cmd.append('-out')
        makeblastdb_cmd.append(many_forward_reads_file_path)

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
            line = p.stdout.readline()
            if not line: break
            self.log(console, line.replace('\n', ''))

        p.stdout.close()
        p.wait()
        self.log(console, 'return code: ' + str(p.returncode))
        if p.returncode != 0:
            raise ValueError('Error running makeblastdb, return code: '+str(p.returncode) + 
                '\n\n'+ '\n'.join(console))

        # Check for db output
        if not os.path.isfile(many_forward_reads_file_path+".psq") and not os.path.isfile(many_forward_reads_file_path+".00.psq"):
            raise ValueError("makeblastdb failed to create DB file '"+many_forward_reads_file_path+".psq'")
        elif not os.path.getsize(many_forward_reads_file_path+".psq") > 0 and not os.path.getsize(many_forward_reads_file_path+".00.psq") > 0:
            raise ValueError("makeblastdb created empty DB file '"+many_forward_reads_file_path+".psq'")


        ### Construct the BLAST command
        #
        # OLD SYNTAX: $blast -q $q -G $G -E $E -m $m -e $e_value -v $limit -b $limit -K $limit -p blastp -i $fasta_file -d $database -o $out_file
        # NEW SYNTAX: blastp -query <queryfile> -db <basename> -out <out_aln_file> -outfmt 0/7 (8 became 7) -evalue <e_value> -dust no (DNA) -seg no (AA) -num_threads <num_cores>
        #
        blast_bin = self.BLASTp
        blast_cmd = [blast_bin]

        # check for necessary files
        if not os.path.isfile(blast_bin):
            raise ValueError("no such file '"+blast_bin+"'")
        if not os.path.isfile(one_forward_reads_file_path):
            raise ValueError("no such file '"+one_forward_reads_file_path+"'")
        elif not os.path.getsize(one_forward_reads_file_path) > 0:
            raise ValueError("empty file '"+one_forward_reads_file_path+"'")
        if not os.path.isfile(many_forward_reads_file_path):
            raise ValueError("no such file '"+many_forward_reads_file_path+"'")
        elif not os.path.getsize(many_forward_reads_file_path):
            raise ValueError("empty file '"+many_forward_reads_file_path+"'")

        # set the output path
        timestamp = int((datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()*1000)
        output_dir = os.path.join(self.scratch,'output.'+str(timestamp))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        output_aln_file_path = os.path.join(output_dir, 'alnout.txt');
        output_filtered_fasta_file_path = os.path.join(output_dir, 'output_filtered.fna');

        # this is command for basic search mode
        blast_cmd.append('-query')
        blast_cmd.append(one_forward_reads_file_path)
        blast_cmd.append('-db')
        blast_cmd.append(many_forward_reads_file_path)
        blast_cmd.append('-out')
        blast_cmd.append(output_aln_file_path)
        blast_cmd.append('-outfmt')
        blast_cmd.append('7')
        blast_cmd.append('-evalue')
        blast_cmd.append(str(params['e_value']))

        # options
        if 'maxaccepts' in params:
            if params['maxaccepts']:
                blast_cmd.append('-max_target_seqs')
                blast_cmd.append(str(params['maxaccepts']))

        # Run BLAST, capture output as it happens
        #
        self.log(console, 'RUNNING BLAST:')
        self.log(console, '    '+' '.join(blast_cmd))
#        report += "\n"+'running BLAST:'+"\n"
#        report += '    '+' '.join(blast_cmd)+"\n"

        p = subprocess.Popen(blast_cmd, \
                             cwd = self.scratch, \
                             stdout = subprocess.PIPE, \
                             stderr = subprocess.STDOUT, \
                             shell = False)

        while True:
            line = p.stdout.readline()
            if not line: break
            self.log(console, line.replace('\n', ''))

        p.stdout.close()
        p.wait()
        self.log(console, 'return code: ' + str(p.returncode))
        if p.returncode != 0:
            raise ValueError('Error running BLAST, return code: '+str(p.returncode) + 
                '\n\n'+ '\n'.join(console))


        # get query_len for filtering later
        #
        query_len = 0
        with open(one_forward_reads_file_path, 'r', 0) as query_file_handle:
            for line in query_file_handle:
                if line.startswith('>'):
                    continue
                query_len += len(re.sub(r" ","", line.rstrip())) 
        

        # Parse the BLAST tabular output and store ids to filter many set to make filtered object to save back to KBase
        #
        self.log(console, 'PARSING BLAST ALIGNMENT OUTPUT')
        if not os.path.isfile(output_aln_file_path):
            raise ValueError("failed to create BLAST output: "+output_aln_file_path)
        elif not os.path.getsize(output_aln_file_path) > 0:
            raise ValueError("created empty file for BLAST output: "+output_aln_file_path)
        hit_seq_ids = dict()
        output_aln_file_handle = open (output_aln_file_path, "r", 0)
        output_aln_buf = output_aln_file_handle.readlines()
        output_aln_file_handle.close()
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
            #self.log(console,'HIT LINE: '+line)  # DEBUG
            hit_info = line.split("\t")
            hit_seq_id     = hit_info[1]
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

            # BLAST SOMETIMES ADDS THIS TO IDs.  NO IDEA WHY, BUT GET RID OF IT!
            if hit_seq_id.startswith('gnl|'):
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

        for hit_seq_id in hit_order:
            hit_buf.append(high_bitscore_line[hit_seq_id])

            #self.log(console,"HIT_SEQ_ID: '"+hit_seq_id+"'")
            if 'ident_thresh' in params and float(params['ident_thresh']) > float(high_bitscore_ident[hit_seq_id]):
                continue
            if 'bitscore' in params and float(params['bitscore']) > float(high_bitscore_score[hit_seq_id]):
                continue
            if 'overlap_fraction' in params and float(params['overlap_fraction']) > float(high_bitscore_alnlen[hit_seq_id])/float(query_len):
                continue
            
            hit_total += 1
            hit_seq_ids[hit_seq_id] = True
            self.log(console, "HIT: '"+hit_seq_id+"'")  # DEBUG
        

        # DEBUG
        for hit_seq_id in hit_seq_ids.keys():
            self.log (console, "HIT_ID: '"+str(hit_seq_id)+"'")

        self.log(console, 'EXTRACTING HITS FROM INPUT')
        self.log(console, 'MANY_TYPE_NAME: '+many_type_name)  # DEBUG


        # SequenceSet input -> SequenceSet output
        #
        if many_type_name == 'SequenceSet':
            seq_total = len(input_many_sequenceSet['sequences'])

            output_sequenceSet = dict()

            if 'sequence_set_id' in input_many_sequenceSet and input_many_sequenceSet['sequence_set_id'] != None:
                output_sequenceSet['sequence_set_id'] = input_many_sequenceSet['sequence_set_id'] + "."+search_tool_name+"_Search_filtered"
            else:
                output_sequenceSet['sequence_set_id'] = search_tool_name+"_Search_filtered"
            if 'description' in input_many_sequenceSet and input_many_sequenceSet['description'] != None:
                output_sequenceSet['description'] = input_many_sequenceSet['description'] + " - "+search_tool_name+"_Search filtered"
            else:
                output_sequenceSet['description'] = search_tool_anme+"_Search filtered"

            self.log(console,"ADDING SEQUENCES TO SEQUENCESET")
            output_sequenceSet['sequences'] = []

            for seq_obj in input_many_sequenceSet['sequences']:
                header_id = seq_obj['sequence_id']
                #header_desc = seq_obj['description']
                #sequence_str = seq_obj['sequence']

                try:
                    #self.log(console,"checking '"+header_id+"'")
                    id_trans = re.sub ('|',':',header_id)
                    in_filtered_set = hit_seq_ids[id_trans]
                    #self.log(console, 'FOUND HIT '+header_id)  # DEBUG
                    output_sequenceSet['sequences'].append(seq_obj)
                except:
                    pass


        # FeatureSet input -> FeatureSet output
        #
        elif many_type_name == 'FeatureSet':
            seq_total = len(input_many_featureSet['elements'].keys())

            output_featureSet = dict()
            if 'description' in input_many_featureSet and input_many_featureSet['description'] != None:
                output_featureSet['description'] = input_many_featureSet['description'] + " - "+search_tool_name+"_Search filtered"
            else:
                output_featureSet['description'] = search_tool_name+"_Search filtered"
            #output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()

            fId_list = input_many_featureSet['elements'].keys()
            self.log(console,"ADDING FEATURES TO FEATURESET")
            for fId in sorted(fId_list):
                for genome_ref in input_many_featureSet['elements'][fId]:
                    try:
                        #self.log(console,"checking '"+fId+"'")
                        #in_filtered_set = hit_seq_ids[fId]
                        id_trans = re.sub ('|',':',genome_ref+genome_id_feature_id_delim+fId)
                        in_filtered_set = hit_seq_ids[id_trans]
                        #self.log(console, 'FOUND HIT '+fId)  # DEBUG
                        try:
                            this_genome_ref_list = output_featureSet['elements'][fId]
                        except:
                            output_featureSet['elements'][fId] = []
                        output_featureSet['elements'][fId].append(genome_ref)
                    except:
                        pass

        # Parse Genome hits into FeatureSet
        #
        elif many_type_name == 'Genome':
            seq_total = 0
            output_featureSet = dict()
#            if 'scientific_name' in input_many_genome and input_many_genome['scientific_name'] != None:
#                output_featureSet['description'] = input_many_genome['scientific_name'] + " - "+search_tool_name+"_Search filtered"
#            else:
#                output_featureSet['description'] = search_tool_name+"_Search filtered"
            output_featureSet['description'] = search_tool_name+"_Search filtered"
            #output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()
            for fid in feature_ids:
                seq_total += 1
                try:
                    id_trans = re.sub ('|',':',fid)  # BLAST seems to make this translation now
                    self.log(console, "FID: '"+str(fid)+"'")  # DEBUG
                    self.log(console, "IDT: '"+str(id_trans)+"'")  # DEBUG

                    in_filtered_set = hit_seq_ids[id_trans]
                    #output_featureSet['element_ordering'].append(fid)
                    output_featureSet['elements'][fid] = [input_many_ref]
                    self.log(console, "GOT ONE: '"+str(fid)+"'")  # DEBUG
                except:
                    pass

        # Parse GenomeSet hits into FeatureSet
        #
        elif many_type_name == 'GenomeSet':
            seq_total = 0

            output_featureSet = dict()
            if 'description' in input_many_genomeSet and input_many_genomeSet['description'] != None:
                output_featureSet['description'] = input_many_genomeSet['description'] + " - "+search_tool_name+"_Search filtered"
            else:
                output_featureSet['description'] = search_tool_name+"_Search filtered"
            #output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()

            self.log(console,"READING HITS FOR GENOMES")  # DEBUG
            for genome_id in feature_ids_by_genome_id.keys():
                self.log(console,"READING HITS FOR GENOME "+genome_id)  # DEBUG
                genome_ref = input_many_genomeSet['elements'][genome_id]['ref']
                for feature_id in feature_ids_by_genome_id[genome_id]:
                    seq_total += 1
                    try:
                        #in_filtered_set = hit_seq_ids[feature['id']]
                        in_filtered_set = hit_seq_ids[genome_ref+genome_id_feature_id_delim+feature_id]
                        #in_filtered_set = hit_seq_ids[feature_id]
                        #self.log(console, 'FOUND HIT: '+feature['id'])  # DEBUG
                        #output_featureSet['element_ordering'].append(feature['id'])
                        try:
                            this_genome_ref_list = output_featureSet['elements'][feature_id]
                        except:
                            output_featureSet['elements'][feature_id] = []
                        output_featureSet['elements'][feature_id].append(genome_ref)
                    except:
                        pass


        # load the method provenance from the context object
        #
        self.log(console,"SETTING PROVENANCE")  # DEBUG
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects'] = []
        provenance[0]['input_ws_objects'].append(input_one_ref)
        provenance[0]['input_ws_objects'].append(input_many_ref)
        provenance[0]['service'] = 'kb_blast'
        provenance[0]['method'] = search_tool_name+'_Search'


        # Upload results
        #
        if len(invalid_msgs) == 0:
            self.log(console,"UPLOADING RESULTS")  # DEBUG

            # input many SequenceSet -> save SequenceSet
            #
            if many_type_name == 'SequenceSet':
                new_obj_info = ws.save_objects({
                            'workspace': params['workspace_name'],
                            'objects':[{
                                    'type': 'KBaseSequences.SequenceSet',
                                    'data': output_sequenceSet,
                                    'name': params['output_filtered_name'],
                                    'meta': {},
                                    'provenance': provenance
                                }]
                        })[0]

            else:  # input FeatureSet, Genome, and GenomeSet -> upload FeatureSet output
                new_obj_info = ws.save_objects({
                            'workspace': params['workspace_name'],
                            'objects':[{
                                    'type': 'KBaseCollections.FeatureSet',
                                    'data': output_featureSet,
                                    'name': params['output_filtered_name'],
                                    'meta': {},
                                    'provenance': provenance
                                }]
                        })[0]

        # build output report object
        #
        self.log(console,"BUILDING REPORT")  # DEBUG
        if len(invalid_msgs) == 0:
            report += 'sequences in many set: '+str(seq_total)+"\n"
            report += 'sequences in hit set:  '+str(hit_total)+"\n"
            report += "\n"
            for line in hit_buf:
                report += line
            reportObj = {
                'objects_created':[{'ref':params['workspace_name']+'/'+params['output_filtered_name'], 'description':search_tool_name+'_Search hits'}],
                'text_message':report
                }
        else:
            report += "FAILURE\n\n"+"\n".join(invalid_msgs)+"\n"
            reportObj = {
                'objects_created':[],
                'text_message':report
                }

        reportName = 'blast_report_'+str(uuid.uuid4())
        report_obj_info = ws.save_objects({
#                'id':info[6],
                'workspace':params['workspace_name'],
                'objects':[
                    {
                        'type':'KBaseReport.Report',
                        'data':reportObj,
                        'name':reportName,
                        'meta':{},
                        'hidden':1,
                        'provenance':provenance
                    }
                ]
            })[0]

        self.log(console,"BUILDING RETURN OBJECT")
#        returnVal = { 'output_report_name': reportName,
#                      'output_report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
#                      'output_filtered_ref': params['workspace_name']+'/'+params['output_filtered_name']
#                      }
        returnVal = { 'report_name': reportName,
                      'report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
                      }
        self.log(console,search_tool_name+"_Search DONE")
        #END BLASTp_Search

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method BLASTp_Search return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def BLASTx_Search(self, ctx, params):
        """
        :param params: instance of type "BLAST_Params" (BLAST Input Params)
           -> structure: parameter "workspace_name" of type "workspace_name"
           (** The workspace object refs are of form: ** **    objects =
           ws.get_objects([{'ref':
           params['workspace_id']+'/'+params['obj_name']}]) ** ** "ref" means
           the entire name combining the workspace id and the object name **
           "id" is a numerical identifier of the workspace or object, and
           should just be used for workspace ** "name" is a string identifier
           of a workspace or object.  This is received from Narrative.),
           parameter "input_one_sequence" of type "sequence", parameter
           "input_one_ref" of type "data_obj_ref", parameter "input_many_ref"
           of type "data_obj_ref", parameter "input_msa_ref" of type
           "data_obj_ref", parameter "output_one_name" of type
           "data_obj_name", parameter "output_filtered_name" of type
           "data_obj_name", parameter "ident_thresh" of Double, parameter
           "e_value" of Double, parameter "bitscore" of Double, parameter
           "overlap_fraction" of Double, parameter "maxaccepts" of Double,
           parameter "rounds" of Double
        :returns: instance of type "BLAST_Output" (BLAST Output) ->
           structure: parameter "report_name" of type "data_obj_name",
           parameter "report_ref" of type "data_obj_ref"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN BLASTx_Search
        console = []
        invalid_msgs = []
        search_tool_name = 'BLASTx'
        self.log(console,'Running '+search_tool_name+'_Search with params=')
        self.log(console, "\n"+pformat(params))
        report = ''
#        report = 'Running '+search_tool_name+'_Search with params='
#        report += "\n"+pformat(params)
        #appropriate_sequence_found_in_one_input = False
        appropriate_sequence_found_in_many_input = False


        #### do some basic checks
        #
        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
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
#        if 'input_one_sequence' not in params:
#            raise ValueError('input_one_sequence parameter is required')
#        if 'input_one_ref' not in params:
#            raise ValueError('input_one_ref parameter is required')
        if 'input_many_ref' not in params:
            raise ValueError('input_many_ref parameter is required')
        if 'output_filtered_name' not in params:
            raise ValueError('output_filtered_name parameter is required')


        # set local names
        input_one_ref  = None
        output_one_ref = None
        input_many_ref = params['input_many_ref']


        # Write the input_one_sequence to a SingleEndLibrary object
        #
        if 'input_one_sequence' in params \
                and params['input_one_sequence'] != None \
                and params['input_one_sequence'] != "Optionally enter DNA sequence...":

            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'])
            ParseFastaStr_retVal = DOTFU.ParseFastaStr ({
                'fasta_str':     params['input_one_sequence'],
                'residue_type': 'NUC',
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
            if 'provenance' in ctx:
                provenance = ctx['provenance']
            # add additional info to provenance here, in this case the input data object reference
            provenance[0]['input_ws_objects'] = []
            provenance[0]['service'] = 'kb_blast'
            provenance[0]['method'] = search_tool_name+'_Search'
                
            # Upload results
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
                ws = workspaceService(self.workspaceURL, token=ctx['token'])
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


        #### Get the input_one object
        ##
        if 'output_one_name' in params and output_one_ref:
            input_one_ref = output_one_ref
        else:
            input_one_ref = params['input_one_ref']

        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            objects = ws.get_objects([{'ref': input_one_ref}])
            input_one_data = objects[0]['data']
            input_one_name = str(objects[0]['info'][1])
            info = objects[0]['info']

            one_type_name = info[2].split('.')[1].split('-')[0]
        except Exception as e:
            raise ValueError('Unable to fetch input_one_ref object from workspace: ' + str(e))
        #to get the full stack trace: traceback.format_exc()

        if 'input_one_sequence' in params \
                and params['input_one_sequence'] != None \
                and params['input_one_sequence'] != "Optionally enter DNA sequence..." \
                and one_type_name != 'SequenceSet':

            self.log(invalid_msgs,"ERROR: Mismatched input type for Query Object: "+input_one_ref+" should be SequenceSet instead of: "+one_type_name)


        # Handle overloading (input_one can be Feature, SingleEndLibrary, or FeatureSet)
        #

        # SequenceSet
        #
        if one_type_name == 'SequenceSet':
            try:
                input_one_sequenceSet = input_one_data
            except Exception as e:
                print(traceback.format_exc())
                raise ValueError('Unable to get sequenceSet object: ' + str(e))

            header_id = input_one_sequenceSet['sequences'][0]['sequence_id']
            sequence_str = input_one_data['sequences'][0]['sequence']

            #PROT_pattern  = re.compile("^[acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWYxX ]+$")
            DNA_pattern   = re.compile("^[acgtuACGTUnryNRY ]+$")
            if not DNA_pattern.match(sequence_str):
                self.log(invalid_msgs,"BAD record for sequence_id: "+header_id+"\n"+sequence_str+"\n")
            else:
                appropriate_sequence_found_in_one_input = True

            one_forward_reads_file_path = os.path.join(self.scratch, header_id+'.fasta')
            one_forward_reads_file_handle = open(one_forward_reads_file_path, 'w', 0)
            self.log(console, 'writing reads file: '+str(one_forward_reads_file_path))
            one_forward_reads_file_handle.write('>'+header_id+"\n")
            one_forward_reads_file_handle.write(sequence_str+"\n")
            one_forward_reads_file_handle.close();
            self.log(console, 'done')        

        # FeatureSet
        #
        elif one_type_name == 'FeatureSet':
            # retrieve sequences for features
            #input_one_featureSet = input_one_data
            genome_id_feature_id_delim = '.f:'
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
                'residue_type':        'nucleotide',
                'feature_type':        'CDS',
                'record_id_pattern':   '%%genome_ref%%'+genome_id_feature_id_delim+'%%feature_id%%',
                'record_desc_pattern': '[%%genome_ref%%]',
                'case':                'upper',
                'linewrap':            50,
                'merge_fasta_files':   'TRUE'
                }

            #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'])
            FeatureSetToFASTA_retVal = DOTFU.FeatureSetToFASTA (FeatureSetToFASTA_params)
            one_forward_reads_file_path = FeatureSetToFASTA_retVal['fasta_file_path']
            if len(FeatureSetToFASTA_retVal['feature_ids_by_genome_ref'].keys()) > 0:
                appropriate_sequence_found_in_one_input = True

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "FeatureSetToFasta() took "+str(end_time-beg_time)+" secs")


        # Feature
        #
        elif one_type_name == 'Feature':
            # export feature to FASTA file
            feature = input_one_data
            one_forward_reads_file_path = os.path.join(self.scratch, input_one_name+".fasta")
            self.log(console, 'writing fasta file: '+one_forward_reads_file_path)
            # BLASTx is nuc-prot
            #if feature['type'] != 'CDS':
            #    self.log(console,params['input_one_ref']+" feature type must be CDS")
            #    self.log(invalid_msgs,params['input_one_ref']+" feature type must be CDS")
            if 'protein_translation' not in feature or feature['protein_translation'] == None:
                #self.log(console,"bad CDS Feature "+params['input_one_name']+": no protein_translation found")
                #raise ValueError ("bad CDS Feature "+params['input_one_name']+": no protein_translation found")
                self.log(console,params['input_one_ref']+" feature type must be CDS")
                self.log(invalid_msgs,params['input_one_ref']+" feature type must be CDS")
            else:
                record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description=genomeRef+"."+feature['id'])
                #record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description=genomeRef+"."+feature['id'])
                SeqIO.write([record], one_forward_reads_file_path, "fasta")
                appropriate_sequence_found_in_one_input = True
        else:
            raise ValueError('Cannot yet handle input_one type of: '+type_name)            


        #### Get the input_many object
        ##
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            objects = ws.get_objects([{'ref': input_many_ref}])
            input_many_data = objects[0]['data']
            info = objects[0]['info']
            input_many_name = str(info[1])
            many_type_name = info[2].split('.')[1].split('-')[0]

        except Exception as e:
            raise ValueError('Unable to fetch input_many_name object from workspace: ' + str(e))
            #to get the full stack trace: traceback.format_exc()

        # Handle overloading (input_many can be FeatureSet, Genome, or GenomeSet)
        #

        # Handle overloading (input_many can be SequenceSet, FeatureSet, Genome, or GenomeSet)
        #

        # SequenceSet
        #
        if many_type_name == 'SequenceSet':
            try:
                input_many_sequenceSet = input_many_data
            except Exception as e:
                print(traceback.format_exc())
                raise ValueError('Unable to get SequenceSet: ' + str(e))

            header_id = input_many_sequenceSet['sequences'][0]['sequence_id']
            many_forward_reads_file_path = os.path.join(self.scratch, header_id+'.fasta')
            many_forward_reads_file_handle = open(many_forward_reads_file_path, 'w', 0)
            self.log(console, 'writing reads file: '+str(many_forward_reads_file_path))

            for seq_obj in input_many_sequenceSet['sequences']:
                header_id = seq_obj['sequence_id']
                sequence_str = seq_obj['sequence']

                PROT_pattern = re.compile("^[acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWYxX ]+$")
                #DNA_pattern = re.compile("^[acgtuACGTUnryNRY ]+$")   
                if not PROT_pattern.match(sequence_str):
                    self.log(invalid_msgs,"BAD record for sequence_id: "+header_id+"\n"+sequence_str+"\n")
                    continue
                appropriate_sequence_found_in_many_input = True
                many_forward_reads_file_handle.write('>'+header_id+"\n")
                many_forward_reads_file_handle.write(sequence_str+"\n")
            many_forward_reads_file_handle.close();
            self.log(console, 'done')


        # FeatureSet
        #
        if many_type_name == 'FeatureSet':
            # retrieve sequences for features
            input_many_featureSet = input_many_data
            genome_id_feature_id_delim = '.f:'
            many_forward_reads_file_dir = self.scratch
            many_forward_reads_file = input_many_name+".fasta"

            # DEBUG
            #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            FeatureSetToFASTA_params = {
                'featureSet_ref':      input_many_ref,
                'file':                many_forward_reads_file,
                'dir':                 many_forward_reads_file_dir,
                'console':             console,
                'invalid_msgs':        invalid_msgs,
                'residue_type':        'nucleotide',
                'feature_type':        'CDS',
                'record_id_pattern':   '%%genome_ref%%'+genome_id_feature_id_delim+'%%feature_id%%',
                'record_desc_pattern': '[%%genome_ref%%]',
                'case':                'upper',
                'linewrap':            50,
                'merge_fasta_files':   'TRUE'
                }

            #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'])
            FeatureSetToFASTA_retVal = DOTFU.FeatureSetToFASTA (FeatureSetToFASTA_params)
            many_forward_reads_file_path = FeatureSetToFASTA_retVal['fasta_file_path']
            feature_ids_by_genome_ref = FeatureSetToFASTA_retVal['feature_ids_by_genome_ref']
            if len(feature_ids_by_genome_ref.keys()) > 0:
                appropriate_sequence_found_in_many_input = True

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "FeatureSetToFasta() took "+str(end_time-beg_time)+" secs")


        # Genome
        #
        elif many_type_name == 'Genome':
            many_forward_reads_file_dir = self.scratch
            many_forward_reads_file = input_many_name+".fasta"

            # DEBUG
            #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            GenomeToFASTA_params = {
                'genome_ref':          input_many_ref,
                'file':                many_forward_reads_file,
                'dir':                 many_forward_reads_file_dir,
                'console':             console,
                'invalid_msgs':        invalid_msgs,
                'residue_type':        'protein',
                'feature_type':        'CDS',
                'record_id_pattern':   '%%feature_id%%',
                'record_desc_pattern': '[%%genome_id%%]',
                'case':                'upper',
                'linewrap':            50
                }

            #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'])
            GenomeToFASTA_retVal = DOTFU.GenomeToFASTA (GenomeToFASTA_params)
            many_forward_reads_file_path = GenomeToFASTA_retVal['fasta_file_path']
            feature_ids = GenomeToFASTA_retVal['feature_ids']
            if len(feature_ids) > 0:
                appropriate_sequence_found_in_many_input = True

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "Genome2Fasta() took "+str(end_time-beg_time)+" secs")


        # GenomeSet
        #
        elif many_type_name == 'GenomeSet':
            input_many_genomeSet = input_many_data
            genome_id_feature_id_delim = '.f:'
            many_forward_reads_file_dir = self.scratch
            many_forward_reads_file = input_many_name+".fasta"

            # DEBUG
            #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            GenomeSetToFASTA_params = {
                'genomeSet_ref':       input_many_ref,
                'file':                many_forward_reads_file,
                'dir':                 many_forward_reads_file_dir,
                'console':             console,
                'invalid_msgs':        invalid_msgs,
                'residue_type':        'protein',
                'feature_type':        'CDS',
                'record_id_pattern':   '%%genome_ref%%'+genome_id_feature_id_delim+'%%feature_id%%',
                'record_desc_pattern': '[%%genome_ref%%]',
                'case':                'upper',
                'linewrap':            50,
                'merge_fasta_files':   'TRUE'
                }

            #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'])
            GenomeSetToFASTA_retVal = DOTFU.GenomeSetToFASTA (GenomeSetToFASTA_params)
            many_forward_reads_file_path = GenomeSetToFASTA_retVal['fasta_file_path_list'][0]
            feature_ids_by_genome_id = GenomeSetToFASTA_retVal['feature_ids_by_genome_id']
            if len(feature_ids_by_genome_id.keys()) > 0:
                appropriate_sequence_found_in_many_input = True

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "FeatureSetToFasta() took "+str(end_time-beg_time)+" secs")


        # Missing proper input_many_type
        #
        else:
            raise ValueError('Cannot yet handle input_many type of: '+type_name)            

        # check for failed input file creation
        #
        if not appropriate_sequence_found_in_one_input:
            self.log(invalid_msgs,"no dna sequences found in '"+input_one_name+"'")
        if not appropriate_sequence_found_in_many_input:
            self.log(invalid_msgs,"no protein sequences found in '"+input_many_name+"'")


        # input data failed validation.  Need to return
        #
        if len(invalid_msgs) > 0:

            # load the method provenance from the context object
            #
            self.log(console,"SETTING PROVENANCE")  # DEBUG
            provenance = [{}]
            if 'provenance' in ctx:
                provenance = ctx['provenance']
            # add additional info to provenance here, in this case the input data object reference
            provenance[0]['input_ws_objects'] = []
            provenance[0]['input_ws_objects'].append(input_one_ref)
            provenance[0]['input_ws_objects'].append(input_many_ref)
            provenance[0]['service'] = 'kb_blast'
            provenance[0]['method'] = search_tool_name+'_Search'


            # build output report object
            #
            self.log(console,"BUILDING REPORT")  # DEBUG
            report += "FAILURE:\n\n"+"\n".join(invalid_msgs)+"\n"
            reportObj = {
                'objects_created':[],
                'text_message':report
                }

            reportName = 'blast_report_'+str(uuid.uuid4())
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            report_obj_info = ws.save_objects({
                    #'id':info[6],
                    'workspace':params['workspace_name'],
                    'objects':[
                        {
                        'type':'KBaseReport.Report',
                        'data':reportObj,
                        'name':reportName,
                        'meta':{},
                        'hidden':1,
                        'provenance':provenance  # DEBUG
                        }
                        ]
                    })[0]

            self.log(console,"BUILDING RETURN OBJECT")
            returnVal = { 'report_name': reportName,
                      'report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
                      }
            self.log(console,search_tool_name+"_Search DONE")
            return [returnVal]


        # FORMAT DB
        #
        # OLD SYNTAX: formatdb -i $database -o T -p F -> $database.nsq or $database.00.nsq
        # NEW SYNTAX: makeblastdb -in $database -parse_seqids -dbtype prot/nucl -out <basename>
        makeblastdb_cmd = [self.Make_BLAST_DB]

        # check for necessary files
        if not os.path.isfile(self.Make_BLAST_DB):
            raise ValueError("no such file '"+self.Make_BLAST_DB+"'")
        if not os.path.isfile(many_forward_reads_file_path):
            raise ValueError("no such file '"+many_forward_reads_file_path+"'")
        elif not os.path.getsize(many_forward_reads_file_path) > 0:
            raise ValueError("empty file '"+many_forward_reads_file_path+"'")

        makeblastdb_cmd.append('-in')
        makeblastdb_cmd.append(many_forward_reads_file_path)
        makeblastdb_cmd.append('-parse_seqids')
        makeblastdb_cmd.append('-dbtype')
        makeblastdb_cmd.append('prot')
        makeblastdb_cmd.append('-out')
        makeblastdb_cmd.append(many_forward_reads_file_path)

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
            line = p.stdout.readline()
            if not line: break
            self.log(console, line.replace('\n', ''))

        p.stdout.close()
        p.wait()
        self.log(console, 'return code: ' + str(p.returncode))
        if p.returncode != 0:
            raise ValueError('Error running makeblastdb, return code: '+str(p.returncode) + 
                '\n\n'+ '\n'.join(console))

        # Check for db output
        if not os.path.isfile(many_forward_reads_file_path+".psq") and not os.path.isfile(many_forward_reads_file_path+".00.psq"):
            raise ValueError("makeblastdb failed to create DB file '"+many_forward_reads_file_path+".psq'")
        elif not os.path.getsize(many_forward_reads_file_path+".psq") > 0 and not os.path.getsize(many_forward_reads_file_path+".00.psq") > 0:
            raise ValueError("makeblastdb created empty DB file '"+many_forward_reads_file_path+".psq'")


        ### Construct the BLAST command
        #
        # OLD SYNTAX: $blast -q $q -G $G -E $E -m $m -e $e_value -v $limit -b $limit -K $limit -p blastx -i $fasta_file -d $database -o $out_file
        # NEW SYNTAX: blastx -query <queryfile> -db <basename> -out <out_aln_file> -outfmt 0/7 (8 became 7) -evalue <e_value> -dust no (DNA) -seg no (AA) -num_threads <num_cores>
        #
        blast_bin = self.BLASTx
        blast_cmd = [blast_bin]

        # check for necessary files
        if not os.path.isfile(blast_bin):
            raise ValueError("no such file '"+blast_bin+"'")
        if not os.path.isfile(one_forward_reads_file_path):
            raise ValueError("no such file '"+one_forward_reads_file_path+"'")
        elif not os.path.getsize(one_forward_reads_file_path) > 0:
            raise ValueError("empty file '"+one_forward_reads_file_path+"'")
        if not os.path.isfile(many_forward_reads_file_path):
            raise ValueError("no such file '"+many_forward_reads_file_path+"'")
        elif not os.path.getsize(many_forward_reads_file_path):
            raise ValueError("empty file '"+many_forward_reads_file_path+"'")

        # set the output path
        timestamp = int((datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()*1000)
        output_dir = os.path.join(self.scratch,'output.'+str(timestamp))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        output_aln_file_path = os.path.join(output_dir, 'alnout.txt');
        output_filtered_fasta_file_path = os.path.join(output_dir, 'output_filtered.fna');

        # this is command for basic search mode
        blast_cmd.append('-query')
        blast_cmd.append(one_forward_reads_file_path)
        blast_cmd.append('-db')
        blast_cmd.append(many_forward_reads_file_path)
        blast_cmd.append('-out')
        blast_cmd.append(output_aln_file_path)
        blast_cmd.append('-outfmt')
        blast_cmd.append('7')
        blast_cmd.append('-evalue')
        blast_cmd.append(str(params['e_value']))

        # options
        if 'maxaccepts' in params:
            if params['maxaccepts']:
                blast_cmd.append('-max_target_seqs')
                blast_cmd.append(str(params['maxaccepts']))

        # Run BLAST, capture output as it happens
        #
        self.log(console, 'RUNNING BLAST:')
        self.log(console, '    '+' '.join(blast_cmd))
#        report += "\n"+'running BLAST:'+"\n"
#        report += '    '+' '.join(blast_cmd)+"\n"

        p = subprocess.Popen(blast_cmd, \
                             cwd = self.scratch, \
                             stdout = subprocess.PIPE, \
                             stderr = subprocess.STDOUT, \
                             shell = False)

        while True:
            line = p.stdout.readline()
            if not line: break
            self.log(console, line.replace('\n', ''))

        p.stdout.close()
        p.wait()
        self.log(console, 'return code: ' + str(p.returncode))
        if p.returncode != 0:
            raise ValueError('Error running BLAST, return code: '+str(p.returncode) + 
                '\n\n'+ '\n'.join(console))


        # get query_len for filtering later
        #
        query_len = 0
        with open(one_forward_reads_file_path, 'r', 0) as query_file_handle:
            for line in query_file_handle:
                if line.startswith('>'):
                    continue
                query_len += len(re.sub(r" ","", line.rstrip())) 
        query_len = query_len/3.0  # BLASTx is nuc-prot

                
        # Parse the BLAST tabular output and store ids to filter many set to make filtered object to save back to KBase
        #
        self.log(console, 'PARSING BLAST ALIGNMENT OUTPUT')
        if not os.path.isfile(output_aln_file_path):
            raise ValueError("failed to create BLAST output: "+output_aln_file_path)
        elif not os.path.getsize(output_aln_file_path) > 0:
            raise ValueError("created empty file for BLAST output: "+output_aln_file_path)
        hit_seq_ids = dict()
        output_aln_file_handle = open (output_aln_file_path, "r", 0)
        output_aln_buf = output_aln_file_handle.readlines()
        output_aln_file_handle.close()
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
            #self.log(console,'HIT LINE: '+line)  # DEBUG
            hit_info = line.split("\t")
            hit_seq_id     = hit_info[1]
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

            # BLAST SOMETIMES ADDS THIS TO IDs.  NO IDEA WHY, BUT GET RID OF IT!
            if hit_seq_id.startswith('gnl|'):
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

        for hit_seq_id in hit_order:
            hit_buf.append(high_bitscore_line[hit_seq_id])

            #self.log(console,"HIT_SEQ_ID: '"+hit_seq_id+"'")
            if 'ident_thresh' in params and float(params['ident_thresh']) > float(high_bitscore_ident[hit_seq_id]):
                continue
            if 'bitscore' in params and float(params['bitscore']) > float(high_bitscore_score[hit_seq_id]):
                continue
            if 'overlap_fraction' in params and float(params['overlap_fraction']) > float(high_bitscore_alnlen[hit_seq_id])/float(query_len):
                continue
            
            hit_total += 1
            hit_seq_ids[hit_seq_id] = True
            self.log(console, "HIT: '"+hit_seq_id+"'")  # DEBUG
        

        self.log(console, 'EXTRACTING HITS FROM INPUT')
        self.log(console, 'MANY_TYPE_NAME: '+many_type_name)  # DEBUG


        # SequenceSet input -> SequenceSet output
        #
        if many_type_name == 'SequenceSet':
            seq_total = len(input_many_sequenceSet['sequences'])

            output_sequenceSet = dict()

            if 'sequence_set_id' in input_many_sequenceSet and input_many_sequenceSet['sequence_set_id'] != None:
                output_sequenceSet['sequence_set_id'] = input_many_sequenceSet['sequence_set_id'] + "."+search_tool_name+"_Search_filtered"
            else:
                output_sequenceSet['sequence_set_id'] = search_tool_name+"_Search_filtered"
            if 'description' in input_many_sequenceSet and input_many_sequenceSet['description'] != None:
                output_sequenceSet['description'] = input_many_sequenceSet['description'] + " - "+search_tool_name+"_Search filtered"
            else:
                output_sequenceSet['description'] = search_tool_name+"_Search filtered"

            self.log(console,"ADDING SEQUENCES TO SEQUENCESET")
            output_sequenceSet['sequences'] = []

            for seq_obj in input_many_sequenceSet['sequences']:
                header_id = seq_obj['sequence_id']
                #header_desc = seq_obj['description']
                #sequence_str = seq_obj['sequence']

                try:
                    #self.log(console,"checking '"+header_id+"'")
                    in_filtered_set = hit_seq_ids[header_id]
                    #self.log(console, 'FOUND HIT '+header_id)  # DEBUG
                    output_sequenceSet['sequences'].append(seq_obj)
                except:
                    pass

        # FeatureSet input -> FeatureSet output
        #
        elif many_type_name == 'FeatureSet':
            seq_total = len(input_many_featureSet['elements'].keys())

            output_featureSet = dict()
            if 'description' in input_many_featureSet and input_many_featureSet['description'] != None:
                output_featureSet['description'] = input_many_featureSet['description'] + " - "+search_tool_name+"_Search filtered"
            else:
                output_featureSet['description'] = search_tool_name+"_Search filtered"
            #output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()

            fId_list = input_many_featureSet['elements'].keys()
            self.log(console,"ADDING FEATURES TO FEATURESET")
            for fId in sorted(fId_list):
                for genome_ref in input_many_featureSet['elements'][fId]:
                    try:
                        #self.log(console,"checking '"+fId+"'")
                        #in_filtered_set = hit_seq_ids[fId]
                        in_filtered_set = hit_seq_ids[genome_ref+genome_id_feature_id_delim+fId]
                        #self.log(console, 'FOUND HIT '+fId)  # DEBUG
                        try:
                            this_genome_ref_list = output_featureSet['elements'][fId]
                        except:
                            output_featureSet['elements'][fId] = []
                        output_featureSet['elements'][fId].append(genome_ref)
                    except:
                        pass

        # Parse Genome hits into FeatureSet
        #
        elif many_type_name == 'Genome':
            seq_total = 0
            output_featureSet = dict()
#            if 'scientific_name' in input_many_genome and input_many_genome['scientific_name'] != None:
#                output_featureSet['description'] = input_many_genome['scientific_name'] + " - "+search_tool_name+"_Search filtered"
#            else:
#                output_featureSet['description'] = search_tool_name+"_Search filtered"
            output_featureSet['description'] = search_tool_name+"_Search filtered"
            #output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()
            for fid in feature_ids:
                seq_total += 1
                try:
                    in_filtered_set = hit_seq_ids[fid]
                    #output_featureSet['element_ordering'].append(fid)
                    output_featureSet['elements'][fid] = [input_many_ref]
                except:
                    pass

        # Parse GenomeSet hits into FeatureSet
        #
        elif many_type_name == 'GenomeSet':
            seq_total = 0

            output_featureSet = dict()
            if 'description' in input_many_genomeSet and input_many_genomeSet['description'] != None:
                output_featureSet['description'] = input_many_genomeSet['description'] + " - "+search_tool_name+"_Search filtered"
            else:
                output_featureSet['description'] = search_tool_name+"_Search filtered"
            #output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()

            self.log(console,"READING HITS FOR GENOMES")  # DEBUG
            for genome_id in feature_ids_by_genome_id.keys():
                self.log(console,"READING HITS FOR GENOME "+genome_id)  # DEBUG
                genome_ref = input_many_genomeSet['elements'][genome_id]['ref']
                for feature_id in feature_ids_by_genome_id[genome_id]:
                    seq_total += 1
                    try:
                        #in_filtered_set = hit_seq_ids[feature['id']]
                        in_filtered_set = hit_seq_ids[genome_ref+genome_id_feature_id_delim+feature_id]
                        #in_filtered_set = hit_seq_ids[feature_id]
                        #self.log(console, 'FOUND HIT: '+feature['id'])  # DEBUG
                        #output_featureSet['element_ordering'].append(feature['id'])
                        try:
                            this_genome_ref_list = output_featureSet['elements'][feature_id]
                        except:
                            output_featureSet['elements'][feature_id] = []
                        output_featureSet['elements'][feature_id].append(genome_ref)
                    except:
                        pass


        # load the method provenance from the context object
        #
        self.log(console,"SETTING PROVENANCE")  # DEBUG
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects'] = []
        provenance[0]['input_ws_objects'].append(input_one_ref)
        provenance[0]['input_ws_objects'].append(input_many_ref)
        provenance[0]['service'] = 'kb_blast'
        provenance[0]['method'] = search_tool_name+'_Search'


        # Upload results
        #
        if len(invalid_msgs) == 0:
            self.log(console,"UPLOADING RESULTS")  # DEBUG

            # input many SequenceSet -> save SequenceSet
            #
            if many_type_name == 'SequenceSet':
                new_obj_info = ws.save_objects({
                            'workspace': params['workspace_name'],
                            'objects':[{
                                    'type': 'KBaseSequences.SequenceSet',
                                    'data': output_sequenceSet,
                                    'name': params['output_filtered_name'],
                                    'meta': {},
                                    'provenance': provenance
                                }]
                        })[0]

            else: # input FeatureSet, Genome, and GenomeSet -> upload FeatureSet output
                new_obj_info = ws.save_objects({
                            'workspace': params['workspace_name'],
                            'objects':[{
                                    'type': 'KBaseCollections.FeatureSet',
                                    'data': output_featureSet,
                                    'name': params['output_filtered_name'],
                                    'meta': {},
                                    'provenance': provenance
                                }]
                        })[0]

        # build output report object
        #
        self.log(console,"BUILDING REPORT")  # DEBUG
        if len(invalid_msgs) == 0:
            report += 'sequences in many set: '+str(seq_total)+"\n"
            report += 'sequences in hit set:  '+str(hit_total)+"\n"
            report += "\n"
            for line in hit_buf:
                report += line
            reportObj = {
                'objects_created':[{'ref':params['workspace_name']+'/'+params['output_filtered_name'], 'description':search_tool_name+'_Search hits'}],
                'text_message':report
                }
        else:
            report += "FAILURE\n\n"+"\n".join(invalid_msgs)+"\n"
            reportObj = {
                'objects_created':[],
                'text_message':report
                }
            
        reportName = 'blast_report_'+str(uuid.uuid4())
        report_obj_info = ws.save_objects({
#                'id':info[6],
                'workspace':params['workspace_name'],
                'objects':[
                    {
                        'type':'KBaseReport.Report',
                        'data':reportObj,
                        'name':reportName,
                        'meta':{},
                        'hidden':1,
                        'provenance':provenance
                    }
                ]
            })[0]

        self.log(console,"BUILDING RETURN OBJECT")
#        returnVal = { 'output_report_name': reportName,
#                      'output_report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
#                      'output_filtered_ref': params['workspace_name']+'/'+params['output_filtered_name']
#                      }
        returnVal = { 'report_name': reportName,
                      'report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
                      }
        self.log(console,search_tool_name+"_Search DONE")
        #END BLASTx_Search

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method BLASTx_Search return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def tBLASTn_Search(self, ctx, params):
        """
        :param params: instance of type "BLAST_Params" (BLAST Input Params)
           -> structure: parameter "workspace_name" of type "workspace_name"
           (** The workspace object refs are of form: ** **    objects =
           ws.get_objects([{'ref':
           params['workspace_id']+'/'+params['obj_name']}]) ** ** "ref" means
           the entire name combining the workspace id and the object name **
           "id" is a numerical identifier of the workspace or object, and
           should just be used for workspace ** "name" is a string identifier
           of a workspace or object.  This is received from Narrative.),
           parameter "input_one_sequence" of type "sequence", parameter
           "input_one_ref" of type "data_obj_ref", parameter "input_many_ref"
           of type "data_obj_ref", parameter "input_msa_ref" of type
           "data_obj_ref", parameter "output_one_name" of type
           "data_obj_name", parameter "output_filtered_name" of type
           "data_obj_name", parameter "ident_thresh" of Double, parameter
           "e_value" of Double, parameter "bitscore" of Double, parameter
           "overlap_fraction" of Double, parameter "maxaccepts" of Double,
           parameter "rounds" of Double
        :returns: instance of type "BLAST_Output" (BLAST Output) ->
           structure: parameter "report_name" of type "data_obj_name",
           parameter "report_ref" of type "data_obj_ref"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN tBLASTn_Search
        console = []
        invalid_msgs = []
        search_tool_name = 'tBLASTn'
        self.log(console,'Running '+search_tool_name+'_Search with params=')
        self.log(console, "\n"+pformat(params))
        report = ''
#        report = 'Running '+search_tool_name+'_Search with params='
#        report += "\n"+pformat(params)
        appropriate_sequence_found_in_one_input = False
        #appropriate_sequence_found_in_many_input = False


        #### do some basic checks
        #
        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
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
#        if 'input_one_sequence' not in params:
#            raise ValueError('input_one_sequence parameter is required')
#        if 'input_one_ref' not in params:
#            raise ValueError('input_one_ref parameter is required')
        if 'input_many_ref' not in params:
            raise ValueError('input_many_ref parameter is required')
        if 'output_filtered_name' not in params:
            raise ValueError('output_filtered_name parameter is required')

        
        # set local names
        input_one_ref  = None
        output_one_ref = None
        input_many_ref = params['input_many_ref']


        # Write the input_one_sequence to file
        #
        if 'input_one_sequence' in params \
                and params['input_one_sequence'] != None \
                and params['input_one_sequence'] != "Optionally enter PROTEIN sequence...":

            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'])
            ParseFastaStr_retVal = DOTFU.ParseFastaStr ({
                'fasta_str':     params['input_one_sequence'],
                'residue_type': 'NUC',
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
            if 'provenance' in ctx:
                provenance = ctx['provenance']
            # add additional info to provenance here, in this case the input data object reference
            provenance[0]['input_ws_objects'] = []
            provenance[0]['service'] = 'kb_blast'
            provenance[0]['method'] = search_tool_name+'_Search'
                
            # Upload results
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
                ws = workspaceService(self.workspaceURL, token=ctx['token'])
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


        #### Get the input_one object
        ##
        if 'output_one_name' in params and output_one_ref:
            input_one_ref = output_one_ref
        else:
            input_one_ref = params['input_one_ref']

        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            objects = ws.get_objects([{'ref': input_one_ref}])
            input_one_data = objects[0]['data']
            input_one_name = str(objects[0]['info'][1])
            info = objects[0]['info']

            one_type_name = info[2].split('.')[1].split('-')[0]
        except Exception as e:
            raise ValueError('Unable to fetch input_one_name object from workspace: ' + str(e))
            #to get the full stack trace: traceback.format_exc()

        if 'input_one_sequence' in params \
                and params['input_one_sequence'] != None \
                and params['input_one_sequence'] != "Optionally enter DNA sequence..." \
                and one_type_name != 'SequenceSet':

            self.log(invalid_msgs,"ERROR: Mismatched input type for Query Object: "+params['input_one_ref']+" should be SequenceSet instead of: "+one_type_name)

        # Handle overloading (input_one can be Feature, or FeatureSet)
        #

        # SequenceSet
        #
        if one_type_name == 'SequenceSet':
            try:
                input_one_sequenceSet = input_one_data
            except Exception as e:
                print(traceback.format_exc())
                raise ValueError('Unable to get sequenceSet object: ' + str(e))

            header_id = input_one_sequenceSet['sequences'][0]['sequence_id']
            sequence_str = input_one_data['sequences'][0]['sequence']

            PROT_pattern = re.compile("^[acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWYxX ]+$")
            #DNA_pattern  = re.compile("^[acgtuACGTUnryNRY ]+$")
            if not PROT_pattern.match(sequence_str):
                self.log(invalid_msgs,"BAD record for sequence_id: "+header_id+"\n"+sequence_str+"\n")
            else:
                appropriate_sequence_found_in_one_input = True

            one_forward_reads_file_path = os.path.join(self.scratch, header_id+'.fasta')
            one_forward_reads_file_handle = open(one_forward_reads_file_path, 'w', 0)
            self.log(console, 'writing reads file: '+str(one_forward_reads_file_path))
            one_forward_reads_file_handle.write('>'+header_id+"\n")
            one_forward_reads_file_handle.write(sequence_str+"\n")
            one_forward_reads_file_handle.close();
            self.log(console, 'done')

        # FeatureSet
        #
        elif one_type_name == 'FeatureSet':
            # retrieve sequences for features
            #input_one_featureSet = input_one_data
            genome_id_feature_id_delim = '.f:'
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
                    'residue_type':        'protein',
                    'feature_type':        'CDS',
                    'record_id_pattern':   '%%genome_ref%%'+genome_id_feature_id_delim+'%%feature_id%%',
                    'record_desc_pattern': '[%%genome_ref%%]',
                    'case':                'upper',
                    'linewrap':            50,
                    'merge_fasta_files':   'TRUE'
                    }

            #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'])
            FeatureSetToFASTA_retVal = DOTFU.FeatureSetToFASTA (FeatureSetToFASTA_params)
            one_forward_reads_file_path = FeatureSetToFASTA_retVal['fasta_file_path']
            if len(FeatureSetToFASTA_retVal['feature_ids_by_genome_ref'].keys()) > 0:
                appropriate_sequence_found_in_one_input = True

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "FeatureSetToFasta() took "+str(end_time-beg_time)+" secs")


        # Feature
        #
        elif one_type_name == 'Feature':
            # export feature to FASTA file
            feature = input_one_data
            one_forward_reads_file_path = os.path.join(self.scratch, input_one_name+".fasta")
            self.log(console, 'writing fasta file: '+one_forward_reads_file_path)
            # tBLASTn is prot-nuc
            #record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description='['+feature['genome_id']+']'+' '+feature['function'])
            #if feature['type'] != 'CDS':
            #    self.log(console,input_one_name+" feature type must be CDS")
            #    self.log(invalid_msgs,input_one_name+" feature type must be CDS")
            if 'protein_translation' not in feature or feature['protein_translation'] == None:
                #self.log(console,"bad CDS Feature "+input_one_name+": no protein_translation found")
                #raise ValueError ("bad CDS Feature "+input_one_name+": no protein_translation found")
                self.log(console,input_one_name+" feature type must be CDS")
                self.log(invalid_msgs,input_one_name+" feature type must be CDS")
            else:
                appropriate_sequence_found_in_one_input = True
                record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description=genomeRef+"."+feature['id'])
                SeqIO.write([record], one_forward_reads_file_path, "fasta")
                appropriate_sequence_found_in_one_input = True                
        else:
            raise ValueError('Cannot yet handle input_one type of: '+type_name)            


        #### Get the input_many object
        ##
        many_forward_reads_file_compression = None
        sequencing_tech = 'N/A'
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            objects = ws.get_objects([{'ref': input_many_ref}])
            input_many_data = objects[0]['data']
            info = objects[0]['info']
            input_many_name = str(info[1])
            many_type_name = info[2].split('.')[1].split('-')[0]

            if many_type_name == 'SingleEndLibrary':
                many_type_namespace = info[2].split('.')[0]
                if many_type_namespace == 'KBaseAssembly':
                    file_name = input_many_data['handle']['file_name']
                elif many_type_namespace == 'KBaseFile':
                    file_name = input_many_data['lib']['file']['file_name']
                else:
                    raise ValueError('bad data type namespace: '+many_type_namespace)
                #self.log(console, 'INPUT_MANY_FILENAME: '+file_name)  # DEBUG
                if file_name[-3:] == ".gz":
                    many_forward_reads_file_compression = 'gz'
                if 'sequencing_tech' in input_many_data:
                    sequencing_tech = input_many_data['sequencing_tech']

        except Exception as e:
            raise ValueError('Unable to fetch input_many_name object from workspace: ' + str(e))
            #to get the full stack trace: traceback.format_exc()


        # Handle overloading (input_many can be SequenceSet, SingleEndLibrary, FeatureSet, Genome, or GenomeSet)
        #
        if many_type_name == 'SequenceSet':
            try:
                input_many_sequenceSet = input_many_data
            except Exception as e:
                print(traceback.format_exc())
                raise ValueError('Unable to get SequenceSet: ' + str(e))

            header_id = input_many_sequenceSet['sequences'][0]['sequence_id']
            many_forward_reads_file_path = os.path.join(self.scratch, header_id+'.fasta')
            many_forward_reads_file_handle = open(many_forward_reads_file_path, 'w', 0)
            self.log(console, 'writing reads file: '+str(many_forward_reads_file_path))

            for seq_obj in input_many_sequenceSet['sequences']:
                header_id = seq_obj['sequence_id']
                sequence_str = seq_obj['sequence']

                #PROT_pattern = re.compile("^[acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWYxX ]+$")
                DNA_pattern = re.compile("^[acgtuACGTUnryNRY ]+$")   
                if not DNA_pattern.match(sequence_str):
                    self.log(invalid_msgs,"BAD record for sequence_id: "+header_id+"\n"+sequence_str+"\n")
                    continue
                appropriate_sequence_found_in_many_input = True
                many_forward_reads_file_handle.write('>'+header_id+"\n")
                many_forward_reads_file_handle.write(sequence_str+"\n")
            many_forward_reads_file_handle.close();
            self.log(console, 'done')

        # SingleEndLibrary
        #
        elif many_type_name == 'SingleEndLibrary':

            # DEBUG
            #for k in data:
            #    self.log(console,"SingleEndLibrary ["+k+"]: "+str(data[k]))

            try:
                if 'lib' in input_many_data:
                    many_forward_reads = input_many_data['lib']['file']
                elif 'handle' in input_many_data:
                    many_forward_reads = input_many_data['handle']
                else:
                    self.log(console,"bad structure for 'many_forward_reads'")
                    raise ValueError("bad structure for 'many_forward_reads'")
                #if 'lib2' in data:
                #    reverse_reads = data['lib2']['file']
                #elif 'handle_2' in data:
                #    reverse_reads = data['handle_2']
                #else:
                #    reverse_reads={}

                ### NOTE: this section is what could be replaced by the transform services
                many_forward_reads_file_path = os.path.join(self.scratch,many_forward_reads['file_name'])
                many_forward_reads_file_handle = open(many_forward_reads_file_path, 'w', 0)
                self.log(console, 'downloading reads file: '+str(many_forward_reads_file_path))
                headers = {'Authorization': 'OAuth '+ctx['token']}
                r = requests.get(many_forward_reads['url']+'/node/'+many_forward_reads['id']+'?download', stream=True, headers=headers)
                for chunk in r.iter_content(1024):
                    appropriate_sequence_found_in_many_input = True
                    many_forward_reads_file_handle.write(chunk)
                many_forward_reads_file_handle.close();
                self.log(console, 'done')
                ### END NOTE


                # remove carriage returns
                new_file_path = many_forward_reads_file_path+"-CRfree"
                new_file_handle = open(new_file_path, 'w', 0)
                many_forward_reads_file_handle = open(many_forward_reads_file_path, 'r', 0)
                for line in many_forward_reads_file_handle:
                    line = re.sub("\r","",line)
                    new_file_handle.write(line)
                many_forward_reads_file_handle.close();
                new_file_handle.close()
                many_forward_reads_file_path = new_file_path


                # convert FASTQ to FASTA (if necessary)
                new_file_path = many_forward_reads_file_path+".fna"
                new_file_handle = open(new_file_path, 'w', 0)
                if many_forward_reads_file_compression == 'gz':
                    many_forward_reads_file_handle = gzip.open(many_forward_reads_file_path, 'r', 0)
                else:
                    many_forward_reads_file_handle = open(many_forward_reads_file_path, 'r', 0)
                header = None
                last_header = None
                last_seq_buf = None
                last_line_was_header = False
                was_fastq = False
                for line in many_forward_reads_file_handle:
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
                many_forward_reads_file_handle.close()
                if was_fastq:
                    many_forward_reads_file_path = new_file_path

            except Exception as e:
                print(traceback.format_exc())
                raise ValueError('Unable to download single-end read library files: ' + str(e))

        # FeatureSet
        #
        elif many_type_name == 'FeatureSet':
            # retrieve sequences for features
            input_many_featureSet = input_many_data
            genome_id_feature_id_delim = '.f:'
            many_forward_reads_file_dir = self.scratch
            many_forward_reads_file = input_many_name+".fasta"

            # DEBUG
            #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            FeatureSetToFASTA_params = {
                'featureSet_ref':      input_many_ref,
                'file':                many_forward_reads_file,
                'dir':                 many_forward_reads_file_dir,
                'console':             console,
                'invalid_msgs':        invalid_msgs,
                'residue_type':        'protein',
                'feature_type':        'CDS',
                'record_id_pattern':   '%%genome_ref%%'+genome_id_feature_id_delim+'%%feature_id%%',
                'record_desc_pattern': '[%%genome_ref%%]',
                'case':                'upper',
                'linewrap':            50,
                'merge_fasta_files':   'TRUE'
                }

            #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'])
            FeatureSetToFASTA_retVal = DOTFU.FeatureSetToFASTA (FeatureSetToFASTA_params)
            many_forward_reads_file_path = FeatureSetToFASTA_retVal['fasta_file_path']
            feature_ids_by_genome_ref = FeatureSetToFASTA_retVal['feature_ids_by_genome_ref']
            if len(feature_ids_by_genome_ref.keys()) > 0:
                appropriate_sequence_found_in_many_input = True

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "FeatureSetToFasta() took "+str(end_time-beg_time)+" secs")


        # Genome
        #
        elif many_type_name == 'Genome':
            many_forward_reads_file_dir = self.scratch
            many_forward_reads_file = input_many_name+".fasta"

            # DEBUG
            #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            GenomeToFASTA_params = {
                'genome_ref':          input_many_ref,
                'file':                many_forward_reads_file,
                'dir':                 many_forward_reads_file_dir,
                'console':             console,
                'invalid_msgs':        invalid_msgs,
                'residue_type':        'nucleotide',
                'feature_type':        'CDS',
                'record_id_pattern':   '%%feature_id%%',
                'record_desc_pattern': '[%%genome_id%%]',
                'case':                'upper',
                'linewrap':            50
                }

            #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'])
            GenomeToFASTA_retVal = DOTFU.GenomeToFASTA (GenomeToFASTA_params)
            many_forward_reads_file_path = GenomeToFASTA_retVal['fasta_file_path']
            feature_ids = GenomeToFASTA_retVal['feature_ids']
            if len(feature_ids) > 0:
                appropriate_sequence_found_in_many_input = True

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "Genome2Fasta() took "+str(end_time-beg_time)+" secs")


        # GenomeSet
        #
        elif many_type_name == 'GenomeSet':
            input_many_genomeSet = input_many_data
            genome_id_feature_id_delim = '.f:'
            many_forward_reads_file_dir = self.scratch
            many_forward_reads_file = input_many_name+".fasta"

            # DEBUG
            #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            GenomeSetToFASTA_params = {
                'genomeSet_ref':       input_many_ref,
                'file':                many_forward_reads_file,
                'dir':                 many_forward_reads_file_dir,
                'console':             console,
                'invalid_msgs':        invalid_msgs,
                'residue_type':        'nucleotide',
                'feature_type':        'CDS',
                'record_id_pattern':   '%%genome_ref%%'+genome_id_feature_id_delim+'%%feature_id%%',
                'record_desc_pattern': '[%%genome_ref%%]',
                'case':                'upper',
                'linewrap':            50,
                'merge_fasta_files':   'TRUE'
                }

            #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'])
            GenomeSetToFASTA_retVal = DOTFU.GenomeSetToFASTA (GenomeSetToFASTA_params)
            many_forward_reads_file_path = GenomeSetToFASTA_retVal['fasta_file_path_list'][0]
            feature_ids_by_genome_id = GenomeSetToFASTA_retVal['feature_ids_by_genome_id']
            if len(feature_ids_by_genome_id.keys()) > 0:
                appropriate_sequence_found_in_many_input = True

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "FeatureSetToFasta() took "+str(end_time-beg_time)+" secs")


        # Missing proper input_many_type
        #
        else:
            raise ValueError('Cannot yet handle input_many type of: '+type_name)            

        # check for failed input file creation
        #
        if not appropriate_sequence_found_in_one_input:
            self.log(invalid_msgs,"no protein sequences found in '"+input_one_name+"'")
        if not appropriate_sequence_found_in_many_input:
            self.log(invalid_msgs,"no dna sequences found in '"+input_many_name+"'")


        # input data failed validation.  Need to return
        #
        if len(invalid_msgs) > 0:

            # load the method provenance from the context object
            #
            self.log(console,"SETTING PROVENANCE")  # DEBUG
            provenance = [{}]
            if 'provenance' in ctx:
                provenance = ctx['provenance']
            # add additional info to provenance here, in this case the input data object reference
            provenance[0]['input_ws_objects'] = []
            provenance[0]['input_ws_objects'].append(input_one_ref)
            provenance[0]['input_ws_objects'].append(input_many_ref)
            provenance[0]['service'] = 'kb_blast'
            provenance[0]['method'] = search_tool_name+'_Search'


            # build output report object
            #
            self.log(console,"BUILDING REPORT")  # DEBUG
            report += "FAILURE:\n\n"+"\n".join(invalid_msgs)+"\n"
            reportObj = {
                'objects_created':[],
                'text_message':report
                }

            reportName = 'blast_report_'+str(uuid.uuid4())
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            report_obj_info = ws.save_objects({
                    #'id':info[6],
                    'workspace':params['workspace_name'],
                    'objects':[
                        {
                        'type':'KBaseReport.Report',
                        'data':reportObj,
                        'name':reportName,
                        'meta':{},
                        'hidden':1,
                        'provenance':provenance  # DEBUG
                        }
                        ]
                    })[0]

            self.log(console,"BUILDING RETURN OBJECT")
            returnVal = { 'report_name': reportName,
                      'report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
                      }
            self.log(console,search_tool_name+"_Search DONE")
            return [returnVal]


        # FORMAT DB
        #
        # OLD SYNTAX: formatdb -i $database -o T -p F -> $database.nsq or $database.00.nsq
        # NEW SYNTAX: makeblastdb -in $database -parse_seqids -dbtype prot/nucl -out <basename>
        makeblastdb_cmd = [self.Make_BLAST_DB]

        # check for necessary files
        if not os.path.isfile(self.Make_BLAST_DB):
            raise ValueError("no such file '"+self.Make_BLAST_DB+"'")
        if not os.path.isfile(many_forward_reads_file_path):
            raise ValueError("no such file '"+many_forward_reads_file_path+"'")
        elif not os.path.getsize(many_forward_reads_file_path) > 0:
            raise ValueError("empty file '"+many_forward_reads_file_path+"'")

        makeblastdb_cmd.append('-in')
        makeblastdb_cmd.append(many_forward_reads_file_path)
        makeblastdb_cmd.append('-parse_seqids')
        makeblastdb_cmd.append('-dbtype')
        makeblastdb_cmd.append('nucl')
        makeblastdb_cmd.append('-out')
        makeblastdb_cmd.append(many_forward_reads_file_path)

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
            line = p.stdout.readline()
            if not line: break
            self.log(console, line.replace('\n', ''))

        p.stdout.close()
        p.wait()
        self.log(console, 'return code: ' + str(p.returncode))
        if p.returncode != 0:
            raise ValueError('Error running makeblastdb, return code: '+str(p.returncode) + 
                '\n\n'+ '\n'.join(console))

        # Check for db output
        if not os.path.isfile(many_forward_reads_file_path+".nsq") and not os.path.isfile(many_forward_reads_file_path+".00.nsq"):
            raise ValueError("makeblastdb failed to create DB file '"+many_forward_reads_file_path+".nsq'")
        elif not os.path.getsize(many_forward_reads_file_path+".nsq") > 0 and not os.path.getsize(many_forward_reads_file_path+".00.nsq") > 0:
            raise ValueError("makeblastdb created empty DB file '"+many_forward_reads_file_path+".nsq'")


        ### Construct the BLAST command
        #
        # OLD SYNTAX: $blast -q $q -G $G -E $E -m $m -e $e_value -v $limit -b $limit -K $limit -p tblastn -i $fasta_file -d $database -o $out_file
        # NEW SYNTAX: tblastn -query <queryfile> -db <basename> -out <out_aln_file> -outfmt 0/7 (8 became 7) -evalue <e_value> -dust no (DNA) -seg no (AA) -num_threads <num_cores>
        #
        blast_bin = self.tBLASTn
        blast_cmd = [blast_bin]

        # check for necessary files
        if not os.path.isfile(blast_bin):
            raise ValueError("no such file '"+blast_bin+"'")
        if not os.path.isfile(one_forward_reads_file_path):
            raise ValueError("no such file '"+one_forward_reads_file_path+"'")
        elif not os.path.getsize(one_forward_reads_file_path) > 0:
            raise ValueError("empty file '"+one_forward_reads_file_path+"'")
        if not os.path.isfile(many_forward_reads_file_path):
            raise ValueError("no such file '"+many_forward_reads_file_path+"'")
        elif not os.path.getsize(many_forward_reads_file_path):
            raise ValueError("empty file '"+many_forward_reads_file_path+"'")

        # set the output path
        timestamp = int((datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()*1000)
        output_dir = os.path.join(self.scratch,'output.'+str(timestamp))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        output_aln_file_path = os.path.join(output_dir, 'alnout.txt');
        output_filtered_fasta_file_path = os.path.join(output_dir, 'output_filtered.fna');

        # this is command for basic search mode
        blast_cmd.append('-query')
        blast_cmd.append(one_forward_reads_file_path)
        blast_cmd.append('-db')
        blast_cmd.append(many_forward_reads_file_path)
        blast_cmd.append('-out')
        blast_cmd.append(output_aln_file_path)
        blast_cmd.append('-outfmt')
        blast_cmd.append('7')
        blast_cmd.append('-evalue')
        blast_cmd.append(str(params['e_value']))

        # options
        if 'maxaccepts' in params:
            if params['maxaccepts']:
                blast_cmd.append('-max_target_seqs')
                blast_cmd.append(str(params['maxaccepts']))

        # Run BLAST, capture output as it happens
        #
        self.log(console, 'RUNNING BLAST:')
        self.log(console, '    '+' '.join(blast_cmd))
#        report += "\n"+'running BLAST:'+"\n"
#        report += '    '+' '.join(blast_cmd)+"\n"

        p = subprocess.Popen(blast_cmd, \
                             cwd = self.scratch, \
                             stdout = subprocess.PIPE, \
                             stderr = subprocess.STDOUT, \
                             shell = False)

        while True:
            line = p.stdout.readline()
            if not line: break
            self.log(console, line.replace('\n', ''))

        p.stdout.close()
        p.wait()
        self.log(console, 'return code: ' + str(p.returncode))
        if p.returncode != 0:
            raise ValueError('Error running BLAST, return code: '+str(p.returncode) + 
                '\n\n'+ '\n'.join(console))


        # get query_len for filtering later
        #
        query_len = 0
        with open(one_forward_reads_file_path, 'r', 0) as query_file_handle:
            for line in query_file_handle:
                if line.startswith('>'):
                    continue
                query_len += len(re.sub(r" ","", line.rstrip())) 
        #query_len = query_len/1.0  # tBLASTn is prot-nuc

                
        # Parse the BLAST tabular output and store ids to filter many set to make filtered object to save back to KBase
        #
        self.log(console, 'PARSING BLAST ALIGNMENT OUTPUT')
        if not os.path.isfile(output_aln_file_path):
            raise ValueError("failed to create BLAST output: "+output_aln_file_path)
        elif not os.path.getsize(output_aln_file_path) > 0:
            raise ValueError("created empty file for BLAST output: "+output_aln_file_path)
        hit_seq_ids = dict()
        output_aln_file_handle = open (output_aln_file_path, "r", 0)
        output_aln_buf = output_aln_file_handle.readlines()
        output_aln_file_handle.close()
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
            #self.log(console,'HIT LINE: '+line)  # DEBUG
            hit_info = line.split("\t")
            hit_seq_id     = hit_info[1]
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

            # BLAST SOMETIMES ADDS THIS TO IDs.  NO IDEA WHY, BUT GET RID OF IT!
            if hit_seq_id.startswith('gnl|'):
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

        for hit_seq_id in hit_order:
            hit_buf.append(high_bitscore_line[hit_seq_id])

            #self.log(console,"HIT_SEQ_ID: '"+hit_seq_id+"'")
            if 'ident_thresh' in params and float(params['ident_thresh']) > float(high_bitscore_ident[hit_seq_id]):
                continue
            if 'bitscore' in params and float(params['bitscore']) > float(high_bitscore_score[hit_seq_id]):
                continue
            if 'overlap_fraction' in params and float(params['overlap_fraction']) > float(high_bitscore_alnlen[hit_seq_id])/float(query_len):
                continue
            
            hit_total += 1
            hit_seq_ids[hit_seq_id] = True
            self.log(console, "HIT: '"+hit_seq_id+"'")  # DEBUG
        

        self.log(console, 'EXTRACTING HITS FROM INPUT')
        self.log(console, 'MANY_TYPE_NAME: '+many_type_name)  # DEBUG


        # SequenceSet input -> SequenceSet output
        #
        if many_type_name == 'SequenceSet':
            seq_total = len(input_many_sequenceSet['sequences'])

            output_sequenceSet = dict()

            if 'sequence_set_id' in input_many_sequenceSet and input_many_sequenceSet['sequence_set_id'] != None:
                output_sequenceSet['sequence_set_id'] = input_many_sequenceSet['sequence_set_id'] + "."+search_tool_name+"_Search_filtered"
            else:
                output_sequenceSet['sequence_set_id'] = search_tool_name+"_Search_filtered"
            if 'description' in input_many_sequenceSet and input_many_sequenceSet['description'] != None:
                output_sequenceSet['description'] = input_many_sequenceSet['description'] + " - "+search_tool_name+"_Search filtered"
            else:
                output_sequenceSet['description'] = search_tool_name+"_Search filtered"

            self.log(console,"ADDING SEQUENCES TO SEQUENCESET")
            output_sequenceSet['sequences'] = []

            for seq_obj in input_many_sequenceSet['sequences']:
                header_id = seq_obj['sequence_id']
                #header_desc = seq_obj['description']
                #sequence_str = seq_obj['sequence']

                try:
                    #self.log(console,"checking '"+header_id+"'")
                    in_filtered_set = hit_seq_ids[header_id]
                    #self.log(console, 'FOUND HIT '+header_id)  # DEBUG
                    output_sequenceSet['sequences'].append(seq_obj)
                except:
                    pass


        # SingleEndLibrary input -> SingleEndLibrary output
        #
        elif many_type_name == 'SingleEndLibrary':

            #  Note: don't use SeqIO.parse because loads everything into memory
            #
#            with open(many_forward_reads_file_path, 'r', -1) as many_forward_reads_file_handle, open(output_filtered_fasta_file_path, 'w', -1) as output_filtered_fasta_file_handle:
            output_filtered_fasta_file_handle = open(output_filtered_fasta_file_path, 'w', -1)
            if many_forward_reads_file_compression == 'gz':
                many_forward_reads_file_handle = gzip.open(many_forward_reads_file_path, 'r', -1)
            else:
                many_forward_reads_file_handle = open(many_forward_reads_file_path, 'r', -1)

            seq_total = 0;
            filtered_seq_total = 0
            last_seq_buf = []
            last_seq_id = None
            last_header = None
            pattern = re.compile('^\S*')
            for line in many_forward_reads_file_handle:
                if line.startswith('>'):
                    #self.log(console, 'LINE: '+line)  # DEBUG
                    seq_total += 1
                    seq_id = line[1:]  # removes '>'
                    seq_id = pattern.findall(seq_id)[0]

                    if last_seq_id != None:
                        #self.log(console, 'ID: '+last_seq_id)  # DEBUG
                        try:
                            in_filtered_set = hit_seq_ids[last_seq_id]
                            #self.log(console, 'FOUND HIT '+last_seq_id)  # DEBUG
                            filtered_seq_total += 1
                            output_filtered_fasta_file_handle.write(last_header)
                            output_filtered_fasta_file_handle.writelines(last_seq_buf)
                        except:
                            pass
                        
                    last_seq_buf = []
                    last_seq_id = seq_id
                    last_header = line
                else:
                    last_seq_buf.append(line)

            if last_seq_id != None:
                #self.log(console, 'ID: '+last_seq_id)  # DEBUG
                try:
                    in_filtered_set = hit_seq_ids[last_seq_id]
                    #self.log(console, 'FOUND HIT: '+last_seq_id)  # DEBUG
                    filtered_seq_total += 1
                    output_filtered_fasta_file_handle.write(last_header)
                    output_filtered_fasta_file_handle.writelines(last_seq_buf)
                except:
                    pass
                
            last_seq_buf = []
            last_seq_id = None
            last_header = None

            many_forward_reads_file_handle.close()
            output_filtered_fasta_file_handle.close()

            if filtered_seq_total != hit_total:
                self.log(console,'hits in BLAST alignment output '+str(hit_total)+' != '+str(filtered_seq_total)+' matched sequences in input file')
                raise ValueError('hits in BLAST alignment output '+str(hit_total)+' != '+str(filtered_seq_total)+' matched sequences in input file')


        # FeatureSet input -> FeatureSet output
        #
        elif many_type_name == 'FeatureSet':
            seq_total = len(input_many_featureSet['elements'].keys())

            output_featureSet = dict()
            if 'description' in input_many_featureSet and input_many_featureSet['description'] != None:
                output_featureSet['description'] = input_many_featureSet['description'] + " - "+search_tool_name+"_Search filtered"
            else:
                output_featureSet['description'] = search_tool_name+"_Search filtered"
            #output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()

            fId_list = input_many_featureSet['elements'].keys()
            self.log(console,"ADDING FEATURES TO FEATURESET")
            for fId in sorted(fId_list):
                for genome_ref in input_many_featureSet['elements'][fId]:
                    try:
                        #self.log(console,"checking '"+fId+"'")
                        #in_filtered_set = hit_seq_ids[fId]
                        in_filtered_set = hit_seq_ids[genome_ref+genome_id_feature_id_delim+fId]
                        #self.log(console, 'FOUND HIT '+fId)  # DEBUG
                        try:
                            this_genome_ref_list = output_featureSet['elements'][fId]
                        except:
                            output_featureSet['elements'][fId] = []
                        output_featureSet['elements'][fId].append(genome_ref)
                    except:
                        pass

        # Parse Genome hits into FeatureSet
        #
        elif many_type_name == 'Genome':
            seq_total = 0
            output_featureSet = dict()
#            if 'scientific_name' in input_many_genome and input_many_genome['scientific_name'] != None:
#                output_featureSet['description'] = input_many_genome['scientific_name'] + " - "+search_tool_name+"_Search filtered"
#            else:
#                output_featureSet['description'] = search_tool_name+"_Search filtered"
            output_featureSet['description'] = search_tool_name+"_Search filtered"
            #output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()
            for fid in feature_ids:
                seq_total += 1
                try:
                    in_filtered_set = hit_seq_ids[fid]
                    #output_featureSet['element_ordering'].append(fid)
                    output_featureSet['elements'][fid] = [input_many_ref]
                except:
                    pass

        # Parse GenomeSet hits into FeatureSet
        #
        elif many_type_name == 'GenomeSet':
            seq_total = 0

            output_featureSet = dict()
            if 'description' in input_many_genomeSet and input_many_genomeSet['description'] != None:
                output_featureSet['description'] = input_many_genomeSet['description'] + " - "+search_tool_name+"_Search filtered"
            else:
                output_featureSet['description'] = search_tool_name+"_Search filtered"
            #output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()

            self.log(console,"READING HITS FOR GENOMES")  # DEBUG
            for genome_id in feature_ids_by_genome_id.keys():
                self.log(console,"READING HITS FOR GENOME "+genome_id)  # DEBUG
                genome_ref = input_many_genomeSet['elements'][genome_id]['ref']
                for feature_id in feature_ids_by_genome_id[genome_id]:
                    seq_total += 1
                    try:
                        #in_filtered_set = hit_seq_ids[feature['id']]
                        in_filtered_set = hit_seq_ids[genome_ref+genome_id_feature_id_delim+feature_id]
                        #in_filtered_set = hit_seq_ids[feature_id]
                        #self.log(console, 'FOUND HIT: '+feature['id'])  # DEBUG
                        #output_featureSet['element_ordering'].append(feature['id'])
                        try:
                            this_genome_ref_list = output_featureSet['elements'][feature_id]
                        except:
                            output_featureSet['elements'][feature_id] = []
                        output_featureSet['elements'][feature_id].append(genome_ref)
                    except:
                        pass


        # load the method provenance from the context object
        #
        self.log(console,"SETTING PROVENANCE")  # DEBUG
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects'] = []
        provenance[0]['input_ws_objects'].append(input_one_ref)
        provenance[0]['input_ws_objects'].append(input_many_ref)
        provenance[0]['service'] = 'kb_blast'
        provenance[0]['method'] = search_tool_name+'_Search'


        # Upload results
        #
        if len(invalid_msgs) == 0:
            self.log(console,"UPLOADING RESULTS")  # DEBUG



            # input SingleEndLibrary -> upload SingleEndLibrary
            #
            if many_type_name == 'SingleEndLibrary':
                
                self.upload_SingleEndLibrary_to_shock_and_ws (ctx,
                                                          console,  # DEBUG
                                                          params['workspace_name'],
                                                          params['output_filtered_name'],
                                                          output_filtered_fasta_file_path,
                                                          provenance,
                                                          sequencing_tech
                                                         )

            # input many SequenceSet -> save SequenceSet
            #
            elif many_type_name == 'SequenceSet':
                new_obj_info = ws.save_objects({
                            'workspace': params['workspace_name'],
                            'objects':[{
                                    'type': 'KBaseSequences.SequenceSet',
                                    'data': output_sequenceSet,
                                    'name': params['output_filtered_name'],
                                    'meta': {},
                                    'provenance': provenance
                                }]
                        })[0]

            else:  # input FeatureSet, Genome, and GenomeSet -> upload FeatureSet output
                new_obj_info = ws.save_objects({
                            'workspace': params['workspace_name'],
                            'objects':[{
                                    'type': 'KBaseCollections.FeatureSet',
                                    'data': output_featureSet,
                                    'name': params['output_filtered_name'],
                                    'meta': {},
                                    'provenance': provenance
                                }]
                        })[0]

        # build output report object
        #
        self.log(console,"BUILDING REPORT")  # DEBUG
        if len(invalid_msgs) == 0:
            report += 'sequences in many set: '+str(seq_total)+"\n"
            report += 'sequences in hit set:  '+str(hit_total)+"\n"
            report += "\n"
            for line in hit_buf:
                report += line
            reportObj = {
                'objects_created':[{'ref':params['workspace_name']+'/'+params['output_filtered_name'], 'description':search_tool_name+'_Search hits'}],
                'text_message':report
                }
        else:
            report += "FAILURE\n\n"+"\n".join(invalid_msgs)+"\n"
            reportObj = {
                'objects_created':[],
                'text_message':report
                }

        reportName = 'blast_report_'+str(uuid.uuid4())
        report_obj_info = ws.save_objects({
#                'id':info[6],
                'workspace':params['workspace_name'],
                'objects':[
                    {
                        'type':'KBaseReport.Report',
                        'data':reportObj,
                        'name':reportName,
                        'meta':{},
                        'hidden':1,
                        'provenance':provenance
                    }
                ]
            })[0]

        self.log(console,"BUILDING RETURN OBJECT")
#        returnVal = { 'output_report_name': reportName,
#                      'output_report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
#                      'output_filtered_ref': params['workspace_name']+'/'+params['output_filtered_name']
#                      }
        returnVal = { 'report_name': reportName,
                      'report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
                      }
        self.log(console,search_tool_name+"_Search DONE")
        #END tBLASTn_Search

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method tBLASTn_Search return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def tBLASTx_Search(self, ctx, params):
        """
        :param params: instance of type "BLAST_Params" (BLAST Input Params)
           -> structure: parameter "workspace_name" of type "workspace_name"
           (** The workspace object refs are of form: ** **    objects =
           ws.get_objects([{'ref':
           params['workspace_id']+'/'+params['obj_name']}]) ** ** "ref" means
           the entire name combining the workspace id and the object name **
           "id" is a numerical identifier of the workspace or object, and
           should just be used for workspace ** "name" is a string identifier
           of a workspace or object.  This is received from Narrative.),
           parameter "input_one_sequence" of type "sequence", parameter
           "input_one_ref" of type "data_obj_ref", parameter "input_many_ref"
           of type "data_obj_ref", parameter "input_msa_ref" of type
           "data_obj_ref", parameter "output_one_name" of type
           "data_obj_name", parameter "output_filtered_name" of type
           "data_obj_name", parameter "ident_thresh" of Double, parameter
           "e_value" of Double, parameter "bitscore" of Double, parameter
           "overlap_fraction" of Double, parameter "maxaccepts" of Double,
           parameter "rounds" of Double
        :returns: instance of type "BLAST_Output" (BLAST Output) ->
           structure: parameter "report_name" of type "data_obj_name",
           parameter "report_ref" of type "data_obj_ref"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN tBLASTx_Search
        console = []
        invalid_msgs = []
        search_tool_name = 'tBLASTx'
        self.log(console,'Running '+search_tool_name+'_Search with params=')
        self.log(console, "\n"+pformat(params))
        report = ''
#        report = 'Running '+search_tool_name+'_Search with params='
#        report += "\n"+pformat(params)


        #### do some basic checks
        #
        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
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
#        if 'input_one_sequence' not in params:
#            raise ValueError('input_one_sequence parameter is required')
#        if 'input_one_ref' not in params:
#            raise ValueError('input_one_ref parameter is required')
        if 'input_many_ref' not in params:
            raise ValueError('input_many_ref parameter is required')
        if 'output_filtered_name' not in params:
            raise ValueError('output_filtered_name parameter is required')


        # set local names
        input_one_ref  = None
        output_one_ref = None
        input_many_ref = params['input_many_ref']


        # Write the input_one_sequence to a SingleEndLibrary object
        #
        if 'input_one_sequence' in params \
                and params['input_one_sequence'] != None \
                and params['input_one_sequence'] != "Optionally enter DNA sequence...":

            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'])
            ParseFastaStr_retVal = DOTFU.ParseFastaStr ({
                'fasta_str':     params['input_one_sequence'],
                'residue_type': 'NUC',
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
            if 'provenance' in ctx:
                provenance = ctx['provenance']
            # add additional info to provenance here, in this case the input data object reference
            provenance[0]['input_ws_objects'] = []
            provenance[0]['service'] = 'kb_blast'
            provenance[0]['method'] = search_tool_name+'_Search'
                
            # Upload results
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
                ws = workspaceService(self.workspaceURL, token=ctx['token'])
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


        #### Get the input_one object
        ##
        if 'output_one_name' in params and output_one_ref:
            input_one_ref = output_one_ref
        else:
            input_one_ref = params['input_one_ref']

        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            objects = ws.get_objects([{'ref': input_one_ref}])
            input_one_data = objects[0]['data']
            input_one_name = str(objects[0]['info'][1])
            info = objects[0]['info']

            one_type_name = info[2].split('.')[1].split('-')[0]
        except Exception as e:
            raise ValueError('Unable to fetch input_one_name object from workspace: ' + str(e))
        #to get the full stack trace: traceback.format_exc()

        if 'input_one_sequence' in params \
                and params['input_one_sequence'] != None \
                and params['input_one_sequence'] != "Optionally enter DNA sequence..." \
                and one_type_name != 'SequenceSet':

            self.log(invalid_msgs,"ERROR: Mismatched input type: input_one_name should be SequenceSet instead of: "+one_type_name)

        # Handle overloading (input_one can be Feature, SingleEndLibrary, or FeatureSet)
        #

        # SequenceSet
        #
        if one_type_name == 'SequenceSet':
            try:
                input_one_sequenceSet = input_one_data
            except Exception as e:
                print(traceback.format_exc())
                raise ValueError('Unable to get sequenceSet object: ' + str(e))

            header_id = input_one_sequenceSet['sequences'][0]['sequence_id']
            sequence_str = input_one_data['sequences'][0]['sequence']

            #PROT_pattern = re.compile("^[acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWYxX ]+$")
            DNA_pattern  = re.compile("^[acgtuACGTUnryNRY ]+$")   
            if not DNA_pattern.match(sequence_str):
                self.log(invalid_msgs,"BAD record for sequence_id: "+header_id+"\n"+sequence_str+"\n")
            else:
                appropriate_sequence_found_in_one_input = True

            one_forward_reads_file_path = os.path.join(self.scratch, header_id+'.fasta')
            one_forward_reads_file_handle = open(one_forward_reads_file_path, 'w', 0)
            self.log(console, 'writing reads file: '+str(one_forward_reads_file_path))
            one_forward_reads_file_handle.write('>'+header_id+"\n")
            one_forward_reads_file_handle.write(sequence_str+"\n")
            one_forward_reads_file_handle.close();
            self.log(console, 'done')

        # FeatureSet
        #
        elif one_type_name == 'FeatureSet':
            # retrieve sequences for features
            #input_one_featureSet = input_one_data
            genome_id_feature_id_delim = '.f:'
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
                'residue_type':        'nucleotide',
                'feature_type':        'CDS',
                'record_id_pattern':   '%%genome_ref%%'+genome_id_feature_id_delim+'%%feature_id%%',
                'record_desc_pattern': '[%%genome_ref%%]',
                'case':                'upper',
                'linewrap':            50,
                'merge_fasta_files':   'TRUE'
                }

            #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'])
            FeatureSetToFASTA_retVal = DOTFU.FeatureSetToFASTA (FeatureSetToFASTA_params)
            one_forward_reads_file_path = FeatureSetToFASTA_retVal['fasta_file_path']
            if len(FeatureSetToFASTA_retVal['feature_ids_by_genome_ref'].keys()) > 0:
                appropriate_sequence_found_in_one_input = True

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "FeatureSetToFasta() took "+str(end_time-beg_time)+" secs")


        # Feature
        #
        elif one_type_name == 'Feature':
            # export feature to FASTA file
            feature = input_one_data
            one_forward_reads_file_path = os.path.join(self.scratch, input_one_name+".fasta")
            self.log(console, 'writing fasta file: '+one_forward_reads_file_path)
            # tBLASTx is nuc-nuc (translated)
            #if feature['type'] != 'CDS':
            #    self.log(console,input_one_name+" feature type must be CDS")
            #    self.log(invalid_msgs,input_one_name+" feature type must be CDS")
            if 'protein_translation' not in feature or feature['protein_translation'] == None:
                #self.log(console,"bad CDS Feature "+params['input_one_name']+": no protein_translation found")
                #raise ValueError ("bad CDS Feature "+params['input_one_name']+": no protein_translation found")
                self.log(console,input_one_name+" feature type must be CDS")
                self.log(invalid_msgs,input_one_name+" feature type must be CDS")
            else:
                record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description=genomeRef+"."+feature['id'])
                #record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description=genomeRef+"."+feature['id'])
                SeqIO.write([record], one_forward_reads_file_path, "fasta")
                appropriate_sequence_found_in_one_input = True
        else:
            raise ValueError('Cannot yet handle input_one type of: '+type_name)            


        #### Get the input_many object
        ##
        many_forward_reads_file_compression = None
        sequencing_tech = 'N/A'
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            objects = ws.get_objects([{'ref': input_many_ref}])
            input_many_data = objects[0]['data']
            info = objects[0]['info']
            input_many_name = str(info[1])
            many_type_name = info[2].split('.')[1].split('-')[0]

            if many_type_name == 'SingleEndLibrary':
                many_type_namespace = info[2].split('.')[0]
                if many_type_namespace == 'KBaseAssembly':
                    file_name = input_many_data['handle']['file_name']
                elif many_type_namespace == 'KBaseFile':
                    file_name = input_many_data['lib']['file']['file_name']
                else:
                    raise ValueError('bad data type namespace: '+many_type_namespace)
                #self.log(console, 'INPUT_MANY_FILENAME: '+file_name)  # DEBUG
                if file_name[-3:] == ".gz":
                    many_forward_reads_file_compression = 'gz'
                if 'sequencing_tech' in input_many_data:
                    sequencing_tech = input_many_data['sequencing_tech']

        except Exception as e:
            raise ValueError('Unable to fetch input_many_name object from workspace: ' + str(e))
            #to get the full stack trace: traceback.format_exc()


        # Handle overloading (input_many can be SequenceSet, SingleEndLibrary, FeatureSet, Genome, or GenomeSet)
        #
        if many_type_name == 'SequenceSet':
            try:
                input_many_sequenceSet = input_many_data
            except Exception as e:
                print(traceback.format_exc())
                raise ValueError('Unable to get SequenceSet: ' + str(e))

            header_id = input_many_sequenceSet['sequences'][0]['sequence_id']
            many_forward_reads_file_path = os.path.join(self.scratch, header_id+'.fasta')
            many_forward_reads_file_handle = open(many_forward_reads_file_path, 'w', 0)
            self.log(console, 'writing reads file: '+str(many_forward_reads_file_path))

            for seq_obj in input_many_sequenceSet['sequences']:
                header_id = seq_obj['sequence_id']
                sequence_str = seq_obj['sequence']

                #PROT_pattern = re.compile("^[acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWYxX ]+$")
                DNA_pattern = re.compile("^[acgtuACGTUnryNRY ]+$")   
                if not DNA_pattern.match(sequence_str):
                    self.log(invalid_msgs,"BAD record for sequence_id: "+header_id+"\n"+sequence_str+"\n")
                    continue
                appropriate_sequence_found_in_many_input = True
                many_forward_reads_file_handle.write('>'+header_id+"\n")
                many_forward_reads_file_handle.write(sequence_str+"\n")
            many_forward_reads_file_handle.close();
            self.log(console, 'done')

        # SingleEndLibrary
        #
        elif many_type_name == 'SingleEndLibrary':

            # DEBUG
            #for k in data:
            #    self.log(console,"SingleEndLibrary ["+k+"]: "+str(data[k]))

            try:
                if 'lib' in input_many_data:
                    many_forward_reads = input_many_data['lib']['file']
                elif 'handle' in input_many_data:
                    many_forward_reads = input_many_data['handle']
                else:
                    self.log(console,"bad structure for 'many_forward_reads'")
                    raise ValueError("bad structure for 'many_forward_reads'")
                #if 'lib2' in data:
                #    reverse_reads = data['lib2']['file']
                #elif 'handle_2' in data:
                #    reverse_reads = data['handle_2']
                #else:
                #    reverse_reads={}

                ### NOTE: this section is what could be replaced by the transform services
                many_forward_reads_file_path = os.path.join(self.scratch,many_forward_reads['file_name'])
                many_forward_reads_file_handle = open(many_forward_reads_file_path, 'w', 0)
                self.log(console, 'downloading reads file: '+str(many_forward_reads_file_path))
                headers = {'Authorization': 'OAuth '+ctx['token']}
                r = requests.get(many_forward_reads['url']+'/node/'+many_forward_reads['id']+'?download', stream=True, headers=headers)
                for chunk in r.iter_content(1024):
                    appropriate_sequence_found_in_many_input = True
                    many_forward_reads_file_handle.write(chunk)
                many_forward_reads_file_handle.close();
                self.log(console, 'done')
                ### END NOTE


                # remove carriage returns
                new_file_path = many_forward_reads_file_path+"-CRfree"
                new_file_handle = open(new_file_path, 'w', 0)
                many_forward_reads_file_handle = open(many_forward_reads_file_path, 'r', 0)
                for line in many_forward_reads_file_handle:
                    line = re.sub("\r","",line)
                    new_file_handle.write(line)
                many_forward_reads_file_handle.close();
                new_file_handle.close()
                many_forward_reads_file_path = new_file_path


                # convert FASTQ to FASTA (if necessary)
                new_file_path = many_forward_reads_file_path+".fna"
                new_file_handle = open(new_file_path, 'w', 0)
                if many_forward_reads_file_compression == 'gz':
                    many_forward_reads_file_handle = gzip.open(many_forward_reads_file_path, 'r', 0)
                else:
                    many_forward_reads_file_handle = open(many_forward_reads_file_path, 'r', 0)
                header = None
                last_header = None
                last_seq_buf = None
                last_line_was_header = False
                was_fastq = False
                for line in many_forward_reads_file_handle:
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
                many_forward_reads_file_handle.close()
                if was_fastq:
                    many_forward_reads_file_path = new_file_path

            except Exception as e:
                print(traceback.format_exc())
                raise ValueError('Unable to download single-end read library files: ' + str(e))

        # FeatureSet
        #
        elif many_type_name == 'FeatureSet':
            # retrieve sequences for features
            input_many_featureSet = input_many_data
            genome_id_feature_id_delim = '.f:'
            many_forward_reads_file_dir = self.scratch
            many_forward_reads_file = input_many_name+".fasta"

            # DEBUG
            #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            FeatureSetToFASTA_params = {
                'featureSet_ref':      input_many_ref,
                'file':                many_forward_reads_file,
                'dir':                 many_forward_reads_file_dir,
                'console':             console,
                'invalid_msgs':        invalid_msgs,
                'residue_type':        'nucleotide',
                'feature_type':        'CDS',
                'record_id_pattern':   '%%genome_ref%%'+genome_id_feature_id_delim+'%%feature_id%%',
                'record_desc_pattern': '[%%genome_ref%%]',
                'case':                'upper',
                'linewrap':            50,
                'merge_fasta_files':   'TRUE'
                }

            #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'])
            FeatureSetToFASTA_retVal = DOTFU.FeatureSetToFASTA (FeatureSetToFASTA_params)
            many_forward_reads_file_path = FeatureSetToFASTA_retVal['fasta_file_path']
            feature_ids_by_genome_ref = FeatureSetToFASTA_retVal['feature_ids_by_genome_ref']
            if len(feature_ids_by_genome_ref.keys()) > 0:
                appropriate_sequence_found_in_many_input = True

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "FeatureSetToFasta() took "+str(end_time-beg_time)+" secs")


        # Genome
        #
        elif many_type_name == 'Genome':
            many_forward_reads_file_dir = self.scratch
            many_forward_reads_file = input_many_name+".fasta"

            # DEBUG
            #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            GenomeToFASTA_params = {
                'genome_ref':          input_many_ref,
                'file':                many_forward_reads_file,
                'dir':                 many_forward_reads_file_dir,
                'console':             console,
                'invalid_msgs':        invalid_msgs,
                'residue_type':        'nucleotide',
                'feature_type':        'CDS',
                'record_id_pattern':   '%%feature_id%%',
                'record_desc_pattern': '[%%genome_id%%]',
                'case':                'upper',
                'linewrap':            50
                }

            #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'])
            GenomeToFASTA_retVal = DOTFU.GenomeToFASTA (GenomeToFASTA_params)
            many_forward_reads_file_path = GenomeToFASTA_retVal['fasta_file_path']
            feature_ids = GenomeToFASTA_retVal['feature_ids']
            if len(feature_ids) > 0:
                appropriate_sequence_found_in_many_input = True

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "Genome2Fasta() took "+str(end_time-beg_time)+" secs")


        # GenomeSet
        #
        elif many_type_name == 'GenomeSet':
            input_many_genomeSet = input_many_data
            genome_id_feature_id_delim = '.f:'
            many_forward_reads_file_dir = self.scratch
            many_forward_reads_file = input_many_name+".fasta"

            # DEBUG
            #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            GenomeSetToFASTA_params = {
                'genomeSet_ref':       input_many_ref,
                'file':                many_forward_reads_file,
                'dir':                 many_forward_reads_file_dir,
                'console':             console,
                'invalid_msgs':        invalid_msgs,
                'residue_type':        'nucleotide',
                'feature_type':        'CDS',
                'record_id_pattern':   '%%genome_ref%%'+genome_id_feature_id_delim+'%%feature_id%%',
                'record_desc_pattern': '[%%genome_ref%%]',
                'case':                'upper',
                'linewrap':            50,
                'merge_fasta_files':   'TRUE'
                }

            #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'])
            GenomeSetToFASTA_retVal = DOTFU.GenomeSetToFASTA (GenomeSetToFASTA_params)
            many_forward_reads_file_path = GenomeSetToFASTA_retVal['fasta_file_path_list'][0]
            feature_ids_by_genome_id = GenomeSetToFASTA_retVal['feature_ids_by_genome_id']
            if len(feature_ids_by_genome_id.keys()) > 0:
                appropriate_sequence_found_in_many_input = True

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "FeatureSetToFasta() took "+str(end_time-beg_time)+" secs")


        # Missing proper input_many_type
        #
        else:
            raise ValueError('Cannot yet handle input_many type of: '+type_name)            


        # check for failed input file creation
        #
        if not appropriate_sequence_found_in_one_input:
            self.log(invalid_msgs,"no dna sequences found in '"+input_one_name+"'")
        if not appropriate_sequence_found_in_many_input:
            self.log(invalid_msgs,"no dna sequences found in '"+input_many_name+"'")


        # input data failed validation.  Need to return
        #
        if len(invalid_msgs) > 0:

            # load the method provenance from the context object
            #
            self.log(console,"SETTING PROVENANCE")  # DEBUG
            provenance = [{}]
            if 'provenance' in ctx:
                provenance = ctx['provenance']
            # add additional info to provenance here, in this case the input data object reference
            provenance[0]['input_ws_objects'] = []
            provenance[0]['input_ws_objects'].append(input_one_ref)
            provenance[0]['input_ws_objects'].append(input_many_ref)
            provenance[0]['service'] = 'kb_blast'
            provenance[0]['method'] = search_tool_name+'_Search'


            # build output report object
            #
            self.log(console,"BUILDING REPORT")  # DEBUG
            report += "FAILURE:\n\n"+"\n".join(invalid_msgs)+"\n"
            reportObj = {
                'objects_created':[],
                'text_message':report
                }

            reportName = 'blast_report_'+str(uuid.uuid4())
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            report_obj_info = ws.save_objects({
                    #'id':info[6],
                    'workspace':params['workspace_name'],
                    'objects':[
                        {
                        'type':'KBaseReport.Report',
                        'data':reportObj,
                        'name':reportName,
                        'meta':{},
                        'hidden':1,
                        'provenance':provenance  # DEBUG
                        }
                        ]
                    })[0]

            self.log(console,"BUILDING RETURN OBJECT")
            returnVal = { 'report_name': reportName,
                      'report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
                      }
            self.log(console,search_tool_name+"_Search DONE")
            return [returnVal]


        # FORMAT DB
        #
        # OLD SYNTAX: formatdb -i $database -o T -p F -> $database.nsq or $database.00.nsq
        # NEW SYNTAX: makeblastdb -in $database -parse_seqids -dbtype prot/nucl -out <basename>
        makeblastdb_cmd = [self.Make_BLAST_DB]

        # check for necessary files
        if not os.path.isfile(self.Make_BLAST_DB):
            raise ValueError("no such file '"+self.Make_BLAST_DB+"'")
        if not os.path.isfile(many_forward_reads_file_path):
            raise ValueError("no such file '"+many_forward_reads_file_path+"'")
        elif not os.path.getsize(many_forward_reads_file_path) > 0:
            raise ValueError("empty file '"+many_forward_reads_file_path+"'")

        makeblastdb_cmd.append('-in')
        makeblastdb_cmd.append(many_forward_reads_file_path)
        makeblastdb_cmd.append('-parse_seqids')
        makeblastdb_cmd.append('-dbtype')
        makeblastdb_cmd.append('nucl')
        makeblastdb_cmd.append('-out')
        makeblastdb_cmd.append(many_forward_reads_file_path)

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
            line = p.stdout.readline()
            if not line: break
            self.log(console, line.replace('\n', ''))

        p.stdout.close()
        p.wait()
        self.log(console, 'return code: ' + str(p.returncode))
        if p.returncode != 0:
            raise ValueError('Error running makeblastdb, return code: '+str(p.returncode) + 
                '\n\n'+ '\n'.join(console))

        # Check for db output
        if not os.path.isfile(many_forward_reads_file_path+".nsq") and not os.path.isfile(many_forward_reads_file_path+".00.nsq"):
            raise ValueError("makeblastdb failed to create DB file '"+many_forward_reads_file_path+".nsq'")
        elif not os.path.getsize(many_forward_reads_file_path+".nsq") > 0 and not os.path.getsize(many_forward_reads_file_path+".00.nsq") > 0:
            raise ValueError("makeblastdb created empty DB file '"+many_forward_reads_file_path+".nsq'")


        ### Construct the BLAST command
        #
        # OLD SYNTAX: $blast -q $q -G $G -E $E -m $m -e $e_value -v $limit -b $limit -K $limit -p tblastx -i $fasta_file -d $database -o $out_file
        # NEW SYNTAX: tblastx -query <queryfile> -db <basename> -out <out_aln_file> -outfmt 0/7 (8 became 7) -evalue <e_value> -dust no (DNA) -seg no (AA) -num_threads <num_cores>
        #
        blast_bin = self.tBLASTx
        blast_cmd = [blast_bin]

        # check for necessary files
        if not os.path.isfile(blast_bin):
            raise ValueError("no such file '"+blast_bin+"'")
        if not os.path.isfile(one_forward_reads_file_path):
            raise ValueError("no such file '"+one_forward_reads_file_path+"'")
        elif not os.path.getsize(one_forward_reads_file_path) > 0:
            raise ValueError("empty file '"+one_forward_reads_file_path+"'")
        if not os.path.isfile(many_forward_reads_file_path):
            raise ValueError("no such file '"+many_forward_reads_file_path+"'")
        elif not os.path.getsize(many_forward_reads_file_path):
            raise ValueError("empty file '"+many_forward_reads_file_path+"'")

        # set the output path
        timestamp = int((datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()*1000)
        output_dir = os.path.join(self.scratch,'output.'+str(timestamp))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        output_aln_file_path = os.path.join(output_dir, 'alnout.txt');
        output_filtered_fasta_file_path = os.path.join(output_dir, 'output_filtered.fna');

        # this is command for basic search mode
        blast_cmd.append('-query')
        blast_cmd.append(one_forward_reads_file_path)
        blast_cmd.append('-db')
        blast_cmd.append(many_forward_reads_file_path)
        blast_cmd.append('-out')
        blast_cmd.append(output_aln_file_path)
        blast_cmd.append('-outfmt')
        blast_cmd.append('7')
        blast_cmd.append('-evalue')
        blast_cmd.append(str(params['e_value']))

        # options
        if 'maxaccepts' in params:
            if params['maxaccepts']:
                blast_cmd.append('-max_target_seqs')
                blast_cmd.append(str(params['maxaccepts']))

        # Run BLAST, capture output as it happens
        #
        self.log(console, 'RUNNING BLAST:')
        self.log(console, '    '+' '.join(blast_cmd))
#        report += "\n"+'running BLAST:'+"\n"
#        report += '    '+' '.join(blast_cmd)+"\n"

        p = subprocess.Popen(blast_cmd, \
                             cwd = self.scratch, \
                             stdout = subprocess.PIPE, \
                             stderr = subprocess.STDOUT, \
                             shell = False)

        while True:
            line = p.stdout.readline()
            if not line: break
            self.log(console, line.replace('\n', ''))

        p.stdout.close()
        p.wait()
        self.log(console, 'return code: ' + str(p.returncode))
        if p.returncode != 0:
            raise ValueError('Error running BLAST, return code: '+str(p.returncode) + 
                '\n\n'+ '\n'.join(console))


        # get query_len for filtering later
        #
        query_len = 0
        with open(one_forward_reads_file_path, 'r', 0) as query_file_handle:
            for line in query_file_handle:
                if line.startswith('>'):
                    continue
                query_len += len(re.sub(r" ","", line.rstrip())) 
        query_len = query_len/3.0  # tBLASTx is nuc-nuc (translated) and reports alnlen in protein length

                
        # Parse the BLAST tabular output and store ids to filter many set to make filtered object to save back to KBase
        #
        self.log(console, 'PARSING BLAST ALIGNMENT OUTPUT')
        if not os.path.isfile(output_aln_file_path):
            raise ValueError("failed to create BLAST output: "+output_aln_file_path)
        elif not os.path.getsize(output_aln_file_path) > 0:
            raise ValueError("created empty file for BLAST output: "+output_aln_file_path)
        hit_seq_ids = dict()
        output_aln_file_handle = open (output_aln_file_path, "r", 0)
        output_aln_buf = output_aln_file_handle.readlines()
        output_aln_file_handle.close()
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
            #self.log(console,'HIT LINE: '+line)  # DEBUG
            hit_info = line.split("\t")
            hit_seq_id     = hit_info[1]
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

            # BLAST SOMETIMES ADDS THIS TO IDs.  NO IDEA WHY, BUT GET RID OF IT!
            if hit_seq_id.startswith('gnl|'):
                hit_seq_id = hit_seq_id[4:]

            try:
                if float(hit_bitscore) > float(high_bitscore_score[hit_seq_id]):
                    self.log(console,"OVERRIDE ID: "+hit_seq_id)  # DEBUG
                    self.log(console,line)  # DEBUG
                    high_bitscore_score[hit_seq_id] = hit_bitscore
                    high_bitscore_ident[hit_seq_id] = hit_ident
                    high_bitscore_alnlen[hit_seq_id] = hit_aln_len
                    high_bitscore_line[hit_seq_id] = line
            except:
                self.log(console,"NEW ID: "+hit_seq_id)  # DEBUG
                self.log(console,line)  # DEBUG
                hit_order.append(hit_seq_id)
                high_bitscore_score[hit_seq_id] = hit_bitscore
                high_bitscore_ident[hit_seq_id] = hit_ident
                high_bitscore_alnlen[hit_seq_id] = hit_aln_len
                high_bitscore_line[hit_seq_id] = line

        for hit_seq_id in hit_order:
            hit_buf.append(high_bitscore_line[hit_seq_id])

            #self.log(console,"HIT_SEQ_ID: '"+hit_seq_id+"'")
            if 'ident_thresh' in params and float(params['ident_thresh']) > float(high_bitscore_ident[hit_seq_id]):
                continue
            if 'bitscore' in params and float(params['bitscore']) > float(high_bitscore_score[hit_seq_id]):
                continue
            if 'overlap_fraction' in params and float(params['overlap_fraction']) > float(high_bitscore_alnlen[hit_seq_id])/float(query_len):
                continue
            
            hit_total += 1
            hit_seq_ids[hit_seq_id] = True
            self.log(console, "HIT: '"+hit_seq_id+"'")  # DEBUG
        

        self.log(console, 'EXTRACTING HITS FROM INPUT')
        self.log(console, 'MANY_TYPE_NAME: '+many_type_name)  # DEBUG


        # SequenceSet input -> SequenceSet output
        #
        if many_type_name == 'SequenceSet':
            seq_total = len(input_many_sequenceSet['sequences'])

            output_sequenceSet = dict()

            if 'sequence_set_id' in input_many_sequenceSet and input_many_sequenceSet['sequence_set_id'] != None:
                output_sequenceSet['sequence_set_id'] = input_many_sequenceSet['sequence_set_id'] + "."+search_tool_name+"_Search_filtered"
            else:
                output_sequenceSet['sequence_set_id'] = search_tool_name+"_Search_filtered"
            if 'description' in input_many_sequenceSet and input_many_sequenceSet['description'] != None:
                output_sequenceSet['description'] = input_many_sequenceSet['description'] + " - "+search_tool_name+"_Search filtered"
            else:
                output_sequenceSet['description'] = search_tool_name+"_Search filtered"

            self.log(console,"ADDING SEQUENCES TO SEQUENCESET")
            output_sequenceSet['sequences'] = []

            for seq_obj in input_many_sequenceSet['sequences']:
                header_id = seq_obj['sequence_id']
                #header_desc = seq_obj['description']
                #sequence_str = seq_obj['sequence']

                try:
                    #self.log(console,"checking '"+header_id+"'")
                    in_filtered_set = hit_seq_ids[header_id]
                    #self.log(console, 'FOUND HIT '+header_id)  # DEBUG
                    output_sequenceSet['sequences'].append(seq_obj)
                except:
                    pass


        # SingleEndLibrary input -> SingleEndLibrary output
        #
        elif many_type_name == 'SingleEndLibrary':

            #  Note: don't use SeqIO.parse because loads everything into memory
            #
#            with open(many_forward_reads_file_path, 'r', -1) as many_forward_reads_file_handle, open(output_filtered_fasta_file_path, 'w', -1) as output_filtered_fasta_file_handle:
            output_filtered_fasta_file_handle = open(output_filtered_fasta_file_path, 'w', -1)
            if many_forward_reads_file_compression == 'gz':
                many_forward_reads_file_handle = gzip.open(many_forward_reads_file_path, 'r', -1)
            else:
                many_forward_reads_file_handle = open(many_forward_reads_file_path, 'r', -1)

            seq_total = 0;
            filtered_seq_total = 0
            last_seq_buf = []
            last_seq_id = None
            last_header = None
            pattern = re.compile('^\S*')
            for line in many_forward_reads_file_handle:
                if line.startswith('>'):
                    #self.log(console, 'LINE: '+line)  # DEBUG
                    seq_total += 1
                    seq_id = line[1:]  # removes '>'
                    seq_id = pattern.findall(seq_id)[0]

                    if last_seq_id != None:
                        #self.log(console, 'ID: '+last_seq_id)  # DEBUG
                        try:
                            in_filtered_set = hit_seq_ids[last_seq_id]
                            #self.log(console, 'FOUND HIT '+last_seq_id)  # DEBUG
                            filtered_seq_total += 1
                            output_filtered_fasta_file_handle.write(last_header)
                            output_filtered_fasta_file_handle.writelines(last_seq_buf)
                        except:
                            pass
                        
                    last_seq_buf = []
                    last_seq_id = seq_id
                    last_header = line
                else:
                    last_seq_buf.append(line)

            if last_seq_id != None:
                #self.log(console, 'ID: '+last_seq_id)  # DEBUG
                try:
                    in_filtered_set = hit_seq_ids[last_seq_id]
                    #self.log(console, 'FOUND HIT: '+last_seq_id)  # DEBUG
                    filtered_seq_total += 1
                    output_filtered_fasta_file_handle.write(last_header)
                    output_filtered_fasta_file_handle.writelines(last_seq_buf)
                except:
                    pass
                
            last_seq_buf = []
            last_seq_id = None
            last_header = None

            many_forward_reads_file_handle.close()
            output_filtered_fasta_file_handle.close()

            if filtered_seq_total != hit_total:
                self.log(console,'hits in BLAST alignment output '+str(hit_total)+' != '+str(filtered_seq_total)+' matched sequences in input file')
                raise ValueError('hits in BLAST alignment output '+str(hit_total)+' != '+str(filtered_seq_total)+' matched sequences in input file')


        # FeatureSet input -> FeatureSet output
        #
        elif many_type_name == 'FeatureSet':
            seq_total = len(input_many_featureSet['elements'].keys())

            output_featureSet = dict()
            if 'description' in input_many_featureSet and input_many_featureSet['description'] != None:
                output_featureSet['description'] = input_many_featureSet['description'] + " - "+search_tool_name+"_Search filtered"
            else:
                output_featureSet['description'] = search_tool_name+"_Search filtered"
            #output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()

            fId_list = input_many_featureSet['elements'].keys()
            self.log(console,"ADDING FEATURES TO FEATURESET")
            for fId in sorted(fId_list):
                for genome_ref in input_many_featureSet['elements'][fId]:
                    try:
                        #self.log(console,"checking '"+fId+"'")
                        #in_filtered_set = hit_seq_ids[fId]
                        in_filtered_set = hit_seq_ids[genome_ref+genome_id_feature_id_delim+fId]
                        #self.log(console, 'FOUND HIT '+fId)  # DEBUG
                        try:
                            this_genome_ref_list = output_featureSet['elements'][fId]
                        except:
                            output_featureSet['elements'][fId] = []
                        output_featureSet['elements'][fId].append(genome_ref)
                    except:
                        pass

        # Parse Genome hits into FeatureSet
        #
        elif many_type_name == 'Genome':
            seq_total = 0
            output_featureSet = dict()
#            if 'scientific_name' in input_many_genome and input_many_genome['scientific_name'] != None:
#                output_featureSet['description'] = input_many_genome['scientific_name'] + " - "+search_tool_name+"_Search filtered"
#            else:
#                output_featureSet['description'] = search_tool_name+"_Search filtered"
            output_featureSet['description'] = search_tool_name+"_Search filtered"
            #output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()
            for fid in feature_ids:
                seq_total += 1
                try:
                    in_filtered_set = hit_seq_ids[fid]
                    #output_featureSet['element_ordering'].append(fid)
                    output_featureSet['elements'][fid] = [input_many_ref]
                except:
                    pass

        # Parse GenomeSet hits into FeatureSet
        #
        elif many_type_name == 'GenomeSet':
            seq_total = 0

            output_featureSet = dict()
            if 'description' in input_many_genomeSet and input_many_genomeSet['description'] != None:
                output_featureSet['description'] = input_many_genomeSet['description'] + " - "+search_tool_name+"_Search filtered"
            else:
                output_featureSet['description'] = search_tool_name+"_Search filtered"
            #output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()

            self.log(console,"READING HITS FOR GENOMES")  # DEBUG
            for genome_id in feature_ids_by_genome_id.keys():
                self.log(console,"READING HITS FOR GENOME "+genome_id)  # DEBUG
                genome_ref = input_many_genomeSet['elements'][genome_id]['ref']
                for feature_id in feature_ids_by_genome_id[genome_id]:
                    seq_total += 1
                    try:
                        #in_filtered_set = hit_seq_ids[feature['id']]
                        in_filtered_set = hit_seq_ids[genome_ref+genome_id_feature_id_delim+feature_id]
                        #in_filtered_set = hit_seq_ids[feature_id]
                        #self.log(console, 'FOUND HIT: '+feature['id'])  # DEBUG
                        #output_featureSet['element_ordering'].append(feature['id'])
                        try:
                            this_genome_ref_list = output_featureSet['elements'][feature_id]
                        except:
                            output_featureSet['elements'][feature_id] = []
                        output_featureSet['elements'][feature_id].append(genome_ref)
                    except:
                        pass


        # load the method provenance from the context object
        #
        self.log(console,"SETTING PROVENANCE")  # DEBUG
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects'] = []
        provenance[0]['input_ws_objects'].append(input_one_ref)
        provenance[0]['input_ws_objects'].append(input_many_ref)
        provenance[0]['service'] = 'kb_blast'
        provenance[0]['method'] = search_tool_name+'_Search'


        # Upload results
        #
        if len(invalid_msgs) == 0:
            self.log(console,"UPLOADING RESULTS")  # DEBUG

            # input SingleEndLibrary -> upload SingleEndLibrary
            #
            if many_type_name == 'SingleEndLibrary':
            
                self.upload_SingleEndLibrary_to_shock_and_ws (ctx,
                                                          console,  # DEBUG
                                                          params['workspace_name'],
                                                          params['output_filtered_name'],
                                                          output_filtered_fasta_file_path,
                                                          provenance,
                                                          sequencing_tech
                                                         )

            # input many SequenceSet -> save SequenceSet
            #
            elif many_type_name == 'SequenceSet':
                new_obj_info = ws.save_objects({
                            'workspace': params['workspace_name'],
                            'objects':[{
                                    'type': 'KBaseSequences.SequenceSet',
                                    'data': output_sequenceSet,
                                    'name': params['output_filtered_name'],
                                    'meta': {},
                                    'provenance': provenance
                                }]
                        })[0]

            else:  # input FeatureSet, Genome, and GenomeSet -> upload FeatureSet output
                new_obj_info = ws.save_objects({
                            'workspace': params['workspace_name'],
                            'objects':[{
                                    'type': 'KBaseCollections.FeatureSet',
                                    'data': output_featureSet,
                                    'name': params['output_filtered_name'],
                                    'meta': {},
                                    'provenance': provenance
                                }]
                        })[0]

        # build output report object
        #
        self.log(console,"BUILDING REPORT")  # DEBUG
        if len(invalid_msgs) == 0:
            report += 'sequences in many set: '+str(seq_total)+"\n"
            report += 'sequences in hit set:  '+str(hit_total)+"\n"
            report += "\n"
            for line in hit_buf:
                report += line
            reportObj = {
                'objects_created':[{'ref':params['workspace_name']+'/'+params['output_filtered_name'], 'description':search_tool_name+'_Search hits'}],
                'text_message':report
                }
        else:
            report += "FAILURE\n\n"+"\n".join(invalid_msgs)+"\n"
            reportObj = {
                'objects_created':[],
                'text_message':report
                }

        reportName = 'blast_report_'+str(uuid.uuid4())
        report_obj_info = ws.save_objects({
#                'id':info[6],
                'workspace':params['workspace_name'],
                'objects':[
                    {
                        'type':'KBaseReport.Report',
                        'data':reportObj,
                        'name':reportName,
                        'meta':{},
                        'hidden':1,
                        'provenance':provenance
                    }
                ]
            })[0]

        self.log(console,"BUILDING RETURN OBJECT")
#        returnVal = { 'output_report_name': reportName,
#                      'output_report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
#                      'output_filtered_ref': params['workspace_name']+'/'+params['output_filtered_name']
#                      }
        returnVal = { 'report_name': reportName,
                      'report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
                      }
        self.log(console,search_tool_name+"_Search DONE")
        #END tBLASTx_Search

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method tBLASTx_Search return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def psiBLAST_msa_start_Search(self, ctx, params):
        """
        :param params: instance of type "BLAST_Params" (BLAST Input Params)
           -> structure: parameter "workspace_name" of type "workspace_name"
           (** The workspace object refs are of form: ** **    objects =
           ws.get_objects([{'ref':
           params['workspace_id']+'/'+params['obj_name']}]) ** ** "ref" means
           the entire name combining the workspace id and the object name **
           "id" is a numerical identifier of the workspace or object, and
           should just be used for workspace ** "name" is a string identifier
           of a workspace or object.  This is received from Narrative.),
           parameter "input_one_sequence" of type "sequence", parameter
           "input_one_ref" of type "data_obj_ref", parameter "input_many_ref"
           of type "data_obj_ref", parameter "input_msa_ref" of type
           "data_obj_ref", parameter "output_one_name" of type
           "data_obj_name", parameter "output_filtered_name" of type
           "data_obj_name", parameter "ident_thresh" of Double, parameter
           "e_value" of Double, parameter "bitscore" of Double, parameter
           "overlap_fraction" of Double, parameter "maxaccepts" of Double,
           parameter "rounds" of Double
        :returns: instance of type "BLAST_Output" (BLAST Output) ->
           structure: parameter "report_name" of type "data_obj_name",
           parameter "report_ref" of type "data_obj_ref"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN psiBLAST_msa_start_Search
        console = []
        invalid_msgs = []
        search_tool_name = 'psiBLAST_msa_start'
        self.log(console,'Running '+search_tool_name+'_Search with params=')
        self.log(console, "\n"+pformat(params))
        report = ''
#        report = 'Running '+search_tool_name+'_Search with params='
#        report += "\n"+pformat(params)
        appropriate_sequence_found_in_one_input = False
        appropriate_sequence_found_in_MSA_input = False
        appropriate_sequence_found_in_many_input = False


        #### do some basic checks
        #
        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
        if 'input_one_ref' not in params:
            raise ValueError('input_one_ref parameter is required')
        if 'input_msa_ref' not in params:
            raise ValueError('input_msa_ref parameter is required')
        if 'input_many_ref' not in params:
            raise ValueError('input_many_ref parameter is required')
        if 'output_filtered_name' not in params:
            raise ValueError('output_filtered_name parameter is required')


        # set local names
        input_one_ref = params['input_one_ref']
        input_msa_ref = params['input_msa_ref']
        input_many_ref = params['input_many_ref']
        

        #### Get the input_one object
        ##
        input_one_id = None
        if 'input_one_ref' in params and params['input_one_ref'] != None:
            try:
                ws = workspaceService(self.workspaceURL, token=ctx['token'])
                objects = ws.get_objects([{'ref': params['input_one_ref']}])
                input_one_data = objects[0]['data']
                input_one_name = str(objects[0]['info'][1])
                info = objects[0]['info']

                one_type_name = info[2].split('.')[1].split('-')[0]
            except Exception as e:
                raise ValueError('Unable to fetch input_one_name object from workspace: ' + str(e))
                #to get the full stack trace: traceback.format_exc()


            # Handle overloading (input_one can be Feature, or FeatureSet)
            #
            if one_type_name == 'FeatureSet':
                # retrieve sequences for features
                #input_one_featureSet = input_one_data
                genome_id_feature_id_delim = '.f:'
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
                    'residue_type':        'protein',
                    'feature_type':        'CDS',
                    'record_id_pattern':   '%%genome_ref%%'+genome_id_feature_id_delim+'%%feature_id%%',
                    'record_desc_pattern': '[%%genome_ref%%]',
                    'case':                'upper',
                    'linewrap':            50,
                    'merge_fasta_files':   'TRUE'
                    }

                #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
                DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'])
                FeatureSetToFASTA_retVal = DOTFU.FeatureSetToFASTA (FeatureSetToFASTA_params)
                one_forward_reads_file_path = FeatureSetToFASTA_retVal['fasta_file_path']
                if len(FeatureSetToFASTA_retVal['feature_ids_by_genome_ref'].keys()) > 0:
                    appropriate_sequence_found_in_one_input = True

                # DEBUG
                #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
                #self.log(console, "FeatureSetToFasta() took "+str(end_time-beg_time)+" secs")


            # Feature
            #
            elif one_type_name == 'Feature':
                # export feature to FASTA file
                feature = input_one_data
                input_one_feature_id = feature['id']
                one_forward_reads_file_path = os.path.join(self.scratch, input_one_name+".fasta")
                self.log(console, 'writing fasta file: '+one_forward_reads_file_path)
                # psiBLAST is prot-prot
                #record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description='['+feature['genome_id']+']'+' '+feature['function'])
                #if feature['type'] != 'CDS':
                #    self.log(console,input_one_name+" feature type must be CDS")
                #    self.log(invalid_msgs,input_one_name+" feature type must be CDS")
                if 'protein_translation' not in feature or feature['protein_translation'] == None:
                    #self.log(console,"bad CDS Feature "+input_one_name+": no protein_translation found")
                    #raise ValueError ("bad CDS Feature "input_one_name+": no protein_translation found")
                    self.log(console,input_one_name+" feature type must be CDS")
                    self.log(invalid_msgs,input_one_name+" feature type must be CDS")
                else:
                    appropriate_sequence_found_in_one_input = True
                    record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description='['+feature['genome_id']+']'+' '+feature['function'])
                    SeqIO.write([record], one_forward_reads_file_path, "fasta")
                    appropriate_sequence_found_in_one_input = True
            else:
                raise ValueError('Cannot yet handle input_one type of: '+type_name)            
        else:
            raise ValueError('Must define either input_one_sequence or input_one_name')


        #### Get the input_msa object
        ##
        if input_one_feature_id == None:
            self.log(invalid_msgs,"input_one_feature_id was not obtained from Query Object: "+input_one_name)
        master_row_idx = 0
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            objects = ws.get_objects([{'ref': input_msa_ref}])
            input_msa_data = objects[0]['data']
            info = objects[0]['info']
            input_msa_name = str(info[1])
            input_msa_type = info[2].split('.')[1].split('-')[0]

        except Exception as e:
            raise ValueError('Unable to fetch input_msa_name object from workspace: ' + str(e))
            #to get the full stack trace: traceback.format_exc()

        if input_msa_type == 'MSA':
            MSA_in = input_msa_data
            row_order = []
            default_row_labels = dict()
            if 'row_order' in MSA_in.keys():
                row_order = MSA_in['row_order']
            else:
                row_order = sorted(MSA_in['alignment'].keys())

            if 'default_row_labels' in MSA_in.keys():
                default_row_labels = MSA_in['default_row_labels']
            else:
                for row_id in row_order:
                    default_row_labels[row_id] = row_id

            # determine row index of query sequence
            for row_id in row_order:
                master_row_idx += 1
                if row_id == input_one_feature_id:
                    break
            if master_row_idx == 0:
                self.log(invalid_msgs,"Failed to find query id "+input_one_feature_id+" from Query Object "+input_one_name+" within MSA: "+input_msa_name)

            
            # export features to Clustal-esque file that PSI-BLAST likes
            input_MSA_file_path = os.path.join(self.scratch, input_msa_name+".fasta")
            self.log(console, 'writing MSA file: '+input_MSA_file_path)
            records = []
            longest_row_id_len = 0
            for row_id in row_order:
                if len(row_id) > longest_row_id_len:
                    longest_row_id_len = len(row_id)
            for row_id in row_order:
                #self.log(console,"row_id: '"+row_id+"'")  # DEBUG
                #self.log(console,"alignment: '"+MSA_in['alignment'][row_id]+"'")  # DEBUG
                # using SeqIO makes multiline sequences.  We want Clustal-esque, but we'll not break them up and hope PSI-BLAST is happy
                #record = SeqRecord(Seq(MSA_in['alignment'][row_id]), id=row_id, description=default_row_labels[row_id])
                #records.append(record)
                #SeqIO.write(records, input_MSA_file_path, "fasta")
                padding = ''
                for i in range(0,longest_row_id_len-len(row_id)):
                    padding += ' '
                records.append(row_id + padding + "\t" +
                               MSA_in['alignment'][row_id]
                               )
            with open(input_MSA_file_path,'w',0) as input_MSA_file_handle:
                input_MSA_file_handle.write("\n".join(records)+"\n")


            # Determine whether nuc or protein sequences
            #
            NUC_MSA_pattern = re.compile("^[\.\-_ACGTUXNRYSWKMBDHVacgtuxnryswkmbdhv \t\n]+$")
            all_seqs_nuc = True
            for row_id in row_order:
                #self.log(console, row_id+": '"+MSA_in['alignment'][row_id]+"'")
                if NUC_MSA_pattern.match(MSA_in['alignment'][row_id]) == None:
                    all_seqs_nuc = False
                    break
                else:
                    appropriate_sequence_found_in_MSA_input = True

        # Missing proper input_type
        #
        else:
            raise ValueError('Cannot yet handle input_name type of: '+type_name)


        #### Get the input_many object
        ##
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            objects = ws.get_objects([{'ref': input_many_ref}])
            input_many_data = objects[0]['data']
            info = objects[0]['info']
            input_many_name = str(info[1])
            many_type_name = info[2].split('.')[1].split('-')[0]

        except Exception as e:
            raise ValueError('Unable to fetch input_many_name object from workspace: ' + str(e))
            #to get the full stack trace: traceback.format_exc()


        # Handle overloading (input_many can be SequenceSet, FeatureSet, Genome, or GenomeSet)
        #
        if many_type_name == 'SequenceSet':
            try:
                input_many_sequenceSet = input_many_data
            except Exception as e:
                print(traceback.format_exc())
                raise ValueError('Unable to get SequenceSet: ' + str(e))

            header_id = input_many_sequenceSet['sequences'][0]['sequence_id']
            many_forward_reads_file_path = os.path.join(self.scratch, header_id+'.fasta')
            many_forward_reads_file_handle = open(many_forward_reads_file_path, 'w', 0)
            self.log(console, 'writing reads file: '+str(many_forward_reads_file_path))

            for seq_obj in input_many_sequenceSet['sequences']:
                header_id = seq_obj['sequence_id']
                sequence_str = seq_obj['sequence']

                PROT_pattern = re.compile("^[acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWYxX ]+$")
                #DNA_pattern = re.compile("^[acgtuACGTUnryNRY ]+$")
                if not PROT_pattern.match(sequence_str):
                    self.log(invalid_msgs,"BAD record for sequence_id: "+header_id+"\n"+sequence_str+"\n")
                    continue
                appropriate_sequence_found_in_many_input = True
                many_forward_reads_file_handle.write('>'+header_id+"\n")
                many_forward_reads_file_handle.write(sequence_str+"\n")
            many_forward_reads_file_handle.close();
            self.log(console, 'done')


        # FeatureSet
        #
        if many_type_name == 'FeatureSet':
            # retrieve sequences for features
            input_many_featureSet = input_many_data
            genome_id_feature_id_delim = '.f:'
            many_forward_reads_file_dir = self.scratch
            many_forward_reads_file = input_many_name+".fasta"

            # DEBUG
            #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            FeatureSetToFASTA_params = {
                'featureSet_ref':      input_many_ref,
                'file':                many_forward_reads_file,
                'dir':                 many_forward_reads_file_dir,
                'console':             console,
                'invalid_msgs':        invalid_msgs,
                'residue_type':        'protein',
                'feature_type':        'CDS',
                'record_id_pattern':   '%%genome_ref%%'+genome_id_feature_id_delim+'%%feature_id%%',
                'record_desc_pattern': '[%%genome_ref%%]',
                'case':                'upper',
                'linewrap':            50,
                'merge_fasta_files':   'TRUE'
                }

            #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'])
            FeatureSetToFASTA_retVal = DOTFU.FeatureSetToFASTA (FeatureSetToFASTA_params)
            many_forward_reads_file_path = FeatureSetToFASTA_retVal['fasta_file_path']
            feature_ids_by_genome_ref = FeatureSetToFASTA_retVal['feature_ids_by_genome_ref']
            if len(feature_ids_by_genome_ref.keys()) > 0:
                appropriate_sequence_found_in_many_input = True

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "FeatureSetToFasta() took "+str(end_time-beg_time)+" secs")


        # Genome
        #
        elif many_type_name == 'Genome':
            many_forward_reads_file_dir = self.scratch
            many_forward_reads_file = input_many_name+".fasta"

            # DEBUG
            #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            GenomeToFASTA_params = {
                'genome_ref':          input_many_ref,
                'file':                many_forward_reads_file,
                'dir':                 many_forward_reads_file_dir,
                'console':             console,
                'invalid_msgs':        invalid_msgs,
                'residue_type':        'protein',
                'feature_type':        'CDS',
                'record_id_pattern':   '%%feature_id%%',
                'record_desc_pattern': '[%%genome_id%%]',
                'case':                'upper',
                'linewrap':            50
                }

            #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'])
            GenomeToFASTA_retVal = DOTFU.GenomeToFASTA (GenomeToFASTA_params)
            many_forward_reads_file_path = GenomeToFASTA_retVal['fasta_file_path']
            feature_ids = GenomeToFASTA_retVal['feature_ids']
            if len(feature_ids) > 0:
                appropriate_sequence_found_in_many_input = True

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "Genome2Fasta() took "+str(end_time-beg_time)+" secs")
            

        # GenomeSet
        #
        elif many_type_name == 'GenomeSet':
            input_many_genomeSet = input_many_data
            genome_id_feature_id_delim = '.f:'
            many_forward_reads_file_dir = self.scratch
            many_forward_reads_file = input_many_name+".fasta"

            # DEBUG
            #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            GenomeSetToFASTA_params = {
                'genomeSet_ref':       input_many_ref,
                'file':                many_forward_reads_file,
                'dir':                 many_forward_reads_file_dir,
                'console':             console,
                'invalid_msgs':        invalid_msgs,
                'residue_type':        'protein',
                'feature_type':        'CDS',
                'record_id_pattern':   '%%genome_ref%%'+genome_id_feature_id_delim+'%%feature_id%%',
                'record_desc_pattern': '[%%genome_ref%%]',
                'case':                'upper',
                'linewrap':            50,
                'merge_fasta_files':   'TRUE'
                }

            #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'])
            GenomeSetToFASTA_retVal = DOTFU.GenomeSetToFASTA (GenomeSetToFASTA_params)
            many_forward_reads_file_path = GenomeSetToFASTA_retVal['fasta_file_path_list'][0]
            feature_ids_by_genome_id = GenomeSetToFASTA_retVal['feature_ids_by_genome_id']
            if len(feature_ids_by_genome_id.keys()) > 0:
                appropriate_sequence_found_in_many_input = True

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "FeatureSetToFasta() took "+str(end_time-beg_time)+" secs")


        # Missing proper input_many_type
        #
        else:
            raise ValueError('Cannot yet handle input_many type of: '+type_name)            


        # check for failed input file creation
        #
        if not appropriate_sequence_found_in_one_input:
            self.log(invalid_msgs,"no protein sequences found in '"+input_one_name+"'")
        if not appropriate_sequence_found_in_MSA_input:
            self.log(invalid_msgs,"no protein sequences found in '"+input_msa_name+"'")
        if not appropriate_sequence_found_in_many_input:
            self.log(invalid_msgs,"no protein sequences found in '"+input_many_name+"'")


        # input data failed validation.  Need to return
        #
        if len(invalid_msgs) > 0:

            # load the method provenance from the context object
            #
            self.log(console,"SETTING PROVENANCE")  # DEBUG
            provenance = [{}]
            if 'provenance' in ctx:
                provenance = ctx['provenance']
            # add additional info to provenance here, in this case the input data object reference
            provenance[0]['input_ws_objects'] = []
            provenance[0]['input_ws_objects'].append(input_one_ref)
            provenance[0]['input_ws_objects'].append(input_msa_ref)
            provenance[0]['input_ws_objects'].append(input_many_ref)
            provenance[0]['service'] = 'kb_blast'
            provenance[0]['method'] = search_tool_name+'_Search'


            # build output report object
            #
            self.log(console,"BUILDING REPORT")  # DEBUG
            report += "FAILURE:\n\n"+"\n".join(invalid_msgs)+"\n"
            reportObj = {
                'objects_created':[],
                'text_message':report
                }

            reportName = 'blast_report_'+str(uuid.uuid4())
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            report_obj_info = ws.save_objects({
                    #'id':info[6],
                    'workspace':params['workspace_name'],
                    'objects':[
                        {
                        'type':'KBaseReport.Report',
                        'data':reportObj,
                        'name':reportName,
                        'meta':{},
                        'hidden':1,
                        'provenance':provenance  # DEBUG
                        }
                        ]
                    })[0]

            self.log(console,"BUILDING RETURN OBJECT")
            returnVal = { 'report_name': reportName,
                      'report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
                      }
            self.log(console,search_tool_name+"_Search DONE")
            return [returnVal]


        # FORMAT DB
        #
        # OLD SYNTAX: formatdb -i $database -o T -p F -> $database.nsq or $database.00.nsq
        # NEW SYNTAX: makeblastdb -in $database -parse_seqids -dbtype prot/nucl -out <basename>
        makeblastdb_cmd = [self.Make_BLAST_DB]

        # check for necessary files
        if not os.path.isfile(self.Make_BLAST_DB):
            raise ValueError("no such file '"+self.Make_BLAST_DB+"'")
        if not os.path.isfile(many_forward_reads_file_path):
            raise ValueError("no such file '"+many_forward_reads_file_path+"'")
        elif not os.path.getsize(many_forward_reads_file_path) > 0:
            raise ValueError("empty file '"+many_forward_reads_file_path+"'")

        makeblastdb_cmd.append('-in')
        makeblastdb_cmd.append(many_forward_reads_file_path)
        makeblastdb_cmd.append('-parse_seqids')
        makeblastdb_cmd.append('-dbtype')
        makeblastdb_cmd.append('prot')
        makeblastdb_cmd.append('-out')
        makeblastdb_cmd.append(many_forward_reads_file_path)

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
            line = p.stdout.readline()
            if not line: break
            self.log(console, line.replace('\n', ''))

        p.stdout.close()
        p.wait()
        self.log(console, 'return code: ' + str(p.returncode))
        if p.returncode != 0:
            raise ValueError('Error running makeblastdb, return code: '+str(p.returncode) + 
                '\n\n'+ '\n'.join(console))

        # Check for db output
        if not os.path.isfile(many_forward_reads_file_path+".psq") and not os.path.isfile(many_forward_reads_file_path+".00.psq"):
            raise ValueError("makeblastdb failed to create DB file '"+many_forward_reads_file_path+".psq'")
        elif not os.path.getsize(many_forward_reads_file_path+".psq") > 0 and not os.path.getsize(many_forward_reads_file_path+".00.psq") > 0:
            raise ValueError("makeblastdb created empty DB file '"+many_forward_reads_file_path+".psq'")


        ### Construct the psiBLAST command
        #
        # OLD SYNTAX: blastpgp -j <rounds> -h <e_value_matrix> -z <database_size:e.g. 1e8> -q $q -G $G -E $E -m $m -e $e_value -v $limit -b $limit -K $limit -i $fasta_file -B <msa_file> -d $database -o $out_file
        # NEW SYNTAX: psiblast -in_msa <msa_queryfile> -msa_master_idx <row_n> -db <basename> -out <out_aln_file> -outfmt 0/7 (8 became 7) -evalue <e_value> -dust no (DNA) -seg no (AA) -num_threads <num_cores>
        #
        blast_bin = self.psiBLAST
        blast_cmd = [blast_bin]

        # check for necessary files
        if not os.path.isfile(blast_bin):
            raise ValueError("no such file '"+blast_bin+"'")
        #if not os.path.isfile(one_forward_reads_file_path):
        #    raise ValueError("no such file '"+one_forward_reads_file_path+"'")
        #elif not os.path.getsize(one_forward_reads_file_path) > 0:
        #    raise ValueError("empty file '"+one_forward_reads_file_path+"'")
        if not os.path.isfile(input_MSA_file_path):
            raise ValueError("no such file '"+input_MSA_file_path+"'")
        elif not os.path.getsize(input_MSA_file_path):
            raise ValueError("empty file '"+input_MSA_file_path+"'")
        if not os.path.isfile(many_forward_reads_file_path):
            raise ValueError("no such file '"+many_forward_reads_file_path+"'")
        elif not os.path.getsize(many_forward_reads_file_path):
            raise ValueError("empty file '"+many_forward_reads_file_path+"'")

        # set the output path
        timestamp = int((datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()*1000)
        output_dir = os.path.join(self.scratch,'output.'+str(timestamp))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        output_aln_file_path = os.path.join(output_dir, 'alnout.txt');
        output_filtered_fasta_file_path = os.path.join(output_dir, 'output_filtered.fna');

        # this is command for basic search mode
#        blast_cmd.append('-query')
#        blast_cmd.append(one_forward_reads_file_path)
        blast_cmd.append('-in_msa')
        blast_cmd.append(input_MSA_file_path)
        blast_cmd.append('-msa_master_idx')
        blast_cmd.append(str(master_row_idx))
        blast_cmd.append('-db')
        blast_cmd.append(many_forward_reads_file_path)
        blast_cmd.append('-out')
        blast_cmd.append(output_aln_file_path)
        blast_cmd.append('-outfmt')
        blast_cmd.append('7')
        blast_cmd.append('-evalue')
        blast_cmd.append(str(params['e_value']))

        # options
        if 'maxaccepts' in params:
            if params['maxaccepts']:
                blast_cmd.append('-max_target_seqs')
                blast_cmd.append(str(params['maxaccepts']))

        # Run BLAST, capture output as it happens
        #
        self.log(console, 'RUNNING BLAST:')
        self.log(console, '    '+' '.join(blast_cmd))
#        report += "\n"+'running BLAST:'+"\n"
#        report += '    '+' '.join(blast_cmd)+"\n"

        p = subprocess.Popen(blast_cmd, \
                             cwd = self.scratch, \
                             stdout = subprocess.PIPE, \
                             stderr = subprocess.STDOUT, \
                             shell = False)

        while True:
            line = p.stdout.readline()
            if not line: break
            self.log(console, line.replace('\n', ''))

        p.stdout.close()
        p.wait()
        self.log(console, 'return code: ' + str(p.returncode))
        if p.returncode != 0:
            raise ValueError('Error running BLAST, return code: '+str(p.returncode) + 
                '\n\n'+ '\n'.join(console))


        # get query_len for filtering later
        #
        query_len = 0
        with open(one_forward_reads_file_path, 'r', 0) as query_file_handle:
            for line in query_file_handle:
                if line.startswith('>'):
                    continue
                query_len += len(re.sub(r" ","", line.rstrip())) 
        

        # Parse the BLAST tabular output and store ids to filter many set to make filtered object to save back to KBase
        #
        self.log(console, 'PARSING BLAST ALIGNMENT OUTPUT')
        if not os.path.isfile(output_aln_file_path):
            raise ValueError("failed to create BLAST output: "+output_aln_file_path)
        elif not os.path.getsize(output_aln_file_path) > 0:
            raise ValueError("created empty file for BLAST output: "+output_aln_file_path)
        hit_seq_ids = dict()
        output_aln_file_handle = open (output_aln_file_path, "r", 0)
        output_aln_buf = output_aln_file_handle.readlines()
        output_aln_file_handle.close()
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
            #self.log(console,'HIT LINE: '+line)  # DEBUG
            hit_info = line.split("\t")
            hit_seq_id     = hit_info[1]
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

            # BLAST SOMETIMES ADDS THIS TO IDs.  NO IDEA WHY, BUT GET RID OF IT!
            if hit_seq_id.startswith('gnl|'):
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

        for hit_seq_id in hit_order:
            hit_buf.append(high_bitscore_line[hit_seq_id])

            #self.log(console,"HIT_SEQ_ID: '"+hit_seq_id+"'")
            if 'ident_thresh' in params and float(params['ident_thresh']) > float(high_bitscore_ident[hit_seq_id]):
                continue
            if 'bitscore' in params and float(params['bitscore']) > float(high_bitscore_score[hit_seq_id]):
                continue
            if 'overlap_fraction' in params and float(params['overlap_fraction']) > float(high_bitscore_alnlen[hit_seq_id])/float(query_len):
                continue
            
            hit_total += 1
            hit_seq_ids[hit_seq_id] = True
            self.log(console, "HIT: '"+hit_seq_id+"'")  # DEBUG
        

        self.log(console, 'EXTRACTING HITS FROM INPUT')
        self.log(console, 'MANY_TYPE_NAME: '+many_type_name)  # DEBUG


        # SequenceSet input -> SequenceSet output
        #
        if many_type_name == 'SequenceSet':
            seq_total = len(input_many_sequenceSet['sequences'])

            output_sequenceSet = dict()

            if 'sequence_set_id' in input_many_sequenceSet and input_many_sequenceSet['sequence_set_id'] != None:
                output_sequenceSet['sequence_set_id'] = input_many_sequenceSet['sequence_set_id'] + "."+search_tool_name+"_Search_filtered"
            else:
                output_sequenceSet['sequence_set_id'] = search_tool_name+"_Search_filtered"
            if 'description' in input_many_sequenceSet and input_many_sequenceSet['description'] != None:
                output_sequenceSet['description'] = input_many_sequenceSet['description'] + " - "+search_tool_name+"_Search filtered"
            else:
                output_sequenceSet['description'] = search_tool_anme+"_Search filtered"

            self.log(console,"ADDING SEQUENCES TO SEQUENCESET")
            output_sequenceSet['sequences'] = []

            for seq_obj in input_many_sequenceSet['sequences']:
                header_id = seq_obj['sequence_id']
                #header_desc = seq_obj['description']
                #sequence_str = seq_obj['sequence']

                try:
                    #self.log(console,"checking '"+header_id+"'")
                    in_filtered_set = hit_seq_ids[header_id]
                    #self.log(console, 'FOUND HIT '+header_id)  # DEBUG
                    output_sequenceSet['sequences'].append(seq_obj)
                except:
                    pass


        # FeatureSet input -> FeatureSet output
        #
        elif many_type_name == 'FeatureSet':
            seq_total = len(input_many_featureSet['elements'].keys())

            output_featureSet = dict()
            if 'description' in input_many_featureSet and input_many_featureSet['description'] != None:
                output_featureSet['description'] = input_many_featureSet['description'] + " - "+search_tool_name+"_Search filtered"
            else:
                output_featureSet['description'] = search_tool_name+"_Search filtered"
            #output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()

            fId_list = input_many_featureSet['elements'].keys()
            self.log(console,"ADDING FEATURES TO FEATURESET")
            for fId in sorted(fId_list):
                for genome_ref in input_many_featureSet['elements'][fId]:
                    try:
                        #self.log(console,"checking '"+fId+"'")
                        #in_filtered_set = hit_seq_ids[fId]
                        in_filtered_set = hit_seq_ids[genome_ref+genome_id_feature_id_delim+fId]
                        #self.log(console, 'FOUND HIT '+fId)  # DEBUG
                        try:
                            this_genome_ref_list = output_featureSet['elements'][fId]
                        except:
                            output_featureSet['elements'][fId] = []
                        output_featureSet['elements'][fId].append(genome_ref)
                    except:
                        pass

        # Parse Genome hits into FeatureSet
        #
        elif many_type_name == 'Genome':
            seq_total = 0
            output_featureSet = dict()
#            if 'scientific_name' in input_many_genome and input_many_genome['scientific_name'] != None:
#                output_featureSet['description'] = input_many_genome['scientific_name'] + " - "+search_tool_name+"_Search filtered"
#            else:
#                output_featureSet['description'] = search_tool_name+"_Search filtered"
            output_featureSet['description'] = search_tool_name+"_Search filtered"
            #output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()
            for fid in feature_ids:
                seq_total += 1
                try:
                    in_filtered_set = hit_seq_ids[fid]
                    #output_featureSet['element_ordering'].append(fid)
                    output_featureSet['elements'][fid] = [input_many_ref]
                except:
                    pass

        # Parse GenomeSet hits into FeatureSet
        #
        elif many_type_name == 'GenomeSet':
            seq_total = 0

            output_featureSet = dict()
            if 'description' in input_many_genomeSet and input_many_genomeSet['description'] != None:
                output_featureSet['description'] = input_many_genomeSet['description'] + " - "+search_tool_name+"_Search filtered"
            else:
                output_featureSet['description'] = search_tool_name+"_Search filtered"
            #output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()

            self.log(console,"READING HITS FOR GENOMES")  # DEBUG
            for genome_id in feature_ids_by_genome_id.keys():
                self.log(console,"READING HITS FOR GENOME "+genome_id)  # DEBUG
                genome_ref = input_many_genomeSet['elements'][genome_id]['ref']
                for feature_id in feature_ids_by_genome_id[genome_id]:
                    seq_total += 1
                    try:
                        #in_filtered_set = hit_seq_ids[feature['id']]
                        in_filtered_set = hit_seq_ids[genome_ref+genome_id_feature_id_delim+feature_id]
                        #in_filtered_set = hit_seq_ids[feature_id]
                        #self.log(console, 'FOUND HIT: '+feature['id'])  # DEBUG
                        #output_featureSet['element_ordering'].append(feature['id'])
                        try:
                            this_genome_ref_list = output_featureSet['elements'][feature_id]
                        except:
                            output_featureSet['elements'][feature_id] = []
                        output_featureSet['elements'][feature_id].append(genome_ref)
                    except:
                        pass


        # load the method provenance from the context object
        #
        self.log(console,"SETTING PROVENANCE")  # DEBUG
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects'] = []
        provenance[0]['input_ws_objects'].append(input_one_ref)
        provenance[0]['input_ws_objects'].append(input_msa_ref)
        provenance[0]['input_ws_objects'].append(input_many_ref)
        provenance[0]['service'] = 'kb_blast'
        provenance[0]['method'] = search_tool_name+'_Search'


        # Upload results
        #
        if len(invalid_msgs) == 0:
            self.log(console,"UPLOADING RESULTS")  # DEBUG

            # input many SequenceSet -> save SequenceSet
            #
            if many_type_name == 'SequenceSet':
                new_obj_info = ws.save_objects({
                            'workspace': params['workspace_name'],
                            'objects':[{
                                    'type': 'KBaseSequences.SequenceSet',
                                    'data': output_sequenceSet,
                                    'name': params['output_filtered_name'],
                                    'meta': {},
                                    'provenance': provenance
                                }]
                        })[0]

            else:  # input FeatureSet, Genome, and GenomeSet -> upload FeatureSet output
                new_obj_info = ws.save_objects({
                            'workspace': params['workspace_name'],
                            'objects':[{
                                    'type': 'KBaseCollections.FeatureSet',
                                    'data': output_featureSet,
                                    'name': params['output_filtered_name'],
                                    'meta': {},
                                    'provenance': provenance
                                }]
                        })[0]

        # build output report object
        #
        self.log(console,"BUILDING REPORT")  # DEBUG
        if len(invalid_msgs) == 0:
            report += 'sequences in many set: '+str(seq_total)+"\n"
            report += 'sequences in hit set:  '+str(hit_total)+"\n"
            report += "\n"
            for line in hit_buf:
                report += line
            reportObj = {
                'objects_created':[{'ref':params['workspace_name']+'/'+params['output_filtered_name'], 'description':search_tool_name+'_Search hits'}],
                'text_message':report
                }
        else:
            report += "FAILURE\n\n"+"\n".join(invalid_msgs)+"\n"
            reportObj = {
                'objects_created':[],
                'text_message':report
                }

        reportName = 'blast_report_'+str(uuid.uuid4())
        report_obj_info = ws.save_objects({
#                'id':info[6],
                'workspace':params['workspace_name'],
                'objects':[
                    {
                        'type':'KBaseReport.Report',
                        'data':reportObj,
                        'name':reportName,
                        'meta':{},
                        'hidden':1,
                        'provenance':provenance
                    }
                ]
            })[0]

        self.log(console,"BUILDING RETURN OBJECT")
#        returnVal = { 'output_report_name': reportName,
#                      'output_report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
#                      'output_filtered_ref': params['workspace_name']+'/'+params['output_filtered_name']
#                      }
        returnVal = { 'report_name': reportName,
                      'report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
                      }
        self.log(console,search_tool_name+"_Search DONE")
        #END psiBLAST_msa_start_Search

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method psiBLAST_msa_start_Search return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK", 'message': "", 'version': self.VERSION, 
                     'git_url': self.GIT_URL, 'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
