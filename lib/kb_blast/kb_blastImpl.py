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

    ######## WARNING FOR GEVENT USERS #######
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    #########################################
    #BEGIN_CLASS_HEADER
    workspaceURL = None
    shockURL = None
    handleURL = None

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

    def get_single_end_read_library(self, ws_data, ws_info, forward):
        pass

    def get_feature_set_seqs(self, ws_data, ws_info):
        pass

    def get_genome_feature_seqs(self, ws_data, ws_info):
        pass

    def get_genome_set_feature_seqs(self, ws_data, ws_info):
        pass

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.workspaceURL = config['workspace-url']
        self.shockURL = config['shock-url']
        self.handleURL = config['handle-service-url']
        self.scratch = os.path.abspath(config['scratch'])
        # HACK!! temporary hack for issue where megahit fails on mac because of silent named pipe error
        #self.host_scratch = self.scratch
        self.scratch = os.path.join('/kb','module','local_scratch')
        # end hack
        if not os.path.exists(self.scratch):
            os.makedirs(self.scratch)

        #END_CONSTRUCTOR
        pass


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
                        })
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
        self.scratch = os.path.abspath(config['scratch'])
        # HACK!! temporary hack for issue where megahit fails on mac because of silent named pipe error
        #self.host_scratch = self.scratch
        self.scratch = os.path.join('/kb','module','local_scratch')
        # end hack
        if not os.path.exists(self.scratch):
            os.makedirs(self.scratch)

        #END_CONSTRUCTOR
        pass

    def BLASTn_Search(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN BLASTn_Search
        console = []
        self.log(console,'Running BLASTn_Search with params=')
        self.log(console, "\n"+pformat(params))
        report = ''
#        report = 'Running BLASTn_Search with params='
#        report += "\n"+pformat(params)


        #### do some basic checks
        #
        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
#        if 'input_one_name' not in params and 'input_one_sequence' not in params:
#            raise ValueError('input_one_sequence or input_one_name parameter is required')
        if 'input_one_name' not in params:
            raise ValueError('input_one_name parameter is required')
        if 'input_many_name' not in params:
            raise ValueError('input_many_name parameter is required')
        if 'output_filtered_name' not in params:
            raise ValueError('output_filtered_name parameter is required')


        # Write the input_one_sequence to a SingleEndLibrary object
        #
        if 'input_one_sequence' in params \
                and params['input_one_sequence'] != None \
                and params['input_one_sequence'] != "Optionally enter DNA sequence...":
            input_one_file_name = params['input_one_name']
            one_forward_reads_file_path = os.path.join(self.scratch,input_one_file_name)
            one_forward_reads_file_handle = open(one_forward_reads_file_path, 'w', 0)
            self.log(console, 'writing query reads file: '+str(one_forward_reads_file_path))

#            input_sequence_buf = params['input_one_sequence'].split("\n")
#            one_forward_reads_file_handle.write('>'+params['input_one_name']+"\n")
#            query_line_seen = False
#            for line in input_sequence_buf:
#                if not line.startswith('>'):
#                    one_forward_reads_file_handle.write(line+"\n")
#                else:
#                    if query_line_seen:
#                        break
#                    query_line_seen = True
#            one_forward_reads_file_handle.close();

            fastq_format = False
            input_sequence_buf = params['input_one_sequence']
            if input_sequence_buf.startswith('@'):
                fastq_format = True
                #self.log(console,"INPUT_SEQ BEFORE: '''\n"+input_sequence_buf+"\n'''")  # DEBUG
            input_sequence_buf = re.sub ('&apos;', "'", input_sequence_buf)
            input_sequence_buf = re.sub ('&quot;', '"', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#39;',  "'", input_sequence_buf)
#        input_sequence_buf = re.sub ('&#34;',  '"', input_sequence_buf)
#        input_sequence_buf = re.sub ('&lt;;',  '<', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#60;',  '<', input_sequence_buf)
#        input_sequence_buf = re.sub ('&gt;',   '>', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#62;',  '>', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#36;',  '$', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#37;',  '%', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#47;',  '/', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#63;',  '?', input_sequence_buf)
##        input_sequence_buf = re.sub ('&#92;',  chr(92), input_sequence_buf)  # FIX LATER
#        input_sequence_buf = re.sub ('&#96;',  '`', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#124;', '|', input_sequence_buf)
#        input_sequence_buf = re.sub ('&amp;', '&', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#38;', '&', input_sequence_buf)
#        self.log(console,"INPUT_SEQ AFTER: '''\n"+input_sequence_buf+"\n'''")  # DEBUG

            DNA_pattern = re.compile("^[acgtuACGTU ]+$")
            space_pattern = re.compile("^[ \t]*$")
            split_input_sequence_buf = input_sequence_buf.split("\n")

            # no header rows, just sequence
            if not input_sequence_buf.startswith('>') and not input_sequence_buf.startswith('@'):
                one_forward_reads_file_handle.write('>'+params['input_one_name']+"\n")
                for line in split_input_sequence_buf:
                    if not space_pattern.match(line):
                        line = re.sub (" ","",line)
                        line = re.sub ("\t","",line)
                        if not DNA_pattern.match(line):
                            raise ValueError ("BAD record:\n"+line+"\n")
                            sys.exit(0)
                        one_forward_reads_file_handle.write(line.lower()+"\n")
                one_forward_reads_file_handle.close()

            else:
                # format checks
                for i,line in enumerate(split_input_sequence_buf):
                    if line.startswith('>') or line.startswith('@'):
                        if not DNA_pattern.match(split_input_sequence_buf[i+1]):
                            if fastq_format:
                                bad_record = "\n".join([split_input_sequence_buf[i],
                                                        split_input_sequence_buf[i+1],
                                                        split_input_sequence_buf[i+2],
                                                        split_input_sequence_buf[i+3]])
                            else:
                                bad_record = "\n".join([split_input_sequence_buf[i],
                                                    split_input_sequence_buf[i+1]])
                            raise ValueError ("BAD record:\n"+bad_record+"\n")
                            sys.exit(0)
                        if fastq_format and line.startswith('@'):
                            format_ok = True
                            seq_len = len(split_input_sequence_buf[i+1])
                            if not seq_len > 0:
                                format_ok = False
                            if not split_input_sequence_buf[i+2].startswith('+'):
                                format_ok = False
                            if not seq_len == len(split_input_sequence_buf[i+3]):
                                format_ok = False
                            if not format_ok:
                                bad_record = "\n".join([split_input_sequence_buf[i],
                                                    split_input_sequence_buf[i+1],
                                                    split_input_sequence_buf[i+2],
                                                    split_input_sequence_buf[i+3]])
                                raise ValueError ("BAD record:\n"+bad_record+"\n")
                                sys.exit(0)

                # write that sucker, removing spaces
                #
                #forward_reads_file_handle.write(input_sequence_buf)        input_sequence_buf = re.sub ('&quot;', '"', input_sequence_buf)
                for i,line in enumerate(split_input_sequence_buf):
                    if line.startswith('>'):
                        record_buf = []
                        record_buf.append(line)
                        for j in range(i+1,len(split_input_sequence_buf)):
                            if split_input_sequence_buf[j].startswith('>'):
                                break
                            seq_line = re.sub (" ","",split_input_sequence_buf[j])
                            seq_line = re.sub ("\t","",seq_line)
                            seq_line = seq_line.lower()
                            record_buf.append(seq_line)
                        record = "\n".join(record_buf)+"\n"
                        one_forward_reads_file_handle.write(record)
                        break  # only want first record
                    elif line.startswith('@'):
                        seq_line = re.sub (" ","",split_input_sequence_buf[i+1])
                        seq_line = re.sub ("\t","",seq_line)
                        seq_line = seq_line.lower()
                        qual_line = re.sub (" ","",split_input_sequence_buf[i+3])
                        qual_line = re.sub ("\t","",qual_line)
                        record = "\n".join([line, seq_line, split_input_sequence_buf[i+2], qual_line])+"\n"
                        one_forward_reads_file_handle.write(record)
                        break  # only want first record

                one_forward_reads_file_handle.close()


            # load the method provenance from the context object
            #
            self.log(console,"SETTING PROVENANCE")  # DEBUG
            provenance = [{}]
            if 'provenance' in ctx:
                provenance = ctx['provenance']
            # add additional info to provenance here, in this case the input data object reference
                provenance[0]['input_ws_objects'] = []
                provenance[0]['service'] = 'kb_blast'
                provenance[0]['method'] = 'BLASTn_Search'

                
                # Upload results
                #
                self.log(console,"UPLOADING QUERY OBJECT")  # DEBUG

                sequencing_tech = 'N/A'
                self.upload_SingleEndLibrary_to_shock_and_ws (ctx,
                                                      console,  # DEBUG
                                                      params['workspace_name'],
                                                      params['input_one_name'],
                                                      one_forward_reads_file_path,
                                                      provenance,
                                                      sequencing_tech
                                                      )

            self.log(console, 'done')

        #### Get the input_one object
        ##
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            objects = ws.get_objects([{'ref': params['workspace_name']+'/'+params['input_one_name']}])
            data = objects[0]['data']
            info = objects[0]['info']
            # Object Info Contents
            # absolute ref = info[6] + '/' + info[0] + '/' + info[4]
            # 0 - obj_id objid
            # 1 - obj_name name
            # 2 - type_string type
            # 3 - timestamp save_date
            # 4 - int version
            # 5 - username saved_by
            # 6 - ws_id wsid
            # 7 - ws_name workspace
            # 8 - string chsum
            # 9 - int size 
            # 10 - usermeta meta
            one_type_name = info[2].split('.')[1].split('-')[0]
        except Exception as e:
            raise ValueError('Unable to fetch input_one_name object from workspace: ' + str(e))
        #to get the full stack trace: traceback.format_exc()

        if 'input_one_sequence' in params \
                and params['input_one_sequence'] != None \
                and params['input_one_sequence'] != "Optionally enter DNA sequence..." \
                and one_type_name != 'SingleEndLibrary':

            raise ValueError("ERROR: Mismatched input type: input_one_name should be SingleEndLibrary instead of: "+one_type_name)
            sys.exit (0)


        # Handle overloading (input_one can be Feature, SingleEndLibrary, or FeatureSet)
        #
        if one_type_name == 'SingleEndLibrary':
            try:
                if 'lib' in data:
                    one_forward_reads = data['lib']['file']
                elif 'handle' in data:
                    one_forward_reads = data['handle']
                else:
                    self.log(console,"bad structure for 'one_forward_reads'")
                    raise ValueError("bad structure for 'one_forward_reads'")

                ### NOTE: this section is what could be replaced by the transform services
                one_forward_reads_file_path = os.path.join(self.scratch,one_forward_reads['file_name'])
                one_forward_reads_file_handle = open(one_forward_reads_file_path, 'w', 0)
                self.log(console, 'downloading reads file: '+str(one_forward_reads_file_path))
                headers = {'Authorization': 'OAuth '+ctx['token']}
                r = requests.get(one_forward_reads['url']+'/node/'+one_forward_reads['id']+'?download', stream=True, headers=headers)
                for chunk in r.iter_content(1024):
                    one_forward_reads_file_handle.write(chunk)
                one_forward_reads_file_handle.close();
                self.log(console, 'done')
                ### END NOTE


                # remove carriage returns
                new_file_path = one_forward_reads_file_path+"-CRfree"
                new_file_handle = open(new_file_path, 'w', 0)
                one_forward_reads_file_handle = open(one_forward_reads_file_path, 'r', 0)
                for line in one_forward_reads_file_handle:
                    line = re.sub("\r","",line)
                    new_file_handle.write(line)
                one_forward_reads_file_handle.close();
                new_file_handle.close()
                one_forward_reads_file_path = new_file_path


                # convert FASTQ to FASTA (if necessary)
                new_file_path = one_forward_reads_file_path+".fna"
                new_file_handle = open(new_file_path, 'w', 0)
                one_forward_reads_file_handle = open(one_forward_reads_file_path, 'r', 0)
                header = None
                last_header = None
                last_seq_buf = None
                last_line_was_header = False
                was_fastq = False
                for line in one_forward_reads_file_handle:
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
                one_forward_reads_file_handle.close()
                if was_fastq:
                    one_forward_reads_file_path = new_file_path

            except Exception as e:
                print(traceback.format_exc())
                raise ValueError('Unable to download single-end read library files: ' + str(e))

        elif one_type_name == 'FeatureSet':
            # retrieve sequences for features
            input_one_featureSet = data
            
            genome2Features = {}
            features = input_one_featureSet['elements']
            for fId in features.keys():
                genomeRef = features[fId][0]
                if genomeRef not in genome2Features:
                    genome2Features[genomeRef] = []
                genome2Features[genomeRef].append(fId)

            # export features to FASTA file
            one_forward_reads_file_path = os.path.join(self.scratch, params['input_one_name']+".fasta")
            self.log(console, 'writing fasta file: '+one_forward_reads_file_path)
            records = []
            for genomeRef in genome2Features:
                genome = ws.get_objects([{'ref':genomeRef}])[0]['data']
                these_genomeFeatureIds = genome2Features[genomeRef]
                for feature in genome['features']:
                    if feature['id'] in these_genomeFeatureIds:
                        # BLASTn is nuc-nuc
                        record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description=genomeRef+"."+feature['id'])
                        #record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description=genomeRef+"."+feature['id'])
                        records.append(record)
            SeqIO.write(records, one_forward_reads_file_path, "fasta")

        elif one_type_name == 'Feature':
            # export feature to FASTA file
            feature = data
            one_forward_reads_file_path = os.path.join(self.scratch, params['input_one_name']+".fasta")
            self.log(console, 'writing fasta file: '+one_forward_reads_file_path)
            # BLASTn is nuc-nuc
            record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description='['+feature['genome_id']+']'+' '+feature['function'])
            #record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description='['+feature['genome_id']+']'+' '+feature['function'])
            SeqIO.write([record], one_forward_reads_file_path, "fasta")

        else:
            raise ValueError('Cannot yet handle input_one type of: '+type_name)            


        #### Get the input_many object
        ##
        many_forward_reads_file_compression = None
        sequencing_tech = 'N/A'
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            objects = ws.get_objects([{'ref': params['workspace_name']+'/'+params['input_many_name']}])
            data = objects[0]['data']
            info = objects[0]['info']
            many_type_name = info[2].split('.')[1].split('-')[0]

            if many_type_name == 'SingleEndLibrary':
                many_type_namespace = info[2].split('.')[0]
                if many_type_namespace == 'KBaseAssembly':
                    file_name = data['handle']['file_name']
                elif many_type_namespace == 'KBaseFile':
                    file_name = data['lib']['file']['file_name']
                else:
                    raise ValueError('bad data type namespace: '+many_type_namespace)
                #self.log(console, 'INPUT_MANY_FILENAME: '+file_name)  # DEBUG
                if file_name[-3:] == ".gz":
                    many_forward_reads_file_compression = 'gz'
                if 'sequencing_tech' in data:
                    sequencing_tech = data['sequencing_tech']

        except Exception as e:
            raise ValueError('Unable to fetch input_many_name object from workspace: ' + str(e))
            #to get the full stack trace: traceback.format_exc()

        # Handle overloading (input_many can be SingleEndLibrary, FeatureSet, Genome, or GenomeSet)
        #
        if many_type_name == 'SingleEndLibrary':

            # DEBUG
            #for k in data:
            #    self.log(console,"SingleEndLibrary ["+k+"]: "+str(data[k]))

            try:
                if 'lib' in data:
                    many_forward_reads = data['lib']['file']
                elif 'handle' in data:
                    many_forward_reads = data['handle']
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
            input_many_featureSet = data

            genome2Features = {}
            features = input_many_featureSet['elements']
            for fId in features.keys():
                genomeRef = features[fId][0]
                if genomeRef not in genome2Features:
                    genome2Features[genomeRef] = []
                genome2Features[genomeRef].append(fId)

            # export features to FASTA file
            many_forward_reads_file_path = os.path.join(self.scratch, params['input_many_name']+".fasta")
            self.log(console, 'writing fasta file: '+many_forward_reads_file_path)
            records = []
            feature_written = dict()
            for genomeRef in genome2Features:
                genome = ws.get_objects([{'ref':genomeRef}])[0]['data']
                these_genomeFeatureIds = genome2Features[genomeRef]
                for feature in genome['features']:
                    if feature['id'] in these_genomeFeatureIds:
                        try:
                            f_written = feature_written[feature['id']]
                        except:
                            feature_written[feature['id']] = True
                            #self.log(console,"kbase_id: '"+feature['id']+"'")  # DEBUG
                            # BLASTn is nuc-nuc
                            record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description=genome['id'])
                            #record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description=genome['id'])
                            records.append(record)
            SeqIO.write(records, many_forward_reads_file_path, "fasta")


        # Genome
        #
        elif many_type_name == 'Genome':
            input_many_genome = data
            input_many_genome_ref = str(info[6])+'/'+str(info[0])+'/'+str(info[4])

            # export features to FASTA file
            many_forward_reads_file_path = os.path.join(self.scratch, params['input_many_name']+".fasta")
            self.log(console, 'writing fasta file: '+many_forward_reads_file_path)
            records = []
            feature_written = dict()
            for feature in input_many_genome['features']:
                try:
                    f_written = feature_written[feature['id']]
                except:
                    feature_written[feature['id']] = True
                    #self.log(console,"kbase_id: '"+feature['id']+"'")  # DEBUG
                    # BLASTn is nuc-nuc
                    record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description=input_many_genome['id'])
                    #record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description=input_many_genome['id'])
                    records.append(record)
            SeqIO.write(records, many_forward_reads_file_path, "fasta")


        # GenomeSet
        #
        elif many_type_name == 'GenomeSet':
            input_many_genomeSet = data

            # export features to FASTA file
            many_forward_reads_file_path = os.path.join(self.scratch, params['input_many_name']+".fasta")
            self.log(console, 'writing fasta file: '+many_forward_reads_file_path)

            records = []
            feature_written = dict()
            for genome_name in input_many_genomeSet['elements'].keys():
                if 'ref' in input_many_genomeSet['elements'][genome_name] and \
                         input_many_genomeSet['elements'][genome_name]['ref'] != None:
                    genome = ws.get_objects([{'ref': input_many_genomeSet['elements'][genome_name]['ref']}])[0]['data']
                    for feature in genome['features']:
                        try:
                            f_written = feature_written[feature['id']]
                        except:
                            feature_written[feature['id']] = True
                            #self.log(console,"kbase_id: '"+feature['id']+"'")  # DEBUG
                            # BLASTn is nuc-nuc
                            record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description=genome['id'])
                            #record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description=genome['id'])
                            records.append(record)

                elif 'data' in input_many_genomeSet['elements'][genome_name] and \
                        input_many_genomeSet['elements'][genome_name]['data'] != None:
                    genome = input_many_genomeSet['elements'][genome_name]['data']
                    for feature in genome['features']:
                        try:
                            f_written = feature_written[feature['id']]
                        except:
                            feature_written[feature['id']] = True
                            #self.log(console,"kbase_id: '"+feature['id']+"'")  # DEBUG
                            # BLASTn is nuc-nuc
                            record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description=genome['id'])
                            #record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description=genome['id'])
                            records.append(record)

                else:
                    raise ValueError('genome '+genome_name+' missing')

            SeqIO.write(records, many_forward_reads_file_path, "fasta")
            
        # Missing proper input_many_type
        #
        else:
            raise ValueError('Cannot yet handle input_many type of: '+type_name)            

        # FORMAT DB
        #
        # OLD SYNTAX: formatdb -i $database -o T -p F -> $database.nsq or $database.00.nsq
        # NEW SYNTAX: makeblastdb -in $database -parse_seqids -dbtype prot/nucl -out <basename>
        makeblastdb_cmd = [self.Make_BLAST_DB]

        # check for necessary files
        if not os.path.isfile(self.Make_BLAST_DB):
            raise ValueError("no such file '"+self.Make_BLAST_DB+"'")
        if not os.path.isfile(one_forward_reads_file_path):
            raise ValueError("no such file '"+one_forward_reads_file_path+"'")
        if not os.path.isfile(many_forward_reads_file_path):
            raise ValueError("no such file '"+many_forward_reads_file_path+"'")

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


        ### Construct the BLAST command
        #
        # OLD SYNTAX: $blast -q $q -G $G -E $E -m $m -e $e_value -v $limit -b $limit -K $limit -p blastn -i $fasta_file -d $database -o $out_file
        # NEW SYNTAX: blastn -query <queryfile> -db <basename> -out <out_aln_file> -outfmt 0/7 (8 became 7) -evalue <e_value> -dust no (DNA) -seg no (DNA) -num_threads <num_cores>
        #
        blast_bin = self.BLASTn
        blast_cmd = [blast_bin]

        # check for necessary files
        if not os.path.isfile(blast_bin):
            raise ValueError("no such file '"+blast_bin+"'")
        if not os.path.isfile(one_forward_reads_file_path):
            raise ValueError("no such file '"+one_forward_reads_file_path+"'")
        if not os.path.isfile(many_forward_reads_file_path):
            raise ValueError("no such file '"+many_forward_reads_file_path+"'")

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


        # Parse the BLAST tabular output and store ids to filter many set to make filtered object to save back to KBase
        #
        self.log(console, 'PARSING BLAST ALIGNMENT OUTPUT')
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
            #self.log(console,"AFTER ident_thresh")
            if 'bitscore' in params and float(params['bitscore']) > float(high_bitscore_score[hit_seq_id]):
                continue
            #self.log(console,"AFTER bitscore")
            # need to fix this by reading query len
            if 'overlap_fraction' in params and float(params['overlap_fraction']) > 1.0:
                continue
            #self.log(console,"AFTER overlap_fraction")
            
            hit_total += 1
            hit_seq_ids[hit_seq_id] = True
            self.log(console, "HIT: '"+hit_seq_id+"'")  # DEBUG
        

        self.log(console, 'EXTRACTING HITS FROM INPUT')
        self.log(console, 'MANY_TYPE_NAME: '+many_type_name)  # DEBUG


        # SingleEndLibrary input -> SingleEndLibrary output
        #
        if many_type_name == 'SingleEndLibrary':

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
                output_featureSet['description'] = input_many_featureSet['description'] + " - BLASTn_Search filtered"
            else:
                output_featureSet['description'] = "BLASTn_Search filtered"
            output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()
            if 'element_ordering' in input_many_featureSet and input_many_featureSet['element_ordering'] != None:
                for fId in input_many_featureSet['element_ordering']:
                    try:
                        in_filtered_set = hit_seq_ids[fId]
                        #self.log(console, 'FOUND HIT '+fId)  # DEBUG
                        output_featureSet['element_ordering'].append(fId)
                        output_featureSet['elements'][fId] = input_many_featureSet['elements'][fId]
                    except:
                        pass
            else:
                fId_list = input_many_featureSet['elements'].keys()
                self.log(console,"ADDING FEATURES TO FEATURESET")
                for fId in sorted(fId_list):
                    try:
                        #self.log(console,"checking '"+fId+"'")
                        in_filtered_set = hit_seq_ids[fId]
                        #self.log(console, 'FOUND HIT '+fId)  # DEBUG
                        output_featureSet['element_ordering'].append(fId)
                        output_featureSet['elements'][fId] = input_many_featureSet['elements'][fId]
                    except:
                        pass

        # Parse Genome hits into FeatureSet
        #
        elif many_type_name == 'Genome':
            seq_total = 0

            output_featureSet = dict()
            if 'scientific_name' in input_many_genome and input_many_genome['scientific_name'] != None:
                output_featureSet['description'] = input_many_genome['scientific_name'] + " - BLASTn_Search filtered"
            else:
                output_featureSet['description'] = "BLASTn_Search filtered"
            output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()
            for feature in input_many_genome['features']:
                seq_total += 1
                try:
                    in_filtered_set = hit_seq_ids[feature['id']]
                    #self.log(console, 'FOUND HIT: '+feature['id'])  # DEBUG
                    output_featureSet['element_ordering'].append(feature['id'])
                    output_featureSet['elements'][feature['id']] = [input_many_genome_ref]
                except:
                    pass

        # Parse GenomeSet hits into FeatureSet
        #
        elif many_type_name == 'GenomeSet':
            seq_total = 0

            output_featureSet = dict()
            if 'description' in input_many_genomeSet and input_many_genomeSet['description'] != None:
                output_featureSet['description'] = input_many_genomeSet['description'] + " - BLASTn_Search filtered"
            else:
                output_featureSet['description'] = "BLASTn_Search filtered"
            output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()

            for genome_name in input_many_genomeSet['elements'].keys():
                if 'ref' in input_many_genomeSet['elements'][genome_name] and \
                        input_many_genomeSet['elements'][genome_name]['ref'] != None:
                    genomeRef = input_many_genomeSet['elements'][genome_name]['ref']
                    genome = ws.get_objects([{'ref':genomeRef}])[0]['data']
                    for feature in genome['features']:
                        seq_total += 1
                        try:
                            in_filtered_set = hit_seq_ids[feature['id']]
                            #self.log(console, 'FOUND HIT: '+feature['id'])  # DEBUG
                            output_featureSet['element_ordering'].append(feature['id'])
                            output_featureSet['elements'][feature['id']] = [genomeRef]
                        except:
                            pass

                elif 'data' in input_many_genomeSet['elements'][genome_name] and \
                        input_many_genomeSet['elements'][genome_name]['data'] != None:
#                    genome = input_many_genomeSet['elements'][genome_name]['data']
#                    for feature in genome['features']:
#                        #self.log(console,"kbase_id: '"+feature['id']+"'")  # DEBUG
#                        seq_total += 1
#                        try:
#                            in_filtered_set = hit_seq_ids[feature['id']]
#                            #self.log(console, 'FOUND HIT: '+feature['id'])  # DEBUG
#                            output_featureSet['element_ordering'].append(feature['id'])
                    raise ValueError ("FAILURE: unable to address genome object that is stored within 'data' field of genomeSet object")
#                            output_featureSet['elements'][feature['id']] = [genomeRef_is_inside_data_within_genomeSet_object_and_that_cant_be_addressed]
#                        except:
#                            pass


        # load the method provenance from the context object
        #
        self.log(console,"SETTING PROVENANCE")  # DEBUG
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects'] = []
        if 'input_one_name' in params and params['input_one_name'] != None:
            provenance[0]['input_ws_objects'].append(params['workspace_name']+'/'+params['input_one_name'])
        provenance[0]['input_ws_objects'].append(params['workspace_name']+'/'+params['input_many_name'])
        provenance[0]['service'] = 'kb_blast'
        provenance[0]['method'] = 'BLASTn_Search'


        # Upload results
        #
        self.log(console,"UPLOADING RESULTS")  # DEBUG

        if many_type_name == 'SingleEndLibrary':
            
            # input SingleEndLibrary -> upload SingleEndLibrary
            #
            self.upload_SingleEndLibrary_to_shock_and_ws (ctx,
                                                          console,  # DEBUG
                                                          params['workspace_name'],
                                                          params['output_filtered_name'],
                                                          output_filtered_fasta_file_path,
                                                          provenance,
                                                          sequencing_tech
                                                         )

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
                        })

        # build output report object
        #
        self.log(console,"BUILDING REPORT")  # DEBUG
        report += 'sequences in many set: '+str(seq_total)+"\n"
        report += 'sequences in hit set:  '+str(hit_total)+"\n"
        report += "\n"
        for line in hit_buf:
            report += line

        reportObj = {
            'objects_created':[{'ref':params['workspace_name']+'/'+params['output_filtered_name'], 'description':'BLASTn_Search hits'}],
            'text_message':report
        }

        reportName = 'blast_report_'+str(hex(uuid.getnode()))
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
        self.log(console,"BLASTn_Search DONE")
        #END BLASTn_Search

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method BLASTn_Search return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]


    def BLASTp_Search(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN BLASTp_Search
        console = []
        self.log(console,'Running BLASTp_Search with params=')
        self.log(console, "\n"+pformat(params))
        report = ''
#        report = 'Running BLASTp_Search with params='
#        report += "\n"+pformat(params)


        #### do some basic checks
        #
        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
#        if 'input_one_name' not in params and 'input_one_sequence' not in params:
#            raise ValueError('input_one_sequence or input_one_name parameter is required')
        if 'input_one_name' not in params:
            raise ValueError('input_one_name parameter is required')
        if 'input_many_name' not in params:
            raise ValueError('input_many_name parameter is required')
        if 'output_filtered_name' not in params:
            raise ValueError('output_filtered_name parameter is required')


        # Write the input_one_sequence to file
        #
        if 'input_one_sequence' in params \
                and params['input_one_sequence'] != None \
                and params['input_one_sequence'] != "Optionally enter PROTEIN sequence...":
            #input_one_file_name = params['input_one_name']
            input_one_name = 'query.faa'
            input_one_file_name = input_one_name
            one_forward_reads_file_path = os.path.join(self.scratch,input_one_file_name)
            one_forward_reads_file_handle = open(one_forward_reads_file_path, 'w', 0)
            self.log(console, 'writing query reads file: '+str(one_forward_reads_file_path))

#            input_sequence_buf = params['input_one_sequence'].split("\n")
#            one_forward_reads_file_handle.write('>'+params['input_one_name']+"\n")
#            query_line_seen = False
#            for line in input_sequence_buf:
#                if not line.startswith('>'):
#                    one_forward_reads_file_handle.write(line+"\n")
#                else:
#                    if query_line_seen:
#                        break
#                    query_line_seen = True
#            one_forward_reads_file_handle.close();

            input_sequence_buf = params['input_one_sequence']
            space_pattern = re.compile("^[ \t]*$")
            split_input_sequence_buf = input_sequence_buf.split("\n")

            # no header rows, just sequence
            if not input_sequence_buf.startswith('>'):
                one_forward_reads_file_handle.write('>'+input_one_name+"\n")
                for line in split_input_sequence_buf:
                    if not space_pattern.match(line):
                        line = re.sub (" ","",line)
                        line = re.sub ("\t","",line)
                        one_forward_reads_file_handle.write(line.upper()+"\n")
                one_forward_reads_file_handle.close()

            else:
                # write that sucker, removing spaces
                #
                #forward_reads_file_handle.write(input_sequence_buf)        input_sequence_buf = re.sub ('&quot;', '"', input_sequence_buf)
                for i,line in enumerate(split_input_sequence_buf):
                    if line.startswith('>'):
                        record_buf = []
                        record_buf.append(line)
                        for j in range(i+1,len(split_input_sequence_buf)):
                            if split_input_sequence_buf[j].startswith('>'):
                                break
                            seq_line = re.sub (" ","",split_input_sequence_buf[j])
                            seq_line = re.sub ("\t","",seq_line)
                            seq_line = seq_line.upper()
                            record_buf.append(seq_line)
                        record = "\n".join(record_buf)+"\n"
                        one_forward_reads_file_handle.write(record)
                        break  # only want first record
                one_forward_reads_file_handle.close()


        #### Get the input_one object
        ##
        elif 'input_one_name' in params and params['input_one_name'] != None:
            try:
                ws = workspaceService(self.workspaceURL, token=ctx['token'])
                objects = ws.get_objects([{'ref': params['workspace_name']+'/'+params['input_one_name']}])
                data = objects[0]['data']
                info = objects[0]['info']
                # Object Info Contents
                # absolute ref = info[6] + '/' + info[0] + '/' + info[4]
                # 0 - obj_id objid
                # 1 - obj_name name
                # 2 - type_string type
                # 3 - timestamp save_date
                # 4 - int version
                # 5 - username saved_by
                # 6 - ws_id wsid
                # 7 - ws_name workspace
                # 8 - string chsum
                # 9 - int size 
                # 10 - usermeta meta
                one_type_name = info[2].split('.')[1].split('-')[0]
            except Exception as e:
                raise ValueError('Unable to fetch input_one_name object from workspace: ' + str(e))
                #to get the full stack trace: traceback.format_exc()


            # Handle overloading (input_one can be Feature, or FeatureSet)
            #
            if one_type_name == 'FeatureSet':
                # retrieve sequences for features
                input_one_featureSet = data
            
                genome2Features = {}
                features = input_one_featureSet['elements']
                for fId in features.keys():
                    genomeRef = features[fId][0]
                    if genomeRef not in genome2Features:
                        genome2Features[genomeRef] = []
                    genome2Features[genomeRef].append(fId)

                # export features to FASTA file
                one_forward_reads_file_path = os.path.join(self.scratch, params['input_one_name']+".fasta")
                self.log(console, 'writing fasta file: '+one_forward_reads_file_path)
                records = []
                for genomeRef in genome2Features:
                    genome = ws.get_objects([{'ref':genomeRef}])[0]['data']
                    these_genomeFeatureIds = genome2Features[genomeRef]
                    for feature in genome['features']:
                        if feature['id'] in these_genomeFeatureIds:
                            # BLASTp is prot-prot
                            #record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description=genomeRef+"."+feature['id'])
                            if feature['type'] != 'CDS':
                                raise ValueError (params['input_one_name']+" feature type must be CDS")
                                sys.exit(0)
                            record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description=genomeRef+"."+feature['id'])
                            records.append(record)
                SeqIO.write(records, one_forward_reads_file_path, "fasta")

            elif one_type_name == 'Feature':
                # export feature to FASTA file
                feature = data
                one_forward_reads_file_path = os.path.join(self.scratch, params['input_one_name']+".fasta")
                self.log(console, 'writing fasta file: '+one_forward_reads_file_path)
                # BLASTp is prot-prot
                #record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description='['+feature['genome_id']+']'+' '+feature['function'])
                if feature['type'] != 'CDS':
                    raise ValueError (params['input_one_name']+" feature type must be CDS")
                    sys.exit(0)
                record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description='['+feature['genome_id']+']'+' '+feature['function'])
                SeqIO.write([record], one_forward_reads_file_path, "fasta")

            else:
                raise ValueError('Cannot yet handle input_one type of: '+type_name)            
        else:
            raise ValueError('Must define either input_one_sequence or input_one_name')
            sys.exit (0)


        #### Get the input_many object
        ##
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            objects = ws.get_objects([{'ref': params['workspace_name']+'/'+params['input_many_name']}])
            data = objects[0]['data']
            info = objects[0]['info']
            many_type_name = info[2].split('.')[1].split('-')[0]

        except Exception as e:
            raise ValueError('Unable to fetch input_many_name object from workspace: ' + str(e))
            #to get the full stack trace: traceback.format_exc()

        # Handle overloading (input_many can be FeatureSet, Genome, or GenomeSet)
        #
        if many_type_name == 'FeatureSet':
            # retrieve sequences for features
            input_many_featureSet = data

            genome2Features = {}
            features = input_many_featureSet['elements']
            for fId in features.keys():
                genomeRef = features[fId][0]
                if genomeRef not in genome2Features:
                    genome2Features[genomeRef] = []
                genome2Features[genomeRef].append(fId)

            # export features to FASTA file
            many_forward_reads_file_path = os.path.join(self.scratch, params['input_many_name']+".fasta")
            self.log(console, 'writing fasta file: '+many_forward_reads_file_path)
            records = []
            feature_written = dict()
            for genomeRef in genome2Features:
                genome = ws.get_objects([{'ref':genomeRef}])[0]['data']
                these_genomeFeatureIds = genome2Features[genomeRef]
                for feature in genome['features']:
                    if feature['id'] in these_genomeFeatureIds:
                        if feature['type'] != 'CDS':
                            raise ValueError (params['input_many_name']+" feature types must all be CDS")
                            sys.exit(0)
                        try:
                            f_written = feature_written[feature['id']]
                        except:
                            feature_written[feature['id']] = True
                            #self.log(console,"kbase_id: '"+feature['id']+"'")  # DEBUG
                            # BLASTp is prot-prot
                            #record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description=genome['id'])
                            record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description=genome['id'])
                            records.append(record)
            SeqIO.write(records, many_forward_reads_file_path, "fasta")


        # Genome
        #
        elif many_type_name == 'Genome':
            input_many_genome = data
            input_many_genome_ref = str(info[6])+'/'+str(info[0])+'/'+str(info[4])

            # export features to FASTA file
            many_forward_reads_file_path = os.path.join(self.scratch, params['input_many_name']+".fasta")
            self.log(console, 'writing fasta file: '+many_forward_reads_file_path)
            records = []
            feature_written = dict()
            for feature in input_many_genome['features']:
                try:
                    f_written = feature_written[feature['id']]
                except:
                    feature_written[feature['id']] = True
                    #self.log(console,"kbase_id: '"+feature['id']+"'")  # DEBUG
                    # BLASTp is prot-prot
                    #record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description=input_many_genome['id'])
                    if feature['type'] != 'CDS':
                        continue
                    record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description=input_many_genome['id'])
                    records.append(record)
            SeqIO.write(records, many_forward_reads_file_path, "fasta")


        # GenomeSet
        #
        elif many_type_name == 'GenomeSet':
            input_many_genomeSet = data

            # export features to FASTA file
            many_forward_reads_file_path = os.path.join(self.scratch, params['input_many_name']+".fasta")
            self.log(console, 'writing fasta file: '+many_forward_reads_file_path)

            records = []
            feature_written = dict()
            for genome_name in input_many_genomeSet['elements'].keys():
                if 'ref' in input_many_genomeSet['elements'][genome_name] and \
                         input_many_genomeSet['elements'][genome_name]['ref'] != None:
                    genome = ws.get_objects([{'ref': input_many_genomeSet['elements'][genome_name]['ref']}])[0]['data']
                    for feature in genome['features']:
                        try:
                            f_written = feature_written[feature['id']]
                        except:
                            feature_written[feature['id']] = True
                            #self.log(console,"kbase_id: '"+feature['id']+"'")  # DEBUG
                            # BLASTp is prot-prot
                            #record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description=genome['id'])
                            if feature['type'] != 'CDS':
                                continue
                            record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description=genome['id'])
                            records.append(record)

                elif 'data' in input_many_genomeSet['elements'][genome_name] and \
                        input_many_genomeSet['elements'][genome_name]['data'] != None:
                    genome = input_many_genomeSet['elements'][genome_name]['data']
                    for feature in genome['features']:
                        try:
                            f_written = feature_written[feature['id']]
                        except:
                            feature_written[feature['id']] = True
                            #self.log(console,"kbase_id: '"+feature['id']+"'")  # DEBUG
                            # BLASTp is prot-prot
                            #record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description=genome['id'])
                            if feature['type'] != 'CDS':
                                continue
                            record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description=genome['id'])
                            records.append(record)

                else:
                    raise ValueError('genome '+genome_name+' missing')

            SeqIO.write(records, many_forward_reads_file_path, "fasta")
            
        # Missing proper input_many_type
        #
        else:
            raise ValueError('Cannot yet handle input_many type of: '+type_name)            

        # FORMAT DB
        #
        # OLD SYNTAX: formatdb -i $database -o T -p F -> $database.nsq or $database.00.nsq
        # NEW SYNTAX: makeblastdb -in $database -parse_seqids -dbtype prot/nucl -out <basename>
        makeblastdb_cmd = [self.Make_BLAST_DB]

        # check for necessary files
        if not os.path.isfile(self.Make_BLAST_DB):
            raise ValueError("no such file '"+self.Make_BLAST_DB+"'")
        if not os.path.isfile(one_forward_reads_file_path):
            raise ValueError("no such file '"+one_forward_reads_file_path+"'")
        if not os.path.isfile(many_forward_reads_file_path):
            raise ValueError("no such file '"+many_forward_reads_file_path+"'")

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


        ### Construct the BLAST command
        #
        # OLD SYNTAX: $blast -q $q -G $G -E $E -m $m -e $e_value -v $limit -b $limit -K $limit -p blastp -i $fasta_file -d $database -o $out_file
        # NEW SYNTAX: blastp -query <queryfile> -db <basename> -out <out_aln_file> -outfmt 0/7 (8 became 7) -evalue <e_value> -dust no (DNA) -seg no (DNA) -num_threads <num_cores>
        #
        blast_bin = self.BLASTp
        blast_cmd = [blast_bin]

        # check for necessary files
        if not os.path.isfile(blast_bin):
            raise ValueError("no such file '"+blast_bin+"'")
        if not os.path.isfile(one_forward_reads_file_path):
            raise ValueError("no such file '"+one_forward_reads_file_path+"'")
        if not os.path.isfile(many_forward_reads_file_path):
            raise ValueError("no such file '"+many_forward_reads_file_path+"'")

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


        # Parse the BLAST tabular output and store ids to filter many set to make filtered object to save back to KBase
        #
        self.log(console, 'PARSING BLAST ALIGNMENT OUTPUT')
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
            #self.log(console,"AFTER ident_thresh")
            if 'bitscore' in params and float(params['bitscore']) > float(high_bitscore_score[hit_seq_id]):
                continue
            #self.log(console,"AFTER bitscore")
            # need to fix this by reading query len
            if 'overlap_fraction' in params and float(params['overlap_fraction']) > 1.0:
                continue
            #self.log(console,"AFTER overlap_fraction")
            
            hit_total += 1
            hit_seq_ids[hit_seq_id] = True
            self.log(console, "HIT: '"+hit_seq_id+"'")  # DEBUG
        

        self.log(console, 'EXTRACTING HITS FROM INPUT')
        self.log(console, 'MANY_TYPE_NAME: '+many_type_name)  # DEBUG


        # FeatureSet input -> FeatureSet output
        #
        if many_type_name == 'FeatureSet':

            seq_total = len(input_many_featureSet['elements'].keys())

            output_featureSet = dict()
            if 'description' in input_many_featureSet and input_many_featureSet['description'] != None:
                output_featureSet['description'] = input_many_featureSet['description'] + " - BLASTp_Search filtered"
            else:
                output_featureSet['description'] = "BLASTp_Search filtered"
            output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()
            if 'element_ordering' in input_many_featureSet and input_many_featureSet['element_ordering'] != None:
                for fId in input_many_featureSet['element_ordering']:
                    try:
                        in_filtered_set = hit_seq_ids[fId]
                        #self.log(console, 'FOUND HIT '+fId)  # DEBUG
                        output_featureSet['element_ordering'].append(fId)
                        output_featureSet['elements'][fId] = input_many_featureSet['elements'][fId]
                    except:
                        pass
            else:
                fId_list = input_many_featureSet['elements'].keys()
                self.log(console,"ADDING FEATURES TO FEATURESET")
                for fId in sorted(fId_list):
                    try:
                        #self.log(console,"checking '"+fId+"'")
                        in_filtered_set = hit_seq_ids[fId]
                        #self.log(console, 'FOUND HIT '+fId)  # DEBUG
                        output_featureSet['element_ordering'].append(fId)
                        output_featureSet['elements'][fId] = input_many_featureSet['elements'][fId]
                    except:
                        pass

        # Parse Genome hits into FeatureSet
        #
        elif many_type_name == 'Genome':
            seq_total = 0

            output_featureSet = dict()
            if 'scientific_name' in input_many_genome and input_many_genome['scientific_name'] != None:
                output_featureSet['description'] = input_many_genome['scientific_name'] + " - BLASTp_Search filtered"
            else:
                output_featureSet['description'] = "BLASTp_Search filtered"
            output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()
            for feature in input_many_genome['features']:
                seq_total += 1
                try:
                    in_filtered_set = hit_seq_ids[feature['id']]
                    #self.log(console, 'FOUND HIT: '+feature['id'])  # DEBUG
                    output_featureSet['element_ordering'].append(feature['id'])
                    output_featureSet['elements'][feature['id']] = [input_many_genome_ref]
                except:
                    pass

        # Parse GenomeSet hits into FeatureSet
        #
        elif many_type_name == 'GenomeSet':
            seq_total = 0

            output_featureSet = dict()
            if 'description' in input_many_genomeSet and input_many_genomeSet['description'] != None:
                output_featureSet['description'] = input_many_genomeSet['description'] + " - BLASTp_Search filtered"
            else:
                output_featureSet['description'] = "BLASTp_Search filtered"
            output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()

            for genome_name in input_many_genomeSet['elements'].keys():
                if 'ref' in input_many_genomeSet['elements'][genome_name] and \
                        input_many_genomeSet['elements'][genome_name]['ref'] != None:
                    genomeRef = input_many_genomeSet['elements'][genome_name]['ref']
                    genome = ws.get_objects([{'ref':genomeRef}])[0]['data']
                    for feature in genome['features']:
                        seq_total += 1
                        try:
                            in_filtered_set = hit_seq_ids[feature['id']]
                            #self.log(console, 'FOUND HIT: '+feature['id'])  # DEBUG
                            output_featureSet['element_ordering'].append(feature['id'])
                            output_featureSet['elements'][feature['id']] = [genomeRef]
                        except:
                            pass

                elif 'data' in input_many_genomeSet['elements'][genome_name] and \
                        input_many_genomeSet['elements'][genome_name]['data'] != None:
#                    genome = input_many_genomeSet['elements'][genome_name]['data']
#                    for feature in genome['features']:
#                        #self.log(console,"kbase_id: '"+feature['id']+"'")  # DEBUG
#                        seq_total += 1
#                        try:
#                            in_filtered_set = hit_seq_ids[feature['id']]
#                            #self.log(console, 'FOUND HIT: '+feature['id'])  # DEBUG
#                            output_featureSet['element_ordering'].append(feature['id'])
                    raise ValueError ("FAILURE: unable to address genome object that is stored within 'data' field of genomeSet object")
#                            output_featureSet['elements'][feature['id']] = [genomeRef_is_inside_data_within_genomeSet_object_and_that_cant_be_addressed]
#                        except:
#                            pass


        # load the method provenance from the context object
        #
        self.log(console,"SETTING PROVENANCE")  # DEBUG
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects'] = []
        if 'input_one_name' in params and params['input_one_name'] != None:
            provenance[0]['input_ws_objects'].append(params['workspace_name']+'/'+params['input_one_name'])
        provenance[0]['input_ws_objects'].append(params['workspace_name']+'/'+params['input_many_name'])
        provenance[0]['service'] = 'kb_blast'
        provenance[0]['method'] = 'BLASTp_Search'


        # Upload results
        #
        self.log(console,"UPLOADING RESULTS")  # DEBUG

        # input FeatureSet, Genome, and GenomeSet -> upload FeatureSet output
        new_obj_info = ws.save_objects({
                            'workspace': params['workspace_name'],
                            'objects':[{
                                    'type': 'KBaseCollections.FeatureSet',
                                    'data': output_featureSet,
                                    'name': params['output_filtered_name'],
                                    'meta': {},
                                    'provenance': provenance
                                }]
                        })

        # build output report object
        #
        self.log(console,"BUILDING REPORT")  # DEBUG
        report += 'sequences in many set: '+str(seq_total)+"\n"
        report += 'sequences in hit set:  '+str(hit_total)+"\n"
        report += "\n"
        for line in hit_buf:
            report += line

        reportObj = {
            'objects_created':[{'ref':params['workspace_name']+'/'+params['output_filtered_name'], 'description':'BLASTp_Search hits'}],
            'text_message':report
        }

        reportName = 'blast_report_'+str(hex(uuid.getnode()))
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
        self.log(console,"BLASTp_Search DONE")
        #END BLASTp_Search

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method BLASTp_Search return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]


    def BLASTx_Search(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN BLASTx_Search
        console = []
        self.log(console,'Running BLASTx_Search with params=')
        self.log(console, "\n"+pformat(params))
        report = ''
#        report = 'Running BLASTx_Search with params='
#        report += "\n"+pformat(params)


        #### do some basic checks
        #
        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
#        if 'input_one_name' not in params and 'input_one_sequence' not in params:
#            raise ValueError('input_one_sequence or input_one_name parameter is required')
        if 'input_one_name' not in params:
            raise ValueError('input_one_name parameter is required')
        if 'input_many_name' not in params:
            raise ValueError('input_many_name parameter is required')
        if 'output_filtered_name' not in params:
            raise ValueError('output_filtered_name parameter is required')


        # Write the input_one_sequence to a SingleEndLibrary object
        #
        if 'input_one_sequence' in params \
                and params['input_one_sequence'] != None \
                and params['input_one_sequence'] != "Optionally enter DNA sequence...":
            input_one_file_name = params['input_one_name']
            one_forward_reads_file_path = os.path.join(self.scratch,input_one_file_name)
            one_forward_reads_file_handle = open(one_forward_reads_file_path, 'w', 0)
            self.log(console, 'writing query reads file: '+str(one_forward_reads_file_path))

#            input_sequence_buf = params['input_one_sequence'].split("\n")
#            one_forward_reads_file_handle.write('>'+params['input_one_name']+"\n")
#            query_line_seen = False
#            for line in input_sequence_buf:
#                if not line.startswith('>'):
#                    one_forward_reads_file_handle.write(line+"\n")
#                else:
#                    if query_line_seen:
#                        break
#                    query_line_seen = True
#            one_forward_reads_file_handle.close();

            fastq_format = False
            input_sequence_buf = params['input_one_sequence']
            if input_sequence_buf.startswith('@'):
                fastq_format = True
                #self.log(console,"INPUT_SEQ BEFORE: '''\n"+input_sequence_buf+"\n'''")  # DEBUG
            input_sequence_buf = re.sub ('&apos;', "'", input_sequence_buf)
            input_sequence_buf = re.sub ('&quot;', '"', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#39;',  "'", input_sequence_buf)
#        input_sequence_buf = re.sub ('&#34;',  '"', input_sequence_buf)
#        input_sequence_buf = re.sub ('&lt;;',  '<', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#60;',  '<', input_sequence_buf)
#        input_sequence_buf = re.sub ('&gt;',   '>', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#62;',  '>', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#36;',  '$', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#37;',  '%', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#47;',  '/', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#63;',  '?', input_sequence_buf)
##        input_sequence_buf = re.sub ('&#92;',  chr(92), input_sequence_buf)  # FIX LATER
#        input_sequence_buf = re.sub ('&#96;',  '`', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#124;', '|', input_sequence_buf)
#        input_sequence_buf = re.sub ('&amp;', '&', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#38;', '&', input_sequence_buf)
#        self.log(console,"INPUT_SEQ AFTER: '''\n"+input_sequence_buf+"\n'''")  # DEBUG

            DNA_pattern = re.compile("^[acgtuACGTU ]+$")
            space_pattern = re.compile("^[ \t]*$")
            split_input_sequence_buf = input_sequence_buf.split("\n")

            # no header rows, just sequence
            if not input_sequence_buf.startswith('>') and not input_sequence_buf.startswith('@'):
                one_forward_reads_file_handle.write('>'+params['input_one_name']+"\n")
                for line in split_input_sequence_buf:
                    if not space_pattern.match(line):
                        line = re.sub (" ","",line)
                        line = re.sub ("\t","",line)
                        if not DNA_pattern.match(line):
                            raise ValueError ("BAD record:\n"+line+"\n")
                            sys.exit(0)
                        one_forward_reads_file_handle.write(line.lower()+"\n")
                one_forward_reads_file_handle.close()

            else:
                # format checks
                for i,line in enumerate(split_input_sequence_buf):
                    if line.startswith('>') or line.startswith('@'):
                        if not DNA_pattern.match(split_input_sequence_buf[i+1]):
                            if fastq_format:
                                bad_record = "\n".join([split_input_sequence_buf[i],
                                                        split_input_sequence_buf[i+1],
                                                        split_input_sequence_buf[i+2],
                                                        split_input_sequence_buf[i+3]])
                            else:
                                bad_record = "\n".join([split_input_sequence_buf[i],
                                                    split_input_sequence_buf[i+1]])
                            raise ValueError ("BAD record:\n"+bad_record+"\n")
                            sys.exit(0)
                        if fastq_format and line.startswith('@'):
                            format_ok = True
                            seq_len = len(split_input_sequence_buf[i+1])
                            if not seq_len > 0:
                                format_ok = False
                            if not split_input_sequence_buf[i+2].startswith('+'):
                                format_ok = False
                            if not seq_len == len(split_input_sequence_buf[i+3]):
                                format_ok = False
                            if not format_ok:
                                bad_record = "\n".join([split_input_sequence_buf[i],
                                                    split_input_sequence_buf[i+1],
                                                    split_input_sequence_buf[i+2],
                                                    split_input_sequence_buf[i+3]])
                                raise ValueError ("BAD record:\n"+bad_record+"\n")
                                sys.exit(0)

                # write that sucker, removing spaces
                #
                #forward_reads_file_handle.write(input_sequence_buf)        input_sequence_buf = re.sub ('&quot;', '"', input_sequence_buf)
                for i,line in enumerate(split_input_sequence_buf):
                    if line.startswith('>'):
                        record_buf = []
                        record_buf.append(line)
                        for j in range(i+1,len(split_input_sequence_buf)):
                            if split_input_sequence_buf[j].startswith('>'):
                                break
                            seq_line = re.sub (" ","",split_input_sequence_buf[j])
                            seq_line = re.sub ("\t","",seq_line)
                            seq_line = seq_line.lower()
                            record_buf.append(seq_line)
                        record = "\n".join(record_buf)+"\n"
                        one_forward_reads_file_handle.write(record)
                        break  # only want first record
                    elif line.startswith('@'):
                        seq_line = re.sub (" ","",split_input_sequence_buf[i+1])
                        seq_line = re.sub ("\t","",seq_line)
                        seq_line = seq_line.lower()
                        qual_line = re.sub (" ","",split_input_sequence_buf[i+3])
                        qual_line = re.sub ("\t","",qual_line)
                        record = "\n".join([line, seq_line, split_input_sequence_buf[i+2], qual_line])+"\n"
                        one_forward_reads_file_handle.write(record)
                        break  # only want first record

                one_forward_reads_file_handle.close()


            # load the method provenance from the context object
            #
            self.log(console,"SETTING PROVENANCE")  # DEBUG
            provenance = [{}]
            if 'provenance' in ctx:
                provenance = ctx['provenance']
            # add additional info to provenance here, in this case the input data object reference
                provenance[0]['input_ws_objects'] = []
                provenance[0]['service'] = 'kb_blast'
                provenance[0]['method'] = 'BLASTx_Search'

                
                # Upload results
                #
                self.log(console,"UPLOADING QUERY OBJECT")  # DEBUG

                sequencing_tech = 'N/A'
                self.upload_SingleEndLibrary_to_shock_and_ws (ctx,
                                                      console,  # DEBUG
                                                      params['workspace_name'],
                                                      params['input_one_name'],
                                                      one_forward_reads_file_path,
                                                      provenance,
                                                      sequencing_tech
                                                      )

            self.log(console, 'done')

        #### Get the input_one object
        ##
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            objects = ws.get_objects([{'ref': params['workspace_name']+'/'+params['input_one_name']}])
            data = objects[0]['data']
            info = objects[0]['info']
            # Object Info Contents
            # absolute ref = info[6] + '/' + info[0] + '/' + info[4]
            # 0 - obj_id objid
            # 1 - obj_name name
            # 2 - type_string type
            # 3 - timestamp save_date
            # 4 - int version
            # 5 - username saved_by
            # 6 - ws_id wsid
            # 7 - ws_name workspace
            # 8 - string chsum
            # 9 - int size 
            # 10 - usermeta meta
            one_type_name = info[2].split('.')[1].split('-')[0]
        except Exception as e:
            raise ValueError('Unable to fetch input_one_name object from workspace: ' + str(e))
        #to get the full stack trace: traceback.format_exc()

        if 'input_one_sequence' in params \
                and params['input_one_sequence'] != None \
                and params['input_one_sequence'] != "Optionally enter DNA sequence..." \
                and one_type_name != 'SingleEndLibrary':

            raise ValueError("ERROR: Mismatched input type: input_one_name should be SingleEndLibrary instead of: "+one_type_name)
            sys.exit (0)


        # Handle overloading (input_one can be Feature, SingleEndLibrary, or FeatureSet)
        #
        if one_type_name == 'SingleEndLibrary':
            try:
                if 'lib' in data:
                    one_forward_reads = data['lib']['file']
                elif 'handle' in data:
                    one_forward_reads = data['handle']
                else:
                    self.log(console,"bad structure for 'one_forward_reads'")
                    raise ValueError("bad structure for 'one_forward_reads'")

                ### NOTE: this section is what could be replaced by the transform services
                one_forward_reads_file_path = os.path.join(self.scratch,one_forward_reads['file_name'])
                one_forward_reads_file_handle = open(one_forward_reads_file_path, 'w', 0)
                self.log(console, 'downloading reads file: '+str(one_forward_reads_file_path))
                headers = {'Authorization': 'OAuth '+ctx['token']}
                r = requests.get(one_forward_reads['url']+'/node/'+one_forward_reads['id']+'?download', stream=True, headers=headers)
                for chunk in r.iter_content(1024):
                    one_forward_reads_file_handle.write(chunk)
                one_forward_reads_file_handle.close();
                self.log(console, 'done')
                ### END NOTE


                # remove carriage returns
                new_file_path = one_forward_reads_file_path+"-CRfree"
                new_file_handle = open(new_file_path, 'w', 0)
                one_forward_reads_file_handle = open(one_forward_reads_file_path, 'r', 0)
                for line in one_forward_reads_file_handle:
                    line = re.sub("\r","",line)
                    new_file_handle.write(line)
                one_forward_reads_file_handle.close();
                new_file_handle.close()
                one_forward_reads_file_path = new_file_path


                # convert FASTQ to FASTA (if necessary)
                new_file_path = one_forward_reads_file_path+".fna"
                new_file_handle = open(new_file_path, 'w', 0)
                one_forward_reads_file_handle = open(one_forward_reads_file_path, 'r', 0)
                header = None
                last_header = None
                last_seq_buf = None
                last_line_was_header = False
                was_fastq = False
                for line in one_forward_reads_file_handle:
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
                one_forward_reads_file_handle.close()
                if was_fastq:
                    one_forward_reads_file_path = new_file_path

            except Exception as e:
                print(traceback.format_exc())
                raise ValueError('Unable to download single-end read library files: ' + str(e))

        elif one_type_name == 'FeatureSet':
            # retrieve sequences for features
            input_one_featureSet = data
            
            genome2Features = {}
            features = input_one_featureSet['elements']
            for fId in features.keys():
                genomeRef = features[fId][0]
                if genomeRef not in genome2Features:
                    genome2Features[genomeRef] = []
                genome2Features[genomeRef].append(fId)

            # export features to FASTA file
            one_forward_reads_file_path = os.path.join(self.scratch, params['input_one_name']+".fasta")
            self.log(console, 'writing fasta file: '+one_forward_reads_file_path)
            records = []
            for genomeRef in genome2Features:
                genome = ws.get_objects([{'ref':genomeRef}])[0]['data']
                these_genomeFeatureIds = genome2Features[genomeRef]
                for feature in genome['features']:
                    if feature['id'] in these_genomeFeatureIds:
                        # BLASTx is nuc-prot
                        if feature['type'] != 'CDS':
                            raise ValueError (params['input_one_name']+" feature type must be CDS")
                            sys.exit(0)
                        record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description=genomeRef+"."+feature['id'])
                        #record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description=genomeRef+"."+feature['id'])
                        records.append(record)
            SeqIO.write(records, one_forward_reads_file_path, "fasta")

        elif one_type_name == 'Feature':
            # export feature to FASTA file
            feature = data
            one_forward_reads_file_path = os.path.join(self.scratch, params['input_one_name']+".fasta")
            self.log(console, 'writing fasta file: '+one_forward_reads_file_path)
            # BLASTx is nuc-prot
            if feature['type'] != 'CDS':
                raise ValueError (params['input_one_name']+" feature type must be CDS")
                sys.exit(0)
            record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description='['+feature['genome_id']+']'+' '+feature['function'])
            #record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description='['+feature['genome_id']+']'+' '+feature['function'])
            SeqIO.write([record], one_forward_reads_file_path, "fasta")

        else:
            raise ValueError('Cannot yet handle input_one type of: '+type_name)            

        #### Get the input_many object
        ##
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            objects = ws.get_objects([{'ref': params['workspace_name']+'/'+params['input_many_name']}])
            data = objects[0]['data']
            info = objects[0]['info']
            many_type_name = info[2].split('.')[1].split('-')[0]

        except Exception as e:
            raise ValueError('Unable to fetch input_many_name object from workspace: ' + str(e))
            #to get the full stack trace: traceback.format_exc()

        # Handle overloading (input_many can be FeatureSet, Genome, or GenomeSet)
        #
        if many_type_name == 'FeatureSet':
            # retrieve sequences for features
            input_many_featureSet = data

            genome2Features = {}
            features = input_many_featureSet['elements']
            for fId in features.keys():
                genomeRef = features[fId][0]
                if genomeRef not in genome2Features:
                    genome2Features[genomeRef] = []
                genome2Features[genomeRef].append(fId)

            # export features to FASTA file
            many_forward_reads_file_path = os.path.join(self.scratch, params['input_many_name']+".fasta")
            self.log(console, 'writing fasta file: '+many_forward_reads_file_path)
            records = []
            feature_written = dict()
            for genomeRef in genome2Features:
                genome = ws.get_objects([{'ref':genomeRef}])[0]['data']
                these_genomeFeatureIds = genome2Features[genomeRef]
                for feature in genome['features']:
                    if feature['id'] in these_genomeFeatureIds:
                        if feature['type'] != 'CDS':
                            raise ValueError (params['input_many_name']+" feature types must all be CDS")
                            sys.exit(0)
                        try:
                            f_written = feature_written[feature['id']]
                        except:
                            feature_written[feature['id']] = True
                            #self.log(console,"kbase_id: '"+feature['id']+"'")  # DEBUG
                            # BLASTx is nuc-prot
                            #record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description=genome['id'])
                            record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description=genome['id'])
                            records.append(record)
            SeqIO.write(records, many_forward_reads_file_path, "fasta")


        # Genome
        #
        elif many_type_name == 'Genome':
            input_many_genome = data
            input_many_genome_ref = str(info[6])+'/'+str(info[0])+'/'+str(info[4])

            # export features to FASTA file
            many_forward_reads_file_path = os.path.join(self.scratch, params['input_many_name']+".fasta")
            self.log(console, 'writing fasta file: '+many_forward_reads_file_path)
            records = []
            feature_written = dict()
            for feature in input_many_genome['features']:
                try:
                    f_written = feature_written[feature['id']]
                except:
                    feature_written[feature['id']] = True
                    #self.log(console,"kbase_id: '"+feature['id']+"'")  # DEBUG
                    # BLASTx is nuc-prot
                    #record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description=input_many_genome['id'])
                    if feature['type'] != 'CDS':
                        continue
                    record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description=input_many_genome['id'])
                    records.append(record)
            SeqIO.write(records, many_forward_reads_file_path, "fasta")


        # GenomeSet
        #
        elif many_type_name == 'GenomeSet':
            input_many_genomeSet = data

            # export features to FASTA file
            many_forward_reads_file_path = os.path.join(self.scratch, params['input_many_name']+".fasta")
            self.log(console, 'writing fasta file: '+many_forward_reads_file_path)

            records = []
            feature_written = dict()
            for genome_name in input_many_genomeSet['elements'].keys():
                if 'ref' in input_many_genomeSet['elements'][genome_name] and \
                         input_many_genomeSet['elements'][genome_name]['ref'] != None:
                    genome = ws.get_objects([{'ref': input_many_genomeSet['elements'][genome_name]['ref']}])[0]['data']
                    for feature in genome['features']:
                        try:
                            f_written = feature_written[feature['id']]
                        except:
                            feature_written[feature['id']] = True
                            #self.log(console,"kbase_id: '"+feature['id']+"'")  # DEBUG
                            # BLASTx is nuc-prot
                            #record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description=genome['id'])
                            if feature['type'] != 'CDS':
                                continue
                            record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description=genome['id'])
                            records.append(record)

                elif 'data' in input_many_genomeSet['elements'][genome_name] and \
                        input_many_genomeSet['elements'][genome_name]['data'] != None:
                    genome = input_many_genomeSet['elements'][genome_name]['data']
                    for feature in genome['features']:
                        try:
                            f_written = feature_written[feature['id']]
                        except:
                            feature_written[feature['id']] = True
                            #self.log(console,"kbase_id: '"+feature['id']+"'")  # DEBUG
                            # BLASTx is nuc-prot
                            #record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description=genome['id'])
                            if feature['type'] != 'CDS':
                                continue
                            record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description=genome['id'])
                            records.append(record)

                else:
                    raise ValueError('genome '+genome_name+' missing')

            SeqIO.write(records, many_forward_reads_file_path, "fasta")
            
        # Missing proper input_many_type
        #
        else:
            raise ValueError('Cannot yet handle input_many type of: '+type_name)            

        # FORMAT DB
        #
        # OLD SYNTAX: formatdb -i $database -o T -p F -> $database.nsq or $database.00.nsq
        # NEW SYNTAX: makeblastdb -in $database -parse_seqids -dbtype prot/nucl -out <basename>
        makeblastdb_cmd = [self.Make_BLAST_DB]

        # check for necessary files
        if not os.path.isfile(self.Make_BLAST_DB):
            raise ValueError("no such file '"+self.Make_BLAST_DB+"'")
        if not os.path.isfile(one_forward_reads_file_path):
            raise ValueError("no such file '"+one_forward_reads_file_path+"'")
        if not os.path.isfile(many_forward_reads_file_path):
            raise ValueError("no such file '"+many_forward_reads_file_path+"'")

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


        ### Construct the BLAST command
        #
        # OLD SYNTAX: $blast -q $q -G $G -E $E -m $m -e $e_value -v $limit -b $limit -K $limit -p blastx -i $fasta_file -d $database -o $out_file
        # NEW SYNTAX: blastx -query <queryfile> -db <basename> -out <out_aln_file> -outfmt 0/7 (8 became 7) -evalue <e_value> -dust no (DNA) -seg no (DNA) -num_threads <num_cores>
        #
        blast_bin = self.BLASTx
        blast_cmd = [blast_bin]

        # check for necessary files
        if not os.path.isfile(blast_bin):
            raise ValueError("no such file '"+blast_bin+"'")
        if not os.path.isfile(one_forward_reads_file_path):
            raise ValueError("no such file '"+one_forward_reads_file_path+"'")
        if not os.path.isfile(many_forward_reads_file_path):
            raise ValueError("no such file '"+many_forward_reads_file_path+"'")

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


        # Parse the BLAST tabular output and store ids to filter many set to make filtered object to save back to KBase
        #
        self.log(console, 'PARSING BLAST ALIGNMENT OUTPUT')
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
            #self.log(console,"AFTER ident_thresh")
            if 'bitscore' in params and float(params['bitscore']) > float(high_bitscore_score[hit_seq_id]):
                continue
            #self.log(console,"AFTER bitscore")
            # need to fix this by reading query len
            if 'overlap_fraction' in params and float(params['overlap_fraction']) > 1.0:
                continue
            #self.log(console,"AFTER overlap_fraction")
            
            hit_total += 1
            hit_seq_ids[hit_seq_id] = True
            self.log(console, "HIT: '"+hit_seq_id+"'")  # DEBUG
        

        self.log(console, 'EXTRACTING HITS FROM INPUT')
        self.log(console, 'MANY_TYPE_NAME: '+many_type_name)  # DEBUG


        # FeatureSet input -> FeatureSet output
        #
        if many_type_name == 'FeatureSet':

            seq_total = len(input_many_featureSet['elements'].keys())

            output_featureSet = dict()
            if 'description' in input_many_featureSet and input_many_featureSet['description'] != None:
                output_featureSet['description'] = input_many_featureSet['description'] + " - BLASTx_Search filtered"
            else:
                output_featureSet['description'] = "BLASTx_Search filtered"
            output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()
            if 'element_ordering' in input_many_featureSet and input_many_featureSet['element_ordering'] != None:
                for fId in input_many_featureSet['element_ordering']:
                    try:
                        in_filtered_set = hit_seq_ids[fId]
                        #self.log(console, 'FOUND HIT '+fId)  # DEBUG
                        output_featureSet['element_ordering'].append(fId)
                        output_featureSet['elements'][fId] = input_many_featureSet['elements'][fId]
                    except:
                        pass
            else:
                fId_list = input_many_featureSet['elements'].keys()
                self.log(console,"ADDING FEATURES TO FEATURESET")
                for fId in sorted(fId_list):
                    try:
                        #self.log(console,"checking '"+fId+"'")
                        in_filtered_set = hit_seq_ids[fId]
                        #self.log(console, 'FOUND HIT '+fId)  # DEBUG
                        output_featureSet['element_ordering'].append(fId)
                        output_featureSet['elements'][fId] = input_many_featureSet['elements'][fId]
                    except:
                        pass

        # Parse Genome hits into FeatureSet
        #
        elif many_type_name == 'Genome':
            seq_total = 0

            output_featureSet = dict()
            if 'scientific_name' in input_many_genome and input_many_genome['scientific_name'] != None:
                output_featureSet['description'] = input_many_genome['scientific_name'] + " - BLASTx_Search filtered"
            else:
                output_featureSet['description'] = "BLASTx_Search filtered"
            output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()
            for feature in input_many_genome['features']:
                seq_total += 1
                try:
                    in_filtered_set = hit_seq_ids[feature['id']]
                    #self.log(console, 'FOUND HIT: '+feature['id'])  # DEBUG
                    output_featureSet['element_ordering'].append(feature['id'])
                    output_featureSet['elements'][feature['id']] = [input_many_genome_ref]
                except:
                    pass

        # Parse GenomeSet hits into FeatureSet
        #
        elif many_type_name == 'GenomeSet':
            seq_total = 0

            output_featureSet = dict()
            if 'description' in input_many_genomeSet and input_many_genomeSet['description'] != None:
                output_featureSet['description'] = input_many_genomeSet['description'] + " - BLASTx_Search filtered"
            else:
                output_featureSet['description'] = "BLASTx_Search filtered"
            output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()

            for genome_name in input_many_genomeSet['elements'].keys():
                if 'ref' in input_many_genomeSet['elements'][genome_name] and \
                        input_many_genomeSet['elements'][genome_name]['ref'] != None:
                    genomeRef = input_many_genomeSet['elements'][genome_name]['ref']
                    genome = ws.get_objects([{'ref':genomeRef}])[0]['data']
                    for feature in genome['features']:
                        seq_total += 1
                        try:
                            in_filtered_set = hit_seq_ids[feature['id']]
                            #self.log(console, 'FOUND HIT: '+feature['id'])  # DEBUG
                            output_featureSet['element_ordering'].append(feature['id'])
                            output_featureSet['elements'][feature['id']] = [genomeRef]
                        except:
                            pass

                elif 'data' in input_many_genomeSet['elements'][genome_name] and \
                        input_many_genomeSet['elements'][genome_name]['data'] != None:
#                    genome = input_many_genomeSet['elements'][genome_name]['data']
#                    for feature in genome['features']:
#                        #self.log(console,"kbase_id: '"+feature['id']+"'")  # DEBUG
#                        seq_total += 1
#                        try:
#                            in_filtered_set = hit_seq_ids[feature['id']]
#                            #self.log(console, 'FOUND HIT: '+feature['id'])  # DEBUG
#                            output_featureSet['element_ordering'].append(feature['id'])
                    raise ValueError ("FAILURE: unable to address genome object that is stored within 'data' field of genomeSet object")
#                            output_featureSet['elements'][feature['id']] = [genomeRef_is_inside_data_within_genomeSet_object_and_that_cant_be_addressed]
#                        except:
#                            pass


        # load the method provenance from the context object
        #
        self.log(console,"SETTING PROVENANCE")  # DEBUG
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects'] = []
        if 'input_one_name' in params and params['input_one_name'] != None:
            provenance[0]['input_ws_objects'].append(params['workspace_name']+'/'+params['input_one_name'])
        provenance[0]['input_ws_objects'].append(params['workspace_name']+'/'+params['input_many_name'])
        provenance[0]['service'] = 'kb_blast'
        provenance[0]['method'] = 'BLASTx_Search'


        # Upload results
        #
        self.log(console,"UPLOADING RESULTS")  # DEBUG

        # input FeatureSet, Genome, and GenomeSet -> upload FeatureSet output
        new_obj_info = ws.save_objects({
                            'workspace': params['workspace_name'],
                            'objects':[{
                                    'type': 'KBaseCollections.FeatureSet',
                                    'data': output_featureSet,
                                    'name': params['output_filtered_name'],
                                    'meta': {},
                                    'provenance': provenance
                                }]
                        })

        # build output report object
        #
        self.log(console,"BUILDING REPORT")  # DEBUG
        report += 'sequences in many set: '+str(seq_total)+"\n"
        report += 'sequences in hit set:  '+str(hit_total)+"\n"
        report += "\n"
        for line in hit_buf:
            report += line

        reportObj = {
            'objects_created':[{'ref':params['workspace_name']+'/'+params['output_filtered_name'], 'description':'BLASTx_Search hits'}],
            'text_message':report
        }

        reportName = 'blast_report_'+str(hex(uuid.getnode()))
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
        self.log(console,"BLASTx_Search DONE")
        #END BLASTx_Search

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method BLASTx_Search return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]


    def tBLASTn_Search(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN tBLASTn_Search
        console = []
        self.log(console,'Running tBLASTn_Search with params=')
        self.log(console, "\n"+pformat(params))
        report = ''
#        report = 'Running tBLASTn_Search with params='
#        report += "\n"+pformat(params)


        #### do some basic checks
        #
        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
#        if 'input_one_name' not in params and 'input_one_sequence' not in params:
#            raise ValueError('input_one_sequence or input_one_name parameter is required')
        if 'input_one_name' not in params:
            raise ValueError('input_one_name parameter is required')
        if 'input_many_name' not in params:
            raise ValueError('input_many_name parameter is required')
        if 'output_filtered_name' not in params:
            raise ValueError('output_filtered_name parameter is required')



        # Write the input_one_sequence to file
        #
        if 'input_one_sequence' in params \
                and params['input_one_sequence'] != None \
                and params['input_one_sequence'] != "Optionally enter PROTEIN sequence...":
            #input_one_file_name = params['input_one_name']
            input_one_name = 'query.faa'
            input_one_file_name = input_one_name
            one_forward_reads_file_path = os.path.join(self.scratch,input_one_file_name)
            one_forward_reads_file_handle = open(one_forward_reads_file_path, 'w', 0)
            self.log(console, 'writing query reads file: '+str(one_forward_reads_file_path))

#            input_sequence_buf = params['input_one_sequence'].split("\n")
#            one_forward_reads_file_handle.write('>'+params['input_one_name']+"\n")
#            query_line_seen = False
#            for line in input_sequence_buf:
#                if not line.startswith('>'):
#                    one_forward_reads_file_handle.write(line+"\n")
#                else:
#                    if query_line_seen:
#                        break
#                    query_line_seen = True
#            one_forward_reads_file_handle.close();

            input_sequence_buf = params['input_one_sequence']
            space_pattern = re.compile("^[ \t]*$")
            split_input_sequence_buf = input_sequence_buf.split("\n")

            # no header rows, just sequence
            if not input_sequence_buf.startswith('>'):
                one_forward_reads_file_handle.write('>'+input_one_name+"\n")
                for line in split_input_sequence_buf:
                    if not space_pattern.match(line):
                        line = re.sub (" ","",line)
                        line = re.sub ("\t","",line)
                        one_forward_reads_file_handle.write(line.upper()+"\n")
                one_forward_reads_file_handle.close()

            else:
                # write that sucker, removing spaces
                #
                #forward_reads_file_handle.write(input_sequence_buf)        input_sequence_buf = re.sub ('&quot;', '"', input_sequence_buf)
                for i,line in enumerate(split_input_sequence_buf):
                    if line.startswith('>'):
                        record_buf = []
                        record_buf.append(line)
                        for j in range(i+1,len(split_input_sequence_buf)):
                            if split_input_sequence_buf[j].startswith('>'):
                                break
                            seq_line = re.sub (" ","",split_input_sequence_buf[j])
                            seq_line = re.sub ("\t","",seq_line)
                            seq_line = seq_line.upper()
                            record_buf.append(seq_line)
                        record = "\n".join(record_buf)+"\n"
                        one_forward_reads_file_handle.write(record)
                        break  # only want first record
                one_forward_reads_file_handle.close()


        #### Get the input_one object
        ##
        elif 'input_one_name' in params and params['input_one_name'] != None:
            try:
                ws = workspaceService(self.workspaceURL, token=ctx['token'])
                objects = ws.get_objects([{'ref': params['workspace_name']+'/'+params['input_one_name']}])
                data = objects[0]['data']
                info = objects[0]['info']
                # Object Info Contents
                # absolute ref = info[6] + '/' + info[0] + '/' + info[4]
                # 0 - obj_id objid
                # 1 - obj_name name
                # 2 - type_string type
                # 3 - timestamp save_date
                # 4 - int version
                # 5 - username saved_by
                # 6 - ws_id wsid
                # 7 - ws_name workspace
                # 8 - string chsum
                # 9 - int size 
                # 10 - usermeta meta
                one_type_name = info[2].split('.')[1].split('-')[0]
            except Exception as e:
                raise ValueError('Unable to fetch input_one_name object from workspace: ' + str(e))
                #to get the full stack trace: traceback.format_exc()


            # Handle overloading (input_one can be Feature, or FeatureSet)
            #
            if one_type_name == 'FeatureSet':
                # retrieve sequences for features
                input_one_featureSet = data
            
                genome2Features = {}
                features = input_one_featureSet['elements']
                for fId in features.keys():
                    genomeRef = features[fId][0]
                    if genomeRef not in genome2Features:
                        genome2Features[genomeRef] = []
                    genome2Features[genomeRef].append(fId)

                # export features to FASTA file
                one_forward_reads_file_path = os.path.join(self.scratch, params['input_one_name']+".fasta")
                self.log(console, 'writing fasta file: '+one_forward_reads_file_path)
                records = []
                for genomeRef in genome2Features:
                    genome = ws.get_objects([{'ref':genomeRef}])[0]['data']
                    these_genomeFeatureIds = genome2Features[genomeRef]
                    for feature in genome['features']:
                        if feature['id'] in these_genomeFeatureIds:
                            # tBLASTn is prot-nuc
                            #record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description=genomeRef+"."+feature['id'])
                            if feature['type'] != 'CDS':
                                raise ValueError (params['input_one_name']+" feature type must be CDS")
                                sys.exit(0)
                            record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description=genomeRef+"."+feature['id'])
                            records.append(record)
                SeqIO.write(records, one_forward_reads_file_path, "fasta")

            elif one_type_name == 'Feature':
                # export feature to FASTA file
                feature = data
                one_forward_reads_file_path = os.path.join(self.scratch, params['input_one_name']+".fasta")
                self.log(console, 'writing fasta file: '+one_forward_reads_file_path)
                # tBLASTn is prot-nuc
                #record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description='['+feature['genome_id']+']'+' '+feature['function'])
                if feature['type'] != 'CDS':
                    raise ValueError (params['input_one_name']+" feature type must be CDS")
                    sys.exit(0)
                record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description='['+feature['genome_id']+']'+' '+feature['function'])
                SeqIO.write([record], one_forward_reads_file_path, "fasta")

            else:
                raise ValueError('Cannot yet handle input_one type of: '+type_name)            
        else:
            raise ValueError('Must define either input_one_sequence or input_one_name')
            sys.exit (0)


        #### Get the input_many object
        ##
        many_forward_reads_file_compression = None
        sequencing_tech = 'N/A'
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            objects = ws.get_objects([{'ref': params['workspace_name']+'/'+params['input_many_name']}])
            data = objects[0]['data']
            info = objects[0]['info']
            many_type_name = info[2].split('.')[1].split('-')[0]

            if many_type_name == 'SingleEndLibrary':
                many_type_namespace = info[2].split('.')[0]
                if many_type_namespace == 'KBaseAssembly':
                    file_name = data['handle']['file_name']
                elif many_type_namespace == 'KBaseFile':
                    file_name = data['lib']['file']['file_name']
                else:
                    raise ValueError('bad data type namespace: '+many_type_namespace)
                #self.log(console, 'INPUT_MANY_FILENAME: '+file_name)  # DEBUG
                if file_name[-3:] == ".gz":
                    many_forward_reads_file_compression = 'gz'
                if 'sequencing_tech' in data:
                    sequencing_tech = data['sequencing_tech']

        except Exception as e:
            raise ValueError('Unable to fetch input_many_name object from workspace: ' + str(e))
            #to get the full stack trace: traceback.format_exc()

        # Handle overloading (input_many can be SingleEndLibrary, FeatureSet, Genome, or GenomeSet)
        #
        if many_type_name == 'SingleEndLibrary':

            # DEBUG
            #for k in data:
            #    self.log(console,"SingleEndLibrary ["+k+"]: "+str(data[k]))

            try:
                if 'lib' in data:
                    many_forward_reads = data['lib']['file']
                elif 'handle' in data:
                    many_forward_reads = data['handle']
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
            input_many_featureSet = data

            genome2Features = {}
            features = input_many_featureSet['elements']
            for fId in features.keys():
                genomeRef = features[fId][0]
                if genomeRef not in genome2Features:
                    genome2Features[genomeRef] = []
                genome2Features[genomeRef].append(fId)

            # export features to FASTA file
            many_forward_reads_file_path = os.path.join(self.scratch, params['input_many_name']+".fasta")
            self.log(console, 'writing fasta file: '+many_forward_reads_file_path)
            records = []
            feature_written = dict()
            for genomeRef in genome2Features:
                genome = ws.get_objects([{'ref':genomeRef}])[0]['data']
                these_genomeFeatureIds = genome2Features[genomeRef]
                for feature in genome['features']:
                    if feature['id'] in these_genomeFeatureIds:
                        try:
                            f_written = feature_written[feature['id']]
                        except:
                            feature_written[feature['id']] = True
                            #self.log(console,"kbase_id: '"+feature['id']+"'")  # DEBUG
                            # tBLASTn is prot-nuc
                            if feature['type'] != 'CDS':
                                raise ValueError (params['input_many_sequence']+" features must all be CDS type")
                                sys.exit(0)
                            record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description=genome['id'])
                            #record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description=genome['id'])
                            records.append(record)
            SeqIO.write(records, many_forward_reads_file_path, "fasta")


        # Genome
        #
        elif many_type_name == 'Genome':
            input_many_genome = data
            input_many_genome_ref = str(info[6])+'/'+str(info[0])+'/'+str(info[4])

            # export features to FASTA file
            many_forward_reads_file_path = os.path.join(self.scratch, params['input_many_name']+".fasta")
            self.log(console, 'writing fasta file: '+many_forward_reads_file_path)
            records = []
            feature_written = dict()
            for feature in input_many_genome['features']:
                try:
                    f_written = feature_written[feature['id']]
                except:
                    feature_written[feature['id']] = True
                    #self.log(console,"kbase_id: '"+feature['id']+"'")  # DEBUG
                    # tBLASTn is prot-nuc
                    if feature['type'] != 'CDS':
                        continue
                    record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description=input_many_genome['id'])
                    #record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description=input_many_genome['id'])
                    records.append(record)
            SeqIO.write(records, many_forward_reads_file_path, "fasta")


        # GenomeSet
        #
        elif many_type_name == 'GenomeSet':
            input_many_genomeSet = data

            # export features to FASTA file
            many_forward_reads_file_path = os.path.join(self.scratch, params['input_many_name']+".fasta")
            self.log(console, 'writing fasta file: '+many_forward_reads_file_path)

            records = []
            feature_written = dict()
            for genome_name in input_many_genomeSet['elements'].keys():
                if 'ref' in input_many_genomeSet['elements'][genome_name] and \
                         input_many_genomeSet['elements'][genome_name]['ref'] != None:
                    genome = ws.get_objects([{'ref': input_many_genomeSet['elements'][genome_name]['ref']}])[0]['data']
                    for feature in genome['features']:
                        try:
                            f_written = feature_written[feature['id']]
                        except:
                            feature_written[feature['id']] = True
                            #self.log(console,"kbase_id: '"+feature['id']+"'")  # DEBUG
                            # tBLASTn is prot-nuc
                            if feature['type'] != 'CDS':
                                continue
                            record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description=genome['id'])
                            #record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description=genome['id'])
                            records.append(record)

                elif 'data' in input_many_genomeSet['elements'][genome_name] and \
                        input_many_genomeSet['elements'][genome_name]['data'] != None:
                    genome = input_many_genomeSet['elements'][genome_name]['data']
                    for feature in genome['features']:
                        try:
                            f_written = feature_written[feature['id']]
                        except:
                            feature_written[feature['id']] = True
                            #self.log(console,"kbase_id: '"+feature['id']+"'")  # DEBUG
                            # tBLASTn is prot-nuc
                            if feature['type'] != 'CDS':
                                continue
                            record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description=genome['id'])
                            #record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description=genome['id'])
                            records.append(record)

                else:
                    raise ValueError('genome '+genome_name+' missing')

            SeqIO.write(records, many_forward_reads_file_path, "fasta")
            
        # Missing proper input_many_type
        #
        else:
            raise ValueError('Cannot yet handle input_many type of: '+type_name)            

        # FORMAT DB
        #
        # OLD SYNTAX: formatdb -i $database -o T -p F -> $database.nsq or $database.00.nsq
        # NEW SYNTAX: makeblastdb -in $database -parse_seqids -dbtype prot/nucl -out <basename>
        makeblastdb_cmd = [self.Make_BLAST_DB]

        # check for necessary files
        if not os.path.isfile(self.Make_BLAST_DB):
            raise ValueError("no such file '"+self.Make_BLAST_DB+"'")
        if not os.path.isfile(one_forward_reads_file_path):
            raise ValueError("no such file '"+one_forward_reads_file_path+"'")
        if not os.path.isfile(many_forward_reads_file_path):
            raise ValueError("no such file '"+many_forward_reads_file_path+"'")

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


        ### Construct the BLAST command
        #
        # OLD SYNTAX: $blast -q $q -G $G -E $E -m $m -e $e_value -v $limit -b $limit -K $limit -p tblastn -i $fasta_file -d $database -o $out_file
        # NEW SYNTAX: tblastn -query <queryfile> -db <basename> -out <out_aln_file> -outfmt 0/7 (8 became 7) -evalue <e_value> -dust no (DNA) -seg no (DNA) -num_threads <num_cores>
        #
        blast_bin = self.tBLASTn
        blast_cmd = [blast_bin]

        # check for necessary files
        if not os.path.isfile(blast_bin):
            raise ValueError("no such file '"+blast_bin+"'")
        if not os.path.isfile(one_forward_reads_file_path):
            raise ValueError("no such file '"+one_forward_reads_file_path+"'")
        if not os.path.isfile(many_forward_reads_file_path):
            raise ValueError("no such file '"+many_forward_reads_file_path+"'")

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


        # Parse the BLAST tabular output and store ids to filter many set to make filtered object to save back to KBase
        #
        self.log(console, 'PARSING BLAST ALIGNMENT OUTPUT')
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
            #self.log(console,"AFTER ident_thresh")
            if 'bitscore' in params and float(params['bitscore']) > float(high_bitscore_score[hit_seq_id]):
                continue
            #self.log(console,"AFTER bitscore")
            # need to fix this by reading query len
            if 'overlap_fraction' in params and float(params['overlap_fraction']) > 1.0:
                continue
            #self.log(console,"AFTER overlap_fraction")
            
            hit_total += 1
            hit_seq_ids[hit_seq_id] = True
            self.log(console, "HIT: '"+hit_seq_id+"'")  # DEBUG
        

        self.log(console, 'EXTRACTING HITS FROM INPUT')
        self.log(console, 'MANY_TYPE_NAME: '+many_type_name)  # DEBUG


        # SingleEndLibrary input -> SingleEndLibrary output
        #
        if many_type_name == 'SingleEndLibrary':

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
                output_featureSet['description'] = input_many_featureSet['description'] + " - tBLASTn_Search filtered"
            else:
                output_featureSet['description'] = "tBLASTn_Search filtered"
            output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()
            if 'element_ordering' in input_many_featureSet and input_many_featureSet['element_ordering'] != None:
                for fId in input_many_featureSet['element_ordering']:
                    try:
                        in_filtered_set = hit_seq_ids[fId]
                        #self.log(console, 'FOUND HIT '+fId)  # DEBUG
                        output_featureSet['element_ordering'].append(fId)
                        output_featureSet['elements'][fId] = input_many_featureSet['elements'][fId]
                    except:
                        pass
            else:
                fId_list = input_many_featureSet['elements'].keys()
                self.log(console,"ADDING FEATURES TO FEATURESET")
                for fId in sorted(fId_list):
                    try:
                        #self.log(console,"checking '"+fId+"'")
                        in_filtered_set = hit_seq_ids[fId]
                        #self.log(console, 'FOUND HIT '+fId)  # DEBUG
                        output_featureSet['element_ordering'].append(fId)
                        output_featureSet['elements'][fId] = input_many_featureSet['elements'][fId]
                    except:
                        pass

        # Parse Genome hits into FeatureSet
        #
        elif many_type_name == 'Genome':
            seq_total = 0

            output_featureSet = dict()
            if 'scientific_name' in input_many_genome and input_many_genome['scientific_name'] != None:
                output_featureSet['description'] = input_many_genome['scientific_name'] + " - tBLASTn_Search filtered"
            else:
                output_featureSet['description'] = "tBLASTn_Search filtered"
            output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()
            for feature in input_many_genome['features']:
                seq_total += 1
                try:
                    in_filtered_set = hit_seq_ids[feature['id']]
                    #self.log(console, 'FOUND HIT: '+feature['id'])  # DEBUG
                    output_featureSet['element_ordering'].append(feature['id'])
                    output_featureSet['elements'][feature['id']] = [input_many_genome_ref]
                except:
                    pass

        # Parse GenomeSet hits into FeatureSet
        #
        elif many_type_name == 'GenomeSet':
            seq_total = 0

            output_featureSet = dict()
            if 'description' in input_many_genomeSet and input_many_genomeSet['description'] != None:
                output_featureSet['description'] = input_many_genomeSet['description'] + " - tBLASTn_Search filtered"
            else:
                output_featureSet['description'] = "tBLASTn_Search filtered"
            output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()

            for genome_name in input_many_genomeSet['elements'].keys():
                if 'ref' in input_many_genomeSet['elements'][genome_name] and \
                        input_many_genomeSet['elements'][genome_name]['ref'] != None:
                    genomeRef = input_many_genomeSet['elements'][genome_name]['ref']
                    genome = ws.get_objects([{'ref':genomeRef}])[0]['data']
                    for feature in genome['features']:
                        seq_total += 1
                        try:
                            in_filtered_set = hit_seq_ids[feature['id']]
                            #self.log(console, 'FOUND HIT: '+feature['id'])  # DEBUG
                            output_featureSet['element_ordering'].append(feature['id'])
                            output_featureSet['elements'][feature['id']] = [genomeRef]
                        except:
                            pass

                elif 'data' in input_many_genomeSet['elements'][genome_name] and \
                        input_many_genomeSet['elements'][genome_name]['data'] != None:
#                    genome = input_many_genomeSet['elements'][genome_name]['data']
#                    for feature in genome['features']:
#                        #self.log(console,"kbase_id: '"+feature['id']+"'")  # DEBUG
#                        seq_total += 1
#                        try:
#                            in_filtered_set = hit_seq_ids[feature['id']]
#                            #self.log(console, 'FOUND HIT: '+feature['id'])  # DEBUG
#                            output_featureSet['element_ordering'].append(feature['id'])
                    raise ValueError ("FAILURE: unable to address genome object that is stored within 'data' field of genomeSet object")
#                            output_featureSet['elements'][feature['id']] = [genomeRef_is_inside_data_within_genomeSet_object_and_that_cant_be_addressed]
#                        except:
#                            pass


        # load the method provenance from the context object
        #
        self.log(console,"SETTING PROVENANCE")  # DEBUG
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects'] = []
        if 'input_one_name' in params and params['input_one_name'] != None:
            provenance[0]['input_ws_objects'].append(params['workspace_name']+'/'+params['input_one_name'])
        provenance[0]['input_ws_objects'].append(params['workspace_name']+'/'+params['input_many_name'])
        provenance[0]['service'] = 'kb_blast'
        provenance[0]['method'] = 'tBLASTn_Search'


        # Upload results
        #
        self.log(console,"UPLOADING RESULTS")  # DEBUG

        if many_type_name == 'SingleEndLibrary':
            
            # input SingleEndLibrary -> upload SingleEndLibrary
            #
            self.upload_SingleEndLibrary_to_shock_and_ws (ctx,
                                                          console,  # DEBUG
                                                          params['workspace_name'],
                                                          params['output_filtered_name'],
                                                          output_filtered_fasta_file_path,
                                                          provenance,
                                                          sequencing_tech
                                                         )

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
                        })

        # build output report object
        #
        self.log(console,"BUILDING REPORT")  # DEBUG
        report += 'sequences in many set: '+str(seq_total)+"\n"
        report += 'sequences in hit set:  '+str(hit_total)+"\n"
        report += "\n"
        for line in hit_buf:
            report += line

        reportObj = {
            'objects_created':[{'ref':params['workspace_name']+'/'+params['output_filtered_name'], 'description':'tBLASTn_Search hits'}],
            'text_message':report
        }

        reportName = 'blast_report_'+str(hex(uuid.getnode()))
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
        self.log(console,"tBLASTn_Search DONE")
        #END tBLASTn_Search

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method tBLASTn_Search return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]


    def tBLASTx_Search(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN tBLASTx_Search
        console = []
        self.log(console,'Running tBLASTx_Search with params=')
        self.log(console, "\n"+pformat(params))
        report = ''
#        report = 'Running tBLASTx_Search with params='
#        report += "\n"+pformat(params)


        #### do some basic checks
        #
        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
#        if 'input_one_name' not in params and 'input_one_sequence' not in params:
#            raise ValueError('input_one_sequence or input_one_name parameter is required')
        if 'input_one_name' not in params:
            raise ValueError('input_one_name parameter is required')
        if 'input_many_name' not in params:
            raise ValueError('input_many_name parameter is required')
        if 'output_filtered_name' not in params:
            raise ValueError('output_filtered_name parameter is required')


        # Write the input_one_sequence to a SingleEndLibrary object
        #
        if 'input_one_sequence' in params \
                and params['input_one_sequence'] != None \
                and params['input_one_sequence'] != "Optionally enter DNA sequence...":
            input_one_file_name = params['input_one_name']
            one_forward_reads_file_path = os.path.join(self.scratch,input_one_file_name)
            one_forward_reads_file_handle = open(one_forward_reads_file_path, 'w', 0)
            self.log(console, 'writing query reads file: '+str(one_forward_reads_file_path))

#            input_sequence_buf = params['input_one_sequence'].split("\n")
#            one_forward_reads_file_handle.write('>'+params['input_one_name']+"\n")
#            query_line_seen = False
#            for line in input_sequence_buf:
#                if not line.startswith('>'):
#                    one_forward_reads_file_handle.write(line+"\n")
#                else:
#                    if query_line_seen:
#                        break
#                    query_line_seen = True
#            one_forward_reads_file_handle.close();

            fastq_format = False
            input_sequence_buf = params['input_one_sequence']
            if input_sequence_buf.startswith('@'):
                fastq_format = True
                #self.log(console,"INPUT_SEQ BEFORE: '''\n"+input_sequence_buf+"\n'''")  # DEBUG
            input_sequence_buf = re.sub ('&apos;', "'", input_sequence_buf)
            input_sequence_buf = re.sub ('&quot;', '"', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#39;',  "'", input_sequence_buf)
#        input_sequence_buf = re.sub ('&#34;',  '"', input_sequence_buf)
#        input_sequence_buf = re.sub ('&lt;;',  '<', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#60;',  '<', input_sequence_buf)
#        input_sequence_buf = re.sub ('&gt;',   '>', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#62;',  '>', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#36;',  '$', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#37;',  '%', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#47;',  '/', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#63;',  '?', input_sequence_buf)
##        input_sequence_buf = re.sub ('&#92;',  chr(92), input_sequence_buf)  # FIX LATER
#        input_sequence_buf = re.sub ('&#96;',  '`', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#124;', '|', input_sequence_buf)
#        input_sequence_buf = re.sub ('&amp;', '&', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#38;', '&', input_sequence_buf)
#        self.log(console,"INPUT_SEQ AFTER: '''\n"+input_sequence_buf+"\n'''")  # DEBUG

            DNA_pattern = re.compile("^[acgtuACGTU ]+$")
            space_pattern = re.compile("^[ \t]*$")
            split_input_sequence_buf = input_sequence_buf.split("\n")

            # no header rows, just sequence
            if not input_sequence_buf.startswith('>') and not input_sequence_buf.startswith('@'):
                one_forward_reads_file_handle.write('>'+params['input_one_name']+"\n")
                for line in split_input_sequence_buf:
                    if not space_pattern.match(line):
                        line = re.sub (" ","",line)
                        line = re.sub ("\t","",line)
                        if not DNA_pattern.match(line):
                            raise ValueError ("BAD record:\n"+line+"\n")
                            sys.exit(0)
                        one_forward_reads_file_handle.write(line.lower()+"\n")
                one_forward_reads_file_handle.close()

            else:
                # format checks
                for i,line in enumerate(split_input_sequence_buf):
                    if line.startswith('>') or line.startswith('@'):
                        if not DNA_pattern.match(split_input_sequence_buf[i+1]):
                            if fastq_format:
                                bad_record = "\n".join([split_input_sequence_buf[i],
                                                        split_input_sequence_buf[i+1],
                                                        split_input_sequence_buf[i+2],
                                                        split_input_sequence_buf[i+3]])
                            else:
                                bad_record = "\n".join([split_input_sequence_buf[i],
                                                    split_input_sequence_buf[i+1]])
                            raise ValueError ("BAD record:\n"+bad_record+"\n")
                            sys.exit(0)
                        if fastq_format and line.startswith('@'):
                            format_ok = True
                            seq_len = len(split_input_sequence_buf[i+1])
                            if not seq_len > 0:
                                format_ok = False
                            if not split_input_sequence_buf[i+2].startswith('+'):
                                format_ok = False
                            if not seq_len == len(split_input_sequence_buf[i+3]):
                                format_ok = False
                            if not format_ok:
                                bad_record = "\n".join([split_input_sequence_buf[i],
                                                    split_input_sequence_buf[i+1],
                                                    split_input_sequence_buf[i+2],
                                                    split_input_sequence_buf[i+3]])
                                raise ValueError ("BAD record:\n"+bad_record+"\n")
                                sys.exit(0)

                # write that sucker, removing spaces
                #
                #forward_reads_file_handle.write(input_sequence_buf)        input_sequence_buf = re.sub ('&quot;', '"', input_sequence_buf)
                for i,line in enumerate(split_input_sequence_buf):
                    if line.startswith('>'):
                        record_buf = []
                        record_buf.append(line)
                        for j in range(i+1,len(split_input_sequence_buf)):
                            if split_input_sequence_buf[j].startswith('>'):
                                break
                            seq_line = re.sub (" ","",split_input_sequence_buf[j])
                            seq_line = re.sub ("\t","",seq_line)
                            seq_line = seq_line.lower()
                            record_buf.append(seq_line)
                        record = "\n".join(record_buf)+"\n"
                        one_forward_reads_file_handle.write(record)
                        break  # only want first record
                    elif line.startswith('@'):
                        seq_line = re.sub (" ","",split_input_sequence_buf[i+1])
                        seq_line = re.sub ("\t","",seq_line)
                        seq_line = seq_line.lower()
                        qual_line = re.sub (" ","",split_input_sequence_buf[i+3])
                        qual_line = re.sub ("\t","",qual_line)
                        record = "\n".join([line, seq_line, split_input_sequence_buf[i+2], qual_line])+"\n"
                        one_forward_reads_file_handle.write(record)
                        break  # only want first record

                one_forward_reads_file_handle.close()


            # load the method provenance from the context object
            #
            self.log(console,"SETTING PROVENANCE")  # DEBUG
            provenance = [{}]
            if 'provenance' in ctx:
                provenance = ctx['provenance']
            # add additional info to provenance here, in this case the input data object reference
                provenance[0]['input_ws_objects'] = []
                provenance[0]['service'] = 'kb_blast'
                provenance[0]['method'] = 'tBLASTx_Search'

                
                # Upload results
                #
                self.log(console,"UPLOADING QUERY OBJECT")  # DEBUG

                sequencing_tech = 'N/A'
                self.upload_SingleEndLibrary_to_shock_and_ws (ctx,
                                                      console,  # DEBUG
                                                      params['workspace_name'],
                                                      params['input_one_name'],
                                                      one_forward_reads_file_path,
                                                      provenance,
                                                      sequencing_tech
                                                      )

            self.log(console, 'done')

        #### Get the input_one object
        ##
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            objects = ws.get_objects([{'ref': params['workspace_name']+'/'+params['input_one_name']}])
            data = objects[0]['data']
            info = objects[0]['info']
            # Object Info Contents
            # absolute ref = info[6] + '/' + info[0] + '/' + info[4]
            # 0 - obj_id objid
            # 1 - obj_name name
            # 2 - type_string type
            # 3 - timestamp save_date
            # 4 - int version
            # 5 - username saved_by
            # 6 - ws_id wsid
            # 7 - ws_name workspace
            # 8 - string chsum
            # 9 - int size 
            # 10 - usermeta meta
            one_type_name = info[2].split('.')[1].split('-')[0]
        except Exception as e:
            raise ValueError('Unable to fetch input_one_name object from workspace: ' + str(e))
        #to get the full stack trace: traceback.format_exc()

        if 'input_one_sequence' in params \
                and params['input_one_sequence'] != None \
                and params['input_one_sequence'] != "Optionally enter DNA sequence..." \
                and one_type_name != 'SingleEndLibrary':

            raise ValueError("ERROR: Mismatched input type: input_one_name should be SingleEndLibrary instead of: "+one_type_name)
            sys.exit (0)


        # Handle overloading (input_one can be Feature, SingleEndLibrary, or FeatureSet)
        #
        if one_type_name == 'SingleEndLibrary':
            try:
                if 'lib' in data:
                    one_forward_reads = data['lib']['file']
                elif 'handle' in data:
                    one_forward_reads = data['handle']
                else:
                    self.log(console,"bad structure for 'one_forward_reads'")
                    raise ValueError("bad structure for 'one_forward_reads'")

                ### NOTE: this section is what could be replaced by the transform services
                one_forward_reads_file_path = os.path.join(self.scratch,one_forward_reads['file_name'])
                one_forward_reads_file_handle = open(one_forward_reads_file_path, 'w', 0)
                self.log(console, 'downloading reads file: '+str(one_forward_reads_file_path))
                headers = {'Authorization': 'OAuth '+ctx['token']}
                r = requests.get(one_forward_reads['url']+'/node/'+one_forward_reads['id']+'?download', stream=True, headers=headers)
                for chunk in r.iter_content(1024):
                    one_forward_reads_file_handle.write(chunk)
                one_forward_reads_file_handle.close();
                self.log(console, 'done')
                ### END NOTE


                # remove carriage returns
                new_file_path = one_forward_reads_file_path+"-CRfree"
                new_file_handle = open(new_file_path, 'w', 0)
                one_forward_reads_file_handle = open(one_forward_reads_file_path, 'r', 0)
                for line in one_forward_reads_file_handle:
                    line = re.sub("\r","",line)
                    new_file_handle.write(line)
                one_forward_reads_file_handle.close();
                new_file_handle.close()
                one_forward_reads_file_path = new_file_path


                # convert FASTQ to FASTA (if necessary)
                new_file_path = one_forward_reads_file_path+".fna"
                new_file_handle = open(new_file_path, 'w', 0)
                one_forward_reads_file_handle = open(one_forward_reads_file_path, 'r', 0)
                header = None
                last_header = None
                last_seq_buf = None
                last_line_was_header = False
                was_fastq = False
                for line in one_forward_reads_file_handle:
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
                one_forward_reads_file_handle.close()
                if was_fastq:
                    one_forward_reads_file_path = new_file_path

            except Exception as e:
                print(traceback.format_exc())
                raise ValueError('Unable to download single-end read library files: ' + str(e))

        elif one_type_name == 'FeatureSet':
            # retrieve sequences for features
            input_one_featureSet = data
            
            genome2Features = {}
            features = input_one_featureSet['elements']
            for fId in features.keys():
                genomeRef = features[fId][0]
                if genomeRef not in genome2Features:
                    genome2Features[genomeRef] = []
                genome2Features[genomeRef].append(fId)

            # export features to FASTA file
            one_forward_reads_file_path = os.path.join(self.scratch, params['input_one_name']+".fasta")
            self.log(console, 'writing fasta file: '+one_forward_reads_file_path)
            records = []
            for genomeRef in genome2Features:
                genome = ws.get_objects([{'ref':genomeRef}])[0]['data']
                these_genomeFeatureIds = genome2Features[genomeRef]
                for feature in genome['features']:
                    if feature['id'] in these_genomeFeatureIds:
                        # tBLASTx is nuc-nuc (translated)
                        if feature['type'] != 'CDS':
                            raise ValueError (params['input_one_name']+" feature type must be CDS")
                            sys.exit(0)
                        record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description=genomeRef+"."+feature['id'])
                        #record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description=genomeRef+"."+feature['id'])
                        records.append(record)
            SeqIO.write(records, one_forward_reads_file_path, "fasta")

        elif one_type_name == 'Feature':
            # export feature to FASTA file
            feature = data
            one_forward_reads_file_path = os.path.join(self.scratch, params['input_one_name']+".fasta")
            self.log(console, 'writing fasta file: '+one_forward_reads_file_path)
            # tBLASTx is nuc-nuc (translated)
            if feature['type'] != 'CDS':
                raise ValueError (params['input_one_name']+" feature type must be CDS")
                sys.exit(0)
            record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description='['+feature['genome_id']+']'+' '+feature['function'])
            #record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description='['+feature['genome_id']+']'+' '+feature['function'])
            SeqIO.write([record], one_forward_reads_file_path, "fasta")

        else:
            raise ValueError('Cannot yet handle input_one type of: '+type_name)            

        #### Get the input_many object
        ##
        many_forward_reads_file_compression = None
        sequencing_tech = 'N/A'
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            objects = ws.get_objects([{'ref': params['workspace_name']+'/'+params['input_many_name']}])
            data = objects[0]['data']
            info = objects[0]['info']
            many_type_name = info[2].split('.')[1].split('-')[0]

            if many_type_name == 'SingleEndLibrary':
                many_type_namespace = info[2].split('.')[0]
                if many_type_namespace == 'KBaseAssembly':
                    file_name = data['handle']['file_name']
                elif many_type_namespace == 'KBaseFile':
                    file_name = data['lib']['file']['file_name']
                else:
                    raise ValueError('bad data type namespace: '+many_type_namespace)
                #self.log(console, 'INPUT_MANY_FILENAME: '+file_name)  # DEBUG
                if file_name[-3:] == ".gz":
                    many_forward_reads_file_compression = 'gz'
                if 'sequencing_tech' in data:
                    sequencing_tech = data['sequencing_tech']

        except Exception as e:
            raise ValueError('Unable to fetch input_many_name object from workspace: ' + str(e))
            #to get the full stack trace: traceback.format_exc()

        # Handle overloading (input_many can be SingleEndLibrary, FeatureSet, Genome, or GenomeSet)
        #
        if many_type_name == 'SingleEndLibrary':

            # DEBUG
            #for k in data:
            #    self.log(console,"SingleEndLibrary ["+k+"]: "+str(data[k]))

            try:
                if 'lib' in data:
                    many_forward_reads = data['lib']['file']
                elif 'handle' in data:
                    many_forward_reads = data['handle']
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
            input_many_featureSet = data

            genome2Features = {}
            features = input_many_featureSet['elements']
            for fId in features.keys():
                genomeRef = features[fId][0]
                if genomeRef not in genome2Features:
                    genome2Features[genomeRef] = []
                genome2Features[genomeRef].append(fId)

            # export features to FASTA file
            many_forward_reads_file_path = os.path.join(self.scratch, params['input_many_name']+".fasta")
            self.log(console, 'writing fasta file: '+many_forward_reads_file_path)
            records = []
            feature_written = dict()
            for genomeRef in genome2Features:
                genome = ws.get_objects([{'ref':genomeRef}])[0]['data']
                these_genomeFeatureIds = genome2Features[genomeRef]
                for feature in genome['features']:
                    if feature['id'] in these_genomeFeatureIds:
                        try:
                            f_written = feature_written[feature['id']]
                        except:
                            feature_written[feature['id']] = True
                            #self.log(console,"kbase_id: '"+feature['id']+"'")  # DEBUG
                            # tBLASTx is nuc-nuc (translated)
                            if feature['type'] != 'CDS':
                                raise ValueError (params['input_many_sequence']+" features must all be CDS type")
                                sys.exit(0)
                            record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description=genome['id'])
                            #record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description=genome['id'])
                            records.append(record)
            SeqIO.write(records, many_forward_reads_file_path, "fasta")


        # Genome
        #
        elif many_type_name == 'Genome':
            input_many_genome = data
            input_many_genome_ref = str(info[6])+'/'+str(info[0])+'/'+str(info[4])

            # export features to FASTA file
            many_forward_reads_file_path = os.path.join(self.scratch, params['input_many_name']+".fasta")
            self.log(console, 'writing fasta file: '+many_forward_reads_file_path)
            records = []
            feature_written = dict()
            for feature in input_many_genome['features']:
                try:
                    f_written = feature_written[feature['id']]
                except:
                    feature_written[feature['id']] = True
                    #self.log(console,"kbase_id: '"+feature['id']+"'")  # DEBUG
                    # tBLASTx is nuc-nuc (translated)
                    if feature['type'] != 'CDS':
                        continue
                    record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description=input_many_genome['id'])
                    #record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description=input_many_genome['id'])
                    records.append(record)
            SeqIO.write(records, many_forward_reads_file_path, "fasta")


        # GenomeSet
        #
        elif many_type_name == 'GenomeSet':
            input_many_genomeSet = data

            # export features to FASTA file
            many_forward_reads_file_path = os.path.join(self.scratch, params['input_many_name']+".fasta")
            self.log(console, 'writing fasta file: '+many_forward_reads_file_path)

            records = []
            feature_written = dict()
            for genome_name in input_many_genomeSet['elements'].keys():
                if 'ref' in input_many_genomeSet['elements'][genome_name] and \
                         input_many_genomeSet['elements'][genome_name]['ref'] != None:
                    genome = ws.get_objects([{'ref': input_many_genomeSet['elements'][genome_name]['ref']}])[0]['data']
                    for feature in genome['features']:
                        try:
                            f_written = feature_written[feature['id']]
                        except:
                            feature_written[feature['id']] = True
                            #self.log(console,"kbase_id: '"+feature['id']+"'")  # DEBUG
                            # tBLASTx is nuc-nuc (translated)
                            if feature['type'] != 'CDS':
                                continue
                            record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description=genome['id'])
                            #record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description=genome['id'])
                            records.append(record)

                elif 'data' in input_many_genomeSet['elements'][genome_name] and \
                        input_many_genomeSet['elements'][genome_name]['data'] != None:
                    genome = input_many_genomeSet['elements'][genome_name]['data']
                    for feature in genome['features']:
                        try:
                            f_written = feature_written[feature['id']]
                        except:
                            feature_written[feature['id']] = True
                            #self.log(console,"kbase_id: '"+feature['id']+"'")  # DEBUG
                            # tBLASTx is nuc-nuc (translated)
                            if feature['type'] != 'CDS':
                                continue
                            record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description=genome['id'])
                            #record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description=genome['id'])
                            records.append(record)

                else:
                    raise ValueError('genome '+genome_name+' missing')

            SeqIO.write(records, many_forward_reads_file_path, "fasta")
            
        # Missing proper input_many_type
        #
        else:
            raise ValueError('Cannot yet handle input_many type of: '+type_name)            

        # FORMAT DB
        #
        # OLD SYNTAX: formatdb -i $database -o T -p F -> $database.nsq or $database.00.nsq
        # NEW SYNTAX: makeblastdb -in $database -parse_seqids -dbtype prot/nucl -out <basename>
        makeblastdb_cmd = [self.Make_BLAST_DB]

        # check for necessary files
        if not os.path.isfile(self.Make_BLAST_DB):
            raise ValueError("no such file '"+self.Make_BLAST_DB+"'")
        if not os.path.isfile(one_forward_reads_file_path):
            raise ValueError("no such file '"+one_forward_reads_file_path+"'")
        if not os.path.isfile(many_forward_reads_file_path):
            raise ValueError("no such file '"+many_forward_reads_file_path+"'")

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


        ### Construct the BLAST command
        #
        # OLD SYNTAX: $blast -q $q -G $G -E $E -m $m -e $e_value -v $limit -b $limit -K $limit -p tblastx -i $fasta_file -d $database -o $out_file
        # NEW SYNTAX: tblastx -query <queryfile> -db <basename> -out <out_aln_file> -outfmt 0/7 (8 became 7) -evalue <e_value> -dust no (DNA) -seg no (DNA) -num_threads <num_cores>
        #
        blast_bin = self.tBLASTx
        blast_cmd = [blast_bin]

        # check for necessary files
        if not os.path.isfile(blast_bin):
            raise ValueError("no such file '"+blast_bin+"'")
        if not os.path.isfile(one_forward_reads_file_path):
            raise ValueError("no such file '"+one_forward_reads_file_path+"'")
        if not os.path.isfile(many_forward_reads_file_path):
            raise ValueError("no such file '"+many_forward_reads_file_path+"'")

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


        # Parse the BLAST tabular output and store ids to filter many set to make filtered object to save back to KBase
        #
        self.log(console, 'PARSING BLAST ALIGNMENT OUTPUT')
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
            #self.log(console,"AFTER ident_thresh")
            if 'bitscore' in params and float(params['bitscore']) > float(high_bitscore_score[hit_seq_id]):
                continue
            #self.log(console,"AFTER bitscore")
            # need to fix this by reading query len
            if 'overlap_fraction' in params and float(params['overlap_fraction']) > 1.0:
                continue
            #self.log(console,"AFTER overlap_fraction")
            
            hit_total += 1
            hit_seq_ids[hit_seq_id] = True
            self.log(console, "HIT: '"+hit_seq_id+"'")  # DEBUG
        

        self.log(console, 'EXTRACTING HITS FROM INPUT')
        self.log(console, 'MANY_TYPE_NAME: '+many_type_name)  # DEBUG


        # SingleEndLibrary input -> SingleEndLibrary output
        #
        if many_type_name == 'SingleEndLibrary':

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
                output_featureSet['description'] = input_many_featureSet['description'] + " - tBLASTx_Search filtered"
            else:
                output_featureSet['description'] = "tBLASTx_Search filtered"
            output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()
            if 'element_ordering' in input_many_featureSet and input_many_featureSet['element_ordering'] != None:
                for fId in input_many_featureSet['element_ordering']:
                    try:
                        in_filtered_set = hit_seq_ids[fId]
                        #self.log(console, 'FOUND HIT '+fId)  # DEBUG
                        output_featureSet['element_ordering'].append(fId)
                        output_featureSet['elements'][fId] = input_many_featureSet['elements'][fId]
                    except:
                        pass
            else:
                fId_list = input_many_featureSet['elements'].keys()
                self.log(console,"ADDING FEATURES TO FEATURESET")
                for fId in sorted(fId_list):
                    try:
                        #self.log(console,"checking '"+fId+"'")
                        in_filtered_set = hit_seq_ids[fId]
                        #self.log(console, 'FOUND HIT '+fId)  # DEBUG
                        output_featureSet['element_ordering'].append(fId)
                        output_featureSet['elements'][fId] = input_many_featureSet['elements'][fId]
                    except:
                        pass

        # Parse Genome hits into FeatureSet
        #
        elif many_type_name == 'Genome':
            seq_total = 0

            output_featureSet = dict()
            if 'scientific_name' in input_many_genome and input_many_genome['scientific_name'] != None:
                output_featureSet['description'] = input_many_genome['scientific_name'] + " - tBLASTx_Search filtered"
            else:
                output_featureSet['description'] = "tBLASTx_Search filtered"
            output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()
            for feature in input_many_genome['features']:
                seq_total += 1
                try:
                    in_filtered_set = hit_seq_ids[feature['id']]
                    #self.log(console, 'FOUND HIT: '+feature['id'])  # DEBUG
                    output_featureSet['element_ordering'].append(feature['id'])
                    output_featureSet['elements'][feature['id']] = [input_many_genome_ref]
                except:
                    pass

        # Parse GenomeSet hits into FeatureSet
        #
        elif many_type_name == 'GenomeSet':
            seq_total = 0

            output_featureSet = dict()
            if 'description' in input_many_genomeSet and input_many_genomeSet['description'] != None:
                output_featureSet['description'] = input_many_genomeSet['description'] + " - tBLASTx_Search filtered"
            else:
                output_featureSet['description'] = "tBLASTx_Search filtered"
            output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()

            for genome_name in input_many_genomeSet['elements'].keys():
                if 'ref' in input_many_genomeSet['elements'][genome_name] and \
                        input_many_genomeSet['elements'][genome_name]['ref'] != None:
                    genomeRef = input_many_genomeSet['elements'][genome_name]['ref']
                    genome = ws.get_objects([{'ref':genomeRef}])[0]['data']
                    for feature in genome['features']:
                        seq_total += 1
                        try:
                            in_filtered_set = hit_seq_ids[feature['id']]
                            #self.log(console, 'FOUND HIT: '+feature['id'])  # DEBUG
                            output_featureSet['element_ordering'].append(feature['id'])
                            output_featureSet['elements'][feature['id']] = [genomeRef]
                        except:
                            pass

                elif 'data' in input_many_genomeSet['elements'][genome_name] and \
                        input_many_genomeSet['elements'][genome_name]['data'] != None:
#                    genome = input_many_genomeSet['elements'][genome_name]['data']
#                    for feature in genome['features']:
#                        #self.log(console,"kbase_id: '"+feature['id']+"'")  # DEBUG
#                        seq_total += 1
#                        try:
#                            in_filtered_set = hit_seq_ids[feature['id']]
#                            #self.log(console, 'FOUND HIT: '+feature['id'])  # DEBUG
#                            output_featureSet['element_ordering'].append(feature['id'])
                    raise ValueError ("FAILURE: unable to address genome object that is stored within 'data' field of genomeSet object")
#                            output_featureSet['elements'][feature['id']] = [genomeRef_is_inside_data_within_genomeSet_object_and_that_cant_be_addressed]
#                        except:
#                            pass


        # load the method provenance from the context object
        #
        self.log(console,"SETTING PROVENANCE")  # DEBUG
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects'] = []
        if 'input_one_name' in params and params['input_one_name'] != None:
            provenance[0]['input_ws_objects'].append(params['workspace_name']+'/'+params['input_one_name'])
        provenance[0]['input_ws_objects'].append(params['workspace_name']+'/'+params['input_many_name'])
        provenance[0]['service'] = 'kb_blast'
        provenance[0]['method'] = 'tBLASTx_Search'


        # Upload results
        #
        self.log(console,"UPLOADING RESULTS")  # DEBUG

        if many_type_name == 'SingleEndLibrary':
            
            # input SingleEndLibrary -> upload SingleEndLibrary
            #
            self.upload_SingleEndLibrary_to_shock_and_ws (ctx,
                                                          console,  # DEBUG
                                                          params['workspace_name'],
                                                          params['output_filtered_name'],
                                                          output_filtered_fasta_file_path,
                                                          provenance,
                                                          sequencing_tech
                                                         )

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
                        })

        # build output report object
        #
        self.log(console,"BUILDING REPORT")  # DEBUG
        report += 'sequences in many set: '+str(seq_total)+"\n"
        report += 'sequences in hit set:  '+str(hit_total)+"\n"
        report += "\n"
        for line in hit_buf:
            report += line

        reportObj = {
            'objects_created':[{'ref':params['workspace_name']+'/'+params['output_filtered_name'], 'description':'tBLASTx_Search hits'}],
            'text_message':report
        }

        reportName = 'blast_report_'+str(hex(uuid.getnode()))
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
        self.log(console,"tBLASTx_Search DONE")
        #END tBLASTx_Search

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method tBLASTx_Search return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]
