# -*- coding: utf-8 -*-
#BEGIN_HEADER
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

# BlastUtil
from kb_blast.Utils.BlastUtil import BlastUtil

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
    VERSION = "1.4.0"
    GIT_URL = "https://github.com/kbaseapps/kb_blast"
    GIT_COMMIT_HASH = "32a1bb326a2a92463136a7f9ad5e9daf3b7acaa2"

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


    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config
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
           "data_obj_ref", parameter "output_filtered_name" of type
           "data_obj_name", parameter "ident_thresh" of Double, parameter
           "e_value" of Double, parameter "bitscore" of Double, parameter
           "overlap_fraction" of Double, parameter "maxaccepts" of Double,
           parameter "output_extra_format" of String, parameter "rounds" of
           Double
        :returns: instance of type "BLAST_Output" (BLAST Output) ->
           structure: parameter "report_name" of type "data_obj_name",
           parameter "report_ref" of type "data_obj_ref"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN BLASTn_Search
        search_tool_name = 'BLASTn'
        bu = BlastUtil(self.config, ctx)
        returnVal = bu.run_BLAST_App (search_tool_name, params)
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
           "data_obj_ref", parameter "output_filtered_name" of type
           "data_obj_name", parameter "ident_thresh" of Double, parameter
           "e_value" of Double, parameter "bitscore" of Double, parameter
           "overlap_fraction" of Double, parameter "maxaccepts" of Double,
           parameter "output_extra_format" of String, parameter "rounds" of
           Double
        :returns: instance of type "BLAST_Output" (BLAST Output) ->
           structure: parameter "report_name" of type "data_obj_name",
           parameter "report_ref" of type "data_obj_ref"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN BLASTp_Search
        search_tool_name = 'BLASTp'
        bu = BlastUtil(self.config, ctx)
        returnVal = bu.run_BLAST_App (search_tool_name, params)
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
           "data_obj_ref", parameter "output_filtered_name" of type
           "data_obj_name", parameter "ident_thresh" of Double, parameter
           "e_value" of Double, parameter "bitscore" of Double, parameter
           "overlap_fraction" of Double, parameter "maxaccepts" of Double,
           parameter "output_extra_format" of String, parameter "rounds" of
           Double
        :returns: instance of type "BLAST_Output" (BLAST Output) ->
           structure: parameter "report_name" of type "data_obj_name",
           parameter "report_ref" of type "data_obj_ref"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN BLASTx_Search
        search_tool_name = 'BLASTx'
        bu = BlastUtil(self.config, ctx)
        returnVal = bu.run_BLAST_App (search_tool_name, params)
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
           "data_obj_ref", parameter "output_filtered_name" of type
           "data_obj_name", parameter "ident_thresh" of Double, parameter
           "e_value" of Double, parameter "bitscore" of Double, parameter
           "overlap_fraction" of Double, parameter "maxaccepts" of Double,
           parameter "output_extra_format" of String, parameter "rounds" of
           Double
        :returns: instance of type "BLAST_Output" (BLAST Output) ->
           structure: parameter "report_name" of type "data_obj_name",
           parameter "report_ref" of type "data_obj_ref"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN tBLASTn_Search
        search_tool_name = 'tBLASTn'
        bu = BlastUtil(self.config, ctx)
        returnVal = bu.run_BLAST_App (search_tool_name, params)
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
           "data_obj_ref", parameter "output_filtered_name" of type
           "data_obj_name", parameter "ident_thresh" of Double, parameter
           "e_value" of Double, parameter "bitscore" of Double, parameter
           "overlap_fraction" of Double, parameter "maxaccepts" of Double,
           parameter "output_extra_format" of String, parameter "rounds" of
           Double
        :returns: instance of type "BLAST_Output" (BLAST Output) ->
           structure: parameter "report_name" of type "data_obj_name",
           parameter "report_ref" of type "data_obj_ref"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN tBLASTx_Search
        search_tool_name = 'tBLASTx'
        bu = BlastUtil(self.config, ctx)
        returnVal = bu.run_BLAST_App (search_tool_name, params)
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
           "data_obj_ref", parameter "output_filtered_name" of type
           "data_obj_name", parameter "ident_thresh" of Double, parameter
           "e_value" of Double, parameter "bitscore" of Double, parameter
           "overlap_fraction" of Double, parameter "maxaccepts" of Double,
           parameter "output_extra_format" of String, parameter "rounds" of
           Double
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
        #appropriate_sequence_found_in_one_input = False
        appropriate_sequence_found_in_MSA_input = False
        appropriate_sequence_found_in_many_input = False
        genome_id_feature_id_delim = '.f:'


        #### do some basic checks
        #
        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
        #if 'input_one_ref' not in params:
        #    raise ValueError('input_one_ref parameter is required')
        if 'input_msa_ref' not in params:
            raise ValueError('input_msa_ref parameter is required')
        if 'input_many_ref' not in params:
            raise ValueError('input_many_ref parameter is required')
        if 'output_filtered_name' not in params:
            raise ValueError('output_filtered_name parameter is required')


        # set local names
        #input_one_ref = params['input_one_ref']
        input_msa_ref = params['input_msa_ref']
        input_many_ref = params['input_many_ref']
        

        #### Get the input_msa object
        ##
#        if input_one_feature_id == None:
#            self.log(invalid_msgs,"input_one_feature_id was not obtained from Query Object: "+input_one_name)
#       master_row_idx = 0
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            #objects = ws.get_objects([{'ref': input_msa_ref}])
            objects = ws.get_objects2({'objects':[{'ref': input_msa_ref}]})['data']
            input_msa_data = objects[0]['data']
            info = objects[0]['info']
            input_msa_name = str(info[1])
            msa_type_name = info[2].split('.')[1].split('-')[0]

        except Exception as e:
            raise ValueError('Unable to fetch input_msa_name object from workspace: ' + str(e))
            #to get the full stack trace: traceback.format_exc()

        if msa_type_name != 'MSA':
            raise ValueError('Cannot yet handle input_msa type of: '+msa_type_name)
        else:
            MSA_in = input_msa_data
            row_order = []
            default_row_labels = dict()
            if 'row_order' in list(MSA_in.keys()):
                row_order = MSA_in['row_order']
            else:
                row_order = sorted(MSA_in['alignment'].keys())

            if 'default_row_labels' in list(MSA_in.keys()):
                default_row_labels = MSA_in['default_row_labels']
            else:
                for row_id in row_order:
                    default_row_labels[row_id] = row_id

            """
            # determine row index of query sequence
            for row_id in row_order:
                master_row_idx += 1
                if row_id == input_one_feature_id:
                    break
            if master_row_idx == 0:
                self.log(invalid_msgs,"Failed to find query id "+input_one_feature_id+" from Query Object "+input_one_name+" within MSA: "+input_msa_name)
            """

            # use longest sequence in MSA to use as the query sequence
            master_row_idx = -1
            longest_seq_len = 0
            longest_seq = ''
            for i,row_id in enumerate(row_order):
                msa_seq = MSA_in['alignment'][row_id].replace('-','')
                if len(msa_seq) > longest_seq_len:
                    master_row_idx = i
                    longest_seq_len = len(msa_seq)
                    longest_seq = msa_seq
            if longest_seq == '':
                raise ValueError ("unable to find longest seq in MSA")
            one_forward_reads_file_path = os.path.join(self.scratch, input_msa_name+"-query.fasta")
            with open (one_forward_reads_file_path, 'w') as input_one_fh:
                input_one_fh.write('>query'+"\n")
                input_one_fh.write(longest_seq+"\n")
            if master_row_idx == -1:
                raise ValueError ("unable to find longest sequence in MSA")
            else:
                master_row_idx += 1  # psiBLAST counts rows starting from 1 instead of 0

            
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
            with open(input_MSA_file_path, 'w') as input_MSA_file_handle:
                input_MSA_file_handle.write("\n".join(records)+"\n")


            # Determine whether nuc or protein sequences
            #
            self.log (console, "CHECKING MSA for PROTEIN seqs...")  # DEBUG                                  
            PROT_MSA_pattern = re.compile("^[\.\-_acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWYxX ]+$")
            DNA_MSA_pattern = re.compile("^[\.\-_ACGTUXNRYSWKMBDHVacgtuxnryswkmbdhv \t\n]+$")               
            appropriate_sequence_found_in_MSA_input = True
            if 'sequence_type' in MSA_in and (MSA_in['sequence_type'] == 'dna' or MSA_in['sequence_type'] == 'DNA'):
                appropriate_sequence_found_in_MSA_input = False
            else:
                for row_id in row_order:
                    #self.log(console, row_id+": '"+MSA_in['alignment'][row_id]+"'")    # DEBUG                   
                    #if not PROT_MSA_pattern.match(MSA_in['alignment'][row_id]):
                    if DNA_MSA_pattern.match(MSA_in['alignment'][row_id]):
                        self.log(invalid_msgs,
                                 "Need protein sequences in MSA. "+
                                 "BAD nucleotide record for MSA row_id: "+row_id+"\n"+MSA_in['alignment'][row_id]+"\n")
                        appropriate_sequence_found_in_MSA_input = False
                        break
                    elif not PROT_MSA_pattern.match(MSA_in['alignment'][row_id]):
                        self.log(invalid_msgs,"BAD record for MSA row_id: "+row_id+"\n"+MSA_in['alignment'][row_id]+"\n")
                        appropriate_sequence_found_in_MSA_input = False
                        break


        #### Get the input_many object
        ##
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            #objects = ws.get_objects([{'ref': input_many_ref}])
            objects = ws.get_objects2({'objects':[{'ref': input_many_ref}]})['data']
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
                print((traceback.format_exc()))
                raise ValueError('Unable to get SequenceSet: ' + str(e))

            header_id = input_many_sequenceSet['sequences'][0]['sequence_id']
            many_forward_reads_file_path = os.path.join(self.scratch, header_id+'.fasta')
            many_forward_reads_file_handle = open(many_forward_reads_file_path, 'w')
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
            if len(list(feature_ids_by_genome_ref.keys())) > 0:
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
            if len(list(feature_ids_by_genome_id.keys())) > 0:
                appropriate_sequence_found_in_many_input = True

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "FeatureSetToFasta() took "+str(end_time-beg_time)+" secs")


        # Missing proper input_many_type
        #
        else:
            raise ValueError('Cannot yet handle input_many type of: '+many_type_name)            


        # check for failed input file creation
        #
        if not appropriate_sequence_found_in_MSA_input or \
           not os.path.isfile(one_forward_reads_file_path) or \
           not os.path.getsize(one_forward_reads_file_path) > 0 or \
           not os.path.isfile(input_MSA_file_path) or \
           not os.path.getsize(input_MSA_file_path):
            self.log(invalid_msgs,"no protein sequences found in MSA'"+input_msa_ref+"'")
        if not appropriate_sequence_found_in_many_input or \
           not os.path.isfile(many_forward_reads_file_path) or \
           not os.path.getsize(many_forward_reads_file_path) > 0:
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
            #provenance[0]['input_ws_objects'].append(input_one_ref)
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
            line = p.stdout.readline().decode()
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

        # check for necessary files
        if not os.path.isfile(blast_bin):
            raise ValueError("no such file '"+blast_bin+"'")
        if not os.path.isfile(one_forward_reads_file_path):
            raise ValueError("no such file '"+one_forward_reads_file_path+"'")
        elif not os.path.getsize(one_forward_reads_file_path) > 0:
            raise ValueError("empty file '"+one_forward_reads_file_path+"'")
        if not os.path.isfile(input_MSA_file_path):
            raise ValueError("no such file '"+input_MSA_file_path+"'")
        elif not os.path.getsize(input_MSA_file_path) > 0:
            raise ValueError("empty file '"+input_MSA_file_path+"'")
        if not os.path.isfile(many_forward_reads_file_path):
            raise ValueError("no such file '"+many_forward_reads_file_path+"'")
        elif not os.path.getsize(many_forward_reads_file_path) > 0:
            raise ValueError("empty file '"+many_forward_reads_file_path+"'")

        # set the output path
        timestamp = int((datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()*1000)
        output_dir = os.path.join(self.scratch,'output.'+str(timestamp))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        output_aln_file_path = os.path.join(output_dir, 'alnout.txt');
        output_extra_file_path = os.path.join(output_dir, 'alnout_extra.txt');
        output_filtered_fasta_file_path = os.path.join(output_dir, 'output_filtered.faa');

        # this is command for extra output
        extra_output = False
        if 'output_extra_format' in params and params['output_extra_format'] != None and params['output_extra_format'] != '' and params['output_extra_format'] != 'none':
            extra_output = True

            blast_cmd = [blast_bin]
            blast_cmd.append('-query')
            blast_cmd.append(one_forward_reads_file_path)
            blast_cmd.append('-db')
            blast_cmd.append(many_forward_reads_file_path)
            blast_cmd.append('-out')
            blast_cmd.append(output_extra_file_path)
            #blast_cmd.append('-html')  # HTML is a flag so doesn't get an arg val
            blast_cmd.append('-outfmt')
            blast_cmd.append(str(params['output_extra_format']))
            blast_cmd.append('-evalue')
            blast_cmd.append(str(params['e_value']))

            # options (not allowed for format 0)
            #if 'maxaccepts' in params:
            #    if params['maxaccepts']:
            #        blast_cmd.append('-max_target_seqs')
            #        blast_cmd.append(str(params['maxaccepts']))

            # Run BLAST, capture output as it happens
            #
            self.log(console, 'RUNNING BLAST (FOR EXTRA OUTPUT):')
            self.log(console, '    '+' '.join(blast_cmd))
            #        report += "\n"+'running BLAST:'+"\n"
            #        report += '    '+' '.join(blast_cmd)+"\n"

            p = subprocess.Popen(blast_cmd, \
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
                raise ValueError('Error running BLAST, return code: '+str(p.returncode) + 
                '\n\n'+ '\n'.join(console))

            # upload BLAST output
            dfu = DFUClient(self.callbackURL)
            try:
                extra_upload_ret = dfu.file_to_shock({'file_path': output_extra_file_path,
# DEBUG
#                                                      'make_handle': 0,
#                                                      'pack': 'zip'})
                                                      'make_handle': 0})
            except:
                raise ValueError ('error loading output_extra file to shock')


        # this is command for basic search mode (with TAB TXT output)
        blast_cmd = [blast_bin]
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
            line = p.stdout.readline().decode()
            if not line: break
            self.log(console, line.replace('\n', ''))

        p.stdout.close()
        p.wait()
        self.log(console, 'return code: ' + str(p.returncode))
        if p.returncode != 0:
            raise ValueError('Error running BLAST, return code: '+str(p.returncode) + 
                '\n\n'+ '\n'.join(console))

        # upload BLAST output
        dfu = DFUClient(self.callbackURL)
        try:
            base_upload_ret = dfu.file_to_shock({'file_path': output_aln_file_path,
# DEBUG
#                                                 'make_handle': 0,
#                                                 'pack': 'zip'})
                                                 'make_handle': 0})
        except:
            raise ValueError ('error loading aln_out file to shock')


        # get query_len for filtering later
        #
        query_len = 0
        with open(one_forward_reads_file_path, 'r') as query_file_handle:
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
        accept_fids = dict()
        output_aln_file_handle = open (output_aln_file_path, 'r')
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

        filtering_fields = dict()
        for hit_seq_id in hit_order:
            hit_buf.append(high_bitscore_line[hit_seq_id])
            filtering_fields[hit_seq_id] = dict()

            #self.log(console,"HIT_SEQ_ID: '"+hit_seq_id+"'")
            filter = False
            #if 'ident_thresh' in params and float(params['ident_thresh']) > float(high_bit#score_ident[hit_seq_id]):
            #    filter = True
            #    filtering_fields[hit_seq_id]['ident_thresh'] = True
            if 'bitscore' in params and float(params['bitscore']) > float(high_bitscore_score[hit_seq_id]):
                filter = True
                filtering_fields[hit_seq_id]['bitscore'] = True
            if 'overlap_fraction' in params and float(params['overlap_fraction']) > 100*float(high_bitscore_alnlen[hit_seq_id])/float(query_len):
                filter = True
                filtering_fields[hit_seq_id]['overlap_fraction'] = True

            if filter:
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

                id_untrans = header_id
                id_trans = re.sub ('\|',':',id_untrans)  # BLAST seems to make this translation now when id format has simple 'kb|blah' format
                if id_trans in hit_seq_ids or id_untrans in hit_seq_ids:
                    #self.log(console, 'FOUND HIT '+header_id)  # DEBUG
                    accept_fids[id_untrans] = True
                    output_sequenceSet['sequences'].append(seq_obj)


        # FeatureSet input -> FeatureSet output
        #
        elif many_type_name == 'FeatureSet':
            seq_total = len(list(input_many_featureSet['elements'].keys()))

            output_featureSet = dict()
            if 'description' in input_many_featureSet and input_many_featureSet['description'] != None:
                output_featureSet['description'] = input_many_featureSet['description'] + " - "+search_tool_name+"_Search filtered"
            else:
                output_featureSet['description'] = search_tool_name+"_Search filtered"
            output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()

            fId_list = list(input_many_featureSet['elements'].keys())
            self.log(console,"ADDING FEATURES TO FEATURESET")
            for fId in sorted(fId_list):
                for genome_ref in input_many_featureSet['elements'][fId]:
                    id_untrans = genome_ref+genome_id_feature_id_delim+fId
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
        elif many_type_name == 'Genome':
            seq_total = 0
            output_featureSet = dict()
#            if 'scientific_name' in input_many_genome and input_many_genome['scientific_name'] != None:
#                output_featureSet['description'] = input_many_genome['scientific_name'] + " - "+search_tool_name+"_Search filtered"
#            else:
#                output_featureSet['description'] = search_tool_name+"_Search filtered"
            output_featureSet['description'] = search_tool_name+"_Search filtered"
            output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()
            for fid in feature_ids:
                seq_total += 1
                id_untrans = fid
                id_trans = re.sub ('\|',':',id_untrans)  # BLAST seems to make this translation now when id format has simple 'kb|blah' format
                if id_trans in hit_seq_ids or id_untrans in hit_seq_ids:
                    #self.log(console, 'FOUND HIT '+fid)  # DEBUG
                    #output_featureSet['element_ordering'].append(fid)
                    accept_fids[id_untrans] = True
                    #fid = input_many_ref+genome_id_feature_id_delim+id_untrans  # don't change fId for output FeatureSet
                    output_featureSet['element_ordering'].append(fid)
                    output_featureSet['elements'][fid] = [input_many_ref]

        # Parse GenomeSet hits into FeatureSet
        #
        elif many_type_name == 'GenomeSet':
            seq_total = 0

            output_featureSet = dict()
            if 'description' in input_many_genomeSet and input_many_genomeSet['description'] != None:
                output_featureSet['description'] = input_many_genomeSet['description'] + " - "+search_tool_name+"_Search filtered"
            else:
                output_featureSet['description'] = search_tool_name+"_Search filtered"
            output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()

            self.log(console,"READING HITS FOR GENOMES")  # DEBUG
            for genome_id in list(feature_ids_by_genome_id.keys()):
                self.log(console,"READING HITS FOR GENOME "+genome_id)  # DEBUG
                genome_ref = input_many_genomeSet['elements'][genome_id]['ref']
                for feature_id in feature_ids_by_genome_id[genome_id]:
                    seq_total += 1
                    id_untrans = genome_ref+genome_id_feature_id_delim+feature_id
                    id_trans = re.sub ('\|',':',id_untrans)  # BLAST seems to make this translation now when id format has simple 'kb|blah' format
                    if id_trans in hit_seq_ids or id_untrans in hit_seq_ids:
                        #self.log(console, 'FOUND HIT: '+feature['id'])  # DEBUG
                        #output_featureSet['element_ordering'].append(feature['id'])
                        accept_fids[id_untrans] = True
                        #feature_id = id_untrans  # don't change fId for output FeatureSet
                        try:
                            this_genome_ref_list = output_featureSet['elements'][feature_id]
                        except:
                            output_featureSet['elements'][feature_id] = []
                            output_featureSet['element_ordering'].append(feature_id)
                        output_featureSet['elements'][feature_id].append(genome_ref)


        # load the method provenance from the context object
        #
        self.log(console,"SETTING PROVENANCE")  # DEBUG
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects'] = []
        #provenance[0]['input_ws_objects'].append(input_one_ref)
        provenance[0]['input_ws_objects'].append(input_msa_ref)
        provenance[0]['input_ws_objects'].append(input_many_ref)
        provenance[0]['service'] = 'kb_blast'
        provenance[0]['method'] = search_tool_name+'_Search'


        # Upload results
        #
        if len(invalid_msgs) == 0 and len(list(hit_seq_ids.keys())) > 0:
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
        if len(invalid_msgs) == 0 and len(hit_order) > 0:

            # text report
            #
            report += 'sequences in search db: '+str(seq_total)+"\n"
            report += 'sequences in hit set: '+str(len(hit_order))+"\n"
            report += 'sequences in accepted hit set: '+str(hit_total)+"\n"
            report += "\n"
            for line in hit_buf:
                report += line
            self.log (console, report)


            # build html report
            if many_type_name == 'Genome':
                feature_id_to_function = GenomeToFASTA_retVal['feature_id_to_function']
                genome_ref_to_sci_name = GenomeToFASTA_retVal['genome_ref_to_sci_name']
            elif many_type_name == 'GenomeSet':
                feature_id_to_function = GenomeSetToFASTA_retVal['feature_id_to_function']
                genome_ref_to_sci_name = GenomeSetToFASTA_retVal['genome_ref_to_sci_name']
            elif many_type_name == 'FeatureSet':
                feature_id_to_function = FeatureSetToFASTA_retVal['feature_id_to_function']
                genome_ref_to_sci_name = FeatureSetToFASTA_retVal['genome_ref_to_sci_name']
                
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

            html_report_lines = []
            html_report_lines += ['<html>']
            html_report_lines += ['<body bgcolor="white">']
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

            for line in hit_buf:
                line = line.strip()
                if line == '' or line.startswith('#'):
                    continue

                [query_id, hit_id, identity, aln_len, mismatches, gap_openings, q_beg, q_end, h_beg, h_end, e_value, bit_score] = line.split("\t")[0:12]

                aln_len_perc = round (100.0*float(aln_len)/float(query_len), 1)
                identity = str(round(float(identity), 1))
                if identity == '100.0':  identity = '100'

                #if many_type_name == 'SingleEndLibrary':
                #    pass
                #elif many_type_name == 'SequenceSet':
                if many_type_name == 'SequenceSet':
                    pass
                elif many_type_name == 'Genome' or \
                        many_type_name == 'GenomeSet' or \
                        many_type_name == 'FeatureSet':

                    if many_type_name != 'Genome':
                        [genome_ref, hit_fid] = hit_id.split(genome_id_feature_id_delim)
                    else:
                        genome_ref = input_many_ref
                        hit_fid = hit_id

                    # can't just use hit_fid because may have pipes translated and can't translate back
                    fid_lookup = None
                    for fid in list(feature_id_to_function[genome_ref].keys()):
                        id_untrans = fid
                        id_trans = re.sub ('\|',':',id_untrans)  # BLAST seems to make this translation now when id format has simple 'kb|blah' format

                        #self.log (console, "SCANNING FIDS.  HIT_FID: '"+str(hit_fid)+"' FID: '"+str(fid)+"' TRANS: '"+str(id_trans)+"'")  # DEBUG

                        if id_untrans == hit_fid or id_trans == hit_fid:
                            #self.log (console, "GOT ONE!")  # DEBUG
                            if many_type_name == 'Genome':
                                accept_id = fid
                            elif many_type_name == 'GenomeSet' or many_type_name == 'FeatureSet':
                                accept_id = genome_ref+genome_id_feature_id_delim+fid
                            if accept_id in accept_fids:
                                row_color = accept_row_color
                            else:
                                row_color = reject_row_color
                            fid_lookup = fid
                            break
                    #self.log (console, "HIT_FID: '"+str(hit_fid)+"' FID_LOOKUP: '"+str(fid_lookup)+"'")  # DEBUG
                    if fid_lookup == None:
                        raise ValueError ("unable to find fid for hit_fid: '"+str(hit_fid))
                    elif fid_lookup not in feature_id_to_function[genome_ref]:
                        raise ValueError ("unable to find function for fid: '"+str(fid_lookup))
                    fid_disp = re.sub (r"^.*\.([^\.]+)\.([^\.]+)$", r"\1.\2", fid_lookup)

                    func_disp = feature_id_to_function[genome_ref][fid_lookup]
                    genome_sci_name = genome_ref_to_sci_name[genome_ref]

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
                    html_report_lines += ['<td style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+genome_sci_name+'</font></td>']
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
                    # bit score
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
            html_file = search_tool_name+'_Search.html'
            html_path = os.path.join (output_dir, html_file)
            with open (html_path, 'w') as html_handle:
                html_handle.write(html_report_str)

            dfu = DFUClient(self.callbackURL)
            try:
                upload_ret = dfu.file_to_shock({'file_path': html_path,
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
            html_buf_lim = 16000  # really 16KB, but whatever
            if len(html_report_str) <= html_buf_lim:
                reportObj['direct_html'] = html_report_str
            else:
                reportObj['direct_html_link_index'] = 0

            reportObj['html_links'] = [{'shock_id': upload_ret['shock_id'],
                                        'name': html_file,
                                        'label': search_tool_name+' Results'}
                                       ]
            reportObj['file_links'] = [{'shock_id': base_upload_ret['shock_id'],
                                        'name': search_tool_name+'_Search-m'+'7'+'.txt',
                                        'label': search_tool_name+' Results: m'+'7'}
                                       ]
            if extra_output:
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
                reportObj['file_links'].append({'shock_id': extra_upload_ret['shock_id'],
                                                'name': search_tool_name+'_Search-m'+str(params['output_extra_format'])+'.'+extension,
                                                'label': search_tool_name+' Results: m'+str(params['output_extra_format'])})
                            
            if hit_total > 0:
                reportObj['objects_created'].append({'ref':str(params['workspace_name'])+'/'+params['output_filtered_name'],'description':search_tool_name+' hits'})
            #reportObj['message'] = report


            # save report object
            #
            SERVICE_VER = 'release'
            reportClient = KBaseReport(self.callbackURL, token=ctx['token'], service_ver=SERVICE_VER)
            #report_info = report.create({'report':reportObj, 'workspace_name':params['workspace_name']})
            report_info = reportClient.create_extended_report(reportObj)

        else:
            if len(hit_order) == 0:  # no hits
                report += "No hits were found\n"
            else:  # data validation error
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
            report_info = dict()
            report_info['name'] = report_obj_info[1]
            report_info['ref'] = str(report_obj_info[6])+'/'+str(report_obj_info[0])+'/'+str(report_obj_info[4])

        self.log(console,"BUILDING RETURN OBJECT")
#        returnVal = { 'output_report_name': reportName,
#                      'output_report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
#                      'output_filtered_ref': params['workspace_name']+'/'+params['output_filtered_name']
#                      }
        returnVal = { 'report_name': report_info['name'],
                      'report_ref': report_info['ref']
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
