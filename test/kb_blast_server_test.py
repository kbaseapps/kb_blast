# -*- coding: utf-8 -*-
import json  # noqa: F401
import os  # noqa: F401
import time
import shutil
from pprint import pprint
import unittest
from configparser import ConfigParser  # py3
from os import environ

from installed_clients.WorkspaceClient import Workspace as workspaceService
from installed_clients.GenomeFileUtilClient import GenomeFileUtil
from kb_blast.authclient import KBaseAuth as _KBaseAuth
from kb_blast.kb_blastImpl import kb_blast
from kb_blast.kb_blastServer import MethodContext


class kb_blastTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = environ.get('KB_AUTH_TOKEN', None)
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('kb_blast'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'kb_blast',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = workspaceService(cls.wsURL)
        cls.serviceImpl = kb_blast(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_kb_blast_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})  # noqa
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    # get obj_ref in form D/D/D
    def get_obj_ref_from_obj_info (self, obj_info):
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = list(range(11))  # object_info tuple

        return '/'.join([str(obj_info[WSID_I]), str(obj_info[OBJID_I]), str(obj_info[VERSION_I])])
    

    # call this method to get the WS object info of a Genome
    #   (will upload the example data if this is the first time the method is called during tests)
    def getGenomeInfo(self, genome_basename, item_i=0):
        if hasattr(self.__class__, 'genomeInfo_list'):
            try:
                info = self.__class__.genomeInfo_list[item_i]
                name = self.__class__.genomeName_list[item_i]
                if info != None:
                    if name != genome_basename:
                        self.__class__.genomeInfo_list[item_i] = None
                        self.__class__.genomeName_list[item_i] = None
                    else:
                        return info
            except:
                pass

        # 1) transform genbank to kbase genome object and upload to ws
        shared_dir = "/kb/module/work/tmp"
        genome_data_file = 'data/genomes/'+genome_basename+'.gbff.gz'
        genome_file = os.path.join(shared_dir, os.path.basename(genome_data_file))
        shutil.copy(genome_data_file, genome_file)

        SERVICE_VER = 'release'
        #SERVICE_VER = 'dev'
        GFU = GenomeFileUtil(os.environ['SDK_CALLBACK_URL'],
                             token=self.getContext()['token'],
                             service_ver=SERVICE_VER
                         )
        print ("UPLOADING genome: "+genome_basename+" to WORKSPACE "+self.getWsName()+" ...")
        genome_upload_result = GFU.genbank_to_genome({'file': {'path': genome_file },
                                                      'workspace_name': self.getWsName(),
                                                      'genome_name': genome_basename
                                                  })
#                                                  })[0]
        pprint(genome_upload_result)
        genome_ref = genome_upload_result['genome_ref']
        new_obj_info = self.getWsClient().get_object_info_new({'objects': [{'ref': genome_ref}]})[0]

        # 2) store it
        if not hasattr(self.__class__, 'genomeInfo_list'):
            self.__class__.genomeInfo_list = []
            self.__class__.genomeName_list = []
        for i in range(item_i+1):
            try:
                assigned = self.__class__.genomeInfo_list[i]
            except:
                self.__class__.genomeInfo_list.append(None)
                self.__class__.genomeName_list.append(None)

        self.__class__.genomeInfo_list[item_i] = new_obj_info
        self.__class__.genomeName_list[item_i] = genome_basename
        return new_obj_info

    # call this method to get the WS object info of an AnnotatedMetagenomeAssembly
    #   (will upload the example data if this is the first time the method is called during tests)
    def getAMAInfo(self, ama_basename, item_i=0):
        if hasattr(self.__class__, 'amaInfo_list'):
            try:
                info = self.__class__.amaInfo_list[item_i]
                name = self.__class__.amaName_list[item_i]
                if info != None:
                    if name != ama_basename:
                        self.__class__.amaInfo_list[item_i] = None
                        self.__class__.amaName_list[item_i] = None
                    else:
                        return info
            except:
                pass

        # 1) transform GFF+FNA to kbase AMA object and upload to ws
        shared_dir = "/kb/module/work/tmp"
        ama_gff_srcfile = 'data/amas/'+ama_basename+'.gff'
        ama_fna_srcfile = 'data/amas/'+ama_basename+'.fa'
        ama_gff_dstfile = os.path.join(shared_dir, os.path.basename(ama_gff_srcfile))
        ama_fna_dstfile = os.path.join(shared_dir, os.path.basename(ama_fna_srcfile))
        shutil.copy(ama_gff_srcfile, ama_gff_dstfile)
        shutil.copy(ama_fna_srcfile, ama_fna_dstfile)

        try:
            SERVICE_VER = 'release'
            #SERVICE_VER = 'dev'
            GFU = GenomeFileUtil(os.environ['SDK_CALLBACK_URL'],
                                 token=self.getContext()['token'],
                                 service_ver=SERVICE_VER
            )
        except:
            raise ValueError ("unable to obtain GenomeFileUtil client")
        print ("UPLOADING AMA: "+ama_basename+" to WORKSPACE "+self.getWsName()+" ...")
        ama_upload_params = {
            "workspace_name": self.getWsName(),
            "genome_name": ama_basename,
            "fasta_file": {"path": ama_fna_dstfile},
            "gff_file": {"path": ama_gff_dstfile},
            "source": "GFF",
            "scientific_name": "TEST AMA",
            "generate_missing_genes": "True"
        }        
        try:
            ama_upload_result = GFU.fasta_gff_to_metagenome(ama_upload_params)
        except:
            raise ValueError("unable to upload test AMA data object")
        print ("AMA UPLOADED")
        pprint(ama_upload_result)

        ama_ref = ama_upload_result['metagenome_ref']
        new_obj_info = self.getWsClient().get_object_info_new({'objects': [{'ref': ama_ref}]})[0]

        # 2) store it
        if not hasattr(self.__class__, 'amaInfo_list'):
            self.__class__.amaInfo_list = []
            self.__class__.amaName_list = []
        for i in range(item_i+1):
            try:
                assigned = self.__class__.amaInfo_list[i]
            except:
                self.__class__.amaInfo_list.append(None)
                self.__class__.amaName_list.append(None)

        self.__class__.amaInfo_list[item_i] = new_obj_info
        self.__class__.amaName_list[item_i] = ama_basename
        return new_obj_info

    #
    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    #


    # Test BLASTn: Single Genome target
    #
    # Uncomment to skip this test
    # HIDE @unittest.skip("skipped test_kb_blast_BLASTn_Search_01")
    def test_kb_blast_BLASTn_Search_01(self):
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = list(range(11))  # object_info tuple

        obj_basename = 'BLASTn'
        obj_out_name = obj_basename+".test_output.FS"
        obj_out_type = "KBaseCollections.FeatureSet"
        expected_hit_cnt = 1
        
        #reference_prok_genomes_WS = 'ReferenceDataManager'  # PROD and CI
        #genome_ref_1 = 'ReferenceDataManager/GCF_001566335.1/1'  # E. coli K-12 MG1655
        genomeInfo_0 = self.getGenomeInfo('GCF_001566335.1_ASM156633v1_genomic', 0)  # E. coli K-12 MG1655
        genome_ref_1 = self.get_obj_ref_from_obj_info(genomeInfo_0)

        # E. coli K-12 MG1655 dnaA
        query_seq_nuc = 'GTGTCACTTTCGCTTTGGCAGCAGTGTCTTGCCCGATTGCAGGATGAGTTACCAGCCACAGAATTCAGTATGTGGATACGCCCATTGCAGGCGGAACTGAGCGATAACACGCTGGCCCTGTACGCGCCAAACCGTTTTGTCCTCGATTGGGTACGGGACAAGTACCTTAATAATATCAATGGACTGCTAACCAGTTTCTGCGGAGCGGATGCCCCACAGCTGCGTTTTGAAGTCGGCACCAAACCGGTGACGCAAACGCCACAAGCGGCAGTGACGAGCAACGTCGCGGCCCCTGCACAGGTGGCGCAAACGCAGCCGCAACGTGCTGCGCCTTCTACGCGCTCAGGTTGGGATAACGTCCCGGCCCCGGCAGAACCGACCTATCGTTCTAACGTAAACGTCAAACACACGTTTGATAACTTCGTTGAAGGTAAATCTAACCAACTGGCGCGCGCGGCGGCTCGCCAGGTGGCGGATAACCCTGGCGGTGCCTATAACCCGTTGTTCCTTTATGGCGGCACGGGTCTGGGTAAAACTCACCTGCTGCATGCGGTGGGTAACGGCATTATGGCGCGCAAGCCGAATGCCAAAGTGGTTTATATGCACTCCGAGCGCTTTGTTCAGGACATGGTTAAAGCCCTGCAAAACAACGCGATCGAAGAGTTTAAACGCTACTACCGTTCCGTAGATGCACTGCTGATCGACGATATTCAGTTTTTTGCTAATAAAGAACGATCTCAGGAAGAGTTTTTCCACACCTTCAACGCCCTGCTGGAAGGTAATCAACAGATCATTCTCACCTCGGATCGCTATCCGAAAGAGATCAACGGCGTTGAGGATCGTTTGAAATCCCGCTTCGGTTGGGGACTGACTGTGGCGATCGAACCGCCAGAGCTGGAAACCCGTGTGGCGATCCTGATGAAAAAGGCCGACGAAAACGACATTCGTTTGCCGGGCGAAGTGGCGTTCTTTATCGCCAAGCGTCTACGATCTAACGTACGTGAGCTGGAAGGGGCGCTGAACCGCGTCATTGCCAATGCCAACTTTACCGGACGGGCGATCACCATCGACTTCGTGCGTGAGGCGCTGCGCGACTTGCTGGCATTGCAGGAAAAACTGGTCACCATCGACAATATTCAGAAGACGGTGGCGGAGTACTACAAGATCAAAGTCGCGGATCTCCTTTCCAAGCGTCGATCCCGCTCGGTGGCGCGTCCGCGCCAGATGGCGATGGCGCTGGCGAAAGAGCTGACTAACCACAGTCTGCCGGAGATTGGCGATGCGTTTGGTGGCCGTGACCACACGACGGTGCTTCATGCCTGCCGTAAGATCGAGCAGTTGCGTGAAGAGAGCCACGATATCAAAGAAGATTTTTCAAATTTAATCAGAACATTGTCATCGTAA'

        parameters = { 'workspace_name': self.getWsName(),
                       'input_one_sequence': query_seq_nuc,
                       #'input_one_ref': "",
                       'output_one_name': obj_basename+'.'+"test_query.SS",
                       'input_many_refs': [genome_ref_1],
                       'output_filtered_name': obj_out_name,
                       'genome_disp_name_config': 'obj_name',
                       'e_value': ".001",
                       'bitscore': "50",
                       'ident_thresh': "97.0",
                       'overlap_fraction': "50.0",
                       'maxaccepts': "1000",
                       'output_extra_format': "none"
                     }

        ret = self.getImpl().BLASTn_Search(self.getContext(), parameters)[0]
        self.assertIsNotNone(ret['report_ref'])

        # check created obj
        #report_obj = self.getWsClient().get_objects2({'objects':[{'ref':ret['report_ref']}]})[0]['data']
        report_obj = self.getWsClient().get_objects([{'ref':ret['report_ref']}])[0]['data']
        self.assertIsNotNone(report_obj['objects_created'][0]['ref'])

        created_obj_0_info = self.getWsClient().get_object_info_new({'objects':[{'ref':report_obj['objects_created'][0]['ref']}]})[0]
        self.assertEqual(created_obj_0_info[NAME_I], obj_out_name)
        self.assertEqual(created_obj_0_info[TYPE_I].split('-')[0], obj_out_type)

        # check number of hits in featureSet output
        featureSet_out_obj = self.getWsClient().get_objects([{'ref':report_obj['objects_created'][0]['ref']}])[0]['data']
        self.assertEqual(expected_hit_cnt, len(featureSet_out_obj['element_ordering']))
        pass


    # Test BLASTn: GenomeSet target
    #
    # Uncomment to skip this test
    # HIDE @unittest.skip("skipped test_kb_blast_BLASTn_Search_01")
    def test_kb_blast_BLASTn_Search_01(self):
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = list(range(11))  # object_info tuple

        obj_basename = 'BLASTn'
        obj_out_name = obj_basename+".test_output.FS"
        obj_out_type = "KBaseCollections.FeatureSet"
        genomeSet_name = 'test_genomeSet.BLASTn.GenomeSet'
        expected_hit_cnt = 3
        
        genomeInfo_0 = self.getGenomeInfo('GCF_001566335.1_ASM156633v1_genomic', 0)  # E. coli K-12 MG1655
        genomeInfo_1 = self.getGenomeInfo('GCF_000021385.1_ASM2138v1_genomic', 1)    # D. vulgaris str. 'Miyazaki F
        genomeInfo_2 = self.getGenomeInfo('GCF_001721825.1_ASM172182v1_genomic', 2)  # Pseudomonas aeruginosa
        genomeInfo_3 = self.getGenomeInfo('GCF_002950035.1_ASM295003v1_genomic', 3)  # Shigella boydii
        genomeInfo_4 = self.getGenomeInfo('GCF_000512125.1_ASM51212v1_genomic', 4)  # Escherichia albertii KF1
        genome_ref_1 = self.get_obj_ref_from_obj_info(genomeInfo_0)
        genome_ref_2 = self.get_obj_ref_from_obj_info(genomeInfo_1)
        genome_ref_3 = self.get_obj_ref_from_obj_info(genomeInfo_2)
        genome_ref_4 = self.get_obj_ref_from_obj_info(genomeInfo_3)
        genome_ref_5 = self.get_obj_ref_from_obj_info(genomeInfo_4)

        genome_ref_list = [genome_ref_1, genome_ref_2, genome_ref_3, genome_ref_4, genome_ref_5]
        genome_scinames = ['FOO', 'BAR', 'FOOBAR', 'BLAH', 'BLAHBLAH']

        # create GenomeSet
        testGS = {
            'description': 'five genomes',
            'elements': dict()
        }
        for genome_i, genome_ref in enumerate(genome_ref_list): 
            testGS['elements'][genome_scinames[genome_i]] = { 'ref': genome_ref }

        obj_info = self.getWsClient().save_objects({'workspace': self.getWsName(),       
                                                    'objects': [
                                                        {
                                                            'type':'KBaseSearch.GenomeSet',
                                                            'data':testGS,
                                                            'name':genomeSet_name,
                                                            'meta':{},
                                                            'provenance':[
                                                                {
                                                                    'service':'kb_blast',
                                                                    'method':'BLASTn_Search'
                                                                }
                                                            ]
                                                        }]
                                                })[0]

        #pprint(obj_info)
        target_genomeSet_ref = "/".join([str(obj_info[WORKSPACE_I]),
                                         str(obj_info[OBJID_I]),
                                         str(obj_info[VERSION_I])])


        # E. coli K-12 MG1655 dnaA
        query_seq_nuc = 'GTGTCACTTTCGCTTTGGCAGCAGTGTCTTGCCCGATTGCAGGATGAGTTACCAGCCACAGAATTCAGTATGTGGATACGCCCATTGCAGGCGGAACTGAGCGATAACACGCTGGCCCTGTACGCGCCAAACCGTTTTGTCCTCGATTGGGTACGGGACAAGTACCTTAATAATATCAATGGACTGCTAACCAGTTTCTGCGGAGCGGATGCCCCACAGCTGCGTTTTGAAGTCGGCACCAAACCGGTGACGCAAACGCCACAAGCGGCAGTGACGAGCAACGTCGCGGCCCCTGCACAGGTGGCGCAAACGCAGCCGCAACGTGCTGCGCCTTCTACGCGCTCAGGTTGGGATAACGTCCCGGCCCCGGCAGAACCGACCTATCGTTCTAACGTAAACGTCAAACACACGTTTGATAACTTCGTTGAAGGTAAATCTAACCAACTGGCGCGCGCGGCGGCTCGCCAGGTGGCGGATAACCCTGGCGGTGCCTATAACCCGTTGTTCCTTTATGGCGGCACGGGTCTGGGTAAAACTCACCTGCTGCATGCGGTGGGTAACGGCATTATGGCGCGCAAGCCGAATGCCAAAGTGGTTTATATGCACTCCGAGCGCTTTGTTCAGGACATGGTTAAAGCCCTGCAAAACAACGCGATCGAAGAGTTTAAACGCTACTACCGTTCCGTAGATGCACTGCTGATCGACGATATTCAGTTTTTTGCTAATAAAGAACGATCTCAGGAAGAGTTTTTCCACACCTTCAACGCCCTGCTGGAAGGTAATCAACAGATCATTCTCACCTCGGATCGCTATCCGAAAGAGATCAACGGCGTTGAGGATCGTTTGAAATCCCGCTTCGGTTGGGGACTGACTGTGGCGATCGAACCGCCAGAGCTGGAAACCCGTGTGGCGATCCTGATGAAAAAGGCCGACGAAAACGACATTCGTTTGCCGGGCGAAGTGGCGTTCTTTATCGCCAAGCGTCTACGATCTAACGTACGTGAGCTGGAAGGGGCGCTGAACCGCGTCATTGCCAATGCCAACTTTACCGGACGGGCGATCACCATCGACTTCGTGCGTGAGGCGCTGCGCGACTTGCTGGCATTGCAGGAAAAACTGGTCACCATCGACAATATTCAGAAGACGGTGGCGGAGTACTACAAGATCAAAGTCGCGGATCTCCTTTCCAAGCGTCGATCCCGCTCGGTGGCGCGTCCGCGCCAGATGGCGATGGCGCTGGCGAAAGAGCTGACTAACCACAGTCTGCCGGAGATTGGCGATGCGTTTGGTGGCCGTGACCACACGACGGTGCTTCATGCCTGCCGTAAGATCGAGCAGTTGCGTGAAGAGAGCCACGATATCAAAGAAGATTTTTCAAATTTAATCAGAACATTGTCATCGTAA'

        parameters = { 'workspace_name': self.getWsName(),
                       'input_one_sequence': query_seq_nuc,
                       #'input_one_ref': "",
                       'output_one_name': obj_basename+'.'+"test_query.SS",
                       'input_many_refs': [target_genomeSet_ref],
                       'output_filtered_name': obj_out_name,
                       'genome_disp_name_config': 'obj_name',
                       #'e_value': ".001",
                       #'bitscore': "50",
                       #'ident_thresh': "10.0",
                       #'overlap_fraction': "50.0",
                       'e_value': ".1",
                       'bitscore': "10",
                       'ident_thresh': "10.0",
                       'overlap_fraction': "10.0",
                       'maxaccepts': "1000",
                       'output_extra_format': "none"
                     }

        ret = self.getImpl().BLASTn_Search(self.getContext(), parameters)[0]
        self.assertIsNotNone(ret['report_ref'])

        # check created obj
        #report_obj = self.getWsClient().get_objects2({'objects':[{'ref':ret['report_ref']}]})[0]['data']
        report_obj = self.getWsClient().get_objects([{'ref':ret['report_ref']}])[0]['data']
        self.assertIsNotNone(report_obj['objects_created'][0]['ref'])

        created_obj_0_info = self.getWsClient().get_object_info_new({'objects':[{'ref':report_obj['objects_created'][0]['ref']}]})[0]
        self.assertEqual(created_obj_0_info[NAME_I], obj_out_name)
        self.assertEqual(created_obj_0_info[TYPE_I].split('-')[0], obj_out_type)

        # check number of hits in featureSet output
        featureSet_out_obj = self.getWsClient().get_objects([{'ref':report_obj['objects_created'][0]['ref']}])[0]['data']
        self.assertEqual(expected_hit_cnt, len(featureSet_out_obj['element_ordering']))
        pass


    # Test BLASTp: Single Genome target
    #
    # Uncomment to skip this test
    # HIDE @unittest.skip("skipped test_kb_blast_BLASTp_Search_01_Genome")
    def test_kb_blast_BLASTp_Search_01_Genome(self):
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = list(range(11))  # object_info tuple

        obj_basename = 'BLASTp_Genome'
        obj_out_name = obj_basename+".test_output.FS"
        obj_out_type = "KBaseCollections.FeatureSet"
        expected_hit_cnt = 1
        
        #reference_prok_genomes_WS = 'ReferenceDataManager'  # PROD and CI
        #genome_ref_1 = 'ReferenceDataManager/GCF_001566335.1/1'  # E. coli K-12 MG1655
        genomeInfo_0 = self.getGenomeInfo('GCF_001566335.1_ASM156633v1_genomic', 0)  # E. coli K-12 MG1655
        genome_ref_1 = self.get_obj_ref_from_obj_info(genomeInfo_0)

        # E. coli K-12 MG1655 dnaA
        query_seq_prot = 'MSLSLWQQCLARLQDELPATEFSMWIRPLQAELSDNTLALYAPNRFVLDWVRDKYLNNINGLLTSFCGADAPQLRFEVGTKPVTQTPQAAVTSNVAAPAQVAQTQPQRAAPSTRSGWDNVPAPAEPTYRSNVNVKHTFDNFVEGKSNQLARAAARQVADNPGGAYNPLFLYGGTGLGKTHLLHAVGNGIMARKPNAKVVYMHSERFVQDMVKALQNNAIEEFKRYYRSVDALLIDDIQFFANKERSQEEFFHTFNALLEGNQQIILTSDRYPKEINGVEDRLKSRFGWGLTVAIEPPELETRVAILMKKADENDIRLPGEVAFFIAKRLRSNVRELEGALNRVIANANFTGRAITIDFVREALRDLLALQEKLVTIDNIQKTVAEYYKIKVADLLSKRRSRSVARPRQMAMALAKELTNHSLPEIGDAFGGRDHTTVLHACRKIEQLREESHDIKEDFSNLIRTLSS'
        
        parameters = { 'workspace_name': self.getWsName(),
                       'input_one_sequence': query_seq_prot,
                       #'input_one_ref': "",
                       'output_one_name': obj_basename+'.'+"test_query.SS",
                       'input_many_refs': [genome_ref_1],
                       'output_filtered_name': obj_out_name,
                       'genome_disp_name_config': 'sci_name',
                       'e_value': ".001",
                       'bitscore': "50",
                       'ident_thresh': "40.0",
                       'overlap_fraction': "50.0",
                       'maxaccepts': "1000",
                       'output_extra_format': "none"
                     }

        ret = self.getImpl().BLASTp_Search(self.getContext(), parameters)[0]
        self.assertIsNotNone(ret['report_ref'])

        # check created obj
        #report_obj = self.getWsClient().get_objects2({'objects':[{'ref':ret['report_ref']}]})[0]['data']
        report_obj = self.getWsClient().get_objects([{'ref':ret['report_ref']}])[0]['data']
        self.assertIsNotNone(report_obj['objects_created'][0]['ref'])

        created_obj_0_info = self.getWsClient().get_object_info_new({'objects':[{'ref':report_obj['objects_created'][0]['ref']}]})[0]
        self.assertEqual(created_obj_0_info[NAME_I], obj_out_name)
        self.assertEqual(created_obj_0_info[TYPE_I].split('-')[0], obj_out_type)

        # check number of hits in featureSet output
        featureSet_out_obj = self.getWsClient().get_objects([{'ref':report_obj['objects_created'][0]['ref']}])[0]['data']
        self.assertEqual(expected_hit_cnt, len(featureSet_out_obj['element_ordering']))
        pass


    # Test BLASTp: GenomeSet target
    #
    # Uncomment to skip this test
    # HIDE @unittest.skip("skipped test_kb_blast_BLASTp_Search_02_GenomeSet")
    def test_kb_blast_BLASTp_Search_02_GenomeSet(self):
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple

        obj_basename = 'BLASTp_GenomeSet'
        obj_out_name = obj_basename+".test_output.FS"
        obj_out_type = "KBaseCollections.FeatureSet"
        genomeSet_name = 'test_genomeSet.BLASTp.GenomeSet'
        expected_hit_cnt = 2

        #reference_prok_genomes_WS = 'ReferenceDataManager'  # PROD and CI
        #genome_ref_1 = 'ReferenceDataManager/GCF_001566335.1/1'  # E. coli K-12 MG1655
        #genome_ref_2 = 'ReferenceDataManager/GCF_002936495.2/1'  # E. coli
        #genome_ref_3 = 'ReferenceDataManager/GCF_002936145.2/1'  # E. coli

        genomeInfo_0 = self.getGenomeInfo('GCF_001566335.1_ASM156633v1_genomic', 0)  # E. coli K-12 MG1655
        genomeInfo_1 = self.getGenomeInfo('GCF_000021385.1_ASM2138v1_genomic', 1)    # D. vulgaris str. 'Miyazaki F
        genomeInfo_2 = self.getGenomeInfo('GCF_001721825.1_ASM172182v1_genomic', 2)  # Pseudomonas aeruginosa
        genome_ref_1 = self.get_obj_ref_from_obj_info(genomeInfo_0)
        genome_ref_2 = self.get_obj_ref_from_obj_info(genomeInfo_1)
        genome_ref_3 = self.get_obj_ref_from_obj_info(genomeInfo_2)

        genome_ref_list = [genome_ref_1, genome_ref_2, genome_ref_3]
        genome_scinames = ['FOO', 'BAR', 'FOOBAR']

        # create GenomeSet
        testGS = {
            'description': 'three genomes',
            'elements': dict()
        }
        for genome_i, genome_ref in enumerate(genome_ref_list): 
            testGS['elements'][genome_scinames[genome_i]] = { 'ref': genome_ref }

        obj_info = self.getWsClient().save_objects({'workspace': self.getWsName(),       
                                                    'objects': [
                                                        {
                                                            'type':'KBaseSearch.GenomeSet',
                                                            'data':testGS,
                                                            'name':genomeSet_name,
                                                            'meta':{},
                                                            'provenance':[
                                                                {
                                                                    'service':'kb_blast',
                                                                    'method':'BLASTp_Search'
                                                                }
                                                            ]
                                                        }]
                                                })[0]

        #pprint(obj_info)
        target_genomeSet_ref = "/".join([str(obj_info[WORKSPACE_I]),
                                         str(obj_info[OBJID_I]),
                                         str(obj_info[VERSION_I])])


        # E. coli K-12 MG1655 dnaA
        query_seq_prot = 'MSLSLWQQCLARLQDELPATEFSMWIRPLQAELSDNTLALYAPNRFVLDWVRDKYLNNINGLLTSFCGADAPQLRFEVGTKPVTQTPQAAVTSNVAAPAQVAQTQPQRAAPSTRSGWDNVPAPAEPTYRSNVNVKHTFDNFVEGKSNQLARAAARQVADNPGGAYNPLFLYGGTGLGKTHLLHAVGNGIMARKPNAKVVYMHSERFVQDMVKALQNNAIEEFKRYYRSVDALLIDDIQFFANKERSQEEFFHTFNALLEGNQQIILTSDRYPKEINGVEDRLKSRFGWGLTVAIEPPELETRVAILMKKADENDIRLPGEVAFFIAKRLRSNVRELEGALNRVIANANFTGRAITIDFVREALRDLLALQEKLVTIDNIQKTVAEYYKIKVADLLSKRRSRSVARPRQMAMALAKELTNHSLPEIGDAFGGRDHTTVLHACRKIEQLREESHDIKEDFSNLIRTLSS'
        
        parameters = { 'workspace_name': self.getWsName(),
                       'input_one_sequence': query_seq_prot,
                       #'input_one_ref': "",
                       'output_one_name': obj_basename+'.'+"test_query.SS",
                       'input_many_refs': [target_genomeSet_ref],
                       'output_filtered_name': obj_out_name,
                       'genome_disp_name_config': 'obj_name_sci_name',
                       'e_value': ".001",
                       'bitscore': "50",
                       'ident_thresh': "40.0",
                       'overlap_fraction': "50.0",
                       'maxaccepts': "1000",
                       'output_extra_format': "none"
                     }

        ret = self.getImpl().BLASTp_Search(self.getContext(), parameters)[0]
        self.assertIsNotNone(ret['report_ref'])

        # check created obj
        #report_obj = self.getWsClient().get_objects2({'objects':[{'ref':ret['report_ref']}]})[0]['data']
        report_obj = self.getWsClient().get_objects([{'ref':ret['report_ref']}])[0]['data']
        self.assertIsNotNone(report_obj['objects_created'][0]['ref'])

        created_obj_0_info = self.getWsClient().get_object_info_new({'objects':[{'ref':report_obj['objects_created'][0]['ref']}]})[0]
        self.assertEqual(created_obj_0_info[NAME_I], obj_out_name)
        self.assertEqual(created_obj_0_info[TYPE_I].split('-')[0], obj_out_type)

        # check number of hits in featureSet output
        featureSet_out_obj = self.getWsClient().get_objects([{'ref':report_obj['objects_created'][0]['ref']}])[0]['data']
        self.assertEqual(expected_hit_cnt, len(featureSet_out_obj['element_ordering']))
        pass


    # Test BLASTp: FeatureSet
    #
    # Uncomment to skip this test
    # HIDE @unittest.skip("skipped test_kb_blast_BLASTp_Search_03_FeatureSet")
    def test_kb_blast_BLASTp_Search_03_FeatureSet(self):
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple

        obj_basename = 'BLASTp_FeatureSet'
        obj_out_name = obj_basename+".test_output.FS"
        obj_out_type = "KBaseCollections.FeatureSet"
        target_1 = obj_basename+'.test_FeatureSet'
        expected_hit_cnt = 2
        
        #reference_prok_genomes_WS = 'ReferenceDataManager'  # PROD and CI
        #genome_ref_1 = 'ReferenceDataManager/GCF_001566335.1/1'  # E. coli K-12 MG1655
        #genome_ref_2 = 'ReferenceDataManager/GCF_002936495.2/1'  # E. coli
        #genome_ref_3 = 'ReferenceDataManager/GCF_002936145.2/1'  # E. coli
        genomeInfo_0 = self.getGenomeInfo('GCF_001566335.1_ASM156633v1_genomic', 0)  # E. coli K-12 MG1655
        genomeInfo_1 = self.getGenomeInfo('GCF_000021385.1_ASM2138v1_genomic', 1)    # D. vulgaris str. 'Miyazaki F
        genomeInfo_2 = self.getGenomeInfo('GCF_001721825.1_ASM172182v1_genomic', 2)  # Pseudomonas aeruginosa
        genome_ref_1 = self.get_obj_ref_from_obj_info(genomeInfo_0)
        genome_ref_2 = self.get_obj_ref_from_obj_info(genomeInfo_1)
        genome_ref_3 = self.get_obj_ref_from_obj_info(genomeInfo_2)

        # build FeatureSet obj
        feature_id_1_0 = 'AWN69_RS07145'  # dnaA 
        feature_id_1_1 = 'AWN69_RS00105'
        feature_id_2_0 = 'DVMF_RS00005'  # dnaA
        feature_id_2_1 = 'DVMF_RS00075'
        feature_id_3_0 = 'A6701_RS00005'  # dnaA
        feature_id_3_1 = 'A6701_RS00105'
        testFS = {
            'description': 'a few features',
            'elements': { feature_id_1_0: [genome_ref_1],
                          feature_id_1_1: [genome_ref_1],
                          feature_id_2_0: [genome_ref_2],
                          feature_id_2_1: [genome_ref_2],
                          feature_id_3_0: [genome_ref_3],
                          feature_id_3_1: [genome_ref_3]
                      }
        }

        obj_info = self.getWsClient().save_objects({'workspace': self.getWsName(),
                                                    'objects': [
                                                        {
                                                            'type':'KBaseCollections.FeatureSet',
                                                            'data':testFS,
                                                            'name':obj_basename+'.test_FeatureSet',
                                                            'meta':{},
                                                            'provenance':[
                                                                {
                                                                    'service':'kb_blast',
                                                                    'method':'BLASTp_Search'
                                                                }
                                                            ]
                                                        }]
                                                })[0]
        #pprint(obj_info)
        target_featureSet_ref = "/".join([str(obj_info[WORKSPACE_I]),
                                          str(obj_info[OBJID_I]),
                                          str(obj_info[VERSION_I])])


        # E. coli K-12 MG1655 dnaA
        query_seq_prot = 'MSLSLWQQCLARLQDELPATEFSMWIRPLQAELSDNTLALYAPNRFVLDWVRDKYLNNINGLLTSFCGADAPQLRFEVGTKPVTQTPQAAVTSNVAAPAQVAQTQPQRAAPSTRSGWDNVPAPAEPTYRSNVNVKHTFDNFVEGKSNQLARAAARQVADNPGGAYNPLFLYGGTGLGKTHLLHAVGNGIMARKPNAKVVYMHSERFVQDMVKALQNNAIEEFKRYYRSVDALLIDDIQFFANKERSQEEFFHTFNALLEGNQQIILTSDRYPKEINGVEDRLKSRFGWGLTVAIEPPELETRVAILMKKADENDIRLPGEVAFFIAKRLRSNVRELEGALNRVIANANFTGRAITIDFVREALRDLLALQEKLVTIDNIQKTVAEYYKIKVADLLSKRRSRSVARPRQMAMALAKELTNHSLPEIGDAFGGRDHTTVLHACRKIEQLREESHDIKEDFSNLIRTLSS'
        
        parameters = { 'workspace_name': self.getWsName(),
                       'input_one_sequence': query_seq_prot,
                       #'input_one_ref': "",
                       'output_one_name': obj_basename+'.'+"test_query.SS",
                       'input_many_refs': [target_featureSet_ref],
                       'output_filtered_name': obj_out_name,
                       'genome_disp_name_config': 'obj_name_ver_sci_name',
                       'e_value': ".001",
                       'bitscore': "50",
                       'ident_thresh': "40.0",
                       'overlap_fraction': "50.0",
                       'maxaccepts': "1000",
                       'output_extra_format': "none"
                     }

        ret = self.getImpl().BLASTp_Search(self.getContext(), parameters)[0]
        self.assertIsNotNone(ret['report_ref'])

        # check created obj
        #report_obj = self.getWsClient().get_objects2({'objects':[{'ref':ret['report_ref']}]})[0]['data']
        report_obj = self.getWsClient().get_objects([{'ref':ret['report_ref']}])[0]['data']
        self.assertIsNotNone(report_obj['objects_created'][0]['ref'])

        created_obj_0_info = self.getWsClient().get_object_info_new({'objects':[{'ref':report_obj['objects_created'][0]['ref']}]})[0]
        self.assertEqual(created_obj_0_info[NAME_I], obj_out_name)
        self.assertEqual(created_obj_0_info[TYPE_I].split('-')[0], obj_out_type)

        # check number of hits in featureSet output
        featureSet_out_obj = self.getWsClient().get_objects([{'ref':report_obj['objects_created'][0]['ref']}])[0]['data']
        self.assertEqual(expected_hit_cnt, len(featureSet_out_obj['element_ordering']))
        pass


    # Test BLASTp: AnnotatedMetagenomeAssembly Target
    #
    # Uncomment to skip this test
    # HIDE @unittest.skip("skipped test_kb_blast_BLASTp_Search_04_AnnotatedMetagenomeAssembly")
    def test_kb_blast_BLASTp_Search_04_AnnotatedMetagenomeAssembly(self):
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = list(range(11))  # object_info tuple

        obj_basename = 'BLASTp_AnnotatedMetagenomeAssembly'
        obj_out_name = obj_basename+".test_output.FS"
        obj_out_type = "KBaseCollections.FeatureSet"
        expected_hit_cnt = 1
        
        # upload test AMA
        amaInfo_0 = self.getAMAInfo("test_ama", 0)
        ama_ref_1 = self.get_obj_ref_from_obj_info(amaInfo_0)

        # gene 5_267 from test_ama.AMA
        query_seq_prot = 'MDRDALTKLVTDLVSIPSVNPLEGPVGNGRGEAELAAFIHSRLTEAGVVCELKEALPGRPNIIARLPGQSEEMIWFDAHMDTVSGEGMAFPPFEPLIEGDRLLGRGSSDNKGSIATMMAALMEVAKSGERPPLTVVFTATADEEYMMRGMLSLFEAGLTAKAGIVAEPTALEIVIAHKGVARFKISTTGKAAHSSRPEEGVNAIYRMGKVLGAIEAYAKRGVGRETHPLLGKGTLSVGIIRGGEYVNVVPDQCEVDVDRRLLPGEDPRRAVSDVRDYLSNALQEEVGLKVSGPTLTVPGLAVSAESPLVQAVAAAVREVTGKAPLTGMQGATHAGQMAAVDIPALVFGPGQMGQAHTATEELDLTQLERAAAVYERLMRTGL'
        
        parameters = { 'workspace_name': self.getWsName(),
                       'input_one_sequence': query_seq_prot,
                       #'input_one_ref': "",
                       'output_one_name': obj_basename+'.'+"test_query.SS",
                       'input_many_refs': [ama_ref_1],
                       'output_filtered_name': obj_out_name,
                       'genome_disp_name_config': 'obj_name_ver_sci_name',
                       'e_value': ".001",
                       'bitscore': "50",
                       'ident_thresh': "40.0",
                       'overlap_fraction': "50.0",
                       'maxaccepts': "1000",
                       'output_extra_format': "none"
                     }

        ret = self.getImpl().BLASTp_Search(self.getContext(), parameters)[0]
        self.assertIsNotNone(ret['report_ref'])

        # check created obj
        #report_obj = self.getWsClient().get_objects2({'objects':[{'ref':ret['report_ref']}]})[0]['data']
        report_obj = self.getWsClient().get_objects([{'ref':ret['report_ref']}])[0]['data']
        self.assertIsNotNone(report_obj['objects_created'][0]['ref'])

        created_obj_0_info = self.getWsClient().get_object_info_new({'objects':[{'ref':report_obj['objects_created'][0]['ref']}]})[0]
        self.assertEqual(created_obj_0_info[NAME_I], obj_out_name)
        self.assertEqual(created_obj_0_info[TYPE_I].split('-')[0], obj_out_type)

        # check number of hits in featureSet output
        featureSet_out_obj = self.getWsClient().get_objects([{'ref':report_obj['objects_created'][0]['ref']}])[0]['data']
        self.assertEqual(expected_hit_cnt, len(featureSet_out_obj['element_ordering']))
        pass


    # Test BLASTp: Multiple targets of different types
    #
    # Uncomment to skip this test
    # HIDE @unittest.skip("skipped test_kb_blast_BLASTp_Search_05_MultipleTargets")
    def test_kb_blast_BLASTp_Search_05_MultipleTargets(self):
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = list(range(11))  # object_info tuple

        obj_basename = 'BLASTp_MultipleTargets'
        obj_out_name = obj_basename+".test_output.FS"
        obj_out_type = "KBaseCollections.FeatureSet"
        expected_hit_cnt = 10
        
        #reference_prok_genomes_WS = 'ReferenceDataManager'  # PROD and CI
        #genome_ref_1 = 'ReferenceDataManager/GCF_001566335.1/1'  # E. coli K-12 MG1655
        #genome_ref_2 = 'ReferenceDataManager/GCF_002936495.2/1'  # E. coli
        #genome_ref_3 = 'ReferenceDataManager/GCF_002936145.2/1'  # E. coli
        genomeInfo_0 = self.getGenomeInfo('GCF_001566335.1_ASM156633v1_genomic', 0)  # E. coli K-12 MG1655
        genomeInfo_1 = self.getGenomeInfo('GCF_000021385.1_ASM2138v1_genomic', 1)    # D. vulgaris str. 'Miyazaki F
        genomeInfo_2 = self.getGenomeInfo('GCF_001721825.1_ASM172182v1_genomic', 2)  # Pseudomonas aeruginosa
        genome_ref_1 = self.get_obj_ref_from_obj_info(genomeInfo_0)
        genome_ref_2 = self.get_obj_ref_from_obj_info(genomeInfo_1)
        genome_ref_3 = self.get_obj_ref_from_obj_info(genomeInfo_2)

        genome_ref_list = [genome_ref_1, genome_ref_2, genome_ref_3]
        genome_scinames = ['FOO', 'BAR', 'FOOBAR']

        # create GenomeSet
        genomeSet_name = 'test_genomeSet.BLASTp_multiple.GenomeSet'
        testGS = {
            'description': 'three genomes',
            'elements': dict()
        }
        for genome_i, genome_ref in enumerate(genome_ref_list): 
            testGS['elements'][genome_scinames[genome_i]] = { 'ref': genome_ref }

        obj_info = self.getWsClient().save_objects({'workspace': self.getWsName(),       
                                                    'objects': [
                                                        {
                                                            'type':'KBaseSearch.GenomeSet',
                                                            'data':testGS,
                                                            'name':genomeSet_name,
                                                            'meta':{},
                                                            'provenance':[
                                                                {
                                                                    'service':'kb_blast',
                                                                    'method':'BLASTp_Search'
                                                                }
                                                            ]
                                                        }]
                                                })[0]

        #pprint(obj_info)
        target_genomeSet_ref = "/".join([str(obj_info[WORKSPACE_I]),
                                         str(obj_info[OBJID_I]),
                                         str(obj_info[VERSION_I])])


        # upload test AMA
        amaInfo_0 = self.getAMAInfo("test_ama", 0)
        ama_ref_1 = self.get_obj_ref_from_obj_info(amaInfo_0)
        ama_name = amaInfo_0[NAME_I]
        
        # gene 5_267 from ama_test.AMA
        query_seq_prot = 'MDRDALTKLVTDLVSIPSVNPLEGPVGNGRGEAELAAFIHSRLTEAGVVCELKEALPGRPNIIARLPGQSEEMIWFDAHMDTVSGEGMAFPPFEPLIEGDRLLGRGSSDNKGSIATMMAALMEVAKSGERPPLTVVFTATADEEYMMRGMLSLFEAGLTAKAGIVAEPTALEIVIAHKGVARFKISTTGKAAHSSRPEEGVNAIYRMGKVLGAIEAYAKRGVGRETHPLLGKGTLSVGIIRGGEYVNVVPDQCEVDVDRRLLPGEDPRRAVSDVRDYLSNALQEEVGLKVSGPTLTVPGLAVSAESPLVQAVAAAVREVTGKAPLTGMQGATHAGQMAAVDIPALVFGPGQMGQAHTATEELDLTQLERAAAVYERLMRTGL'
        
        parameters = { 'workspace_name': self.getWsName(),
                       'input_one_sequence': query_seq_prot,
                       #'input_one_ref': "",
                       'output_one_name': obj_basename+'.'+"test_query.SS",
                       'input_many_refs': [ama_ref_1, target_genomeSet_ref],
                       'output_filtered_name': obj_out_name,
                       'genome_disp_name_config': 'obj_name_ver_sci_name',
                       'e_value': ".001",
                       'bitscore': "50",
                       'ident_thresh': "10.0",
                       'overlap_fraction': "25.0",
                       'maxaccepts': "1000",
                       'output_extra_format': "none"
                     }

        ret = self.getImpl().BLASTp_Search(self.getContext(), parameters)[0]
        self.assertIsNotNone(ret['report_ref'])

        # check created obj
        #report_obj = self.getWsClient().get_objects2({'objects':[{'ref':ret['report_ref']}]})[0]['data']
        report_obj = self.getWsClient().get_objects([{'ref':ret['report_ref']}])[0]['data']
        self.assertIsNotNone(report_obj['objects_created'][0]['ref'])

        created_obj_0_info = self.getWsClient().get_object_info_new({'objects':[{'ref':report_obj['objects_created'][0]['ref']}]})[0]
        self.assertEqual(created_obj_0_info[NAME_I], obj_out_name)
        self.assertEqual(created_obj_0_info[TYPE_I].split('-')[0], obj_out_type)
        created_obj_1_info = self.getWsClient().get_object_info_new({'objects':[{'ref':report_obj['objects_created'][1]['ref']}]})[0]
        self.assertEqual(created_obj_1_info[NAME_I], obj_out_name+'-'+genomeSet_name)
        self.assertEqual(created_obj_1_info[TYPE_I].split('-')[0], obj_out_type)
        created_obj_2_info = self.getWsClient().get_object_info_new({'objects':[{'ref':report_obj['objects_created'][2]['ref']}]})[0]
        self.assertEqual(created_obj_2_info[NAME_I], obj_out_name+'-'+ama_name)
        self.assertEqual(created_obj_2_info[TYPE_I].split('-')[0], obj_out_type)

        # check number of hits in featureSet output
        featureSet_out_obj = self.getWsClient().get_objects([{'ref':report_obj['objects_created'][0]['ref']}])[0]['data']
        self.assertEqual(expected_hit_cnt, len(featureSet_out_obj['element_ordering']))
        pass


    # Test BLASTx: Single Genome target
    #
    # Uncomment to skip this test
    # HIDE @unittest.skip("skipped test_kb_blast_BLASTx_Search_01")
    def test_kb_blast_BLASTx_Search_01(self):
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = list(range(11))  # object_info tuple

        obj_basename = 'BLASTx'
        obj_out_name = obj_basename+'.'+"test_output.FS"
        obj_out_type = "KBaseCollections.FeatureSet"
        expected_hit_cnt = 1
        
        #reference_prok_genomes_WS = 'ReferenceDataManager'  # PROD and CI
        #genome_ref_1 = 'ReferenceDataManager/GCF_001566335.1/1'  # E. coli K-12 MG1655
        genomeInfo_0 = self.getGenomeInfo('GCF_001566335.1_ASM156633v1_genomic', 0)  # E. coli K-12 MG1655
        genome_ref_1 = self.get_obj_ref_from_obj_info(genomeInfo_0)

        # E. coli K-12 MG1655 dnaA
        query_seq_nuc = 'GTGTCACTTTCGCTTTGGCAGCAGTGTCTTGCCCGATTGCAGGATGAGTTACCAGCCACAGAATTCAGTATGTGGATACGCCCATTGCAGGCGGAACTGAGCGATAACACGCTGGCCCTGTACGCGCCAAACCGTTTTGTCCTCGATTGGGTACGGGACAAGTACCTTAATAATATCAATGGACTGCTAACCAGTTTCTGCGGAGCGGATGCCCCACAGCTGCGTTTTGAAGTCGGCACCAAACCGGTGACGCAAACGCCACAAGCGGCAGTGACGAGCAACGTCGCGGCCCCTGCACAGGTGGCGCAAACGCAGCCGCAACGTGCTGCGCCTTCTACGCGCTCAGGTTGGGATAACGTCCCGGCCCCGGCAGAACCGACCTATCGTTCTAACGTAAACGTCAAACACACGTTTGATAACTTCGTTGAAGGTAAATCTAACCAACTGGCGCGCGCGGCGGCTCGCCAGGTGGCGGATAACCCTGGCGGTGCCTATAACCCGTTGTTCCTTTATGGCGGCACGGGTCTGGGTAAAACTCACCTGCTGCATGCGGTGGGTAACGGCATTATGGCGCGCAAGCCGAATGCCAAAGTGGTTTATATGCACTCCGAGCGCTTTGTTCAGGACATGGTTAAAGCCCTGCAAAACAACGCGATCGAAGAGTTTAAACGCTACTACCGTTCCGTAGATGCACTGCTGATCGACGATATTCAGTTTTTTGCTAATAAAGAACGATCTCAGGAAGAGTTTTTCCACACCTTCAACGCCCTGCTGGAAGGTAATCAACAGATCATTCTCACCTCGGATCGCTATCCGAAAGAGATCAACGGCGTTGAGGATCGTTTGAAATCCCGCTTCGGTTGGGGACTGACTGTGGCGATCGAACCGCCAGAGCTGGAAACCCGTGTGGCGATCCTGATGAAAAAGGCCGACGAAAACGACATTCGTTTGCCGGGCGAAGTGGCGTTCTTTATCGCCAAGCGTCTACGATCTAACGTACGTGAGCTGGAAGGGGCGCTGAACCGCGTCATTGCCAATGCCAACTTTACCGGACGGGCGATCACCATCGACTTCGTGCGTGAGGCGCTGCGCGACTTGCTGGCATTGCAGGAAAAACTGGTCACCATCGACAATATTCAGAAGACGGTGGCGGAGTACTACAAGATCAAAGTCGCGGATCTCCTTTCCAAGCGTCGATCCCGCTCGGTGGCGCGTCCGCGCCAGATGGCGATGGCGCTGGCGAAAGAGCTGACTAACCACAGTCTGCCGGAGATTGGCGATGCGTTTGGTGGCCGTGACCACACGACGGTGCTTCATGCCTGCCGTAAGATCGAGCAGTTGCGTGAAGAGAGCCACGATATCAAAGAAGATTTTTCAAATTTAATCAGAACATTGTCATCGTAA'

        parameters = { 'workspace_name': self.getWsName(),
                       'input_one_sequence': query_seq_nuc,
                       #'input_one_ref': "",
                       'output_one_name': obj_basename+'.'+"test_query.SS",
                       'input_many_refs': [genome_ref_1],
                       'output_filtered_name': obj_out_name,
                       'genome_disp_name_config': 'obj_name_ver',
                       'e_value': ".001",
                       'bitscore': "50",
                       'ident_thresh': "40.0",
                       'overlap_fraction': "25.0",   # 50.0 is too high for this query
                       'maxaccepts': "1000",
                       'output_extra_format': "none"
                     }

        ret = self.getImpl().BLASTx_Search(self.getContext(), parameters)[0]
        self.assertIsNotNone(ret['report_ref'])

        # check created obj
        #report_obj = self.getWsClient().get_objects2({'objects':[{'ref':ret['report_ref']}]})[0]['data']
        report_obj = self.getWsClient().get_objects([{'ref':ret['report_ref']}])[0]['data']
        self.assertIsNotNone(report_obj['objects_created'][0]['ref'])

        created_obj_0_info = self.getWsClient().get_object_info_new({'objects':[{'ref':report_obj['objects_created'][0]['ref']}]})[0]
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = list(range(11))  # object_info tuple
        self.assertEqual(created_obj_0_info[NAME_I], obj_out_name)
        self.assertEqual(created_obj_0_info[TYPE_I].split('-')[0], obj_out_type)

        # check number of hits in featureSet output
        featureSet_out_obj = self.getWsClient().get_objects([{'ref':report_obj['objects_created'][0]['ref']}])[0]['data']
        self.assertEqual(expected_hit_cnt, len(featureSet_out_obj['element_ordering']))
        pass


    # Test tBLASTx
    #
    # SKIPPING tBLASTx test because App disabled
    @unittest.skip("skipped test_kb_blast_tBLASTx_Search_01")
    def test_kb_blast_tBLASTx_Search_01(self):
        obj_basename = 'tBLASTx'
        obj_out_name = obj_basename+'.'+"test_output.FS"
        obj_out_type = "KBaseCollections.FeatureSet"

        reference_prok_genomes_WS = 'ReferenceDataManager'  # PROD and CI
        genome_ref_1 = 'ReferenceDataManager/GCF_001566335.1/1'  # E. coli K-12 MG1655

        # E. coli K-12 MG1655 dnaA
        query_seq_nuc = 'GTGTCACTTTCGCTTTGGCAGCAGTGTCTTGCCCGATTGCAGGATGAGTTACCAGCCACAGAATTCAGTATGTGGATACGCCCATTGCAGGCGGAACTGAGCGATAACACGCTGGCCCTGTACGCGCCAAACCGTTTTGTCCTCGATTGGGTACGGGACAAGTACCTTAATAATATCAATGGACTGCTAACCAGTTTCTGCGGAGCGGATGCCCCACAGCTGCGTTTTGAAGTCGGCACCAAACCGGTGACGCAAACGCCACAAGCGGCAGTGACGAGCAACGTCGCGGCCCCTGCACAGGTGGCGCAAACGCAGCCGCAACGTGCTGCGCCTTCTACGCGCTCAGGTTGGGATAACGTCCCGGCCCCGGCAGAACCGACCTATCGTTCTAACGTAAACGTCAAACACACGTTTGATAACTTCGTTGAAGGTAAATCTAACCAACTGGCGCGCGCGGCGGCTCGCCAGGTGGCGGATAACCCTGGCGGTGCCTATAACCCGTTGTTCCTTTATGGCGGCACGGGTCTGGGTAAAACTCACCTGCTGCATGCGGTGGGTAACGGCATTATGGCGCGCAAGCCGAATGCCAAAGTGGTTTATATGCACTCCGAGCGCTTTGTTCAGGACATGGTTAAAGCCCTGCAAAACAACGCGATCGAAGAGTTTAAACGCTACTACCGTTCCGTAGATGCACTGCTGATCGACGATATTCAGTTTTTTGCTAATAAAGAACGATCTCAGGAAGAGTTTTTCCACACCTTCAACGCCCTGCTGGAAGGTAATCAACAGATCATTCTCACCTCGGATCGCTATCCGAAAGAGATCAACGGCGTTGAGGATCGTTTGAAATCCCGCTTCGGTTGGGGACTGACTGTGGCGATCGAACCGCCAGAGCTGGAAACCCGTGTGGCGATCCTGATGAAAAAGGCCGACGAAAACGACATTCGTTTGCCGGGCGAAGTGGCGTTCTTTATCGCCAAGCGTCTACGATCTAACGTACGTGAGCTGGAAGGGGCGCTGAACCGCGTCATTGCCAATGCCAACTTTACCGGACGGGCGATCACCATCGACTTCGTGCGTGAGGCGCTGCGCGACTTGCTGGCATTGCAGGAAAAACTGGTCACCATCGACAATATTCAGAAGACGGTGGCGGAGTACTACAAGATCAAAGTCGCGGATCTCCTTTCCAAGCGTCGATCCCGCTCGGTGGCGCGTCCGCGCCAGATGGCGATGGCGCTGGCGAAAGAGCTGACTAACCACAGTCTGCCGGAGATTGGCGATGCGTTTGGTGGCCGTGACCACACGACGGTGCTTCATGCCTGCCGTAAGATCGAGCAGTTGCGTGAAGAGAGCCACGATATCAAAGAAGATTTTTCAAATTTAATCAGAACATTGTCATCGTAA'

        parameters = { 'workspace_name': self.getWsName(),
                       'input_one_sequence': query_seq_nuc,
                       #'input_one_ref': "",
                       'output_one_name': obj_basename+'.'+"test_query.SS",
                       'input_many_refs': [genome_ref_1],
                       'output_filtered_name': obj_out_name,
                       'e_value': ".001",
                       'bitscore': "50",
                       'ident_thresh': "40",
                       'overlap_fraction': "50",
                       'maxaccepts': "1000",
                       'output_extra_format': "none"
                     }

        ret = self.getImpl().tBLASTx_Search(self.getContext(), parameters)[0]
        self.assertIsNotNone(ret['report_ref'])

        # check created obj
        #report_obj = self.getWsClient().get_objects2({'objects':[{'ref':ret['report_ref']}]})[0]['data']
        report_obj = self.getWsClient().get_objects([{'ref':ret['report_ref']}])[0]['data']
        self.assertIsNotNone(report_obj['objects_created'][0]['ref'])

        created_obj_0_info = self.getWsClient().get_object_info_new({'objects':[{'ref':report_obj['objects_created'][0]['ref']}]})[0]
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = list(range(11))  # object_info tuple
        self.assertEqual(created_obj_0_info[NAME_I], obj_out_name)
        self.assertEqual(created_obj_0_info[TYPE_I].split('-')[0], obj_out_type)
        pass


    # Test tBLASTn
    #
    # SKIPPING tBLASTn test because App disabled
    @unittest.skip("skipped test_kb_blast_tBLASTn_Search_01")
    def test_kb_blast_tBLASTn_Search_01(self):
        obj_basename = 'tBLASTn'
        obj_out_name = obj_basename+'.'+"test_output.FS"
        obj_out_type = "KBaseCollections.FeatureSet"

        reference_prok_genomes_WS = 'ReferenceDataManager'  # PROD and CI
        genome_ref_1 = 'ReferenceDataManager/GCF_001566335.1/1'  # E. coli K-12 MG1655

        # E. coli K-12 MG1655 dnaA
        query_seq_prot = 'MSLSLWQQCLARLQDELPATEFSMWIRPLQAELSDNTLALYAPNRFVLDWVRDKYLNNINGLLTSFCGADAPQLRFEVGTKPVTQTPQAAVTSNVAAPAQVAQTQPQRAAPSTRSGWDNVPAPAEPTYRSNVNVKHTFDNFVEGKSNQLARAAARQVADNPGGAYNPLFLYGGTGLGKTHLLHAVGNGIMARKPNAKVVYMHSERFVQDMVKALQNNAIEEFKRYYRSVDALLIDDIQFFANKERSQEEFFHTFNALLEGNQQIILTSDRYPKEINGVEDRLKSRFGWGLTVAIEPPELETRVAILMKKADENDIRLPGEVAFFIAKRLRSNVRELEGALNRVIANANFTGRAITIDFVREALRDLLALQEKLVTIDNIQKTVAEYYKIKVADLLSKRRSRSVARPRQMAMALAKELTNHSLPEIGDAFGGRDHTTVLHACRKIEQLREESHDIKEDFSNLIRTLSS'
        
        parameters = { 'workspace_name': self.getWsName(),
                       'input_one_sequence': query_seq_prot,
                       #'input_one_ref': "",
                       'output_one_name': obj_basename+'.'+"test_query.SS",
                       'input_many_refs': [genome_ref_1],
                       'output_filtered_name': obj_out_name,
                       'e_value': ".001",
                       'bitscore': "50",
                       'ident_thresh': "40",
                       'overlap_fraction': "50",
                       'maxaccepts': "1000",
                       'output_extra_format': "none"
                     }

        ret = self.getImpl().tBLASTn_Search(self.getContext(), parameters)[0]
        self.assertIsNotNone(ret['report_ref'])

        # check created obj
        #report_obj = self.getWsClient().get_objects2({'objects':[{'ref':ret['report_ref']}]})[0]['data']
        report_obj = self.getWsClient().get_objects([{'ref':ret['report_ref']}])[0]['data']
        self.assertIsNotNone(report_obj['objects_created'][0]['ref'])

        created_obj_0_info = self.getWsClient().get_object_info_new({'objects':[{'ref':report_obj['objects_created'][0]['ref']}]})[0]
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = list(range(11))  # object_info tuple
        self.assertEqual(created_obj_0_info[NAME_I], obj_out_name)
        self.assertEqual(created_obj_0_info[TYPE_I].split('-')[0], obj_out_type)
        pass


    # Test psiBLAST
    #
    # SKIPPING psiBLAST_msa_start test because App disabled
    @unittest.skip("skipped test_kb_blast_psiBLAST_msa_start_Search_01")
    def test_kb_blast_psiBLAST_msa_start_Search_01(self):
        obj_basename = 'psiBLAST_msa_start'
        obj_out_name = obj_basename+'.'+"test_output.FS"
        obj_out_type = "KBaseCollections.FeatureSet"

        reference_prok_genomes_WS = 'ReferenceDataManager'  # PROD and CI
        genome_ref_1 = 'ReferenceDataManager/GCF_000021385.1/1'  # D. vulgaris str. 'Miyazaki F'

        # MSA
        MSA_json_file = os.path.join('data', 'DsrA.MSA.json')
        with open (MSA_json_file, 'r') as MSA_json_fh:
            MSA_obj = json.load(MSA_json_fh)

        provenance = [{}]
        MSA_info = self.getWsClient().save_objects({
            'workspace': self.getWsName(), 
            'objects': [
                {
                    'type': 'KBaseTrees.MSA',
                    'data': MSA_obj,
                    'name': 'test_MSA',
                    'meta': {},
                    'provenance': provenance
                }
            ]})[0]

        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = list(range(11))  # object_info tuple
        MSA_ref = str(MSA_info[WSID_I])+'/'+str(MSA_info[OBJID_I])+'/'+str(MSA_info[VERSION_I])

        
        parameters = { 'workspace_name': self.getWsName(),
                       #'input_one_sequence': "",
                       #'input_one_ref': "",
                       'input_msa_ref': MSA_ref,
                       #'output_one_name': obj_basename+'.'+"test_query.SS",
                       'input_many_refs': [genome_ref_1],
                       'output_filtered_name': obj_out_name,
                       'e_value': ".001",
                       'bitscore': "50",
                       'ident_thresh': "10",
                       'overlap_fraction': "50",
                       'maxaccepts': "1000",
                       'output_extra_format': "none"
                     }

        ret = self.getImpl().psiBLAST_msa_start_Search(self.getContext(), parameters)[0]
        self.assertIsNotNone(ret['report_ref'])

        # check created obj
        #report_obj = self.getWsClient().get_objects2({'objects':[{'ref':ret['report_ref']}]})[0]['data']
        report_obj = self.getWsClient().get_objects([{'ref':ret['report_ref']}])[0]['data']
        self.assertIsNotNone(report_obj['objects_created'][0]['ref'])

        created_obj_0_info = self.getWsClient().get_object_info_new({'objects':[{'ref':report_obj['objects_created'][0]['ref']}]})[0]
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = list(range(11))  # object_info tuple
        self.assertEqual(created_obj_0_info[NAME_I], obj_out_name)
        self.assertEqual(created_obj_0_info[TYPE_I].split('-')[0], obj_out_type)
        pass


    # SKIPPING psiBLAST tests because App disabled
    @unittest.skip("skipped test_kb_blast_psiBLAST_msa_start_Search_02_nuc_MSA")
    def test_kb_blast_psiBLAST_msa_start_Search_02_nuc_MSA(self):
        obj_basename = 'psiBLAST_msa_start'
        obj_out_name = obj_basename+'.'+"test_output.FS"
        obj_out_type = "KBaseCollections.FeatureSet"

        reference_prok_genomes_WS = 'ReferenceDataManager'  # PROD and CI
        genome_ref_1 = 'ReferenceDataManager/GCF_000021385.1/1'  # D. vulgaris str. 'Miyazaki F'

        # MSA
        MSA_json_file = os.path.join('data', 'ExbD_nuc.MSA.json')
        with open (MSA_json_file, 'r') as MSA_json_fh:
            MSA_obj = json.load(MSA_json_fh)

        provenance = [{}]
        MSA_info = self.getWsClient().save_objects({
            'workspace': self.getWsName(), 
            'objects': [
                {
                    'type': 'KBaseTrees.MSA',
                    'data': MSA_obj,
                    'name': 'test_MSA_nuc',
                    'meta': {},
                    'provenance': provenance
                }
            ]})[0]

        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = list(range(11))  # object_info tuple
        MSA_ref = str(MSA_info[WSID_I])+'/'+str(MSA_info[OBJID_I])+'/'+str(MSA_info[VERSION_I])

        
        parameters = { 'workspace_name': self.getWsName(),
                       #'input_one_sequence': "",
                       #'input_one_ref': "",
                       'input_msa_ref': MSA_ref,
                       #'output_one_name': obj_basename+'.'+"test_query.SS",
                       'input_many_refs': [genome_ref_1],
                       'output_filtered_name': obj_out_name,
                       'e_value': ".001",
                       'bitscore': "50",
                       'ident_thresh': "10",
                       'overlap_fraction': "50",
                       'maxaccepts': "1000",
                       'output_extra_format': "none"
                     }

        ret = self.getImpl().psiBLAST_msa_start_Search(self.getContext(), parameters)[0]
        self.assertIsNotNone(ret['report_ref'])

        # check created obj
        #report_obj = self.getWsClient().get_objects2({'objects':[{'ref':ret['report_ref']}]})[0]['data']
        report_obj = self.getWsClient().get_objects([{'ref':ret['report_ref']}])[0]['data']
        self.assertEqual(report_obj['text_message'][0:7],"FAILURE")
        pass
