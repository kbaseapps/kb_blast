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

    # retrieve stored obj info
    def _get_stored_obj_info (self, obj_type, obj_name, item_i=0):
        infoAttr = obj_type + 'Info_list' # e.g. 'ama' or 'genome'
        nameAttr = obj_type + 'Name_list'
        if hasattr(self.__class__, infoAttr):
            try:
                info_list = getattr(self.__class__, infoAttr)
                name_list = getattr(self.__class__, nameAttr)
                info      = info_list[item_i]
                name      = name_list[item_i]
                if info != None:
                    if name != obj_name:
                        info_list[item_i] = None
                        name_list[item_i] = None
                        setattr (self.__class__, infoAttr, info_list)
                        setattr (self.__class__, nameAttr, name_list)
                    else:
                        return info
            except:
                pass
        return None

    # save stored obj info
    def _save_stored_obj_info (self, obj_type, obj_info, obj_name, item_i=0):
        infoAttr = obj_type + 'Info_list' # e.g. 'ama' or 'genome'
        nameAttr = obj_type + 'Name_list'
        if not hasattr(self.__class__, infoAttr):
            setattr (self.__class__, infoAttr, [])
            setattr (self.__class__, nameAttr, [])

        info_list = getattr(self.__class__, infoAttr)
        name_list = getattr(self.__class__, nameAttr)
        for i in range(item_i+1):
            try:
                assigned = info_list[i]
            except:
                info_list.append(None)
                name_list.append(None)
        info_list[item_i] = obj_info
        name_list[item_i] = obj_name
        setattr (self.__class__, infoAttr, info_list)
        setattr (self.__class__, nameAttr, name_list)
        return
        
    # call this method to get the WS object info of a Genome
    #   (will upload the example data if this is the first time the method is called during tests)
    def getGenomeInfo(self, genome_basename, item_i=0):
        info = self._get_stored_obj_info ('genome', genome_basename, item_i)
        if info != None:
            return info

        # 1) transform genbank to kbase genome object and upload to ws
        shared_dir = "/kb/module/work/tmp"
        genome_data_file = 'data/genomes/'+genome_basename+'.gbff.gz'
        genome_file = os.path.join(shared_dir, os.path.basename(genome_data_file))
        shutil.copy(genome_data_file, genome_file)

        SERVICE_VER = 'release'
        GFU = GenomeFileUtil(os.environ['SDK_CALLBACK_URL'],
                             token=self.getContext()['token'],
                             service_ver=SERVICE_VER
                         )
        print ("UPLOADING genome: "+genome_basename+" to WORKSPACE "+self.getWsName()+" ...")
        genome_upload_result = GFU.genbank_to_genome({'file': {'path': genome_file },
                                                      'workspace_name': self.getWsName(),
                                                      'genome_name': genome_basename
                                                  })
        pprint(genome_upload_result)
        genome_ref = genome_upload_result['genome_ref']
        new_obj_info = self.getWsClient().get_object_info_new({'objects': [{'ref': genome_ref}]})[0]

        # 2) store it
        self._save_stored_obj_info ('genome', new_obj_info, genome_basename, item_i)
        return new_obj_info

    # call this method to get the WS object info of an AnnotatedMetagenomeAssembly
    #   (will upload the example data if this is the first time the method is called during tests)
    def getAMAInfo(self, ama_basename, item_i=0):
        info = self._get_stored_obj_info ('ama', ama_basename, item_i)
        if info != None:
            return info

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
        self._save_stored_obj_info ('ama', new_obj_info, ama_basename, item_i)
        return new_obj_info

    #
    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    #


    # Test BLASTn: Single Genome target
    #
    # Uncomment to skip this test
    # HIDE @unittest.skip("skipped test_kb_blast_BLASTn_Search_01")
    def test_kb_blast_BLASTn_Search_01_Genome(self):
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = list(range(11))  # object_info tuple

        obj_basename = 'BLASTn'
        obj_out_name = obj_basename+".test_output.FS"
        obj_out_type = "KBaseCollections.FeatureSet"
        expected_hit_cnt = 1
        
        genomeInfo_0 = self.getGenomeInfo('GCF_001566335.1_ASM156633v1_genomic', 0)  # E. coli K-12 MG1655
        genome_ref_0 = self.get_obj_ref_from_obj_info(genomeInfo_0)

        # E. coli K-12 MG1655 dnaA
        query_seq_nuc = 'GTGTCACTTTCGCTTTGGCAGCAGTGTCTTGCCCGATTGCAGGATGAGTTACCAGCCACAGAATTCAGTATGTGGATACGCCCATTGCAGGCGGAACTGAGCGATAACACGCTGGCCCTGTACGCGCCAAACCGTTTTGTCCTCGATTGGGTACGGGACAAGTACCTTAATAATATCAATGGACTGCTAACCAGTTTCTGCGGAGCGGATGCCCCACAGCTGCGTTTTGAAGTCGGCACCAAACCGGTGACGCAAACGCCACAAGCGGCAGTGACGAGCAACGTCGCGGCCCCTGCACAGGTGGCGCAAACGCAGCCGCAACGTGCTGCGCCTTCTACGCGCTCAGGTTGGGATAACGTCCCGGCCCCGGCAGAACCGACCTATCGTTCTAACGTAAACGTCAAACACACGTTTGATAACTTCGTTGAAGGTAAATCTAACCAACTGGCGCGCGCGGCGGCTCGCCAGGTGGCGGATAACCCTGGCGGTGCCTATAACCCGTTGTTCCTTTATGGCGGCACGGGTCTGGGTAAAACTCACCTGCTGCATGCGGTGGGTAACGGCATTATGGCGCGCAAGCCGAATGCCAAAGTGGTTTATATGCACTCCGAGCGCTTTGTTCAGGACATGGTTAAAGCCCTGCAAAACAACGCGATCGAAGAGTTTAAACGCTACTACCGTTCCGTAGATGCACTGCTGATCGACGATATTCAGTTTTTTGCTAATAAAGAACGATCTCAGGAAGAGTTTTTCCACACCTTCAACGCCCTGCTGGAAGGTAATCAACAGATCATTCTCACCTCGGATCGCTATCCGAAAGAGATCAACGGCGTTGAGGATCGTTTGAAATCCCGCTTCGGTTGGGGACTGACTGTGGCGATCGAACCGCCAGAGCTGGAAACCCGTGTGGCGATCCTGATGAAAAAGGCCGACGAAAACGACATTCGTTTGCCGGGCGAAGTGGCGTTCTTTATCGCCAAGCGTCTACGATCTAACGTACGTGAGCTGGAAGGGGCGCTGAACCGCGTCATTGCCAATGCCAACTTTACCGGACGGGCGATCACCATCGACTTCGTGCGTGAGGCGCTGCGCGACTTGCTGGCATTGCAGGAAAAACTGGTCACCATCGACAATATTCAGAAGACGGTGGCGGAGTACTACAAGATCAAAGTCGCGGATCTCCTTTCCAAGCGTCGATCCCGCTCGGTGGCGCGTCCGCGCCAGATGGCGATGGCGCTGGCGAAAGAGCTGACTAACCACAGTCTGCCGGAGATTGGCGATGCGTTTGGTGGCCGTGACCACACGACGGTGCTTCATGCCTGCCGTAAGATCGAGCAGTTGCGTGAAGAGAGCCACGATATCAAAGAAGATTTTTCAAATTTAATCAGAACATTGTCATCGTAA'

        parameters = { 'workspace_name': self.getWsName(),
                       'input_one_sequence': query_seq_nuc,
                       #'input_one_ref': "",
                       'output_one_name': obj_basename+'.'+"test_query.SS",
                       'input_many_refs': [genome_ref_0],
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
    # HIDE @unittest.skip("skipped test_kb_blast_BLASTn_Search_01_GenomeSet")
    def test_kb_blast_BLASTn_Search_02_GenomeSet(self):
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = list(range(11))  # object_info tuple

        obj_basename = 'BLASTn'
        obj_out_name = obj_basename+".test_output.FS"
        obj_out_type = "KBaseCollections.FeatureSet"
        genomeSet_name = 'test_genomeSet.BLASTn.GenomeSet'
        expected_hit_cnt = 3
        
        load_genomes = [
            { 'file': 'GCF_001566335.1_ASM156633v1_genomic',
              'sciname': 'E. coli K-12 MG1655'
            },
            { 'file': 'GCF_000021385.1_ASM2138v1_genomic',
              'sciname': 'D. vulgaris str. Miyazaki F'
            },
            { 'file': 'GCF_001721825.1_ASM172182v1_genomic',
              'sciname': 'Pseudomonas aeruginosa'
            },
            { 'file': 'GCF_002950035.1_ASM295003v1_genomic',
              'sciname': 'Shigella boydii'
            },
            { 'file': 'GCF_000512125.1_ASM51212v1_genomic',
              'sciname': 'Escherichia albertii KF1'
            }
        ]
        for genome_i,genome in enumerate(load_genomes):
            load_genomes[genome_i]['ref'] = self.get_obj_ref_from_obj_info(self.getGenomeInfo(genome['file'], genome_i))

        # create GenomeSet
        testGS = {
            'description': 'five genomes',
            'elements': dict()
        }
        for genome_i,genome in enumerate(load_genomes): 
            testGS['elements'][genome['sciname']] = { 'ref': genome['ref'] }

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
        target_genomeSet_ref = self.get_obj_ref_from_obj_info(obj_info)

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
    def test_kb_blast_BLASTp_Search_03_Genome(self):
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = list(range(11))  # object_info tuple

        obj_basename = 'BLASTp_Genome'
        obj_out_name = obj_basename+".test_output.FS"
        obj_out_type = "KBaseCollections.FeatureSet"
        expected_hit_cnt = 1
        
        genomeInfo_0 = self.getGenomeInfo('GCF_001566335.1_ASM156633v1_genomic', 0)  # E. coli K-12 MG1655
        genome_ref_0 = self.get_obj_ref_from_obj_info(genomeInfo_0)

        # E. coli K-12 MG1655 dnaA
        query_seq_prot = 'MSLSLWQQCLARLQDELPATEFSMWIRPLQAELSDNTLALYAPNRFVLDWVRDKYLNNINGLLTSFCGADAPQLRFEVGTKPVTQTPQAAVTSNVAAPAQVAQTQPQRAAPSTRSGWDNVPAPAEPTYRSNVNVKHTFDNFVEGKSNQLARAAARQVADNPGGAYNPLFLYGGTGLGKTHLLHAVGNGIMARKPNAKVVYMHSERFVQDMVKALQNNAIEEFKRYYRSVDALLIDDIQFFANKERSQEEFFHTFNALLEGNQQIILTSDRYPKEINGVEDRLKSRFGWGLTVAIEPPELETRVAILMKKADENDIRLPGEVAFFIAKRLRSNVRELEGALNRVIANANFTGRAITIDFVREALRDLLALQEKLVTIDNIQKTVAEYYKIKVADLLSKRRSRSVARPRQMAMALAKELTNHSLPEIGDAFGGRDHTTVLHACRKIEQLREESHDIKEDFSNLIRTLSS'
        
        parameters = { 'workspace_name': self.getWsName(),
                       'input_one_sequence': query_seq_prot,
                       #'input_one_ref': "",
                       'output_one_name': obj_basename+'.'+"test_query.SS",
                       'input_many_refs': [genome_ref_0],
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
    def test_kb_blast_BLASTp_Search_04_GenomeSet(self):
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple

        obj_basename = 'BLASTp_GenomeSet'
        obj_out_name = obj_basename+".test_output.FS"
        obj_out_type = "KBaseCollections.FeatureSet"
        genomeSet_name = 'test_genomeSet.BLASTp.GenomeSet'
        expected_hit_cnt = 2

        load_genomes = [
            { 'file': 'GCF_001566335.1_ASM156633v1_genomic',
              'sciname': 'E. coli K-12 MG1655'
            },
            { 'file': 'GCF_000021385.1_ASM2138v1_genomic',
              'sciname': 'D. vulgaris str. Miyazaki F'
            },
            { 'file': 'GCF_001721825.1_ASM172182v1_genomic',
              'sciname': 'Pseudomonas aeruginosa'
            },
        ]
        for genome_i,genome in enumerate(load_genomes):
            load_genomes[genome_i]['ref'] = self.get_obj_ref_from_obj_info(self.getGenomeInfo(genome['file'], genome_i))

        # create GenomeSet
        testGS = {
            'description': 'three genomes',
            'elements': dict()
        }
        for genome_i,genome in enumerate(load_genomes): 
            testGS['elements'][genome['sciname']] = { 'ref': genome['ref'] }

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
        target_genomeSet_ref = self.get_obj_ref_from_obj_info(obj_info)


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
    def test_kb_blast_BLASTp_Search_05_FeatureSet(self):
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple

        obj_basename = 'BLASTp_FeatureSet'
        obj_out_name = obj_basename+".test_output.FS"
        obj_out_type = "KBaseCollections.FeatureSet"
        target_1 = obj_basename+'.test_FeatureSet'
        expected_hit_cnt = 2
        
        load_genomes = [
            { 'file': 'GCF_001566335.1_ASM156633v1_genomic',
              'sciname': 'E. coli K-12 MG1655'
            },
            { 'file': 'GCF_000021385.1_ASM2138v1_genomic',
              'sciname': 'D. vulgaris str. Miyazaki F'
            },
            { 'file': 'GCF_001721825.1_ASM172182v1_genomic',
              'sciname': 'Pseudomonas aeruginosa'
            },
        ]
        for genome_i,genome in enumerate(load_genomes):
            load_genomes[genome_i]['ref'] = self.get_obj_ref_from_obj_info(self.getGenomeInfo(genome['file'], genome_i))

        # build FeatureSet obj
        feature_ids = [ [ 'AWN69_RS07145',  # dnaA
                          'AWN69_RS00105'
                        ],
                        [ 'DVMF_RS00005',  # dnaA
                          'DVMF_RS00075',
                        ],
                        [ 'A6701_RS00005',  # dnaA
                          'A6701_RS00105'
                        ]
                      ]
        testFS = {
            'description': 'a few features',
            'elements': { feature_ids[0][0]: [load_genomes[0]['ref']],
                          feature_ids[0][1]: [load_genomes[0]['ref']],
                          feature_ids[1][0]: [load_genomes[1]['ref']],
                          feature_ids[1][1]: [load_genomes[1]['ref']],
                          feature_ids[2][0]: [load_genomes[2]['ref']],
                          feature_ids[2][1]: [load_genomes[2]['ref']]
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
        target_featureSet_ref = self.get_obj_ref_from_obj_info(obj_info)

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
    def test_kb_blast_BLASTp_Search_06_AnnotatedMetagenomeAssembly(self):
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
    def test_kb_blast_BLASTp_Search_07_MultipleTargets(self):
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = list(range(11))  # object_info tuple

        obj_basename = 'BLASTp_MultipleTargets'
        obj_out_name = obj_basename+".test_output.FS"
        obj_out_type = "KBaseCollections.FeatureSet"
        genomeSet_name = 'test_genomeSet_multiple.BLASTp.GenomeSet'
        expected_hit_cnt = 10
        
        load_genomes = [
            { 'file': 'GCF_001566335.1_ASM156633v1_genomic',
              'sciname': 'E. coli K-12 MG1655'
            },
            { 'file': 'GCF_000021385.1_ASM2138v1_genomic',
              'sciname': 'D. vulgaris str. Miyazaki F'
            },
            { 'file': 'GCF_001721825.1_ASM172182v1_genomic',
              'sciname': 'Pseudomonas aeruginosa'
            },
        ]
        for genome_i,genome in enumerate(load_genomes):
            load_genomes[genome_i]['ref'] = self.get_obj_ref_from_obj_info(self.getGenomeInfo(genome['file'], genome_i))

        # create GenomeSet
        testGS = {
            'description': 'three genomes',
            'elements': dict()
        }
        for genome_i,genome in enumerate(load_genomes): 
            testGS['elements'][genome['sciname']] = { 'ref': genome['ref'] }

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
        target_genomeSet_ref = self.get_obj_ref_from_obj_info(obj_info)

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
    def test_kb_blast_BLASTx_Search_08_Genome(self):
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = list(range(11))  # object_info tuple

        obj_basename = 'BLASTx'
        obj_out_name = obj_basename+'.'+"test_output.FS"
        obj_out_type = "KBaseCollections.FeatureSet"
        expected_hit_cnt = 1
        
        genomeInfo_0 = self.getGenomeInfo('GCF_001566335.1_ASM156633v1_genomic', 0)  # E. coli K-12 MG1655
        genome_ref_0 = self.get_obj_ref_from_obj_info(genomeInfo_0)

        # E. coli K-12 MG1655 dnaA
        query_seq_nuc = 'GTGTCACTTTCGCTTTGGCAGCAGTGTCTTGCCCGATTGCAGGATGAGTTACCAGCCACAGAATTCAGTATGTGGATACGCCCATTGCAGGCGGAACTGAGCGATAACACGCTGGCCCTGTACGCGCCAAACCGTTTTGTCCTCGATTGGGTACGGGACAAGTACCTTAATAATATCAATGGACTGCTAACCAGTTTCTGCGGAGCGGATGCCCCACAGCTGCGTTTTGAAGTCGGCACCAAACCGGTGACGCAAACGCCACAAGCGGCAGTGACGAGCAACGTCGCGGCCCCTGCACAGGTGGCGCAAACGCAGCCGCAACGTGCTGCGCCTTCTACGCGCTCAGGTTGGGATAACGTCCCGGCCCCGGCAGAACCGACCTATCGTTCTAACGTAAACGTCAAACACACGTTTGATAACTTCGTTGAAGGTAAATCTAACCAACTGGCGCGCGCGGCGGCTCGCCAGGTGGCGGATAACCCTGGCGGTGCCTATAACCCGTTGTTCCTTTATGGCGGCACGGGTCTGGGTAAAACTCACCTGCTGCATGCGGTGGGTAACGGCATTATGGCGCGCAAGCCGAATGCCAAAGTGGTTTATATGCACTCCGAGCGCTTTGTTCAGGACATGGTTAAAGCCCTGCAAAACAACGCGATCGAAGAGTTTAAACGCTACTACCGTTCCGTAGATGCACTGCTGATCGACGATATTCAGTTTTTTGCTAATAAAGAACGATCTCAGGAAGAGTTTTTCCACACCTTCAACGCCCTGCTGGAAGGTAATCAACAGATCATTCTCACCTCGGATCGCTATCCGAAAGAGATCAACGGCGTTGAGGATCGTTTGAAATCCCGCTTCGGTTGGGGACTGACTGTGGCGATCGAACCGCCAGAGCTGGAAACCCGTGTGGCGATCCTGATGAAAAAGGCCGACGAAAACGACATTCGTTTGCCGGGCGAAGTGGCGTTCTTTATCGCCAAGCGTCTACGATCTAACGTACGTGAGCTGGAAGGGGCGCTGAACCGCGTCATTGCCAATGCCAACTTTACCGGACGGGCGATCACCATCGACTTCGTGCGTGAGGCGCTGCGCGACTTGCTGGCATTGCAGGAAAAACTGGTCACCATCGACAATATTCAGAAGACGGTGGCGGAGTACTACAAGATCAAAGTCGCGGATCTCCTTTCCAAGCGTCGATCCCGCTCGGTGGCGCGTCCGCGCCAGATGGCGATGGCGCTGGCGAAAGAGCTGACTAACCACAGTCTGCCGGAGATTGGCGATGCGTTTGGTGGCCGTGACCACACGACGGTGCTTCATGCCTGCCGTAAGATCGAGCAGTTGCGTGAAGAGAGCCACGATATCAAAGAAGATTTTTCAAATTTAATCAGAACATTGTCATCGTAA'

        parameters = { 'workspace_name': self.getWsName(),
                       'input_one_sequence': query_seq_nuc,
                       #'input_one_ref': "",
                       'output_one_name': obj_basename+'.'+"test_query.SS",
                       'input_many_refs': [genome_ref_0],
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
