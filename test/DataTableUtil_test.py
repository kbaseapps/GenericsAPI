# -*- coding: utf-8 -*-
import inspect
import os
import unittest
import time
from mock import patch
import shutil

from configparser import ConfigParser

from GenericsAPI.Utils.DataTableUtil import DataTableUtil
from GenericsAPI.GenericsAPIImpl import GenericsAPI
from GenericsAPI.GenericsAPIServer import MethodContext
from GenericsAPI.authclient import KBaseAuth as _KBaseAuth
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace as workspaceService


class DataTableTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('GenericsAPI'):
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
                            {'service': 'GenericsAPI',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = workspaceService(cls.wsURL)
        cls.serviceImpl = GenericsAPI(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        cls.dfu = DataFileUtil(cls.callback_url)
        cls.data_table_util = DataTableUtil(cls.cfg)

        suffix = int(time.time() * 1000)
        cls.wsName = "test_view_matrix_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})
        cls.wsId = ret[0]

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        return self.__class__.wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    def getDataTableUtil(self):
        return self.__class__.data_table_util

    def mock_download_staging_file(params):
        print('Mocking DataFileUtilClient.download_staging_file')
        print(params)

        file_path = params.get('staging_file_subdir_path')

        return {'copy_file_path': file_path}

    # @patch.object(DataFileUtil, "download_staging_file", side_effect=mock_download_staging_file)
    # def loadAmpliconMatrix(self, download_staging_file):
    #     if hasattr(self.__class__, 'amplicon_matrix_ref'):
    #         return self.__class__.amplicon_matrix_ref

    #     col_attribute = {'attributes': [{'attribute': 'BODY_SITE', 'source': 'upload'},
    #                                     {'attribute': 'BarcodeSequence', 'source': 'upload'},
    #                                     {'attribute': 'Description', 'source': 'upload'},
    #                                     {'attribute': 'LinkerPrimerSequence', 'source': 'upload'}],
    #                      'instances': {'Sample1': ['gut', 'CGCTTATCGAGA', 'human gut', 'CATGCTGCCTCCCGTAGGAGT'],
    #                                    'Sample2': ['gut', 'CATACCAGTAGC', 'human gut', 'CATGCTGCCTCCCGTAGGAGT'],
    #                                    'Sample3': ['gut', 'CTCTCTACCTGT', 'human gut', 'CATGCTGCCTCCCGTAGGAGT'],
    #                                    'Sample4': ['skin', 'CTCTCGGCCTGT', 'human skin', 'CATGCTGCCTCCCGTAGGAGT'],
    #                                    'Sample5': ['skin', 'CTCTCTACCAAT', 'human skin', 'CATGCTGCCTCCCGTAGGAGT'],
    #                                    'Sample6': ['skin', 'CTAACTACCAAT', 'human skin', 'CATGCTGCCTCCCGTAGGAGT']},
    #                      'ontology_mapping_method': 'BIOM file'}

    #     info = self.dfu.save_objects({
    #                         'id': self.wsId,
    #                         'objects': [{
    #                             'type': 'KBaseExperiments.AttributeMapping',
    #                             'data': col_attribute,
    #                             'name': 'test_AmpliconMatrix_col_attributes'
    #                         }]
    #                     })[0]

    #     col_attributemapping_ref = "%s/%s/%s" % (info[6], info[0], info[4])

    #     row_attribute = {'attributes': [{'attribute': 'taxonomy', 'source': 'upload'}],
    #                      'instances': {'GG_OTU_1': ["['k__Bacteria', 'p__Proteobacteria', 'c__Gammaproteobacteria', 'o__Enterobacteriales', 'f__Enterobacteriaceae', 'g__Escherichia', 's__']"],
    #                                    'GG_OTU_2': ["['k__Bacteria', 'p__Cyanobacteria', 'c__Nostocophycideae', 'o__Nostocales', 'f__Nostocaceae', 'g__Dolichospermum', 's__']"],
    #                                    'GG_OTU_3': ["['k__Archaea', 'p__Euryarchaeota', 'c__Methanomicrobia', 'o__Methanosarcinales', 'f__Methanosarcinaceae', 'g__Methanosarcina', 's__']"],
    #                                    'GG_OTU_4': ["['k__Bacteria', 'p__Firmicutes', 'c__Clostridia', 'o__Halanaerobiales', 'f__Halanaerobiaceae', 'g__Halanaerobium', 's__Halanaerobiumsaccharolyticum']"],
    #                                    'GG_OTU_5': ["['k__Bacteria', 'p__Proteobacteria', 'c__Gammaproteobacteria', 'o__Enterobacteriales', 'f__Enterobacteriaceae', 'g__Escherichia', 's__']"]},
    #                      'ontology_mapping_method': 'BIOM file'}

    #     info = self.dfu.save_objects({
    #                         'id': self.wsId,
    #                         'objects': [{
    #                             'type': 'KBaseExperiments.AttributeMapping',
    #                             'data': row_attribute,
    #                             'name': 'test_AmpliconMatrix_row_attributes'
    #                         }]
    #                     })[0]

    #     row_attributemapping_ref = "%s/%s/%s" % (info[6], info[0], info[4])

    #     matrix_data = {'amplicon_type': '16S',
    #                    'attributes': {'generated_by': 'QIIME revision XYZ'},
    #                    'clustering_cutoff': 0.3,
    #                    'clustering_method': 'clustering_method',
    #                    'col_attributemapping_ref': col_attributemapping_ref,
    #                    'col_mapping': {'Sample1': 'Sample1',
    #                                    'Sample2': 'Sample2',
    #                                    'Sample3': 'Sample3',
    #                                    'Sample4': 'Sample4',
    #                                    'Sample5': 'Sample5',
    #                                    'Sample6': 'Sample6'},
    #                    'data': {'col_ids': ['Sample1', 'Sample2', 'Sample3', 'Sample4', 'Sample5', 'Sample6'],
    #                             'row_ids': ['GG_OTU_1', 'GG_OTU_2', 'GG_OTU_3', 'GG_OTU_4', 'GG_OTU_5'],
    #                             'values': [[0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
    #                                        [5.0, 1.0, 0.0, 2.0, 3.0, 1.0],
    #                                        [0.0, 0.0, 1.0, 4.0, 2.0, 0.0],
    #                                        [2.0, 1.0, 1.0, 0.0, 0.0, 1.0],
    #                                        [0.0, 1.0, 1.0, 0.0, 0.0, 0.0]]},
    #                    'description': 'OTU data',
    #                    'forward_primer_sequence': 'forward_primer_sequence',
    #                    'reverse_primer_sequence': 'reverse_primer_sequence',
    #                    'row_attributemapping_ref': row_attributemapping_ref,
    #                    'row_mapping': {'GG_OTU_1': 'GG_OTU_1',
    #                                    'GG_OTU_2': 'GG_OTU_2',
    #                                    'GG_OTU_3': 'GG_OTU_3',
    #                                    'GG_OTU_4': 'GG_OTU_4',
    #                                    'GG_OTU_5': 'GG_OTU_5'},
    #                    'scale': 'raw',
    #                    'search_attributes': ['generated_by|QIIME revision XYZ'],
    #                    'sequencing_platform': 'Illumina',
    #                    'sequencing_quality_filter_cutoff': 'sequencing_quality_filter_cutoff',
    #                    'target_gene_region': 'V1'}

    #     info = self.dfu.save_objects({
    #                         'id': self.wsId,
    #                         'objects': [{
    #                             'type': 'KBaseMatrices.AmpliconMatrix',
    #                             'data': matrix_data,
    #                             'name': 'test_AmpliconMatrix'
    #                         }]
    #                     })[0]

    #     amplicon_matrix_ref = "%s/%s/%s" % (info[6], info[0], info[4])

    #     self.__class__.amplicon_matrix_ref = amplicon_matrix_ref
    #     print('Loaded AmpliconMatrix: ' + amplicon_matrix_ref)
    #     return amplicon_matrix_ref

    @patch.object(DataFileUtil, "download_staging_file", side_effect=mock_download_staging_file)
    def loadAmpliconMatrix(self, download_staging_file):
        if hasattr(self.__class__, 'amplicon_matrix_ref'):
            return self.__class__.amplicon_matrix_ref

        taxonomic_abundance_tsv = os.path.join(self.scratch, 'amplicon_test.tsv')
        shutil.copy(os.path.join('data', 'amplicon_test.tsv'), taxonomic_abundance_tsv)

        taxonomic_fasta = os.path.join(self.scratch, 'phyloseq_test.fa')
        shutil.copy(os.path.join('data', 'phyloseq_test.fa'), taxonomic_fasta)

        params = {'obj_type': 'AmpliconMatrix',
                  'matrix_name': 'test_AmpliconMatrix',
                  'workspace_id': self.wsId,
                  'taxonomic_abundance_tsv': taxonomic_abundance_tsv,
                  'taxonomic_fasta': taxonomic_fasta,
                  'metadata_keys': 'taxonomy_id, taxonomy, taxonomy_source, consensus_sequence',
                  'scale': 'raw',
                  'description': "OTU data",
                  'amplicon_type': '16S',
                  'target_gene_region': 'V1',
                  'forward_primer_sequence': 'forward_primer_sequence',
                  'reverse_primer_sequence': 'reverse_primer_sequence',
                  'sequencing_platform': 'Illumina',
                  'sequencing_quality_filter_cutoff': 'sequencing_quality_filter_cutoff',
                  'clustering_cutoff': 0.3,
                  'clustering_method': 'clustering_method',
                  'input_local_file': True
                  }

        returnVal = self.getImpl().import_matrix_from_biom(self.ctx, params)[0]

        amplicon_matrix_ref = returnVal['matrix_obj_ref']

        self.__class__.amplicon_matrix_ref = amplicon_matrix_ref
        print('Loaded AmpliconMatrix: ' + amplicon_matrix_ref)
        return amplicon_matrix_ref

    def fail_view_matrix_as_table(self, params, error, exception=ValueError, contains=False):
        with self.assertRaises(exception) as context:
            self.getImpl().view_matrix(self.ctx, params)
        if contains:
            self.assertIn(error, str(context.exception.args[0]))
        else:
            self.assertEqual(error, str(context.exception.args[0]))

    def start_test(self):
        testname = inspect.stack()[1][3]
        print('\n*** starting test: ' + testname + ' **')

    def test_view_matrix_as_table_ok(self):
        self.start_test()
        matrix_ref = self.loadAmpliconMatrix()

        params = {'input_matrix_ref': matrix_ref,
                  'workspace_name': self.wsName,
                  'with_attribute_info': 1}

        ret = self.getImpl().view_matrix(self.ctx, params)[0]

        self.assertIn('report_name', ret)
        self.assertIn('report_ref', ret)

    def test_init_ok(self):
        self.start_test()
        class_attri = ['scratch', 'token', 'callback_url', 'dfu', 'matrix_types']

        data_table_util = self.getDataTableUtil()
        self.assertTrue(set(class_attri) <= set(data_table_util.__dict__.keys()))
        self.assertEqual(data_table_util.scratch, self.cfg.get('scratch'))
