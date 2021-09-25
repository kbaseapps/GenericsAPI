# -*- coding: utf-8 -*-
import inspect
import json
import os
import time
import unittest
from configparser import ConfigParser
from os import environ
from mock import patch
import shutil
from Bio import SeqIO
import requests

from installed_clients.DataFileUtilClient import DataFileUtil
from GenericsAPI.GenericsAPIImpl import GenericsAPI
from GenericsAPI.GenericsAPIServer import MethodContext
from GenericsAPI.authclient import KBaseAuth as _KBaseAuth
from installed_clients.WorkspaceClient import Workspace as workspaceService
from installed_clients.sample_uploaderClient import sample_uploader
from installed_clients.FakeObjectsForTestsClient import FakeObjectsForTests
from installed_clients.SampleServiceClient import SampleService
from installed_clients.AbstractHandleClient import AbstractHandle as HandleService


class BioMultiTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.token = environ.get('KB_AUTH_TOKEN', None)
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('GenericsAPI'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(cls.token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': cls.token,
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
        cls.shockURL = cls.cfg['shock-url']
        cls.dfu = DataFileUtil(cls.callback_url)
        cls.sample_uploader = sample_uploader(cls.callback_url, service_ver="dev")
        cls.sample_url = cls.cfg.get('kbase-endpoint') + '/sampleservice'
        cls.sample_ser = SampleService(cls.sample_url)
        cls.hs = HandleService(url=cls.cfg['handle-service-url'],
                               token=cls.token)

        suffix = int(time.time() * 1000)
        cls.wsName = "test_GenericsAPI_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})
        cls.wsId = ret[0]

        small_file = os.path.join(cls.scratch, 'test.txt')
        with open(small_file, "w") as f:
            f.write("empty content")
        cls.test_shock = cls.dfu.file_to_shock({'file_path': small_file, 'make_handle': True})
        cls.handles_to_delete = []
        cls.nodes_to_delete = []
        cls.handles_to_delete.append(cls.test_shock['handle']['hid'])
        cls.nodes_to_delete.append(cls.test_shock['shock_id'])

        cls.prepare_data()

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')
        if hasattr(cls, 'nodes_to_delete'):
            for node in cls.nodes_to_delete:
                cls.delete_shock_node(node)
        if hasattr(cls, 'handles_to_delete'):
            cls.hs.delete_handles(cls.hs.hids_to_handles(cls.handles_to_delete))
            print('Deleted handles ' + str(cls.handles_to_delete))

    @classmethod
    def delete_shock_node(cls, node_id):
        header = {'Authorization': 'Oauth {0}'.format(cls.token)}
        requests.delete(cls.shockURL + '/node/' + node_id, headers=header,
                        allow_redirects=True)
        print('Deleted shock node ' + node_id)

    def createAnObject(self):
        if hasattr(self.__class__, 'fake_object_ref'):
            return self.__class__.fake_object_ref

        obj_name = 'test_obj.1'
        foft = FakeObjectsForTests(self.callback_url)
        info = foft.create_any_objects({'ws_name': self.wsName, 'obj_names': [obj_name]})[0]

        fake_object_ref = "%s/%s/%s" % (info[6], info[0], info[4])

        self.__class__.fake_object_ref = fake_object_ref
        print('Loaded Fake Object: ' + fake_object_ref)
        return fake_object_ref

    @classmethod
    def prepare_data(cls):
        workspace_id = cls.dfu.ws_name_to_id(cls.wsName)
        object_type = 'KBaseExperiments.AttributeMapping'
        attribute_mapping_object_name = 'test_attribute_mapping'
        attribute_mapping_data = json.load(open('data/biom_am.json'))
        save_object_params = {
            'id': workspace_id,
            'objects': [{'type': object_type,
                         'data': attribute_mapping_data,
                         'name': attribute_mapping_object_name}]
        }

        dfu_oi = cls.dfu.save_objects(save_object_params)[0]
        cls.attribute_mapping_ref = str(dfu_oi[6]) + '/' + str(dfu_oi[0]) + '/' + str(dfu_oi[4])

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        return self.__class__.wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    def start_test(self):
        testname = inspect.stack()[1][3]
        print('\n*** starting test: ' + testname + ' **')

    def mock_download_staging_file(params):
        print('Mocking DataFileUtilClient.download_staging_file')
        print(params)

        file_path = params.get('staging_file_subdir_path')

        return {'copy_file_path': file_path}

    def mock_file_to_shock(params):
        print('Mocking DataFileUtilClient.file_to_shock')
        print(params)

        return BioMultiTest().test_shock

    def loadSampleSet(self):
        if hasattr(self.__class__, 'sample_set_ref'):
            return self.__class__.sample_set_ref

        sample_set_file_name = 'sample_set_test.xls'
        sample_set_file_path = os.path.join(self.scratch, sample_set_file_name)
        shutil.copy(os.path.join('data', sample_set_file_name), sample_set_file_path)

        params = {
            'workspace_name': self.wsName,
            'workspace_id': self.wsId,
            'sample_file': sample_set_file_path,
            'file_format': "SESAR",
            'set_name': 'test1',
            'description': "this is a test sample set."
        }
        import_samples_rec = self.sample_uploader.import_samples(params)

        report_data = self.dfu.get_objects(
                    {"object_refs": [import_samples_rec['report_ref']]})['data'][0]['data']

        sample_set_ref = report_data['objects_created'][0]['ref']

        self.__class__.sample_set_ref = sample_set_ref
        print('Loaded SampleSet: ' + sample_set_ref)
        return sample_set_ref

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

    @patch.object(DataFileUtil, "download_staging_file", side_effect=mock_download_staging_file)
    @patch.object(DataFileUtil, "file_to_shock", side_effect=mock_file_to_shock)
    def test_import_matrix_from_biom_1_0_tsv_fasta(self, download_staging_file, file_to_shock):
        self.start_test()

        params = {'obj_type': 'AmpliconMatrix',
                  'matrix_name': 'test_AmpliconMatrix',
                  'workspace_id': self.wsId,
                  "tsv_fasta": {
                    "tsv_file_tsv_fasta": os.path.join('data', 'amplicon_test.tsv'),
                    "fasta_file_tsv_fasta": os.path.join('data', 'phyloseq_test.fa'),
                    'metadata_keys_tsv_fasta': 'taxonomy_id, taxonomy, taxonomy_source, consensus_sequence',
                    },
                  'scale': 'raw',
                  'description': "OTU data",
                  'amplicon_type': '16S',
                  'target_gene_region': 'V1',
                  'forward_primer_sequence': 'forward_primer_sequence',
                  'reverse_primer_sequence': 'reverse_primer_sequence',
                  'sequencing_platform': 'Illumina',
                  'sequencing_quality_filter_cutoff': 'sequencing_quality_filter_cutoff',
                  'clustering_cutoff': 0.3,
                  'clustering_method': 'clustering_method'
                  }
        returnVal = self.getImpl().import_matrix_from_biom(self.ctx, params)[0]
        self.assertIn('matrix_obj_ref', returnVal)
        self.assertIn('report_name', returnVal)
        self.assertIn('report_ref', returnVal)
        obj = self.dfu.get_objects(
            {'object_refs': [returnVal['matrix_obj_ref']]}
        )['data'][0]['data']
        self.assertIn('description', obj)
        self.assertEqual(obj['description'], 'OTU data')
        self.assertIn('row_attributemapping_ref', obj)
        self.assertNotIn('col_attributemapping_ref', obj)

    # @unittest.skip("narrative UI no longer support this option")
    @patch.object(DataFileUtil, "download_staging_file", side_effect=mock_download_staging_file)
    @patch.object(DataFileUtil, "file_to_shock", side_effect=mock_file_to_shock)
    def test_import_matrix_from_biom_1_0_biom_fasta(self, download_staging_file, file_to_shock):
        self.start_test()

        params = {'obj_type': 'AmpliconMatrix',
                  'matrix_name': 'test_AmpliconMatrix',
                  'workspace_id': self.wsId,
                  "biom_fasta": {
                        "biom_file_biom_fasta": os.path.join('data', 'phyloseq_test.biom'),
                        "fasta_file_biom_fasta": os.path.join('data', 'phyloseq_test.fa')
                        },
                  'scale': 'raw',
                  'description': "OTU data",
                  'amplicon_type': '16S',
                  'target_gene_region': 'V1',
                  'forward_primer_sequence': 'forward_primer_sequence',
                  'reverse_primer_sequence': 'reverse_primer_sequence',
                  'sequencing_platform': 'Illumina',
                  'sequencing_quality_filter_cutoff': 'sequencing_quality_filter_cutoff',
                  'clustering_cutoff': 0.3,
                  'clustering_method': 'clustering_method'
                  }
        returnVal = self.getImpl().import_matrix_from_biom(self.ctx, params)[0]
        self.assertIn('matrix_obj_ref', returnVal)
        self.assertIn('report_name', returnVal)
        self.assertIn('report_ref', returnVal)
        obj = self.dfu.get_objects(
            {'object_refs': [returnVal['matrix_obj_ref']]}
        )['data'][0]['data']
        self.assertIn('description', obj)
        self.assertEqual(obj['description'], 'OTU data')
        self.assertIn('attributes', obj)
        self.assertEqual(obj['attributes'], {'generated_by': 'QIIME revision XYZ'})
        self.assertIn('row_attributemapping_ref', obj)
        self.assertIn('col_attributemapping_ref', obj)

    @patch.object(DataFileUtil, "file_to_shock", side_effect=mock_file_to_shock)
    def test_import_matrix_from_tsv_fasta(self, file_to_shock):
        self.start_test()

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
        self.assertIn('matrix_obj_ref', returnVal)
        self.assertIn('report_name', returnVal)
        self.assertIn('report_ref', returnVal)
        obj = self.dfu.get_objects(
            {'object_refs': [returnVal['matrix_obj_ref']]}
        )['data'][0]['data']
        self.assertIn('description', obj)
        self.assertEqual(obj['description'], 'OTU data')
        self.assertIn('row_attributemapping_ref', obj)
        self.assertNotIn('col_attributemapping_ref', obj)

    @patch.object(DataFileUtil, "download_staging_file", side_effect=mock_download_staging_file)
    @patch.object(DataFileUtil, "file_to_shock", side_effect=mock_file_to_shock)
    def test_import_matrix_from_biom_1_0_tsv_fasta_with_sample_set(self, download_staging_file,
                                                                   file_to_shock):
        self.start_test()

        sample_set_ref = self.loadSampleSet()

        params = {'obj_type': 'AmpliconMatrix',
                  'matrix_name': 'test_AmpliconMatrix',
                  'workspace_id': self.wsId,
                  "tsv_fasta": {
                        "tsv_file_tsv_fasta": os.path.join('data', 'amplicon_test.tsv'),
                        "fasta_file_tsv_fasta": os.path.join('data', 'phyloseq_test.fa'),
                        'metadata_keys_tsv_fasta': 'taxonomy_id, taxonomy, taxonomy_source, consensus_sequence',
                        },
                  'scale': 'raw',
                  'description': "OTU data",
                  'sample_set_ref': sample_set_ref,
                  'amplicon_type': '16S',
                  'target_gene_region': 'V1',
                  'forward_primer_sequence': 'forward_primer_sequence',
                  'reverse_primer_sequence': 'reverse_primer_sequence',
                  'sequencing_platform': 'Illumina',
                  'sequencing_quality_filter_cutoff': 'sequencing_quality_filter_cutoff',
                  'clustering_cutoff': 0.3,
                  'clustering_method': 'clustering_method'
                  }
        returnVal = self.getImpl().import_matrix_from_biom(self.ctx, params)[0]
        self.assertIn('matrix_obj_ref', returnVal)
        self.assertIn('report_name', returnVal)
        self.assertIn('report_ref', returnVal)
        obj = self.dfu.get_objects(
            {'object_refs': [returnVal['matrix_obj_ref']]}
        )['data'][0]['data']
        self.assertIn('description', obj)
        self.assertEqual(obj['description'], 'OTU data')
        self.assertIn('row_attributemapping_ref', obj)
        self.assertIn('col_attributemapping_ref', obj)

    @patch.object(DataFileUtil, "file_to_shock", side_effect=mock_file_to_shock)
    def test_import_chemical_abundance_matrix_with_sample_set(self, file_to_shock):
        self.start_test()

        sample_set_ref = self.loadSampleSet()

        obj_type = 'ChemicalAbundanceMatrix'
        params = {'obj_type': obj_type,
                  'matrix_name': 'test_ChemicalAbundanceMatrix',
                  'workspace_name': self.wsName,
                  'input_file_path': os.path.join('data', 'metabolite_data.xlsx'),
                  'scale': 'log2',
                  'biochemistry_ref': 'kbase/default',
                  'description': "a biochem matrix",
                  'sample_set_ref': sample_set_ref,
                  }
        returnVal = self.serviceImpl.import_matrix_from_excel(self.ctx, params)[0]
        matrix_obj_ref = returnVal['matrix_obj_ref']

        self._check_matrix_to_sample_links(matrix_obj_ref, sample_set_ref)

    @patch.object(DataFileUtil, "download_staging_file", side_effect=mock_download_staging_file)
    @patch.object(DataFileUtil, "file_to_shock", side_effect=mock_file_to_shock)
    def test_import_amplicon_matrix_with_sample_set(self, download_staging_file, file_to_shock):
        self.start_test()

        sample_set_ref = self.loadSampleSet()

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
                  'input_local_file': True,
                  'sample_set_ref': sample_set_ref,
                  }
        returnVal = self.getImpl().import_matrix_from_biom(self.ctx, params)[0]
        matrix_obj_ref = returnVal['matrix_obj_ref']

        self._check_matrix_to_sample_links(matrix_obj_ref, sample_set_ref)

    def _check_matrix_to_sample_links(self, matrix_obj_ref, sample_set_ref):
        sample_set_obj = self.dfu.get_objects(
            {'object_refs': [sample_set_ref]}
        )['data'][0]['data']
        sample_infos = [{'id': d['id'], 'version': d['version']} for d in sample_set_obj['samples']]

        links_data = self.sample_ser.get_data_links_from_data({'upa': matrix_obj_ref})['links']
        links_sample = [
            link
            for sample_info in sample_infos
            for link in self.sample_ser.get_data_links_from_sample(sample_info)['links']
            if link['upa'] == matrix_obj_ref
        ]

        sorted_dicts = lambda l: sorted(l, key=lambda x: str(x))

        assert sorted_dicts(links_data) == sorted_dicts(links_sample), '%s vs %s' % (links_data, links_sample)
        assert len(links_data) == len(sample_set_obj['samples'])

    @patch.object(DataFileUtil, "download_staging_file", side_effect=mock_download_staging_file)
    @patch.object(DataFileUtil, "file_to_shock", side_effect=mock_file_to_shock)
    def test_import_with_external_am(self, download_staging_file, file_to_shock):
        self.start_test()

        biom_file_biom_fasta = os.path.join(self.scratch, 'phyloseq_test.biom')
        shutil.copy(os.path.join('data', 'phyloseq_test.biom'), biom_file_biom_fasta)

        fasta_file_biom_fasta = os.path.join(self.scratch, 'phyloseq_test.fa')
        shutil.copy(os.path.join('data', 'phyloseq_test.fa'), fasta_file_biom_fasta)

        params = {'obj_type': 'AmpliconMatrix',
                  'matrix_name': 'test_AmpliconMatrix',
                  'workspace_id': self.wsId,
                  "biom_fasta": {
                        "biom_file_biom_fasta": biom_file_biom_fasta,
                        "fasta_file_biom_fasta": fasta_file_biom_fasta
                        },
                  'scale': 'raw',
                  'description': "OTU data",
                  'col_attributemapping_ref': self.attribute_mapping_ref,
                  'amplicon_type': '16S',
                  'target_gene_region': 'V1',
                  'forward_primer_sequence': 'forward_primer_sequence',
                  'reverse_primer_sequence': 'reverse_primer_sequence',
                  'sequencing_platform': 'Illumina',
                  'sequencing_quality_filter_cutoff': 'sequencing_quality_filter_cutoff',
                  'clustering_cutoff': 0.3,
                  'clustering_method': 'clustering_method'
                  }

        returnVal = self.getImpl().import_matrix_from_biom(self.ctx, params)[0]
        self.assertIn('matrix_obj_ref', returnVal)
        self.assertIn('report_name', returnVal)
        self.assertIn('report_ref', returnVal)
        obj = self.dfu.get_objects(
            {'object_refs': [returnVal['matrix_obj_ref']]}
        )['data'][0]['data']
        self.assertIn('description', obj)
        self.assertEqual(obj['description'], 'OTU data')
        self.assertIn('attributes', obj)
        self.assertEqual(obj['attributes'], {'generated_by': 'QIIME revision XYZ'})
        self.assertIn('row_attributemapping_ref', obj)
        self.assertIn('col_attributemapping_ref', obj)
        self.assertEqual(obj['col_attributemapping_ref'], self.attribute_mapping_ref)

    @patch.object(DataFileUtil, "download_staging_file", side_effect=mock_download_staging_file)
    def test_bad_import_matrix_params(self, download_staging_file):
        self.start_test()

        with self.assertRaisesRegex(ValueError, "parameter is required, but missing"):
            params = {'obj_type': 'AmpliconMatrix',
                      'matrix_name': 'test_AmpliconMatrix',
                      'workspace_id': self.wsId,
                      "biom_fasta": {
                            "biom_file_biom_fasta": os.path.join('data', 'phyloseq_test.biom'),
                            "fasta_file_biom_fasta": os.path.join('data', 'phyloseq_test.fa')
                            },
                      'amplicon_type': '16S',
                      'target_gene_region': 'V1',
                      'forward_primer_sequence': 'forward_primer_sequence',
                      'reverse_primer_sequence': 'reverse_primer_sequence',
                      'sequencing_platform': 'Illumina',
                      'sequencing_quality_filter_cutoff': 'sequencing_quality_filter_cutoff',
                      'clustering_cutoff': 0.3,
                      'clustering_method': 'clustering_method'
                      }
            self.getImpl().import_matrix_from_biom(self.ctx, params)[0]

        with self.assertRaisesRegex(ValueError, "Unknown matrix object type"):
            params = {'obj_type': 'foo',
                      'matrix_name': 'test_AmpliconMatrix',
                      'workspace_id': self.wsId,
                      "biom_fasta": {
                            "biom_file_biom_fasta": os.path.join('data', 'phyloseq_test.biom'),
                            "fasta_file_biom_fasta": os.path.join('data', 'phyloseq_test.fa')
                            },
                      'scale': 'log2',
                      'amplicon_type': '16S',
                      'target_gene_region': 'V1',
                      'forward_primer_sequence': 'forward_primer_sequence',
                      'reverse_primer_sequence': 'reverse_primer_sequence',
                      'sequencing_platform': 'Illumina',
                      'sequencing_quality_filter_cutoff': 'sequencing_quality_filter_cutoff',
                      'clustering_cutoff': 0.3,
                      'clustering_method': 'clustering_method'
                      }
            self.getImpl().import_matrix_from_biom(self.ctx, params)[0]

        with self.assertRaisesRegex(ValueError, "Unknown scale type"):
            params = {'obj_type': 'AmpliconMatrix',
                      'matrix_name': 'test_AmpliconMatrix',
                      'workspace_id': self.wsId,
                      "biom_fasta": {
                            "biom_file_biom_fasta": os.path.join('data', 'phyloseq_test.biom'),
                            "fasta_file_biom_fasta": os.path.join('data', 'phyloseq_test.fa')
                            },
                      'scale': 'foo',
                      'amplicon_type': '16S',
                      'target_gene_region': 'V1',
                      'forward_primer_sequence': 'forward_primer_sequence',
                      'reverse_primer_sequence': 'reverse_primer_sequence',
                      'sequencing_platform': 'Illumina',
                      'sequencing_quality_filter_cutoff': 'sequencing_quality_filter_cutoff',
                      'clustering_cutoff': 0.3,
                      'clustering_method': 'clustering_method'
                      }
            self.getImpl().import_matrix_from_biom(self.ctx, params)[0]

        with self.assertRaisesRegex(ValueError, "IDs from the uploaded matrix do not match"):
            params = {'obj_type': 'AmpliconMatrix',
                      'matrix_name': 'test_AmpliconMatrix',
                      'workspace_id': self.wsId,
                      "biom_fasta": {
                            "biom_file_biom_fasta": os.path.join('data', 'phyloseq_test.biom'),
                            "fasta_file_biom_fasta": os.path.join('data', 'phyloseq_test.fa')
                            },
                      'scale': 'raw',
                      'description': "OTU data",
                      'row_attributemapping_ref': self.attribute_mapping_ref,
                      'amplicon_type': '16S',
                      'target_gene_region': 'V1',
                      'forward_primer_sequence': 'forward_primer_sequence',
                      'reverse_primer_sequence': 'reverse_primer_sequence',
                      'sequencing_platform': 'Illumina',
                      'sequencing_quality_filter_cutoff': 'sequencing_quality_filter_cutoff',
                      'clustering_cutoff': 0.3,
                      'clustering_method': 'clustering_method'
                      }
            self.getImpl().import_matrix_from_biom(self.ctx, params)[0]

    def test_fetch_sequence(self):
        self.start_test()

        fake_object_ref = self.createAnObject()
        amplicon_matrix_ref = self.loadAmpliconMatrix()

        with self.assertRaisesRegex(ValueError, "Unexpected data type"):
            self.getImpl().fetch_sequence(self.ctx, fake_object_ref)[0]

        fasta_file_path = self.getImpl().fetch_sequence(self.ctx, amplicon_matrix_ref)[0]
        fastq_dict = SeqIO.index(fasta_file_path, "fasta")

        expected_seq_ids = ['GG_OTU_1', 'GG_OTU_2', 'GG_OTU_3', 'GG_OTU_4', 'GG_OTU_5', 'GG_OTU_6']
        self.assertCountEqual(expected_seq_ids, list(fastq_dict.keys()))

        expected_seq = 'ATCGATCGATCGTACGATCG'
        for seq_id in expected_seq_ids:
            self.assertEqual(expected_seq, str(fastq_dict.get(seq_id).seq))
