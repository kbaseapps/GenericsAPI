# -*- coding: utf-8 -*-
import inspect
import json
import os
import time
import unittest
from unittest.mock import patch
from configparser import ConfigParser
from os import environ
from mock import patch
import shutil
import sys

from installed_clients.DataFileUtilClient import DataFileUtil
from GenericsAPI.GenericsAPIImpl import GenericsAPI
from GenericsAPI.GenericsAPIServer import MethodContext
from GenericsAPI.authclient import KBaseAuth as _KBaseAuth
from installed_clients.GenomeFileUtilClient import GenomeFileUtil
from installed_clients.WorkspaceClient import Workspace as workspaceService


class MatrixUtilTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = environ.get('KB_AUTH_TOKEN', None)
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
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

        cls.gfu = GenomeFileUtil(cls.callback_url)
        cls.dfu = DataFileUtil(cls.callback_url)

        suffix = int(time.time() * 1000)
        cls.wsName = "test_GenericsAPI_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})
        cls.wsId = ret[0]
        #cls.prepare_data()

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        return self.__class__.wsName

    def getWsId(self):
        return self.__class__.wsId

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    def mock_download_staging_file(params):
        print('Mocking DataFileUtilClient.download_staging_file')
        print(params)

        file_path = params.get('staging_file_subdir_path')

        return {'copy_file_path': file_path}

    @patch.object(DataFileUtil, "download_staging_file", side_effect=mock_download_staging_file)
    def loadAmpliconMatrix(self, download_staging_file):
        if hasattr(self.__class__, 'amplicon_matrix_ref'):
            return self.__class__.amplicon_matrix_ref

        params = {'obj_type': 'AmpliconMatrix',
                  'matrix_name': 'test_AmpliconMatrix',
                  'workspace_name': self.wsName,
                  "biom_fasta": {
                        "biom_file_biom_fasta": os.path.join('data', 'phyloseq_test.biom'),
                        "fasta_file_biom_fasta": os.path.join('data', 'phyloseq_test.fa')
                        },
                  'scale': 'raw',
                  'description': "OTU data",
                  'amplicon_set_name': 'test_AmpliconSet',
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

        amplicon_matrix_ref = returnVal['matrix_obj_ref']

        self.__class__.amplicon_matrix_ref = amplicon_matrix_ref
        print('Loaded AmpliconMatrix: ' + amplicon_matrix_ref)
        return amplicon_matrix_ref


    def test_perform_rarefy_default_size(self):
        # all samples should be rarefied
        # in logging, look for 'Start generating html report in xxx' to find html report
        ampmat_ref = self.loadAmpliconMatrix()

        ret = self.getImpl().perform_rarefy(
            self.ctx, {
                'workspace_id': self.getWsId(),
                'input_matrix_ref': ampmat_ref,
                'dimension': 'col',
                'seed_number': 7,
                'subsample_size': None,
                'new_matrix_name': 'new_matrix_name',
            })

        report_ref = ret[0]['report_ref']
        report_obj = self.dfu.get_objects({'object_refs': [report_ref]})['data'][0]['data']

        warnings = report_obj['warnings']
        self.assertTrue(len(warnings) == 0, 'length is %d' % len(warnings))

    def test_perform_rarefy_med_large_size(self):
        # some samples should not be rarefied
        # Sample names are ['Sample1', 'Sample2', 'Sample3', 'Sample4', 'Sample5', 'Sample6']
        # Sample sizes are [7, 3, 4, 6, 5, 2]
        ampmat_ref = self.loadAmpliconMatrix()

        ret = self.getImpl().perform_rarefy(
            self.ctx, {
                'workspace_id': self.getWsId(),
                'input_matrix_ref': ampmat_ref,
                'dimension': 'col',
                'seed_number': 7,
                'subsample_size': 5,
                'new_matrix_name': 'new_matrix_name',
            })

        report_ref = ret[0]['report_ref']
        report_obj = self.dfu.get_objects({'object_refs': [report_ref]})['data'][0]['data']

        warnings = report_obj['warnings']
        self.assertTrue(len(warnings) == 1, 'length is %d' % len(warnings))
        msg = (
            "At subsampling size 5, samples ['Sample2', 'Sample3', 'Sample6'] are too small "
            "and will not be rarefied. Smallest sample size is 2")
        self.assertTrue(warnings[0] == msg, 'warnings[0] is\n`%s`, msg is\n`%s`' % (warnings[0], msg))

    def test_perform_rarefy_very_large_size(self):
        # all samples should not be rarefied
        ampmat_ref = self.loadAmpliconMatrix()

        ret = self.getImpl().perform_rarefy(
            self.ctx, {
                'workspace_id': self.getWsId(),
                'input_matrix_ref': ampmat_ref,
                'dimension': 'col',
                'seed_number': 7,
                'subsample_size': 1000,
                'new_matrix_name': 'new_matrix_name',
            })

        report_ref = ret[0]['report_ref']
        report_obj = self.dfu.get_objects({'object_refs': [report_ref]})['data'][0]['data']

        warnings = report_obj['warnings']
        self.assertTrue(len(warnings) == 1, 'length is %d' % len(warnings))
        msg = (
            "At subsampling size 1000, samples ['Sample1', 'Sample2', 'Sample3', 'Sample4', "
            "'Sample5', 'Sample6'] are too small "
            "and will not be rarefied. Smallest sample size is 2")
        self.assertTrue(warnings[0] == msg, 'warnings[0] is\n`%s`, msg is\n`%s`' % (warnings[0], msg))


    def test_perform_rarefy_correctness(self):
        # TODO test correctness on larger matrix
        pass
