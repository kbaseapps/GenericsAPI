# -*- coding: utf-8 -*-
import os
import time
import unittest
from unittest.mock import patch
from configparser import ConfigParser
from os import environ
import functools
import numpy as np

from installed_clients.DataFileUtilClient import DataFileUtil
from GenericsAPI.GenericsAPIImpl import GenericsAPI
from GenericsAPI.GenericsAPIServer import MethodContext
from GenericsAPI.authclient import KBaseAuth as _KBaseAuth
from installed_clients.WorkspaceClient import Workspace as workspaceService

dec = '###' * 5
skipped_tests = []


def tag_kb_env(e):
    @functools.wraps(f)  # preserve wrapped's function name
    def decorator(f):
        f.kb_env = e
        return f
    return decorator


def skip_cond(select_run=None, kb_env=None):
    def decorator(f):
        @functools.wraps(f)  # preserve wrapped's function name
        def f_new(self, *a, **kw):
            if kb_env is not None and not hasattr(f, "kb_env"):
                raise Exception("Tag function (e.g., @tag_kb_env('ci')) with kb_env first to skip with this feature")

            if kb_env is not None and kb_env != f.kb_env:
                skipped_tests.append(f.__name__)
                self.skipTest("Test does not operate in this KBase environment")
            if select_run is not None and f.__name__ not in select_run:
                skipped_tests.append(f.__name__)
                print(dec, 'Skipping test %s because not in select_run' % f.__name__, dec)
                self.skipTest("Test is not in list of select_run")
            print(dec, 'Running test %s' % f.__name__, dec)
            f(self, *a, **kw)
        return f_new
    return decorator


select_run = None  # ['test_rarefy_defaultSubsample']


class MatrixUtilTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        print('In setUpClass')
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

        cls.dfu = DataFileUtil(cls.callback_url)

        suffix = int(time.time() * 1000)
        cls.wsName = "test_GenericsAPI_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})
        cls.wsId = ret[0]

        cls.prepare_data()

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')
        ran_tests = list(set(all_tests) - set(skipped_tests))
        print('All %d tests: %s' % (len(all_tests), all_tests))
        print('Skipped %d tests: %s' % (len(skipped_tests), skipped_tests))
        print('Ran %d tests: %s' % (len(ran_tests), ran_tests))

    def shortDescription(self):
        '''Override unittest using test*() docstrings in lieu of test*() method name in output summary'''
        return None

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

    @classmethod
    @patch.object(DataFileUtil, "download_staging_file", side_effect=mock_download_staging_file)
    def loadAmpliconMatrix(cls, download_staging_file):
        if hasattr(cls, 'amplicon_matrix_ref'):
            return

        print('Executing loadAmpliconMatrix')

        params = {'obj_type': 'AmpliconMatrix',
                  'matrix_name': 'test_AmpliconMatrix',
                  'workspace_name': cls.wsName,
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
        returnVal = cls.serviceImpl.import_matrix_from_biom(cls.ctx, params)[0]

        cls.amplicon_matrix_ref = returnVal['matrix_obj_ref']

        print('Loaded AmpliconMatrix: ' + cls.amplicon_matrix_ref)

    @classmethod
    def prepare_data(cls):
        cls.loadAmpliconMatrix()

        # the toy matrix loaded with patched self.loadAmpliconMatrix
        # sample names are ['Sample1', 'Sample2', 'Sample3', 'Sample4', 'Sample5', 'Sample6']
        # sample sizes are [7, 3, 4, 6, 5, 2]
        cls.matrix = [
            [0,0,1,0,0,0],
            [5,1,0,2,3,1],
            [0,0,1,4,2,0],
            [2,1,1,0,0,1],
            [0,1,1,0,0,0]
        ]

        cls.matrix_subsample5 = [ # rarefying with seed 7, subsample 5
            [0,0,1,0,0,0],        # confirm with `set.seed(7); t(rrarefy(t(m), 5))
            [4,1,0,2,3,1],
            [0,0,1,3,2,0],
            [1,1,1,0,0,1],
            [0,1,1,0,0,0]
        ]

        cls.matrix_subsample2 = [ # rarefying with seed 7, subsample 2
            [0,0,0,0,0,0],        # confirm with `set.seed(7); t(rrarefy(t(m), 2))
            [2,0,0,1,2,1],
            [0,0,1,1,0,0],
            [0,1,1,0,0,1],
            [0,1,0,0,0,0]
        ]

        cls.matrix_bootstrap9Reps_median = [ # rarefying with seed 7, subsample 5
            [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
            [4.0, 1.0, 0.0, 2.0, 3.0, 1.0],
            [0.0, 0.0, 1.0, 3.0, 2.0, 0.0],
            [1.0, 1.0, 1.0, 0.0, 0.0, 1.0],
            [0.0, 1.0, 1.0, 0.0, 0.0, 0.0]]

        cls.matrix_bootstrap9Reps_mean = np.array([ # rarefying with seed 7, subsample 5
            [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
            [3.7777777777777777, 1.0, 0.0, 1.5555555555555556, 3.0, 1.0],
            [0.0, 0.0, 1.0, 3.4444444444444446, 2.0, 0.0],
            [1.2222222222222223, 1.0, 1.0, 0.0, 0.0, 1.0],
            [0.0, 1.0, 1.0, 0.0, 0.0, 0.0]
        ])

        ## In logging, look for 'Start generating html report in xxx' to find html report ##
        ## Count "In setUpClass" and "Executing loadAmpliconMatrix" -> once each

    def get_out_data(self, ret, matrix_out=True):
        report_obj = self.dfu.get_objects({'object_refs': [ret[0]['report_ref']]})['data'][0]['data']
        warnings = report_obj['warnings']
        matrix_out = self.dfu.get_objects({
            'object_refs': [report_obj['objects_created'][0]['ref']]
        })['data'][0]['data']['data']['values']

        return warnings, matrix_out

    def assert_matrices_equal(self, m1, m2):
        self.assertTrue(
            m1 == m2,
            'm1 is\n`%s`\nm2 is\n`%s`' % (m1, m2)
        )

    ##########
    ##########
    @skip_cond(select_run=select_run)
    def test_rarefy_defaultSubsample(self):
        '''
        All samples should be rarefied to least sample size, 2
        No warnings about anything too big to rarefy
        Check output matrix against test data matrix
        '''

        ret = self.getImpl().perform_rarefy(
            self.ctx, {
                'workspace_id': self.getWsId(),
                'input_matrix_ref': self.amplicon_matrix_ref,
                'dimension': 'col',
                'seed_number': 7,
                'subsample_size': None,
                'new_matrix_name': 'new_matrix_name',
            })

        warnings, matrix_out = self.get_out_data(ret)

        self.assertTrue(len(warnings) == 0, 'length is %d' % len(warnings))
        self.assert_matrices_equal(matrix_out, self.matrix_subsample2)

    ##########
    ##########
    @skip_cond(select_run=select_run)
    def test_rarefy_medSubsample(self):
        '''
        Some samples should and should not be rarefied
        That means warnings for the ones too big to rarefy
        Check output matrix against test data matrix
        '''

        ret = self.getImpl().perform_rarefy(
            self.ctx, {
                'workspace_id': self.getWsId(),
                'input_matrix_ref': self.amplicon_matrix_ref,
                'dimension': 'col',
                'seed_number': 7,
                'subsample_size': 5,
                'new_matrix_name': 'new_matrix_name',
            })

        warnings, matrix_out = self.get_out_data(ret)

        self.assertTrue(len(warnings) == 1, 'length is %d' % len(warnings))
        msg = (
            "At subsampling size 5, samples ['Sample2', 'Sample3', 'Sample6'] are too small "
            "and will not be rarefied. Smallest sample size is 2")
        self.assertTrue(warnings[0] == msg, 'warnings[0] is\n`%s`, msg is\n`%s`' % (warnings[0], msg))

        self.assert_matrices_equal(matrix_out, self.matrix_subsample5)

    ##########
    ##########
    @skip_cond(select_run=select_run)
    def test_rarefy_medSubsample_bootstrap9Reps(self):
        '''
        At subsample 5 and 9 bootstrap reps
        Check output matrices against test data matrices
        '''
        ret = self.getImpl().perform_rarefy(
            self.ctx, {
                'workspace_id': self.getWsId(),
                'input_matrix_ref': self.amplicon_matrix_ref,
                'dimension': 'col',
                'seed_number': 7,
                'subsample_size': 5,
                'bootstrap': {
                    'num_rare_reps': 9,
                    'central_tendency': 'median',
                },
                'new_matrix_name': 'new_matrix_name',
            })

        _, matrix_out_median = self.get_out_data(ret)

        ret = self.getImpl().perform_rarefy(
            self.ctx, {
                'workspace_id': self.getWsId(),
                'input_matrix_ref': self.amplicon_matrix_ref,
                'dimension': 'col',
                'seed_number': 7,
                'subsample_size': 5,
                'bootstrap': {
                    'num_rare_reps': 9,
                    'central_tendency': 'mean',
                },
                'new_matrix_name': 'new_matrix_name',
            })

        _, matrix_out_mean = self.get_out_data(ret)

        self.assert_matrices_equal(matrix_out_median, self.matrix_bootstrap9Reps_median)
        self.assertTrue(
            np.all(np.abs(np.array(matrix_out_mean) - self.matrix_bootstrap9Reps_mean) < 1e-10)  # allow some error for mean
        )

    ##########
    ##########
    @skip_cond(select_run=select_run)
    def test_rarefy_largeSubsample(self):
        '''
        All samples, bootstrapped or not, should not be rarefied and should be returned as is
        All resulting matrices in this test should be same
        '''

        # run perform_rarefy with no bootstrap
        ret = self.getImpl().perform_rarefy(
            self.ctx, {
                'workspace_id': self.getWsId(),
                'input_matrix_ref': self.amplicon_matrix_ref,
                'dimension': 'col',
                'seed_number': 7,
                'subsample_size': 1000,
                'new_matrix_name': 'new_matrix_name',
            })

        warnings, matrix_out_noBootstrap = self.get_out_data(ret)

        self.assertTrue(len(warnings) == 1, 'length is %d' % len(warnings))
        msg = (
            "At subsampling size 1000, samples ['Sample1', 'Sample2', 'Sample3', 'Sample4', "
            "'Sample5', 'Sample6'] are too small "
            "and will not be rarefied. Smallest sample size is 2")
        self.assertTrue(warnings[0] == msg, 'warnings[0] is\n`%s`, msg is\n`%s`' % (warnings[0], msg))

        # run perform_rarefy with bootstrap with rep=1 and central_tendency=median
        ret = self.getImpl().perform_rarefy(
            self.ctx, {
                'workspace_id': self.getWsId(),
                'input_matrix_ref': self.amplicon_matrix_ref,
                'dimension': 'col',
                'seed_number': 7,
                'subsample_size': 1000,
                'bootstrap': {
                    'num_rare_reps': 1,
                    'central_tendency': 'median',
                },
                'new_matrix_name': 'new_matrix_name',
            })

        _, matrix_out_median = self.get_out_data(ret)

        # run perform_rarefy with bootstrap with rep=1 and central_tendency=mean
        ret = self.getImpl().perform_rarefy(
            self.ctx, {
                'workspace_id': self.getWsId(),
                'input_matrix_ref': self.amplicon_matrix_ref,
                'dimension': 'col',
                'seed_number': 7,
                'subsample_size': 1000,
                'bootstrap': {
                    'num_rare_reps': 1,
                    'central_tendency': 'mean',
                },
                'new_matrix_name': 'new_matrix_name',
            })

        _, matrix_out_mean = self.get_out_data(ret)

        # these resulting 3 matrices should match
        # since seeded same
        # and ran just 1nce
        self.assertTrue(matrix_out_noBootstrap, matrix_out_median)
        self.assertTrue(matrix_out_median, matrix_out_mean)


all_tests = [k for k, v in MatrixUtilTest.__dict__.items() if k.startswith('test') and callable(v)]
