# -*- coding: utf-8 -*-
import inspect
import json
import os
import time
import unittest
from unittest import TestCase
from unittest.mock import patch, create_autospec, call
from configparser import ConfigParser
from os import environ
from mock import patch
import shutil
import sys
import functools
import re

import pandas as pd
import numpy as np

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.SampleServiceClient import SampleService
from GenericsAPI.GenericsAPIImpl import GenericsAPI
from GenericsAPI.Utils.MatrixUtil import MatrixUtil
from GenericsAPI.Utils.MatrixValidation import MatrixValidationException
from GenericsAPI.Utils import MatrixValidation as vd
from GenericsAPI.GenericsAPIServer import MethodContext
from GenericsAPI.authclient import KBaseAuth as _KBaseAuth
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

        #cls.gfu = GenomeFileUtil(cls.callback_url)
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

    def shortDescription(self):
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
        if hasattr(cls, 'amplicon_matrix_ref'): return

        print('Executing loadAmpliconMatrix')

        biom_file_biom_fasta = os.path.join(cls.scratch, 'phyloseq_test.biom')
        shutil.copy(os.path.join('data', 'phyloseq_test.biom'), biom_file_biom_fasta)

        fasta_file_biom_fasta = os.path.join(cls.scratch, 'phyloseq_test.fa')
        shutil.copy(os.path.join('data', 'phyloseq_test.fa'), fasta_file_biom_fasta)

        params = {'obj_type': 'AmpliconMatrix',
                  'matrix_name': 'test_AmpliconMatrix',
                  'workspace_id': cls.wsId,
                  "biom_fasta": {
                        "biom_file_biom_fasta": biom_file_biom_fasta,
                        "fasta_file_biom_fasta": fasta_file_biom_fasta
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
        cls.loadAmpliconMatrix() # TODO replace with mocking in tests

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

        cls.matrix_bootstrap9Reps_mean = [ # rarefying with seed 7, subsample 5
            [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
            [3.7777777777777777, 1.0, 0.0, 1.5555555555555556, 3.0, 1.0],
            [0.0, 0.0, 1.0, 3.4444444444444446, 2.0, 0.0],
            [1.2222222222222223, 1.0, 1.0, 0.0, 0.0, 1.0],
            [0.0, 1.0, 1.0, 0.0, 0.0, 0.0]
        ]


    def get_out_data(self, ret, matrix_out=True, attri_out=False):
        report_obj = self.dfu.get_objects({'object_refs': [ret[0]['report_ref']]})['data'][0]['data']
        warnings = report_obj['warnings']

        if matrix_out or attri_out:
            obj_data = self.dfu.get_objects({
                'object_refs': [report_obj['objects_created'][0]['ref']]
            })['data'][0]['data']

            if attri_out:
                row_attri = obj_data.get('row_attributemapping_ref')
                col_attri = obj_data.get('col_attributemapping_ref')

                row_attri_size = 0
                if row_attri:
                    row_attri_obj = self.dfu.get_objects({'object_refs': [row_attri]})['data'][0]['data']
                    row_attri_size = len(row_attri_obj['instances'])

                col_attri_size = 0
                if col_attri:
                    col_attri_obj = self.dfu.get_objects({'object_refs': [col_attri]})['data'][0]['data']
                    col_attri_size = len(col_attri_obj['instances'])

            if matrix_out:
                matrix_out = obj_data['data']['values']

                if attri_out:
                    return warnings, matrix_out, row_attri_size, col_attri_size
                else:
                    return warnings, matrix_out

        return warnings

    def assert_matrices_equal(self, m1, m2):
        if not isinstance(m1, np.ndarray): m1 = np.array(m1)
        if not isinstance(m2, np.ndarray): m2 = np.array(m2)

        self.assertTrue(
            np.allclose(m1, m2) and np.allclose(m2, m1),
            'm1 is\n`%s`\nm2 is\n`%s`' % (m1, m2)
        )

    '''
    {
        "input_matrix_ref": "test_amplicon_matrix",
        "operations": [
            "abundance_filtering",
            "relative_abundance",
            "standardization",
            "ratio_transformation",
            "logit",
            "sqrt",
            "log"
        ],
        "abundance_filtering_params": {
            "abundance_filtering_row_threshold": 19,
            "abundance_filtering_columns_threshold": 0,
            "abundance_filtering_row_sum_threshold": 125,
            "abundance_filtering_columns_sum_threshold": 3500
        },
        "perform_relative_abundance": 1,
        "standardization_params": {
            "standardization_with_mean": 0,
            "standardization_with_std": 1,
            "standardization_dimension": "col"
        },
        "ratio_transformation_params": {
            "ratio_transformation_method": "clr",
            "ratio_transformation_dimension": "col"
        },
        "log_params": {
            "base": 10,
            "offset": 1e-10
        }
        "new_matrix_name": "test_amplicon_transformed_matrix"
    }
    '''

    def test_transform_pipeline(self):

        with self.subTest(): # TODO subTest not catching?
            '''
            '''
            ret = self.getImpl().transform_matrix(
                self.ctx, {
                    'workspace_id': self.getWsId(),
                    'input_matrix_ref': self.amplicon_matrix_ref,
                    'operations': [
                        'abundance_filtering',
                        "log",
                    ],
                    "abundance_filtering_params": {
                        "abundance_filtering_row_threshold": 0,
                        "abundance_filtering_columns_threshold": 0,
                        "abundance_filtering_row_sum_threshold": 0,
                        "abundance_filtering_columns_sum_threshold": 3
                    },
                    "log_params": {
                        "log_base": 2.718281828459045,
                        "log_offset": 1e-10
                    }
                })

            out1 = [
                [-2.30258509e+01,  1.00000008e-10, -2.30258509e+01, -2.30258509e+01],
                [ 1.60943791e+00, -2.30258509e+01,  6.93147181e-01,  1.09861229e+00],
                [-2.30258509e+01,  1.00000008e-10,  1.38629436e+00,  6.93147181e-01],
                [ 6.93147181e-01,  1.00000008e-10, -2.30258509e+01, -2.30258509e+01],
                [-2.30258509e+01,  1.00000008e-10, -2.30258509e+01, -2.30258509e+01]
            ]

            _, out2, row_attri_size, col_attri_size = self.get_out_data(ret, matrix_out=True, attri_out=True)

            self.assert_matrices_equal(out1, out2)
            self.assertEqual(row_attri_size, 5)
            self.assertEqual(col_attri_size, 4)

        with self.subTest():
            '''
            Sqrt
            '''
            ret = self.getImpl().transform_matrix(
                self.ctx, {
                    'workspace_id': self.getWsId(),
                    'input_matrix_ref': self.amplicon_matrix_ref,
                    'operations': [
                        'abundance_filtering',
                        "relative_abundance", # fixed axis='col' currently
                        "sqrt",
                    ],
                    "abundance_filtering_params": {
                        "abundance_filtering_row_threshold": 0,
                        "abundance_filtering_columns_threshold": 0,
                        "abundance_filtering_row_sum_threshold": 2,
                        "abundance_filtering_columns_sum_threshold": 0
                    },
                })

            out1 = [
                [ 0.845154,  0.707107,         0,   0.57735,  0.774597,  0.707107],
                [        0,         0,  0.707107,  0.816497,  0.632456,         0],
                [ 0.534522,  0.707107,  0.707107,         0,         0,  0.707107]
            ]

            _, out2, row_attri_size, col_attri_size = self.get_out_data(ret, matrix_out=True, attri_out=True)

            self.assert_matrices_equal(out1, out2)
            self.assertEqual(row_attri_size, 3)
            self.assertEqual(col_attri_size, 6)


        with self.subTest():
            '''
            Ratio transform
            '''
            ret = self.getImpl().transform_matrix(
                self.ctx, {
                    'workspace_id': self.getWsId(),
                    'input_matrix_ref': self.amplicon_matrix_ref,
                    'operations': [
                        'abundance_filtering',
                        "relative_abundance",
                        'ratio_transformation',
                    ],
                    "abundance_filtering_params": {
                        "abundance_filtering_row_threshold": 0,
                        "abundance_filtering_columns_threshold": 0,
                        "abundance_filtering_row_sum_threshold": 0,
                        "abundance_filtering_columns_sum_threshold": 0
                    },
                    "ratio_transformation_params": {
                        "ratio_transformation_method": "clr",
                        "ratio_transformation_dimension": "col"
                    }
                })

            # nan/inf


        with self.subTest():
            '''
            '''
            ret = self.getImpl().transform_matrix(
                self.ctx, {
                    'workspace_id': self.getWsId(),
                    'input_matrix_ref': self.amplicon_matrix_ref,
                    'operations': [
                        'abundance_filtering',
                        "relative_abundance",
                        'log',
                    ],
                    "abundance_filtering_params": {
                        "abundance_filtering_row_threshold": 0,
                        "abundance_filtering_columns_threshold": 0,
                        "abundance_filtering_row_sum_threshold": 0,
                        "abundance_filtering_columns_sum_threshold": 0
                    },
                    "log_params": {
                        "log_base": 2,
                        "log_offset": 1
                    },
                    "standardization_params": {
                        "standardization_with_mean": 1,
                        "standardization_with_std": 1,
                        "standardization_dimension": "col"
                    },
                })

            out1 = [
                [0.0, 0.0, 0.32192809488736235, 0.0, 0.0, 0.0],
                [0.7776075786635522, 0.4150374992788437, 0.0, 0.4150374992788437, 0.6780719051126378, 0.5849625007211562],
                [0.0, 0.0, 0.32192809488736235, 0.7369655941662061, 0.4854268271702417, 0.0],
                [0.36257007938470814, 0.4150374992788437, 0.32192809488736235, 0.0, 0.0, 0.5849625007211562],
                [0.0, 0.4150374992788437, 0.32192809488736235, 0.0, 0.0, 0.0]
            ]

            _, out2 = self.get_out_data(ret)
            self.assert_matrices_equal(out1, out2)

        with self.subTest():
            '''
            Logit
            '''
            ret = self.getImpl().transform_matrix(
                self.ctx, {
                    'workspace_id': self.getWsId(),
                    'input_matrix_ref': self.amplicon_matrix_ref,
                    'operations': [
                        'log', # make positive
                        "relative_abundance", # scale
                        'logit', # requires domain (0,1)
                    ],
                    "abundance_filtering_params": {
                        "abundance_filtering_row_threshold": 0,
                        "abundance_filtering_columns_threshold": 0,
                        "abundance_filtering_row_sum_threshold": 0,
                        "abundance_filtering_columns_sum_threshold": 0
                    },
                    "log_params": {
                        "log_base": 2,
                        "log_offset": 3
                    },
                    "standardization_params": {
                        "standardization_with_mean": 1,
                        "standardization_with_std": 1,
                        "standardization_dimension": "col"
                    },
                    "ratio_transformation_params": {
                        "ratio_transformation_method": "clr",
                        "ratio_transformation_dimension": "col"
                    },
                })

            out1 = [
                [-1.6785465 , -1.56560692, -1.33302049, -1.65559934, -1.62843694, -1.50933445],
                [-0.85821174, -1.27674801, -1.61888079, -1.18076985, -1.00711303, -1.21711914],
                [-1.6785465 , -1.56560692, -1.33302049, -0.9245813 , -1.15092049, -1.50933445],
                [-1.20592537, -1.27674801, -1.33302049, -1.65559934, -1.62843694, -1.21711914],
                [-1.6785465 , -1.27674801, -1.33302049, -1.65559934, -1.62843694, -1.50933445]
            ]

            _, out2 = self.get_out_data(ret)

            self.assert_matrices_equal(out1, out2)


        with self.subTest():
            '''
            Random
            '''
            ret = self.getImpl().transform_matrix(
                self.ctx, {
                    'workspace_id': self.getWsId(),
                    'input_matrix_ref': self.amplicon_matrix_ref,
                    'operations': [
                        'abundance_filtering',
                        'relative_abundance',
                        'ratio_transformation',
                    ],
                    'abundance_filtering_params': {
                        'abundance_filtering_row_threshold': -1,
                        'abundance_filtering_columns_threshold': -1,
                        'abundance_filtering_row_sum_threshold': -1,
                        'abundance_filtering_columns_sum_threshold': -1,
                    }
                })

            # numericized inf


    def test_transform_pipeline_throws(self):

        with self.assertRaises(MatrixValidationException) as cm:
            '''
            Logit, out of domain
            '''
            ret = self.getImpl().transform_matrix(
                self.ctx, {
                    'workspace_id': self.getWsId(),
                    'input_matrix_ref': self.amplicon_matrix_ref,
                    'operations': [
                        'abundance_filtering',
                        "relative_abundance",
                        'logit', # throws, requires domain (0,1)
                    ],
                    "abundance_filtering_params": {
                        "abundance_filtering_row_threshold": 0,
                        "abundance_filtering_columns_threshold": 0,
                        "abundance_filtering_row_sum_threshold": 0,
                        "abundance_filtering_columns_sum_threshold": 0
                    },
                })
        print(cm.exception)

        with self.assertRaises(Exception) as cm:
            '''
            singular operation
            '''
            ret = self.getImpl().transform_matrix(
                self.ctx, {
                    'workspace_id': self.getWsId(),
                    'input_matrix_ref': self.amplicon_matrix_ref,
                    'operations': [
                        'ratio_transformation',
                        'standardization',
                    ],
                })

        with self.assertRaises(Exception) as cm:
            '''
            Unknown op
            '''
            ret = self.getImpl().transform_matrix(
                self.ctx, {
                    'workspace_id': self.getWsId(),
                    'input_matrix_ref': self.amplicon_matrix_ref,
                    'operations': [
                        'abundance_filtering',
                        "relative_abundance",
                        'fake_op' #
                    ],
                    "abundance_filtering_params": {
                        "abundance_filtering_row_threshold": 0,
                        "abundance_filtering_columns_threshold": 0,
                        "abundance_filtering_row_sum_threshold": 0,
                        "abundance_filtering_columns_sum_threshold": 0
                    },
                })
        print(cm.exception)

        with self.assertRaises(Exception):
            '''
            Try to save matrix with inf
            '''
            ret = self.getImpl().transform_matrix(
                self.ctx, {
                    'workspace_id': self.getWsId(),
                    'input_matrix_ref': self.amplicon_matrix_ref,
                    'operations': [
                        'log' #
                    ],
                    "log_params": {
                        "log_base": 10,
                        "log_offset": 0
                    },
                })


    def test_transform_unit_op(self):
        '''
        '''
        mu = MatrixUtil
        m = np.array(self.matrix)
        df = pd.DataFrame(m)


        ## Logit ##

        out1 = [
            [-2.19722458, -2.19722458, -1.38629436, -2.19722458, -2.19722458, -2.19722458],
            [ 0.40546511, -1.38629436, -2.19722458, -0.84729786, -0.40546511, -1.38629436],
            [-2.19722458, -2.19722458, -1.38629436,  0.        , -0.84729786, -2.19722458],
            [-0.84729786, -1.38629436, -1.38629436, -2.19722458, -2.19722458, -1.38629436],
            [-2.19722458, -1.38629436, -1.38629436, -2.19722458, -2.19722458, -2.19722458]
        ]
        out2 = mu._logit((df+1) / 10)

        self.assert_matrices_equal(out1, out2)


        ## Sqrt ##

        out1 = [
            [0.        , 0.        , 1.        , 0.        , 0.        , 0.        ],
            [2.23606798, 1.        , 0.        , 1.41421356, 1.73205081, 1.        ],
            [0.        , 0.        , 1.        , 2.        , 1.41421356, 0.        ],
            [1.41421356, 1.        , 1.        , 0.        , 0.        , 1.        ],
            [0.        , 1.        , 1.        , 0.        , 0.        , 0.        ]
        ]
        out2 = mu._sqrt(df)

        self.assert_matrices_equal(out1, out2)

        ## Log ##

        out1 = [
            [0.        , 0.        , 0.30103   , 0.        , 0.        , 0.        ],
            [0.77815125, 0.30103   , 0.        , 0.47712125, 0.60205999, 0.30103   ],
            [0.        , 0.        , 0.30103   , 0.69897   , 0.47712125, 0.        ],
            [0.47712125, 0.30103   , 0.30103   , 0.        , 0.        , 0.30103   ],
            [0.        , 0.30103   , 0.30103   , 0.        , 0.        , 0.        ]
        ]
        out2 = mu._log(df, base=10, a=1)

        self.assert_matrices_equal(out1, out2)

        out1 = [
            [-2.30258509e+01, -2.30258509e+01,  1.00000008e-10, -2.30258509e+01, -2.30258509e+01, -2.30258509e+01],
            [ 1.60943791e+00,  1.00000008e-10, -2.30258509e+01,  6.93147181e-01,  1.09861229e+00,  1.00000008e-10],
            [-2.30258509e+01, -2.30258509e+01,  1.00000008e-10,  1.38629436e+00,  6.93147181e-01, -2.30258509e+01],
            [ 6.93147181e-01,  1.00000008e-10,  1.00000008e-10, -2.30258509e+01, -2.30258509e+01,  1.00000008e-10],
            [-2.30258509e+01,  1.00000008e-10,  1.00000008e-10, -2.30258509e+01, -2.30258509e+01, -2.30258509e+01]
        ]
        out2 = mu._log(df, base=2.718281828459045, a=1e-10)

        self.assert_matrices_equal(out1, out2)


    def test_transform_op_validation(self):
        mu = MatrixUtil
        m = np.array(self.matrix)
        m0 = np.zeros((0,0)) # ?
        m1 = np.zeros(m.shape); m1[0,0] = 1

        DF = pd.DataFrame

        ## Logit ##
        ## In range (0,1) ##
        with self.subTest():
            valid = [(m+1)/10, m1/2+.1]
            invalid = [m, m*0, m*0+1, m-2*m, m1, -m1]

            for m_v in valid:
                mu._logit(DF(m_v))
            for m_inv in invalid:
                with self.assertRaises(MatrixValidationException) as cm:
                    mu._logit(DF(m_inv))
                print(cm.exception)

        ## Log ##
        ## Nonnegative after offset ##
        with self.subTest():
            valid = [m, m-1, m1, m1-1, m1*0, m1-1e-6]
            invalid = [m-2, m*-2, m1-2, m1*-2]

            for m_v in valid:
                mu._log(DF(m_v), base=10, a=1)
            for m_inv in invalid:
                with self.assertRaises(MatrixValidationException) as cm:
                    mu._log(DF(m_inv), base=10, a=1)
                    mu._log(DF(m_inv), base=2.718281828459045, a=1e-10)
                print(cm.exception)

        ## Sqrt ##
        ## Nonnegative ##
        with self.subTest():
            valid = [m, m1, m1*0]
            invalid = [m-1, m-1e-6, m*-1]

            for m_v in valid:
                mu._sqrt(DF(m_v))
            for m_inv in invalid:
                with self.assertRaises(MatrixValidationException) as cm:
                    mu._sqrt(DF(m_inv))
                print(cm.exception)


    def test_transform_identity(self):
        '''
        Test identity operations
        '''
        ret = self.getImpl().transform_matrix(
            self.ctx, {
                'workspace_id': self.getWsId(),
                'input_matrix_ref': self.amplicon_matrix_ref,
                'operations': [
                    'abundance_filtering',
                    'standardization',
                    'abundance_filtering',
                    'standardization',
                    'abundance_filtering',
                    'standardization',
                ],
                'abundance_filtering_params': {
                    'abundance_filtering_row_threshold': -1,
                    'abundance_filtering_columns_threshold': -1,
                    'abundance_filtering_row_sum_threshold': -1,
                    'abundance_filtering_columns_sum_threshold': -1,
                },
                'standardization_params': {
                    'standardization_with_mean': 0,
                    'standardization_with_std': 0,
                },
            })

        _, matrix_out = self.get_out_data(ret)

        self.assert_matrices_equal(matrix_out, self.matrix)


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
        self.assert_matrices_equal(matrix_out_mean, self.matrix_bootstrap9Reps_mean)


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


    def test_link_matrix_to_samples(self):
        matrix_obj = {
            'data': {
                'col_ids': [
                    'Sample0',
                    'Sample1',
                    'Sample2',
                    'Sample3',
                ]
            }
        }
        sample_set_obj = {
            'samples': [
                {
                    'name': 'Sample0',
                    'id': '0',
                    'version': 100,
                }, {
                    'name': 'Sample1',
                    'id': '1',
                    'version': 100,
                }, {
                    'name': 'Sample2',
                    'id': '2',
                    'version': 100,
                }, {
                    'name': 'Sample4',
                    'id': '4',
                    'version': 100
                }
            ]
        }
        get_sample_rets = [
            {'node_tree': [{'id': 'Sample0'}]},
            {'node_tree': [{'id': 'Sample1'}]},
            {'node_tree': [{'id': 'Sample2'}]},
        ]

        mock_dfu = create_autospec(DataFileUtil, instance=True, spec_set=True)
        mock_dfu.get_objects.return_value = {'data': [{'data': sample_set_obj}]}
        mock_ss = create_autospec(SampleService, instance=True, spec_set=True)
        mock_ss.get_sample.side_effect = get_sample_rets
        
        serviceImpl = GenericsAPI(self.cfg)

        serviceImpl.matrix_util.dfu = mock_dfu
        serviceImpl.matrix_util.sample_ser = mock_ss

        serviceImpl.matrix_util._link_matrix_to_samples(
            'dummy/matrix/ref',
            matrix_obj,
            'dummy/ss/ref',
        )

        mock_ss.create_data_link.assert_has_calls([
            call({
                'upa': 'dummy/matrix/ref',
                'dataid': 'Sample0',
                'id': '0',
                'version': 100,
                'node': 'Sample0',
                'update': 1,
            }),
            call({
                'upa': 'dummy/matrix/ref',
                'dataid': 'Sample1',
                'id': '1',
                'version': 100,
                'node': 'Sample1',
                'update': 1,
            }),
            call({
                'upa': 'dummy/matrix/ref',
                'dataid': 'Sample2',
                'id': '2',
                'version': 100,
                'node': 'Sample2',
                'update': 1,
            })
        ])

