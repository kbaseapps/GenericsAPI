# -*- coding: utf-8 -*-
import inspect
import os
import unittest
import time
import uuid
from mock import patch
import shutil
import requests
import pandas as pd
import numpy as np
from configparser import ConfigParser

from GenericsAPI.Utils.CorrelationUtil import CorrelationUtil
from GenericsAPI.GenericsAPIImpl import GenericsAPI
from GenericsAPI.GenericsAPIServer import MethodContext
from GenericsAPI.authclient import KBaseAuth as _KBaseAuth
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace as workspaceService
from installed_clients.AbstractHandleClient import AbstractHandle as HandleService


class CorrUtilTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
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
        cls.dfu = DataFileUtil(cls.callback_url)
        cls.corr_util = CorrelationUtil(cls.cfg)
        cls.hs = HandleService(url=cls.cfg['handle-service-url'],
                               token=cls.token)

        suffix = int(time.time() * 1000)
        cls.wsName = "test_corr_util_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})
        cls.wsId = ret[0]
        cls.shockURL = cls.cfg['shock-url']

        small_file = os.path.join(cls.scratch, 'test.txt')
        with open(small_file, "w") as f:
            f.write("empty content")
        cls.test_shock = cls.dfu.file_to_shock({'file_path': small_file, 'make_handle': True})
        cls.handles_to_delete = []
        cls.nodes_to_delete = []
        cls.handles_to_delete.append(cls.test_shock['handle']['hid'])
        cls.nodes_to_delete.append(cls.test_shock['shock_id'])

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

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        return self.__class__.wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    def getCorrUtil(self):
        return self.__class__.corr_util

    def mock_generate_corr_report(corr_matrix_obj_ref, workspace_name, corr_matrix_plot_path,
                                  scatter_plot_path=None, df1=None, df2=None):
        print('Mocking CorrelationUtil._generate_corr_report')

        return {'report_name': 'fake_report_name', 'report_ref': 'fake_report_ref'}

    def mock_file_to_shock(params):
        print('Mocking DataFileUtilClient.file_to_shock')
        print(params)

        return CorrUtilTest().test_shock

    def loadDF(self):

        if hasattr(self.__class__, 'test_random_df'):
            return self.__class__.test_random_df

        # Build test DF
        index = ['idx_{}'.format(idx) for idx in range(0, 6)]
        columns = ['col_{}'.format(idx) for idx in range(0, 10)]

        np.random.seed(0)
        test_random_df = pd.DataFrame(np.random.uniform(low=0, high=100,
                                                        size=(len(index), len(columns))),
                                      index=index, columns=columns)

        self.__class__.test_random_df = test_random_df
        print('Created Random DF:\n' + test_random_df.to_string())
        return test_random_df

    def loadExpressionMatrix(self):
        if hasattr(self.__class__, 'expr_matrix_ref'):
            return self.__class__.expr_matrix_ref

        data = {'data': {'col_ids': ['instance_1', 'instance_2', 'instance_3', 'instance_4'],
                         'row_ids': ['WRI_RS00010_CDS_1', 'WRI_RS00015_CDS_1', 'WRI_RS00025_CDS_1'],
                         'values': [[0.1, 0.2, 0.3, 0.4],
                                    [0.5, 0.6, 0.7, 0.8],
                                    [None, None, 1.1, 1.2]]},
                'scale': 'log2'}

        info = self.dfu.save_objects({
            'id': self.wsId,
            'objects': [{
                'type': 'KBaseMatrices.ExpressionMatrix',
                'data': data,
                'name': 'test_ExpressionMatrix'
            }]
        })[0]

        expr_matrix_ref = "%s/%s/%s" % (info[6], info[0], info[4])

        self.__class__.expr_matrix_ref = expr_matrix_ref
        print('Loaded ExpressionMatrix: ' + expr_matrix_ref)
        return expr_matrix_ref

    def loadExpressionMatrix2(self):
        if hasattr(self.__class__, 'expr_matrix_ref_2'):
            return self.__class__.expr_matrix_ref_2

        data = {'data': {'col_ids': ['instance_1', 'instance_4', 'instance_3'],
                         'row_ids': ['gene_1', 'gene_2', 'gene_3'],
                         'values': [[0.1, 0.4, 0.3], [0.5, None, 0.7], [None, 1.2, 1.1]]},
                'scale': 'log2'}

        info = self.dfu.save_objects({
            'id': self.wsId,
            'objects': [{
                'type': 'KBaseMatrices.ExpressionMatrix',
                'data': data,
                'name': 'test_ExpressionMatrix_2'
            }]
        })[0]

        expr_matrix_ref_2 = "%s/%s/%s" % (info[6], info[0], info[4])

        self.__class__.expr_matrix_ref_2 = expr_matrix_ref_2
        print('Loaded ExpressionMatrix: ' + expr_matrix_ref_2)
        return expr_matrix_ref_2

    def mock_download_staging_file(params):
        print('Mocking DataFileUtilClient.download_staging_file')
        print(params)

        file_path = params.get('staging_file_subdir_path')

        return {'copy_file_path': file_path}

    @patch.object(DataFileUtil, "download_staging_file", side_effect=mock_download_staging_file)
    @patch.object(DataFileUtil, "file_to_shock", side_effect=mock_file_to_shock)
    def loadAmpliconMatrix(self, download_staging_file, file_to_shock):
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
                  'target_gene': '16S',
                  'target_subfragment': ['V1'],
                  'sequencing_technology': 'Illumina',
                  'sequencing_instrument': 'Illumina Genome Analyzer',
                  'taxon_calling': {'taxon_calling_method': ['clustering'],
                                    'clustering_cutoff': 0.3,
                                    'clustering_method': 'clustering_method'},
                  'input_local_file': True
                  }

        returnVal = self.getImpl().import_matrix_from_biom(self.ctx, params)[0]

        amplicon_matrix_ref = returnVal['matrix_obj_ref']

        self.__class__.amplicon_matrix_ref = amplicon_matrix_ref
        print('Loaded AmpliconMatrix: ' + amplicon_matrix_ref)
        return amplicon_matrix_ref

    def loadCorrData(self):
        if hasattr(self.__class__, 'corr_data'):
            return self.__class__.corr_data

        corr_data = {'row_ids': ['WRI_RS00010_CDS_1',
                                 'WRI_RS00015_CDS_1',
                                 'WRI_RS00025_CDS_1'],
                     'values': [[1.0, 0.99, 0.91],
                                [0.99, 1.0, 0.91],
                                [0.91, 0.91, 1.0]],
                     'col_ids': ['WRI_RS00010_CDS_1',
                                 'WRI_RS00015_CDS_1',
                                 'WRI_RS00025_CDS_1']}

        self.__class__.corr_data = corr_data
        print('Loaded Correlation Data:\n{}\n'.format(corr_data))
        return corr_data

    def loadCorrMatrix(self):

        if hasattr(self.__class__, 'corr_matrix_ref'):
            return self.__class__.corr_matrix_ref

        object_type = 'KBaseExperiments.CorrelationMatrix'
        corr_matrix_object_name = 'test_corr_matrix'
        corr_matrix_data = {'correlation_parameters': {'method': 'pearson'},
                            'coefficient_data': self.loadCorrData(),
                            'significance_data': {'row_ids': ['WRI_RS00010_CDS_1',
                                                              'WRI_RS00015_CDS_1',
                                                              'WRI_RS00025_CDS_1'],
                                                  'values': [[0.0, 0.0, 0.0879],
                                                             [0.0, 0.0, 0.0879],
                                                             [0.0879, 0.0879, 0.0]],
                                                  'col_ids': ['WRI_RS00010_CDS_1',
                                                              'WRI_RS00015_CDS_1',
                                                              'WRI_RS00025_CDS_1']}}

        save_object_params = {
            'id': self.wsId,
            'objects': [{'type': object_type,
                         'data': corr_matrix_data,
                         'name': corr_matrix_object_name}]
        }

        dfu_oi = self.dfu.save_objects(save_object_params)[0]
        corr_matrix_ref = str(dfu_oi[6]) + '/' + str(dfu_oi[0]) + '/' + str(dfu_oi[4])

        self.__class__.corr_matrix_ref = corr_matrix_ref
        print('Loaded Correlation Matrix: ' + corr_matrix_ref)
        return corr_matrix_ref

    def fail_df_to_corr(self, df, error, method='pearson', dimension='col',
                        exception=ValueError, contains=False):
        with self.assertRaises(exception) as context:
            self.getCorrUtil().df_to_corr(df, method=method, dimension=dimension)
        if contains:
            self.assertIn(error, str(context.exception.args))
        else:
            self.assertEqual(error, str(context.exception.args[0]))

    def fail_plot_scatter_matrix(self, df, error, alpha=0.2, dimension='col',
                                 exception=ValueError, contains=False):
        with self.assertRaises(exception) as context:
            self.getCorrUtil().plot_scatter_matrix(df, alpha=alpha, dimension=dimension)
        if contains:
            self.assertIn(error, str(context.exception.args))
        else:
            self.assertEqual(error, str(context.exception.args[0]))

    def fail_compute_correlation_matrix(self, params, error, exception=ValueError,
                                        contains=False):
        with self.assertRaises(exception) as context:
            self.getImpl().compute_correlation_matrix(self.ctx, params)
        if contains:
            self.assertIn(error, str(context.exception.args[0]))
        else:
            self.assertEqual(error, str(context.exception.args[0]))

    def fail_compute_correlation_across_matrices(self, params, error, exception=ValueError,
                                                 contains=False):
        with self.assertRaises(exception) as context:
            self.getImpl().compute_correlation_across_matrices(self.ctx, params)
        if contains:
            self.assertIn(error, str(context.exception.args[0]))
        else:
            self.assertEqual(error, str(context.exception.args[0]))

    def start_test(self):
        testname = inspect.stack()[1][3]
        print('\n*** starting test: ' + testname + ' **')

    def test_export_corr_matrix_excel(self):
        self.start_test()

        corr_matrix_ref = self.loadCorrMatrix()

        params = {'input_ref': corr_matrix_ref}

        ret = self.getImpl().export_corr_matrix_excel(self.ctx, params)[0]

        assert ret and ('shock_id' in ret)

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        os.makedirs(output_directory)

        self.dfu.shock_to_file({'shock_id': ret['shock_id'],
                                'file_path': output_directory,
                                'unpack': 'unpack'})

        xl_files = [file for file in os.listdir(output_directory) if file.endswith('.xlsx')]
        self.assertEqual(len(xl_files), 1)

        xl = pd.ExcelFile(os.path.join(output_directory, xl_files[0]))
        expected_sheet_names = ['coefficient_data', 'significance_data']
        self.assertCountEqual(xl.sheet_names, expected_sheet_names)

        df = pd.read_excel(os.path.join(output_directory, xl_files[0]), index_col=0,
                           sheet_name="coefficient_data")
        expected_index = ['WRI_RS00010_CDS_1', 'WRI_RS00015_CDS_1', 'WRI_RS00025_CDS_1']
        expected_col = ['WRI_RS00010_CDS_1', 'WRI_RS00015_CDS_1', 'WRI_RS00025_CDS_1']
        self.assertCountEqual(df.index.tolist(), expected_index)
        self.assertCountEqual(df.columns.tolist(), expected_col)

        df = pd.read_excel(os.path.join(output_directory, xl_files[0]), index_col=0,
                           sheet_name="significance_data")
        expected_index = ['WRI_RS00010_CDS_1', 'WRI_RS00015_CDS_1', 'WRI_RS00025_CDS_1']
        expected_col = ['WRI_RS00010_CDS_1', 'WRI_RS00015_CDS_1', 'WRI_RS00025_CDS_1']
        self.assertCountEqual(df.index.tolist(), expected_index)
        self.assertCountEqual(df.columns.tolist(), expected_col)

    def test_compute_correlation_across_matrices_fail(self):
        self.start_test()

        invalidate_params = {'missing_matrix_ref_1': 'matrix_ref_1',
                             'matrix_ref_2': 'matrix_ref_2',
                             'workspace_name': 'workspace_name',
                             'corr_matrix_name': 'corr_matrix_name'}
        error_msg = '"matrix_ref_1" parameter is required, but missing'
        self.fail_compute_correlation_across_matrices(invalidate_params, error_msg)

        invalidate_params = {'matrix_ref_1': 'matrix_ref_1',
                             'matrix_ref_2': 'matrix_ref_2',
                             'missing_workspace_name': 'workspace_name',
                             'corr_matrix_name': 'corr_matrix_name'}
        error_msg = '"workspace_name" parameter is required, but missing'
        self.fail_compute_correlation_across_matrices(invalidate_params, error_msg)

    @unittest.skip("redundant and takes long time to run")
    def test_compute_correlation_across_matrices_amplicon_matrix_ok(self):
        self.start_test()
        expr_matrix_ref = self.loadAmpliconMatrix()
        expr_matrix_ref_2 = self.loadAmpliconMatrix()

        params = {'matrix_ref_1': expr_matrix_ref,
                  'matrix_ref_2': expr_matrix_ref_2,
                  'workspace_name': self.wsName,
                  'corr_matrix_name': 'test_corr_matrix',
                  'plot_corr_matrix': True,
                  'compute_significance': True}

        ret = self.getImpl().compute_correlation_across_matrices(self.ctx, params)[0]

        self.assertIn('corr_matrix_obj_ref', ret)
        corr_matrix_obj_ref = ret.get('corr_matrix_obj_ref')

        res = self.dfu.get_objects({'object_refs': [corr_matrix_obj_ref]})['data'][0]
        obj_data = res['data']

        self.assertEqual(obj_data.get('correlation_parameters').get('method'), 'pearson')

        expected_index = ['GG_OTU_1', 'GG_OTU_2',
                          'GG_OTU_3', 'GG_OTU_4',
                          'GG_OTU_5']

        self.assertCountEqual(obj_data.get('coefficient_data').get('row_ids'), expected_index)
        self.assertCountEqual(obj_data.get('coefficient_data').get('col_ids'), expected_index)

        self.assertCountEqual(obj_data.get('significance_data').get('row_ids'), expected_index)
        self.assertCountEqual(obj_data.get('significance_data').get('col_ids'), expected_index)

    @patch.object(DataFileUtil, "file_to_shock", side_effect=mock_file_to_shock)
    def test_compute_correlation_across_matrices_ok(self, file_to_shock):
        self.start_test()
        expr_matrix_ref = self.loadExpressionMatrix()
        expr_matrix_ref_2 = self.loadExpressionMatrix2()

        params = {'matrix_ref_1': expr_matrix_ref,
                  'matrix_ref_2': expr_matrix_ref_2,
                  'workspace_name': self.wsName,
                  'corr_matrix_name': 'test_corr_matrix',
                  'plot_corr_matrix': True,
                  'compute_significance': True}

        ret = self.getImpl().compute_correlation_across_matrices(self.ctx, params)[0]

        self.assertIn('corr_matrix_obj_ref', ret)
        corr_matrix_obj_ref = ret.get('corr_matrix_obj_ref')

        res = self.dfu.get_objects({'object_refs': [corr_matrix_obj_ref]})['data'][0]
        obj_data = res['data']

        self.assertEqual(obj_data.get('correlation_parameters').get('method'), 'pearson')

        expected_index = ['WRI_RS00010_CDS_1', 'WRI_RS00015_CDS_1', 'WRI_RS00025_CDS_1']
        expected_col = ['gene_1', 'gene_2', 'gene_3']
        self.assertCountEqual(obj_data.get('coefficient_data').get('row_ids'), expected_index)
        self.assertCountEqual(obj_data.get('coefficient_data').get('col_ids'), expected_col)

        self.assertCountEqual(obj_data.get('significance_data').get('row_ids'), expected_index)
        self.assertCountEqual(obj_data.get('significance_data').get('col_ids'), expected_col)

    def test_compute_correlation_matrix_fail(self):
        self.start_test()

        invalidate_params = {'missing_input_obj_ref': 'input_obj_ref',
                             'workspace_name': 'workspace_name'}
        error_msg = '"input_obj_ref" parameter is required, but missing'
        self.fail_compute_correlation_matrix(invalidate_params, error_msg)

        invalidate_params = {'input_obj_ref': 'input_obj_ref',
                             'missing_workspace_name': 'workspace_name'}
        error_msg = '"workspace_name" parameter is required, but missing'
        self.fail_compute_correlation_matrix(invalidate_params, error_msg)

    @patch.object(DataFileUtil, "file_to_shock", side_effect=mock_file_to_shock)
    def test_comp_corr_matrix_ok(self, file_to_shock):
        self.start_test()
        expr_matrix_ref = self.loadExpressionMatrix()

        params = {'input_obj_ref': expr_matrix_ref,
                  'workspace_name': self.wsName,
                  'corr_matrix_name': 'test_corr_matrix',
                  'plot_corr_matrix': True,
                  'plot_scatter_matrix': True,
                  'compute_significance': True}

        ret = self.getImpl().compute_correlation_matrix(self.ctx, params)[0]

        self.assertIn('corr_matrix_obj_ref', ret)
        corr_matrix_obj_ref = ret.get('corr_matrix_obj_ref')

        res = self.dfu.get_objects({'object_refs': [corr_matrix_obj_ref]})['data'][0]
        obj_data = res['data']

        self.assertEqual(obj_data.get('correlation_parameters').get('method'), 'pearson')

        corr_items = ['WRI_RS00010_CDS_1', 'WRI_RS00015_CDS_1', 'WRI_RS00025_CDS_1']
        self.assertCountEqual(obj_data.get('coefficient_data').get('row_ids'), corr_items)
        self.assertCountEqual(obj_data.get('coefficient_data').get('col_ids'), corr_items)

        self.assertCountEqual(obj_data.get('significance_data').get('row_ids'), corr_items)
        self.assertCountEqual(obj_data.get('significance_data').get('col_ids'), corr_items)

    def test_init_ok(self):
        self.start_test()
        class_attri = ['scratch', 'token', 'callback_url', 'ws_url']

        corr_util = self.getCorrUtil()
        self.assertTrue(set(class_attri) <= set(corr_util.__dict__.keys()))
        self.assertEqual(corr_util.scratch, self.cfg.get('scratch'))

    def test_df_to_corr_fail(self):
        self.start_test()
        df = self.loadDF()

        # unavailable correlation method
        self.fail_df_to_corr(df, 'Input correlation method', method='fake_method', contains=True)

        # unavailable dimension
        self.fail_df_to_corr(df, 'Input dimension', dimension='fake_dimension', contains=True)

    def test_df_to_corr_ok(self):
        self.start_test()
        df = self.loadDF()

        # test default
        corr_df = self.getCorrUtil().df_to_corr(df)
        self.assertCountEqual(corr_df.index.tolist(), corr_df.columns.tolist())
        self.assertCountEqual(corr_df.index.tolist(), df.columns.tolist())

        # test correlation on rows
        corr_df = self.getCorrUtil().df_to_corr(df, dimension='row')
        self.assertCountEqual(corr_df.index.tolist(), corr_df.columns.tolist())
        self.assertCountEqual(corr_df.index.tolist(), df.index.tolist())

    def test_plot_scatter_matrix_fail(self):
        self.start_test()
        df = self.loadDF()

        # unavailable dimension
        self.fail_plot_scatter_matrix(df, 'Input dimension', dimension='fake_dim', contains=True)

        # fail pd.plotting.scatter_matrix
        self.fail_plot_scatter_matrix(df, 'Running scatter_matrix returned an error', alpha=10, contains=True)

    def test_plot_scatter_matrix_ok(self):
        self.start_test()
        df = self.loadDF()

        # default
        scatter_plot_path = self.getCorrUtil().plot_scatter_matrix(df)
        self.assertTrue(os.path.isfile(scatter_plot_path))
        self.assertGreater(os.path.getsize(scatter_plot_path), 1024)  # file size greate than 1KB

        # scatter plot on rows
        scatter_plot_path = self.getCorrUtil().plot_scatter_matrix(df, dimension='row')
        self.assertTrue(os.path.isfile(scatter_plot_path))
        self.assertGreater(os.path.getsize(scatter_plot_path), 1024)  # file size greate than 1KB

    def test_plot_corr_matrix_ok(self):
        self.start_test()
        df = self.loadDF()
        corr_df = self.getCorrUtil().df_to_corr(df, dimension='col')

        # default
        corr_matrix_plot_path = self.getCorrUtil().plot_corr_matrix(corr_df)
        self.assertTrue(os.path.isfile(corr_matrix_plot_path))
        self.assertGreater(os.path.getsize(corr_matrix_plot_path), 1024)  # file size greate than 1KB
