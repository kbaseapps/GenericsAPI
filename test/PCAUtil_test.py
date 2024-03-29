# -*- coding: utf-8 -*-
import inspect
import os  # noqa: F401
import unittest
import time
from configparser import ConfigParser
import uuid
import pandas as pd
from mock import patch
import requests

from GenericsAPI.Utils.PCAUtil import PCAUtil
from GenericsAPI.GenericsAPIImpl import GenericsAPI
from GenericsAPI.GenericsAPIServer import MethodContext
from GenericsAPI.authclient import KBaseAuth as _KBaseAuth
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace as workspaceService
from installed_clients.AbstractHandleClient import AbstractHandle as HandleService


class PCAUtilTest(unittest.TestCase):

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
        cls.shockURL = cls.cfg['shock-url']
        cls.serviceImpl = GenericsAPI(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        cls.dfu = DataFileUtil(cls.callback_url)
        cls.pca_util = PCAUtil(cls.cfg)
        cls.hs = HandleService(url=cls.cfg['handle-service-url'],
                               token=cls.token)

        suffix = int(time.time() * 1000)
        cls.wsName = "test_pca_util_" + str(suffix)
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

    def getPCAUtil(self):
        return self.__class__.pca_util

    def mock_generate_pca_report(pca_ref, score_plots, loading_plots, bi_plots,
                                 workspace_name, n_components):
        print('Mocking PCAUtil._generate_pca_report')

        return {'report_name': 'fake_report_name', 'report_ref': 'fake_report_ref'}

    def mock_file_to_shock(params):
        print('Mocking DataFileUtilClient.file_to_shock')
        print(params)

        return PCAUtilTest().test_shock

    def loadExpressionMatrix(self):
        if hasattr(self.__class__, 'expr_matrix_ref') and hasattr(self.__class__,
                                                                  'asso_matrix_ref'):
            return self.__class__.expr_matrix_ref, self.__class__.asso_matrix_ref

        # matrix_file_name = 'test_import.xlsx'
        col_attribute = {'attributes': [{'attribute': 'test_attribute_1',
                                         'attribute_ont_id': 'OBI_0500020',
                                         'source': 'upload',
                                         'unit': 'Hour',
                                         'unit_ont_id': 'UO_0000032'},
                                        {'attribute': 'test_attribute_2',
                                         'attribute_ont_id': 'CHEBI:9168',
                                         'source': 'upload',
                                         'unit': 'nanogram per milliliter',
                                         'unit_ont_id': 'UO_0000275'},
                                        {'attribute': 'test_attribute_3',
                                         'attribute_ont_id': 'CHEBI:9168',
                                         'source': 'upload',
                                         'unit': 'nanogram per milliliter',
                                         'unit_ont_id': 'UO_0000275'}],
                         'instances': {'test_col_instance_1': ['1', '5', '9'],
                                       'test_col_instance_2': ['2', '6', '10'],
                                       'test_col_instance_3': ['3', '7', '11'],
                                       'test_col_instance_4': ['4', '8', '12']},
                         'ontology_mapping_method': 'User Curation'}

        info = self.dfu.save_objects({
                            'id': self.wsId,
                            'objects': [{
                                'type': 'KBaseExperiments.AttributeMapping',
                                'data': col_attribute,
                                'name': 'test_ExpressionMatrix_col_attribute_mapping'
                            }]
                        })[0]

        col_attributemapping_ref = "%s/%s/%s" % (info[6], info[0], info[4])

        row_attribute = {'attributes': [{'attribute': 'test_attribute_1',
                                         'attribute_ont_id': 'OBI_0500020',
                                         'source': 'upload',
                                         'unit': 'Hour',
                                         'unit_ont_id': 'UO_0000032'},
                                        {'attribute': 'test_attribute_2',
                                         'attribute_ont_id': 'CHEBI:9168',
                                         'source': 'upload',
                                         'unit': 'nanogram per milliliter',
                                         'unit_ont_id': 'UO_0000275'},
                                        {'attribute': 'test_attribute_3',
                                         'attribute_ont_id': 'CHEBI:9168',
                                         'source': 'upload',
                                         'unit': 'nanogram per milliliter',
                                         'unit_ont_id': 'UO_0000275'}],
                         'instances': {'test_row_instance_1': ['1', '4', '7'],
                                       'test_row_instance_2': ['3', '4', '8'],
                                       'test_row_instance_3': ['3', '6', '7']},
                         'ontology_mapping_method': 'User Curation'}

        info = self.dfu.save_objects({
                            'id': self.wsId,
                            'objects': [{
                                'type': 'KBaseExperiments.AttributeMapping',
                                'data': row_attribute,
                                'name': 'test_ExpressionMatrix_row_attribute_mapping'
                            }]
                        })[0]

        row_attributemapping_ref = "%s/%s/%s" % (info[6], info[0], info[4])

        matrix_data = {'attributes': {'Instrument': 'Old Faithful',
                                      'Scientist': 'Marie Currie'},
                       'col_attributemapping_ref': col_attributemapping_ref,
                       'col_mapping': {'instance_1': 'test_col_instance_1',
                                       'instance_2': 'test_col_instance_2',
                                       'instance_3': 'test_col_instance_3',
                                       'instance_4': 'test_col_instance_4'},
                       'col_normalization': 'test_col_normalization',
                       'data': {'col_ids': ['instance_1', 'instance_2', 'instance_3', 'instance_4'],
                                'row_ids': ['WRI_RS00010_CDS_1', 'WRI_RS00015_CDS_1', 'WRI_RS00025_CDS_1'],
                                'values': [[0.1, 0.2, 0.3, 0.4],
                                           [0.5, 0.6, 0.7, 0.8],
                                           [None, None, 1.1, 1.2]]},
                       'description': 'test_desc',
                       'row_attributemapping_ref': row_attributemapping_ref,
                       'row_mapping': {'WRI_RS00010_CDS_1': 'test_row_instance_1',
                                       'WRI_RS00015_CDS_1': 'test_row_instance_2',
                                       'WRI_RS00025_CDS_1': 'test_row_instance_3'},
                       'row_normalization': 'test_row_normalization',
                       'scale': 'log2',
                       'search_attributes': ['Scientist | Marie Currie',
                                             'Instrument | Old Faithful']}

        info = self.dfu.save_objects({
                            'id': self.wsId,
                            'objects': [{
                                'type': 'KBaseMatrices.ExpressionMatrix',
                                'data': matrix_data,
                                'name': 'test_ExpressionMatrix'
                            }]
                        })[0]

        expr_matrix_ref = "%s/%s/%s" % (info[6], info[0], info[4])

        self.__class__.expr_matrix_ref = expr_matrix_ref
        print('Loaded ExpressionMatrix: ' + expr_matrix_ref)

        # load associated matrix
        matrix_data = {'attributes': {'Instrument': 'Old Faithful',
                                      'Scientist': 'Marie Currie'},
                       'col_attributemapping_ref': col_attributemapping_ref,
                       'col_mapping': {'instance_1': 'instance_1',
                                       'instance_2': 'instance_2',
                                       'instance_3': 'instance_3',
                                       'instance_4': 'instance_4'},
                       'col_normalization': 'test_col_normalization',
                       'data': {'col_ids': ['instance_1', 'instance_2', 'instance_3',
                                            'instance_4'],
                                'row_ids': ['WRI_RS00010_CDS_1', 'WRI_RS00015_CDS_1',
                                            'WRI_RS00025_CDS_1', 'WRI_RS00030_CDS_1',
                                            'WRI_RS00035_CDS_1'],
                                'values': [[0.1, 0.2, 0.3, 0.4],
                                           [0.5, 0.6, 0.7, 0.8],
                                           [0.9, 1, 1.1, 1.2],
                                           [0.9, 1, 1.1, 1.2],
                                           [0.9, 1, 1.1, 1.2]]},
                       'description': 'test_desc',
                       'row_attributemapping_ref': row_attributemapping_ref,
                       'row_mapping': {'WRI_RS00010_CDS_1': 'WRI_RS00010_CDS_1',
                                       'WRI_RS00015_CDS_1': 'WRI_RS00015_CDS_1',
                                       'WRI_RS00025_CDS_1': 'WRI_RS00025_CDS_1',
                                       'WRI_RS00030_CDS_1': 'WRI_RS00030_CDS_1',
                                       'WRI_RS00035_CDS_1': 'WRI_RS00035_CDS_1'},
                       'row_normalization': 'test_row_normalization',
                       'scale': 'log2',
                       'search_attributes': ['Scientist | Marie Currie',
                                             'Instrument | Old Faithful']}

        info = self.dfu.save_objects({
            'id': self.wsId,
            'objects': [{'type': 'KBaseMatrices.ExpressionMatrix',
                         'data': matrix_data,
                         'name': 'test_associated_ExpressionMatrix'}]})[0]

        asso_matrix_ref = "%s/%s/%s" % (info[6], info[0], info[4])

        self.__class__.asso_matrix_ref = asso_matrix_ref
        print('Loaded Associated ExpressionMatrix: ' + asso_matrix_ref)

        return expr_matrix_ref, asso_matrix_ref

    def loadPCAMatrix(self):

        if hasattr(self.__class__, 'pca_matrix_ref'):
            return self.__class__.pca_matrix_ref

        object_type = 'KBaseExperiments.PCAMatrix'
        pca_matrix_object_name = 'test_PCA_matrix'
        pca_matrix_data = {'explained_variance_ratio': [0.628769688409428, 0.371230311590572],
                           'explained_variance': [0.628769688409428, 0.371230311590572],
                           'pca_parameters': {'dimension': 'row', 'n_components': '2'},
                           'rotation_matrix': {'col_ids': ['principal_component_1',
                                                           'principal_component_2'],
                                               'row_ids': ['WRI_RS00010_CDS_1',
                                                           'WRI_RS00015_CDS_1',
                                                           'WRI_RS00025_CDS_1'],
                                               'values': [[-0.45, 1.06],
                                                          [-0.69, -0.92],
                                                          [1.14, -0.13]]}}

        save_object_params = {
            'id': self.wsId,
            'objects': [{'type': object_type,
                         'data': pca_matrix_data,
                         'name': pca_matrix_object_name}]
        }

        dfu_oi = self.dfu.save_objects(save_object_params)[0]
        pca_matrix_ref = str(dfu_oi[6]) + '/' + str(dfu_oi[0]) + '/' + str(dfu_oi[4])

        self.__class__.pca_matrix_ref = pca_matrix_ref
        print('Loaded Correlation Matrix: ' + pca_matrix_ref)
        return pca_matrix_ref

    def fail_run_pca(self, params, error, exception=ValueError, contains=False):
        with self.assertRaises(exception) as context:
            self.getImpl().run_pca(self.ctx, params)
        if contains:
            self.assertIn(error, str(context.exception.args))
        else:
            self.assertEqual(error, str(context.exception.args[0]))

    def start_test(self):
        testname = inspect.stack()[1][3]
        print('\n*** starting test: ' + testname + ' **')

    def test_run_pca_fail(self):
        self.start_test()

        invalidate_params = {'missing_input_obj_ref': 'input_obj_ref',
                             'workspace_name': 'workspace_name'}
        error_msg = '"input_obj_ref" parameter is required, but missing'
        self.fail_run_pca(invalidate_params, error_msg)

        invalidate_params = {'input_obj_ref': 'input_obj_ref',
                             'missing_workspace_name': 'workspace_name'}
        error_msg = '"workspace_name" parameter is required, but missing'
        self.fail_run_pca(invalidate_params, error_msg)

    @patch.object(DataFileUtil, "file_to_shock", side_effect=mock_file_to_shock)
    def test_run_pca_ok(self, file_to_shock):
        self.start_test()

        expr_matrix_ref, asso_matrix_ref = self.loadExpressionMatrix()

        params = {'input_obj_ref': expr_matrix_ref,
                  'workspace_name': self.wsName,
                  'pca_matrix_name': 'test_pca_matrix_obj',
                  'scale_size_by': {"attribute_size": ["test_attribute_1"]},
                  'color_marker_by': {"attribute_color": ["test_attribute_2"]},
                  'n_components': 3,
                  'dimension': 'row'}

        ret = self.getImpl().run_pca(self.ctx, params)[0]

        self.assertTrue('report_name' in ret)
        self.assertTrue('report_ref' in ret)
        self.assertTrue('pca_ref' in ret)

        pca_matrix_ref = ret.get('pca_ref')

        pca_matrix_data = self.dfu.get_objects(
                    {"object_refs": [pca_matrix_ref]})['data'][0]['data']

        self.assertTrue('explained_variance_ratio' in pca_matrix_data)
        self.assertTrue('rotation_matrix' in pca_matrix_data)
        self.assertEqual(len(pca_matrix_data.get('explained_variance_ratio')), 3)

        expected_row_ids = ['WRI_RS00010_CDS_1', 'WRI_RS00015_CDS_1', 'WRI_RS00025_CDS_1']
        expected_col_ids = ['principal_component_1', 'principal_component_2',
                            'principal_component_3']
        self.assertCountEqual(pca_matrix_data['rotation_matrix']['row_ids'], expected_row_ids)
        self.assertCountEqual(pca_matrix_data['rotation_matrix']['col_ids'], expected_col_ids)

    @patch.object(DataFileUtil, "file_to_shock", side_effect=mock_file_to_shock)
    def test_run_pca_size_by_asso_matrix_ok(self, file_to_shock):
        self.start_test()

        expr_matrix_ref, asso_matrix_ref = self.loadExpressionMatrix()

        params = {'input_obj_ref': expr_matrix_ref,
                  'workspace_name': self.wsName,
                  'pca_matrix_name': 'test_pca_matrix_obj',
                  'associated_matrix_obj_ref': asso_matrix_ref,
                  'scale_size_by': {'row_size': ['WRI_RS00010_CDS_1']},
                  'color_marker_by': {"attribute_color": ["test_attribute_2"]},
                  'n_components': 3,
                  'dimension': 'col'}

        ret = self.getImpl().run_pca(self.ctx, params)[0]

        self.assertTrue('report_name' in ret)
        self.assertTrue('report_ref' in ret)
        self.assertTrue('pca_ref' in ret)

        pca_matrix_ref = ret.get('pca_ref')

        pca_matrix_data = self.dfu.get_objects(
                    {"object_refs": [pca_matrix_ref]})['data'][0]['data']

        self.assertTrue('explained_variance_ratio' in pca_matrix_data)
        self.assertTrue('rotation_matrix' in pca_matrix_data)
        self.assertEqual(len(pca_matrix_data.get('explained_variance_ratio')), 3)

        expected_row_ids = ['instance_1', 'instance_2', 'instance_3', 'instance_4']
        expected_col_ids = ['principal_component_1', 'principal_component_2',
                            'principal_component_3']
        self.assertCountEqual(pca_matrix_data['rotation_matrix']['row_ids'], expected_row_ids)
        self.assertCountEqual(pca_matrix_data['rotation_matrix']['col_ids'], expected_col_ids)

    def test_export_pca_matrix_excel_ok(self):
        self.start_test()

        pca_ref = self.loadPCAMatrix()

        params = {'input_ref': pca_ref}

        ret = self.getImpl().export_pca_matrix_excel(self.ctx, params)[0]

        assert ret and ('shock_id' in ret)

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        os.makedirs(output_directory)

        self.dfu.shock_to_file({'shock_id': ret['shock_id'],
                                'file_path': output_directory,
                                'unpack': 'unpack'})

        xl_files = [file for file in os.listdir(output_directory) if file.endswith('.xlsx')]
        self.assertEqual(len(xl_files), 1)

        xl = pd.ExcelFile(os.path.join(output_directory, xl_files[0]))
        expected_sheet_names = ['principal_component_matrix']
        self.assertCountEqual(xl.sheet_names, expected_sheet_names)

        df = pd.read_excel(os.path.join(output_directory, xl_files[0]), index_col=0,
                           sheet_name="principal_component_matrix")
        expected_index = ['WRI_RS00010_CDS_1', 'WRI_RS00015_CDS_1', 'WRI_RS00025_CDS_1']
        expected_col = ['principal_component_1', 'principal_component_2']
        self.assertCountEqual(df.index.tolist(), expected_index)
        self.assertCountEqual(df.columns.tolist(), expected_col)

    def test_init_ok(self):
        self.start_test()
        class_attri = ['scratch', 'token', 'callback_url', 'ws_url']

        network_util = self.getPCAUtil()
        self.assertTrue(set(class_attri) <= set(network_util.__dict__.keys()))
        self.assertEqual(network_util.scratch, self.cfg.get('scratch'))
