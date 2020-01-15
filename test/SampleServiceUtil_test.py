# -*- coding: utf-8 -*-
import os  # noqa: F401
import unittest
import time
from configparser import ConfigParser
import shutil

from GenericsAPI.Utils.SampleServiceUtil import SampleServiceUtil
from GenericsAPI.GenericsAPIImpl import GenericsAPI
from GenericsAPI.GenericsAPIServer import MethodContext
from GenericsAPI.authclient import KBaseAuth as _KBaseAuth
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace as workspaceService
from installed_clients.sample_uploaderClient import sample_uploader


class SampleServiceTest(unittest.TestCase):

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

        cls.sampleservice_util = SampleServiceUtil(cls.cfg)

        suffix = int(time.time() * 1000)
        cls.wsName = "test_pca_util_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})
        cls.wsId = ret[0]

        cls.sample_id = '80d16006-62ac-4a36-99fe-f5861c4cc8c8'  # pre-saved sample

        cls.sample_uploader = sample_uploader(cls.callback_url, service_ver="dev")

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

    def getSampleServiceUtil(self):
        return self.__class__.sampleservice_util

    def loadSampleSet(self):
        if hasattr(self.__class__, 'sample_set_ref'):
            return self.__class__.sample_set_ref

        sample_set_file_name = 'ANLPW_JulySamples_IGSN_v2-forKB.csv'
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

    def test_get_sample_service_url(self):
        sampleservice_util = self.getSampleServiceUtil()

        ss_url = sampleservice_util.get_sample_service_url()

        print('Getting sample_service URL: {}'.format(ss_url))

        self.assertTrue('SampleService' in ss_url)

    def test_get_sample(self):
        sampleservice_util = self.getSampleServiceUtil()

        sample_rec = sampleservice_util.get_sample(self.sample_id)

        self.assertEqual(sample_rec['id'], self.sample_id)
        self.assertEqual(sample_rec['name'], 'OC-1')

    def test_sample_set_to_attribute_mapping(self):
        sampleservice_util = self.getSampleServiceUtil()
        sample_set_ref = self.loadSampleSet()

        am_data = sampleservice_util.sample_set_to_attribute_mapping(sample_set_ref)

        attributes = am_data['attributes']
        attri_names = [attribute['attribute'] for attribute in attributes]

        attri_names_expected = ['id', 'type', 'parent', 'Elevation start', 'Latitude',
                                'Collection date', 'Longitude', 'Name of physiographic feature',
                                'Locality Description', 'Primary physiographic feature', 'Purpose',
                                'Current archive contact', 'Relation Type', 'Field program/cruise',
                                'Collector/Chief Scientist', 'Navigation type',
                                'Coordinate Precision?', 'Material', 'Location Description',
                                'Related Identifiers', 'Collection method', 'Current archive']
        self.assertCountEqual(attri_names, attri_names_expected)

        instances = am_data['instances']
        sample_names = instances.keys()
        sample_names_expected = ['PB-Low-5', 'PB-High-5', 'PB-Low-6', 'PB-High-6', 'PB-Low-7',
                                 'PB-High-7', 'PB-Low-8', 'PB-High-8']

        self.assertCountEqual(sample_names, sample_names_expected)
