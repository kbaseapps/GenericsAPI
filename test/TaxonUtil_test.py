# -*- coding: utf-8 -*-
import inspect
import os  # noqa: F401
import unittest
import time
from configparser import ConfigParser

from GenericsAPI.Utils.TaxonUtil import TaxonUtil
from GenericsAPI.GenericsAPIImpl import GenericsAPI
from GenericsAPI.GenericsAPIServer import MethodContext
from GenericsAPI.authclient import KBaseAuth as _KBaseAuth
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace as workspaceService


class TaxonUtilTest(unittest.TestCase):

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
        cls.taxon_util = TaxonUtil(cls.cfg)

        suffix = int(time.time() * 1000)
        cls.wsName = "test_taxon_util_" + str(suffix)
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

    def getTaxonUtil(self):
        return self.__class__.taxon_util

    def start_test(self):
        testname = inspect.stack()[1][3]
        print('\n*** starting test: ' + testname + ' **')

    def test_process_taxonomic_str_ok(self):
        self.start_test()

        taxonomic_str = 'd:Bacteria,p:Firmicutes,c:Bacilli,o:Lactobacillales,f:Streptococcaceae,g:Streptococcus;'
        processed_taxonomic_str = self.getTaxonUtil().process_taxonomic_str(taxonomic_str)
        expect_processed_taxonomic_str = 'Bacteria;Firmicutes;Bacilli;Lactobacillales;Streptococcaceae;Streptococcus;'
        self.assertEqual(processed_taxonomic_str, expect_processed_taxonomic_str)

        taxonomic_str = 'd:Bacteria,p:Firmicutes,c:Bacilli,f:Streptococcaceae,g:Streptococcus;'
        processed_taxonomic_str = self.getTaxonUtil().process_taxonomic_str(taxonomic_str)
        expect_processed_taxonomic_str = 'Bacteria;Firmicutes;Bacilli;;Streptococcaceae;Streptococcus;'
        self.assertEqual(processed_taxonomic_str, expect_processed_taxonomic_str)
