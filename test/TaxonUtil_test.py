# -*- coding: utf-8 -*-
import inspect
import unittest

from GenericsAPI.Utils.TaxonUtil import TaxonUtil


class TaxonUtilTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cfg = {}
        cls.taxon_util = TaxonUtil(cls.cfg)

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
