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

        taxonomic_str = 'Archaea;Nanoarchaeota;Woesearchaeia;UBA12501;UBA11576;UBA11576;(UBA11576)'
        processed_taxonomic_str = self.getTaxonUtil().process_taxonomic_str(taxonomic_str)
        expect_processed_taxonomic_str = 'Archaea;Nanoarchaeota;Woesearchaeia;UBA12501;UBA11576;UBA11576;(UBA11576);'
        self.assertEqual(processed_taxonomic_str, expect_processed_taxonomic_str)

        taxonomic_str = 'Bacteria;Proteobacteria;Gammaproteobacteria;Xanthomonadales;Sinobacteraceae;Nevskia;'
        processed_taxonomic_str = self.getTaxonUtil().process_taxonomic_str(taxonomic_str)
        expect_processed_taxonomic_str = 'Bacteria;Proteobacteria;Gammaproteobacteria;Xanthomonadales;Sinobacteraceae;Nevskia;'
        self.assertEqual(processed_taxonomic_str, expect_processed_taxonomic_str)

        taxonomic_str = 'A,B,C,D,E,F'
        processed_taxonomic_str = self.getTaxonUtil().process_taxonomic_str(taxonomic_str)
        expect_processed_taxonomic_str = 'A;B;C;D;E;F;'
        self.assertEqual(processed_taxonomic_str, expect_processed_taxonomic_str)

        taxonomic_str = 'A,B,C,D;E,(F)'
        processed_taxonomic_str = self.getTaxonUtil().process_taxonomic_str(taxonomic_str)
        expect_processed_taxonomic_str = 'A;B;C;D;E;(F);'
        self.assertEqual(processed_taxonomic_str, expect_processed_taxonomic_str)

        taxonomic_str = 'A;B;C;D;E;F'
        processed_taxonomic_str = self.getTaxonUtil().process_taxonomic_str(taxonomic_str)
        expect_processed_taxonomic_str = 'A;B;C;D;E;F;'
        self.assertEqual(processed_taxonomic_str, expect_processed_taxonomic_str)

        taxonomic_str = 'd:A,p:B,c:C,o:D,f:E,g:F'
        processed_taxonomic_str = self.getTaxonUtil().process_taxonomic_str(taxonomic_str)
        expect_processed_taxonomic_str = 'A;B;C;D;E;F;'
        self.assertEqual(processed_taxonomic_str, expect_processed_taxonomic_str)

        taxonomic_str = 'd:A,p:B,c:C,f:D,g:F;'
        processed_taxonomic_str = self.getTaxonUtil().process_taxonomic_str(taxonomic_str)
        expect_processed_taxonomic_str = 'A;B;C;;D;F;'
        self.assertEqual(processed_taxonomic_str, expect_processed_taxonomic_str)

        taxonomic_str = 'k__A; p__B; c__C; o__D; f__E; g__F; s__G'
        processed_taxonomic_str = self.getTaxonUtil().process_taxonomic_str(taxonomic_str)
        expect_processed_taxonomic_str = 'A;B;C;D;E;F;G;'
        self.assertEqual(processed_taxonomic_str, expect_processed_taxonomic_str)

        taxonomic_str = 'k__A; p__B; f__E; g__F; s__G'
        processed_taxonomic_str = self.getTaxonUtil().process_taxonomic_str(taxonomic_str)
        expect_processed_taxonomic_str = 'A;B;;;E;F;G;'
        self.assertEqual(processed_taxonomic_str, expect_processed_taxonomic_str)

        taxonomic_str = 'k__A; p__B; c__; o__; f__; g__; s__'
        processed_taxonomic_str = self.getTaxonUtil().process_taxonomic_str(taxonomic_str)
        expect_processed_taxonomic_str = 'A;B;;;;;;'
        self.assertEqual(processed_taxonomic_str, expect_processed_taxonomic_str)
