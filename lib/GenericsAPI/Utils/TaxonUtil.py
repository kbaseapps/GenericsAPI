
import logging
import csv
import re


class TaxonUtil:

    def _fetch_taxon_level(self, taxon_char):

        taxon_level_mapping = {'l': 'Life', 'd': 'Domain', 'k': 'Kingdom', 'p': 'Phylum',
                               'c': 'Class', 'o': 'Order', 'f': 'Family', 'g': 'Genus',
                               's': 'Species'}

        return taxon_level_mapping.get(taxon_char.lower(), taxon_char)

    def _convert_taxonomic_str(self, lineage, taxon_level_delimiter):
        logging.info('start converting lineage: {} with taxon level delimiter [{}]'.format(
                                                                            lineage,
                                                                            taxon_level_delimiter))
        taxon_levels = 'pcofgs'

        processed_taxonomic_str = ''
        taxon_str = lineage[0].rstrip(';')
        taxon_items = taxon_str.split(taxon_level_delimiter)
        scientific_name = taxon_items[-1]
        pre_taxon_level_char = taxon_items[0]
        processed_taxonomic_str += scientific_name + ';'

        for taxon_str in lineage[1:]:
            taxon_str = taxon_str.rstrip(';')
            taxon_items = taxon_str.split(taxon_level_delimiter)

            scientific_name = taxon_items[-1]
            taxon_level_char = taxon_items[0]

            taxon_level_diff = 1
            if pre_taxon_level_char in taxon_levels and taxon_level_char in taxon_levels:
                taxon_level_diff = max(taxon_levels.index(taxon_level_char) -
                                       taxon_levels.index(pre_taxon_level_char), 1)

            pre_taxon_level_char = taxon_level_char

            processed_taxonomic_str += (taxon_level_diff - 1) * ';'
            processed_taxonomic_str += scientific_name + ';'

        logging.info('converted taxonomic string: {}'.format(processed_taxonomic_str))

        return processed_taxonomic_str

    def _convert_taxonomic_str2(self, lineage, taxon_level_delimiter):
        logging.info('start converting lineage: {} with taxon level delimiter [{}]'.format(
                                                                            lineage,
                                                                            taxon_level_delimiter))

        processed_taxonomic_str = ''

        for taxon_str in lineage:
            taxon_str = taxon_str.rstrip(';')
            taxon_items = taxon_str.split(taxon_level_delimiter)

            scientific_name = taxon_items[-1]

            processed_taxonomic_str += scientific_name + ';'

        return processed_taxonomic_str

    def __init__(self, config):
        self.taxon_wsname = config.get('taxon-workspace-name')

    def process_taxonomic_str(self, taxonomic_str):
        logging.info('start processing taxonomic string: {}'.format(taxonomic_str))

        try:
            if not isinstance(taxonomic_str, str):
                raise ValueError('input taxonomic string is not a str type')

            # remove whitespaces
            taxonomic_str = taxonomic_str.replace(' ', '').replace('\t', '')

            if taxonomic_str.isalpha():
                return taxonomic_str

            # count non-alphanumeric characters
            delimiters = re.sub(r'\w+', '', taxonomic_str)
            delimiters = ''.join(set(delimiters))

            if len(delimiters) == 1:
                return taxonomic_str

            delimiter = csv.Sniffer().sniff(taxonomic_str).delimiter
            lineage = [x.strip() for x in taxonomic_str.split(delimiter)]
            logging.info('identified delimiter [{}] from the original taxonomic string'.format(
                                                                                        delimiter))

            if all(['__' in x for x in lineage]):
                taxon_level_delimiter = '__'
            elif all([':' in x for x in lineage]):
                taxon_level_delimiter = ':'
            else:
                taxon_level_delimiter = None

            if taxon_level_delimiter is None:
                return taxonomic_str

            processed_taxonomic_str = self._convert_taxonomic_str(lineage, taxon_level_delimiter)

        except Exception:
            logging.warning('failed to process taxonomic string')
            processed_taxonomic_str = taxonomic_str

        return processed_taxonomic_str
