
import logging
import csv
import re


class TaxonUtil:

    def _fetch_taxon_level(self, taxon_char):

        taxon_level_mapping = {'l': 'Life', 'd': 'Domain', 'k': 'Kingdom', 'p': 'Phylum',
                               'c': 'Class', 'o': 'Order', 'f': 'Family', 'g': 'Genus',
                               's': 'Species'}

        return taxon_level_mapping.get(taxon_char.lower(), taxon_char)

    def __init__(self, config):
        self.taxon_wsname = config['taxon-workspace-name']

    def process_taxonomic_str(self, taxonmic_str):
        logging.info('start processing taxonmic string: {}'.format(taxonmic_str))

        try:
            if not isinstance(taxonmic_str, str):
                raise ValueError('input taxonomic string is not a str type')

            # remove whitespaces
            taxonmic_str = taxonmic_str.replace(' ', '').replace('\t', '')

            if taxonmic_str.isalpha():
                return taxonmic_str

            # count non-alphanumeric characters
            delimiters = re.sub(r'\w+', '', taxonmic_str)
            delimiters = ''.join(set(delimiters))

            if len(delimiters) == 1:
                return taxonmic_str

            delimiter = csv.Sniffer().sniff(taxonmic_str).delimiter
            lineage = [x.strip() for x in taxonmic_str.split(delimiter)]

            if all(['__' in x for x in lineage]):
                taxon_level_delimiter = '__'
            elif all([':' in x for x in lineage]):
                taxon_level_delimiter = ':'
            else:
                taxon_level_delimiter = None

            if taxon_level_delimiter is None:
                return taxonmic_str

            processed_taxonomic_str = ''
            for taxon_str in lineage:
                taxon_str = taxon_str.rstrip(';')
                taxon_items = taxon_str.split(taxon_level_delimiter)

                scientific_name = taxon_items[-1]
                taxon_level_char = taxon_items[0]

                processed_taxonomic_str += scientific_name + ';'

        except Exception:
            logging.warning('failed to process taxonomic string')
            processed_taxonomic_str = taxonmic_str

        return processed_taxonomic_str
