
import logging
import csv
import re


class TaxonUtil:

    ROOT_CANDIDATES = ['bacteria', 'archaea', 'viridae', 'eukaryota', 'virus', 'eukarya']

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

        # remove empty trailing slots
        processed_taxonomic_str = processed_taxonomic_str.rstrip(';') + ';'

        # remove unclassified scientific names
        if 'unclassified' in processed_taxonomic_str.lower():
            processed_taxonomic_str = self._remove_unclassified(processed_taxonomic_str)

        # prepend genus name to species epithet if missing
        if processed_taxonomic_str.count(';') == 7:
            processed_taxonomic_str = self._prepend_genus_name(processed_taxonomic_str)

        processed_taxonomic_str = self._remove_root(processed_taxonomic_str)

        logging.info('converted taxonomic string: {}'.format(processed_taxonomic_str))

        return processed_taxonomic_str

    def _remove_unclassified(self, taxonomic_str, delimiter=';'):
        # remove unclassified scientific names
        scientific_names = taxonomic_str.split(delimiter)
        for idx, scientific_name in enumerate(scientific_names):
            if 'unclassified' in scientific_name.lower():
                logging.info('removing unclassified scientific name: {}'.format(
                                                                            scientific_names[idx]))
                scientific_names[idx] = ''
        taxonomic_str = delimiter.join(scientific_names)

        return taxonomic_str

    def _prepend_genus_name(self, taxonomic_str, delimiter=';'):
        # prepend genus name to species epithet if missing
        scientific_names = taxonomic_str.split(delimiter)
        genus_name = scientific_names[-3]
        species_name = scientific_names[-2]

        if not species_name.lower().startswith(genus_name.lower()):
            logging.info("prepend genus name [{}] to species epithet [{}]".format(genus_name,
                                                                                  species_name))
            species_name = genus_name + ' ' + species_name
            scientific_names[-2] = species_name

        taxonomic_str = delimiter.join(scientific_names)

        return taxonomic_str

    def _remove_root(self, taxonomic_str, delimiter=';'):
        # make sure taxonomic str starting with one of ROOT_CANDIDATES
        remove_root = False

        for root_candidate in self.ROOT_CANDIDATES:
            if root_candidate in taxonomic_str.lower():
                remove_root = True
                break

        if remove_root:
            scientific_names = taxonomic_str.split(delimiter)
            starting_name = ''
            for scientific_name in scientific_names:
                if scientific_name.lower() in self.ROOT_CANDIDATES:
                    starting_name = scientific_name
                    break
            if starting_name:
                idx = scientific_names.index(starting_name)
                if idx != 0:
                    msg = 'Removing root scientific name(s)\n'
                    msg += 'New starting scientific name is [{}]'.format(starting_name)
                    logging.info(msg)
                    scientific_names = scientific_names[idx:]
            taxonomic_str = delimiter.join(scientific_names)

        return taxonomic_str

    def _process_taxonomic_str(self, taxonomic_str, delimiter):

        if taxonomic_str.endswith(delimiter):
            taxonomic_str = taxonomic_str.replace(delimiter, ';')
        else:
            taxonomic_str = taxonomic_str.replace(delimiter, ';') + ';'
        # remove unclassified scientific names
        if 'unclassified' in taxonomic_str.lower():
            taxonomic_str = self._remove_unclassified(taxonomic_str)
        # prepend genus name to species epithet if missing
        if taxonomic_str.count(';') == 7:
            taxonomic_str = self._prepend_genus_name(taxonomic_str)

        taxonomic_str = self._remove_root(taxonomic_str)

        return taxonomic_str

    def __init__(self, config):
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        self.taxon_wsname = config.get('taxon-workspace-name')

    def process_taxonomic_str(self, taxonomic_str):
        '''
        parse user's input taxonomic string into a standardized syntax (ideally, a 7 slot string)
        general rules:
            1. return original string if we cannot find a pattern to parse
            2. return original string if user's taxonomic string has 8 or more than 8 slots
               (string with replaced delimiters (e.g. commas replaced with semicolons))
            3. append trailing semicolon (if missing) if user's taxonomic string has 6 slots
            4. prepend genus name to species epithet IF missing (when s is specified or length = 7)
            5. start parsed taxonomic string with one of:
               Bacteria, Archaea, Viridae, Eukaryota; Virus; Eukarya; Viruses
        '''

        logging.info('start processing taxonomic string: {}'.format(taxonomic_str))

        try:
            if not isinstance(taxonomic_str, str):
                raise ValueError('input taxonomic string is not a str type')

            # remove whitespaces
            taxonomic_str = taxonomic_str.strip()

            if taxonomic_str.isalpha():
                return taxonomic_str + ';'

            # count non-alphanumeric characters
            delimiters = re.sub(r'[a-zA-Z0-9]+', '', taxonomic_str)
            delimiters = ''.join(set(delimiters))

            if len(delimiters) == 1:
                taxonomic_str = self._process_taxonomic_str(taxonomic_str, delimiters)
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
                taxonomic_str = self._process_taxonomic_str(taxonomic_str, delimiter)
                return taxonomic_str

            processed_taxonomic_str = self._convert_taxonomic_str(lineage, taxon_level_delimiter)

        except Exception:
            logging.warning('failed to process taxonomic string')
            processed_taxonomic_str = taxonomic_str

        return processed_taxonomic_str
