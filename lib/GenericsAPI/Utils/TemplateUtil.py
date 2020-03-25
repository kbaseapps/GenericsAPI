import errno
import os
import uuid
import logging
import pandas as pd

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport
from GenericsAPI.Utils.SampleServiceUtil import SampleServiceUtil


class TemplateUtil:

    def _mkdir_p(self, path):
        """
        _mkdir_p: make directory for given path
        """
        if not path:
            return
        try:
            os.makedirs(path)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise

    def _generate_report(self, template_file, workspace_id):

        file_links = [{'path': template_file,
                       'name': os.path.basename(template_file),
                       'label': 'Chemical Abundance Matrix template file',
                       'description': 'use this file for Chemical Abundance Matrix uploader'}]

        message = 'Successfully created a template for Chemical Abundance Matrix uploader'
        report_params = {'message': message,
                         'workspace_id': workspace_id,
                         'file_links': file_links,
                         'report_object_name': 'chem_abun_template_' + str(uuid.uuid4())}

        kbase_report_client = KBaseReport(self.callback_url, token=self.token)
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output

    def _fetch_chemical_datas(self, chemical_data_included, chemical_ids_included):
        chemical_datas = list()

        rename_map = {'aggregate_mz': 'Aggregate M/Z',
                      'name': 'Compound Name',
                      'formula': 'Predicted Formula',
                      'smiles': 'Predicted Structure (smiles)',
                      'inchi': 'Predicted Structure (inchi)',
                      'inchikey': 'Predicted Structure (inchi-key)',
                      'mass': 'Theoretical Mass',
                      'retention_time': 'Retention Time',
                      'polarity': 'Polarity',
                      'kegg': 'KEGG',
                      'chembi': 'ChemBi',
                      'modelseed': 'ModelSEED',
                      }

        if chemical_data_included:
            items = [rename_map[key] for key in chemical_data_included if chemical_data_included[key]]
            chemical_datas.extend(items)

        if chemical_ids_included:
            items = [rename_map[key] for key in chemical_ids_included if chemical_ids_included[key]]
            chemical_datas.extend(items)

        if not chemical_datas:
            raise ValueError('Please provide at least one of chemical data or chemical ID')

        return chemical_datas

    def _create_template_file(self, chemical_datas, sample_names):

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        logging.info('Start building tempalte file in dir: {}'.format(output_directory))
        self._mkdir_p(output_directory)

        template_file = os.path.join(output_directory, 'chemical_abundance_matrix_template.xlsx')

        file_df = pd.DataFrame.from_dict([], orient='index', columns=chemical_datas + sample_names)
        file_df.index.name = 'ID (unique value)'

        file_df.to_excel(template_file)

        return template_file

    def __init__(self, config):
        self.callback_url = config['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.token = config['KB_AUTH_TOKEN']

        self.dfu = DataFileUtil(self.callback_url)
        self.sampleservice_util = SampleServiceUtil(config)

        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)

    def build_chemical_abundance_template(self, params):
        logging.info('Start building tempalte for Chemical Abundance:\n{}\n'.format(params))
        workspace_id = params.get('workspace_id')
        sample_set_ref = params.get('sample_set_ref')
        chemical_data_included = params.get('chemical_data_included')
        chemical_ids_included = params.get('chemical_ids_included')

        chemical_datas = self._fetch_chemical_datas(chemical_data_included, chemical_ids_included)

        sample_names = list()
        if not sample_set_ref:
            raise ValueError('Please provide a Sample Set object')
        else:
            sample_set = self.dfu.get_objects(
                    {"object_refs": [sample_set_ref]})['data'][0]['data']

            samples = sample_set['samples']

            for sample in samples:
                sample_id = sample['id']
                sample_data = self.sampleservice_util.get_sample(sample_id)

                sample_names.append(sample_data['name'])

        template_file = self._create_template_file(chemical_datas, sample_names)

        report_output = self._generate_report(template_file, workspace_id)

        return report_output
