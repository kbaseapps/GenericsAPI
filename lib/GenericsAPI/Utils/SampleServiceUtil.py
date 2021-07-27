
import json
import requests
import uuid
import logging

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.SampleServiceClient import SampleService


class SampleServiceUtil:

    def _fetch_attri_from_meta(self, meta_data):
        attributes = list()

        for key, value in meta_data.items():
            attri = {'attribute': key, 'source': 'SampleService'}
            if 'units' in value:
                attri['unit'] = value.get('units')

            attributes.append(attri)

        return attributes

    def __init__(self, config):
        self.callback_url = config['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.token = config['KB_AUTH_TOKEN']
        self.srv_wiz_url = config['srv-wiz-url']
        self.sample_url = config.get('kbase-endpoint') + '/sampleservice'
        self.dfu = DataFileUtil(self.callback_url)
        self.sample_ser = SampleService(self.sample_url)

        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)

    def get_sample_service_url(self):

        # payload = {
        #     "method": "ServiceWizard.get_service_status",
        #     "id": '',
        #     "params": [{"module_name": "SampleService", "version": "dev"}],  # TODO: change to beta/release
        #     "version": "1.1"
        # }

        # sw_resp = requests.post(url=self.srv_wiz_url, data=json.dumps(payload))
        # wiz_resp = sw_resp.json()
        # if wiz_resp.get('error'):
        #     raise RuntimeError("ServiceWizard Error - " + str(wiz_resp['error']))

        # return wiz_resp['result'][0]['url']

        return self.sample_url

    def get_sample(self, sample_id, version=None):

        sample_url = self.get_sample_service_url()

        headers = {
            "Authorization": self.token,
            "Content-Type": "application/json"
        }
        params = {
            "id": sample_id,
            "version": version
        }
        payload = {
            "method": "SampleService.get_sample",
            "id": str(uuid.uuid4()),
            "params": [params],
            "version": "1.1"
        }
        resp = requests.post(url=sample_url, headers=headers, data=json.dumps(payload))
        resp_json = resp.json()
        if resp_json.get('error'):
            raise RuntimeError(f"Error from SampleService - {resp_json['error']}")
        sample = resp_json['result'][0]

        # sample = self.sample_ser.get_sample(params)[0]

        return sample

    def save_sample(self, sample):
        """
        copied from sample_uploader
        (https://github.com/kbaseapps/sample_uploader/blob/master/lib/sample_uploader/utils/sample_utils.py#L79)
        """
        logging.info('start saving sample')

        sample_url = self.get_sample_service_url()
        headers = {"Authorization": self.token}
        params = {
            "sample": sample,
            "prior_version": None,
        }
        payload = {
            "method": "SampleService.create_sample",
            "id": str(uuid.uuid4()),
            "params": [params],
            "version": "1.1"
        }
        resp = requests.post(url=sample_url, headers=headers, data=json.dumps(payload, default=str))
        if not resp.ok:
            raise RuntimeError(f'Error from SampleService - {resp.text}')
        resp_json = resp.json()
        if resp_json.get('error'):
            raise RuntimeError(f"Error from SampleService - {resp_json['error']}")
        sample_id = resp_json['result'][0]['id']
        return sample_id

    def upsert_sample(self, upsert_sample, sample_id):
        """
        upsert sample data into existing sample
        """
        logging.info('start upserting sample')

        old_sample = self.get_sample(sample_id)

        try:
            old_sample.get('node_tree')[0].get('meta_user').update(upsert_sample)
            sample_id = self.save_sample(old_sample)
        except Exception:
            raise ValueError('failed to upsert sample')

        return sample_id

    def create_data_link(self, sample_id, upa, version, dataid=None, node='root'):
        """
        'params': [{'upa': '48710/1/1',
                    'dataid': 'column2',
                    'id': id_,
                    'version': 1,
                    'node': 'root',}]
        """
        logging.info('start creating data link')

        sample_url = self.get_sample_service_url()
        headers = {"Authorization": self.token}
        params = {
            "id": sample_id,
            "upa": upa,
            "version": version,
            "dataid": dataid,
            "node": node,
        }
        payload = {
            "method": "SampleService.create_data_link",
            "id": str(uuid.uuid4()),
            "params": [params],
            "version": "1.1"
        }
        resp = requests.post(url=sample_url, headers=headers, data=json.dumps(payload, default=str))
        if not resp.ok:
            raise RuntimeError(f'Error from SampleService - {resp.text}')
        resp_json = resp.json()
        if resp_json.get('error'):
            raise RuntimeError(f"Error from SampleService - {resp_json['error']}")

    def sample_set_to_attribute_mapping(self, sample_set_ref):

        logging.info('start converting sample set to attribute mapping')

        sample_set = self.dfu.get_objects(
                    {"object_refs": [sample_set_ref]})['data'][0]['data']

        samples = sample_set['samples']

        am_data = {'ontology_mapping_method': "SampleService", 'instances': {}}

        attributes = list()

        sample_datas = []
        for sample in samples:
            sample_id = sample.get('id')
            version = sample.get('version')

            sample_data = self.get_sample(sample_id, version=version)

            node_tree = sample_data.get('node_tree', [{}])

            for node in node_tree:
                meta_controlled = node.get('meta_controlled')
                meta_user = node.get('meta_user')

                if meta_controlled:
                    attributes_controlled = self._fetch_attri_from_meta(meta_controlled)
                    attributes += attributes_controlled

                if meta_user:
                    attributes_user = self._fetch_attri_from_meta(meta_user)
                    attributes += attributes_user

            sample_data['sample_id'] = sample_id
            sample_datas.append(sample_data)

        attributes = [i for n, i in enumerate(attributes) if i not in attributes[n + 1:]]
        attributes = [{'attribute': 'sample_id', 'source': 'SampleService'},
                      {'attribute': 'version', 'source': 'SampleService'},
                      {'attribute': 'type', 'source': 'SampleService'},
                      {'attribute': 'parent', 'source': 'SampleService'}] + attributes

        am_data['attributes'] = attributes

        instances = am_data['instances']

        for sample_data in sample_datas:
            instance = list()
            node_tree = sample_data.get('node_tree', [{}])

            for node in node_tree:

                meta_controlled = node.get('meta_controlled', {})
                meta_user = node.get('meta_user', {})

                for attribute in attributes:
                    attri_name = attribute['attribute']

                    if attri_name in ['sample_id', 'version']:
                        instance.append(sample_data.get(attri_name))
                    elif attri_name in ['type', 'parent']:
                        instance.append(node.get(attri_name))
                    else:
                        attri_value = meta_user.get(attri_name, {}).get('value') or meta_controlled.get(attri_name, {}).get('value')
                        instance.append(attri_value)

            instance = [str(i) for i in instance]
            instances[sample_data['name']] = instance

        return am_data
