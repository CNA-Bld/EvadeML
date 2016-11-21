import hashlib
import json
import logging

import requests

logger = logging.getLogger('pdfrate_online')


def sha1(filepath):
    with open(filepath, 'rb') as f:
        return hashlib.sha1(f.read()).hexdigest()


def load_data():
    try:
        with open('pdfrate_result.json', 'r') as f:
            return json.load(f)
    except:
        return {}


def save_data(result):
    with open('pdfrate_result.json', 'w') as f:
        json.dump(result, f)


try:
    with open('pdfrate_fp.json', 'r') as f:
        fp_data = json.load(f)
except:
    fp_data = {}


def add_fp_data(new_data):
    fp_data[new_data['fileinfo']['sha1']] = new_data
    with open('pdfrate_fp.json', 'w') as f:
        json.dump(fp_data, f)


def request_online(filename, retry=0):
    try:
        files = {'filesubmission': open(filename, 'rb')}
        r = requests.post('https://csmutz.com/pdfrate/submitfile', allow_redirects=False, files=files)
        data = r.json()
        _ = data['results'].values()[-1]
        return data
    except:
        if retry < 10:
            return request_online(filename, retry + 1)
        else:
            raise


def get_fitness_func(field='contagio_bm', threshold=50):
    def fitness_func(filenames):
        result = load_data()
        for filename in filenames:
            sha = sha1(filename)
            if sha in fp_data:
                result[filename] = fp_data[sha]
            else:
                result[filename] = request_online(filename)
                add_fp_data(result[filename])
            logger.info("%s: %s (%s)" % (filename, result[filename]['results'].values()[-1][field], field))
        save_data(result)
        return [float(result[filename]['results'].values()[-1][field]) > threshold for filename in filenames]

    return fitness_func
