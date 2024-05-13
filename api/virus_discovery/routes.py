import os
import re
import uuid
from json import dumps

from email_validator import validate_email, EmailNotValidError
from flask import Blueprint, request, send_file

from shared import config
from virus_discovery.model import VirusDiscoveryJob

virus_discovery_blueprint = Blueprint('virus_discovery', __name__)


@virus_discovery_blueprint.route('/virus-discovery', methods=['POST'])
def submit_virus_discovery_job():
    job_args = parse_form()
    if job_args is not None:
        if VirusDiscoveryJob.create(job_args):
            return dumps({'status': 'success'})
        else:
            return dumps({'status': 'failed', 'error': 'Could not save to DB'}), 500
    else:
        return dumps({'status': 'failed', 'error': 'Please fill in all the fields with the appropriate values'}), 400


@virus_discovery_blueprint.route('/virus-discovery/<int:job_id>/<string:job_hash>', methods=['GET'])
def virus_discovery_job_results(job_id=0, job_hash=None):
    if job_id != 0 and job_hash is not None:
        job = VirusDiscoveryJob.get(job_id)
        if job is not None:
            if job.verify_hash(job_hash):
                return dumps({'status': 'success', 'results': job.get_results()})
            else:
                return dumps({'status': 'failed', 'error': 'Unauthorized access'}), 400
        else:
            return dumps({'status': 'failed', 'error': 'Could not find a matching result'}), 400
    else:
        return dumps({'status': 'failed', 'error': 'Please provide all the necessary information'}), 400


@virus_discovery_blueprint.route('/virus-discovery/<int:job_id>/<string:job_hash>/download', methods=['GET'])
def virus_discovery_job_results_zipped(job_id=0, job_hash=None):
    if job_id != 0 and job_hash is not None:
        job = VirusDiscoveryJob.get(job_id)
        if job is not None:
            if job.verify_hash(job_hash):
                return send_file(job.get_zipped_results(), mimetype='application/zip', as_attachment=True)
            else:
                return dumps({'status': 'failed', 'error': 'Unauthorized access'}), 400
        else:
            return dumps({'status': 'failed', 'error': 'Could not find a matching result'}), 400
    else:
        return dumps({'status': 'failed', 'error': 'Please provide all the necessary information'}), 400


def parse_form():
    job_args = {
        'email': parse_email(),
        'sample_name': parse_sample_name(),
        'genome': parse_file('reference_genome'),
        'adapter': parse_file('adapter', config['defaults']['adapter']),
        'window': parse_sliding_window(config['defaults']['sliding_window']),
        'min_len': parse_int('min_length', 1, 1000, config['defaults']['min_length'])
    }
    if 'sequencing_technology' in request.form:
        job_args['sequencing_technology'] = request.form['sequencing_technology']
        if request.form['sequencing_technology'] == 'single':
            job_args['single_file'] = parse_file('single_file')
        elif request.form['sequencing_technology'] == 'paired':
            job_args['forward_file'] = parse_file('forward_file')
            job_args['reverse_file'] = parse_file('reverse_file')
        else:
            job_args['sequencing_technology'] = None
    else:
        job_args['sequencing_technology'] = None
    if None in job_args.values():
        return None
    else:
        return job_args


def parse_email():
    if 'email' not in request.form:
        return None
    else:
        try:
            validate_email(request.form['email'])
        except EmailNotValidError:
            return None
        return request.form['email']


def parse_sample_name():
    if 'sample_name' not in request.form:
        return None
    else:
        if re.search('^\w+$', request.form['sample_name']) is None:
            return None
        return request.form['sample_name']


def parse_sliding_window(default=None):
    if 'sliding_window' not in request.form:
        return default
    else:
        if len(request.form['sliding_window']) == 0:
            return default
        if re.search('^\d+:\d+$', request.form['sliding_window']) is None:
            return None
        return request.form['sliding_window']


def parse_int(field, min_w=None, max_w=None, default=None):
    if field not in request.form:
        return default
    else:
        if len(request.form['sliding_window']) == 0:
            return default
        try:
            x = int(request.form[field])
            if min_w is not None:
                if x < min_w:
                    return None
            if max_w is not None:
                if x > max_w:
                    return None
        except ValueError:
            return None
        return x


def parse_file(name, default=None):
    if name not in request.files:
        return default
    else:
        filepath = None
        file = request.files[name]
        if file.filename != '':
            if file and allowed_file(file.filename):
                filename = '{}'.format(str(uuid.uuid4()))
                filepath = os.path.join(config['app']['uploads_path'], filename)
                file.save(filepath)
        return filepath


def allowed_file(filename):
    filetypes = ['txt', 'fa', 'fq', 'fasta', 'fastq', 'gz']
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in filetypes
