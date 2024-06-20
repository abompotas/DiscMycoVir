import os
import re
import uuid
from json import dumps

from email_validator import validate_email, EmailNotValidError
from flask import Blueprint, request, send_file

from shared import config
from virus_discovery.model import VirusDiscoveryJob

virus_discovery_blueprint = Blueprint('virus_discovery', __name__)


@virus_discovery_blueprint.route('/job', methods=['POST'])
def submit_job():
    job_args = parse_job_form()
    if job_args is not None:
        if VirusDiscoveryJob.create(job_args):
            return dumps({'status': 'success'})
        else:
            return dumps({'status': 'failed', 'error': 'Could not save to DB'}), 500
    else:
        return dumps({'status': 'failed', 'error': 'Please fill in all the fields with the appropriate values'}), 400


@virus_discovery_blueprint.route('/analysis/<int:job_id>/<string:job_hash>', methods=['GET'])
def analysis_reports(job_id=0, job_hash=None):
    if job_id != 0 and job_hash is not None:
        job = VirusDiscoveryJob.get(job_id)
        if job is not None:
            if job.verify_hash(job_hash, 'analysis'):
                return dumps({'status': 'success', 'results': job.get_analysis_reports()})
            else:
                return dumps({'status': 'failed', 'error': 'Unauthorized access'}), 403
        else:
            return dumps({'status': 'failed', 'error': 'Could not find a matching result'}), 404
    else:
        return dumps({'status': 'failed', 'error': 'Please provide all the necessary information'}), 400


@virus_discovery_blueprint.route('/analysis/<int:job_id>/<string:job_hash>/download', methods=['GET'])
def analysis_reports_zipped(job_id=0, job_hash=None):
    if job_id != 0 and job_hash is not None:
        job = VirusDiscoveryJob.get(job_id)
        if job is not None:
            if job.verify_hash(job_hash, 'analysis'):
                zip_path = job.get_analysis_reports_zipped()
                if zip_path is not None:
                    return send_file(zip_path, mimetype='application/zip', as_attachment=True)
                else:
                    return dumps({'status': 'failed', 'error': 'Could not create the archive'}), 500
            else:
                return dumps({'status': 'failed', 'error': 'Unauthorized access'}), 403
        else:
            return dumps({'status': 'failed', 'error': 'Could not find a matching result'}), 404
    else:
        return dumps({'status': 'failed', 'error': 'Please provide all the necessary information'}), 400


@virus_discovery_blueprint.route('/trimming/<int:job_id>/<string:job_hash>', methods=['PUT'])
def trimming(job_id=0, job_hash=None):
    if job_id != 0 and job_hash is not None:
        job = VirusDiscoveryJob.get(job_id)
        if job is not None:
            if job.verify_hash(job_hash, 'analysis'):
                job_args = parse_trimming_form()
                if job_args is not None:
                    if job.schedule_for_trimming(job_args):
                        return dumps({'status': 'success'})
                    else:
                        return dumps({'status': 'failed', 'error': 'This job has already been submitted'}), 400
                else:
                    return dumps({'status': 'failed',
                                  'error': 'Please fill in all the fields with the appropriate values'}), 400
            else:
                return dumps({'status': 'failed', 'error': 'Unauthorized access'}), 403
        else:
            return dumps({'status': 'failed', 'error': 'Could not find a matching result'}), 404
    else:
        return dumps({'status': 'failed', 'error': 'Please provide all the necessary information'}), 400


@virus_discovery_blueprint.route('/discovery/<int:job_id>/<string:job_hash>', methods=['PUT'])
def discovery(job_id=0, job_hash=None):
    if job_id != 0 and job_hash is not None:
        job = VirusDiscoveryJob.get(job_id)
        if job is not None:
            if job.verify_hash(job_hash, 'analysis'):
                if job.schedule_for_discovery():
                    return dumps({'status': 'success'})
                else:
                    return dumps({'status': 'failed', 'error': 'This job has already been submitted'}), 400
            else:
                return dumps({'status': 'failed', 'error': 'Unauthorized access'}), 403
        else:
            return dumps({'status': 'failed', 'error': 'Could not find a matching result'}), 404
    else:
        return dumps({'status': 'failed', 'error': 'Please provide all the necessary information'}), 400


@virus_discovery_blueprint.route('/results/<int:job_id>/<string:job_hash>', methods=['GET'])
def results(job_id=0, job_hash=None):
    if job_id != 0 and job_hash is not None:
        job = VirusDiscoveryJob.get(job_id)
        if job is not None:
            if job.verify_hash(job_hash):
                final_results = job.get_final_results()
                if final_results is not None:
                    return dumps({'status': 'success', 'results': final_results})
                else:
                    return dumps({'status': 'failed', 'error': 'Could not retrieve results'}), 500
            else:
                return dumps({'status': 'failed', 'error': 'Unauthorized access'}), 403
        else:
            return dumps({'status': 'failed', 'error': 'Could not find a matching result'}), 404
    else:
        return dumps({'status': 'failed', 'error': 'Please provide all the necessary information'}), 400


@virus_discovery_blueprint.route('/results/<int:job_id>/<string:job_hash>/download', methods=['GET'])
def results_zipped(job_id=0, job_hash=None):
    if job_id != 0 and job_hash is not None:
        job = VirusDiscoveryJob.get(job_id)
        if job is not None:
            if job.verify_hash(job_hash, 'results'):
                zip_path = job.get_all_results_zipped()
                if zip_path is not None:
                    return send_file(zip_path, mimetype='application/zip', as_attachment=True)
                else:
                    return dumps({'status': 'failed', 'error': 'Could not create the archive'}), 500
            else:
                return dumps({'status': 'failed', 'error': 'Unauthorized access'}), 403
        else:
            return dumps({'status': 'failed', 'error': 'Could not find a matching result'}), 404
    else:
        return dumps({'status': 'failed', 'error': 'Please provide all the necessary information'}), 400


def parse_job_form():
    job_args = {
        'email': parse_email(),
        'sample_name': parse_sample_name(),
        'genome': parse_file('reference_genome')
    }
    if 'sequencing_technology' in request.form:
        job_args['sequencing_technology'] = request.form['sequencing_technology']
        if request.form['sequencing_technology'] == 'single':
            job_args['single_file'] = parse_file('single_file')
        elif request.form['sequencing_technology'] == 'paired':
            job_args['forward_file'] = parse_file('forward_file')
            job_args['reverse_file'] = parse_file('reverse_file')
        else:
            return None
    else:
        return None
    if None in job_args.values():
        return None
    else:
        return job_args


def parse_trimming_form():
    job_args = {
        'adapter': parse_file('adapter', config['defaults']['adapter']),
        'window': parse_sliding_window(config['defaults']['sliding_window']),
        'min_len': parse_int('min_length', 1, 1000, config['defaults']['min_length'])
    }
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
        filename = ''
        file = request.files[name]
        if file.filename != '':
            if file and allowed_file(file.filename):
                filename = str(uuid.uuid4())
                filepath = os.path.join(config['app']['uploads_path'], filename)
                file.save(filepath)
        return filename


def allowed_file(filename):
    filetypes = ['txt', 'fa', 'fq', 'fasta', 'fastq', 'gz']
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in filetypes
