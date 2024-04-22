import os
import uuid

from flask import Blueprint, json, request, send_file
from email_validator import validate_email, EmailNotValidError

from virus_discovery.model import VirusDiscoveryJob
from shared import config

virus_discovery_blueprint = Blueprint('virus_discovery', __name__)


@virus_discovery_blueprint.route('/virus_discovery', methods=['POST'])
def submit_virus_discovery_job():
    job_args = validate_form(request.form)
    if job_args is not None:
        filenames = save_protein_files(request.files)
        if 'pdb' in filenames:
            job_args['pdb-file'] = filenames['pdb']
            if 'xtc' in filenames:
                job_args['xtc-file'] = filenames['xtc']
            if PocketomeJob.create(job_args):
                return json.dumps({'status': 'success'})
            else:
                return json.dumps({'status': 'failed', 'error': 'Could not save to DB'}), 500
        else:
            return json.dumps({'status': 'failed', 'error': 'Please provide a valid PDB file'}), 400
    else:
        return json.dumps({'status': 'failed', 'error': 'Please fill in all the fields with the appropriate values'}), 400


@virus_discovery_blueprint.route('/virus_discovery/<int:job_id>/<string:job_hash>', methods=['GET'])
def virus_discovery_job_results(job_id=0, job_hash=None):
    if job_id != 0 and job_hash is not None:
        job = PocketomeJob.get(job_id)
        if job is not None:
            if job.verify_hash(job_hash):
                return json.dumps({'status': 'success', 'results': job.get_results()})
            else:
                return json.dumps({'status': 'failed', 'error': 'Unauthorized access'}), 400
        else:
            return json.dumps({'status': 'failed', 'error': 'Could not find a matching result'}), 400
    else:
        return json.dumps({'status': 'failed', 'error': 'Please provide all the necessary information'}), 400


@virus_discovery_blueprint.route('/virus_discovery/<int:job_id>/<string:job_hash>/download', methods=['GET'])
def virus_discovery_job_results_zipped(job_id=0, job_hash=None):
    if job_id != 0 and job_hash is not None:
        job = PocketomeJob.get(job_id)
        if job is not None:
            if job.verify_hash(job_hash):
                return send_file(job.get_zipped_results(), mimetype='application/zip', as_attachment=True)
            else:
                return json.dumps({'status': 'failed', 'error': 'Unauthorized access'}), 400
        else:
            return json.dumps({'status': 'failed', 'error': 'Could not find a matching result'}), 400
    else:
        return json.dumps({'status': 'failed', 'error': 'Please provide all the necessary information'}), 400


def validate_form(form):
    job_args = {}
    fields = [{'name': 'email', 'required': True, 'type': 'email'},
              {'name': 'protein-chain', 'required': True, 'type': 'str', 'len': 1},
              {'name': 'ligand-chain', 'required': True, 'type': 'str', 'len': 1},
              {'name': 'kflag', 'required': True, 'type': 'int', 'min': 0},
              {'name': 'dist', 'required': True, 'type': 'float', 'min': 0},
              {'name': 'sasa-threshold', 'required': True, 'type': 'float', 'min': 0, 'max': 1},
              {'name': 'dock-threshold', 'required': True, 'type': 'float', 'min': 0, 'max': 1}]
    for f in fields:
        if validate_field(form, f):
            job_args[f['name']] = form[f['name']]
        else:
            return None
    return job_args


def validate_field(form, field):
    if 'name' not in field:
        return False
    if 'type' not in field:
        return False
    if field['name'] in form:
        if field['type'] == 'str':
            if not validate_str(form, field):
                return False
        elif field['type'] == 'email':
            try:
                validate_email(form[field['name']])
            except EmailNotValidError:
                return False
        elif field['type'] == 'int':
            if not validate_num(form, field):
                return False
        elif field['type'] == 'float':
            if not validate_num(form, field, True):
                return False
        return True
    return not field['required']


def validate_str(form, field):
    if 'len' in field:
        if len(form[field['name']]) != field['len']:
            return False
    return True


def validate_num(form, field, isfloat=False):
    try:
        if isfloat:
            x = float(form[field['name']])
        else:
            x = int(form[field['name']])
        if 'min' in field:
            if x < field['min']:
                return False
        if 'max' in field:
            if x > field['max']:
                return False
    except ValueError:
        return False
    return True


def allowed_file(filename, filetype):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ['txt', filetype]


def save_protein_files(files):
    filetypes = ['pdb', 'xtc']
    filenames = {}
    for ft in filetypes:
        protein_file = ft + '-file'
        if protein_file in files:
            file = files[protein_file]
            if file.filename != '':
                if file and allowed_file(file.filename, ft):
                    filenames[ft] = str(uuid.uuid4()) + '.' + ft
                    filepath = os.path.join(config['app']['uploads_path'], filenames[ft])
                    file.save(filepath)
    return filenames
