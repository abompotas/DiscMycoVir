import os
import base64
import shutil

from datetime import datetime
from hashlib import sha256

from shared import config, db


def get_text_output(filepath, array=True):
    with open(filepath, 'r') as txt:
        if array:
            contents = txt.readlines()
        else:
            contents = txt.read()
    return contents


def get_image_output(filepath):
    with open(filepath, 'rb') as img:
        contents = img.read()
    return base64.b64encode(contents).decode("utf-8")


class VirusDiscoveryJob(db.Model):
    __tablename__ = 'virus_discovery_jobs'

    id = db.Column('id', db.Integer, primary_key=True)
    submitted = db.Column('submitted', db.TIMESTAMP, nullable=False, default=datetime.now())
    started = db.Column('started', db.TIMESTAMP)
    completed = db.Column('completed', db.TIMESTAMP)
    user = db.Column('user', db.String(255), nullable=False)
    sample_name = db.Column('sample_name', db.String(255), nullable=False)
    genome = db.Column('genome', db.String(255), nullable=False)
    paired = db.Column('paired', db.Boolean, nullable=False, default=False)
    adapter = db.Column('adapter', db.String(255), nullable=False)
    min_len = db.Column('min_len', db.Integer, nullable=False, default=50)
    window = db.Column('window', db.String(255), nullable=False)
    forward_file = db.Column('forward_file', db.String(255), nullable=False)
    reverse_file = db.Column('reverse_file', db.String(255))

    @staticmethod
    def get(job_id):
        return VirusDiscoveryJob.query.get(job_id)

    @staticmethod
    def create(job_args):
        new_job = None
        if job_args['sequencing_technology'] == 'single':
            new_job = VirusDiscoveryJob(user=job_args['email'], sample_name=job_args['sample_name'],
                                        genome=job_args['genome'], paired=False,
                                        forward_file=job_args['single_file'], reverse_file=None,
                                        adapter=job_args['adapter'], min_len=job_args['min_len'],
                                        window=job_args['window'])
        elif job_args['sequencing_technology'] == 'paired':
            new_job = VirusDiscoveryJob(user=job_args['email'], sample_name=job_args['sample_name'],
                                        genome=job_args['genome'], paired=True,
                                        forward_file=job_args['forward_file'], reverse_file=job_args['reverse_file'],
                                        adapter=job_args['adapter'], min_len=job_args['min_len'],
                                        window=job_args['window'])
        if new_job:
            db.session.add(new_job)
            db.session.commit()
            return True
        return False

    def verify_hash(self, job_hash):
        calc_hash = sha256('{}%%{}'.format(self.id, self.user).encode('utf-8')).hexdigest()
        if calc_hash == job_hash:
            return True
        return False

    def get_results(self):
        root_path = os.path.join(config['app']['output_path'], str(self.id))
        visualisation_path = root_path + 'visualisation' + os.sep
        centroids = []
        for f in sorted(os.listdir(visualisation_path)):
            if f.endswith('.pdb'):
                centroids.append(get_text_output(visualisation_path + f, False))
        return {
            'results': get_text_output(root_path + 'Results.txt'),
            'graphs': {
                'filtering': 'data:image/png;base64,'
                             + get_image_output(root_path + 'filtering.png'),
            },
            'centroids': centroids
        }

    def get_zipped_results(self):
        root_path = os.path.join(config['app']['output_path'], str(self.id))
        results_path = os.path.join(root_path, 'results')
        zip_path = os.path.join(root_path, 'virus_discovery-{}'.format(str(self.id)))
        zip_file = '{}.zip'.format(zip_path)
        if not os.path.exists(zip_file):
            shutil.make_archive(zip_path, 'zip', results_path)
        return zip_file
