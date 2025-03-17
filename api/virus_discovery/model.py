import os
import shutil
import uuid
from datetime import datetime
from hashlib import sha256

from shared import config, db
from .helpers import get_text_output, parse_blast_xml


class VirusDiscoveryJob(db.Model):
    __tablename__ = 'virus_discovery_jobs'

    id = db.Column('id', db.Integer, primary_key=True)
    stage = db.Column('stage', db.Integer, nullable=False, default=False)
    submitted = db.Column('submitted', db.TIMESTAMP, nullable=False, default=datetime.now())
    started_analysis = db.Column('started_analysis', db.TIMESTAMP)
    completed_analysis = db.Column('completed_analysis', db.TIMESTAMP)
    started_trimming = db.Column('started_trimming', db.TIMESTAMP)
    completed_trimming = db.Column('completed_trimming', db.TIMESTAMP)
    started_discovery = db.Column('started_discovery', db.TIMESTAMP)
    completed_assembly = db.Column('completed_assembly', db.TIMESTAMP)
    started_blast = db.Column('started_blast', db.TIMESTAMP)
    completed_discovery = db.Column('completed_discovery', db.TIMESTAMP)
    user = db.Column('user', db.String(255), nullable=False)
    sample_name = db.Column('sample_name', db.String(255), nullable=False)
    genome = db.Column('genome', db.String(255), nullable=False)
    fasta = db.Column('fasta', db.Boolean, nullable=False, default=False)
    paired = db.Column('paired', db.Boolean, nullable=False, default=False)
    adapter = db.Column('adapter', db.String(255), nullable=True)
    min_len = db.Column('min_len', db.Integer, nullable=True)
    window = db.Column('window', db.String(255), nullable=True)
    forward_file = db.Column('forward_file', db.String(255), nullable=False)
    reverse_file = db.Column('reverse_file', db.String(255))

    @staticmethod
    def get(job_id):
        return VirusDiscoveryJob.query.get(job_id)

    @staticmethod
    def create(job_args):
        new_job = None
        fa = True
        if job_args['input_format'] == 'fq':
            fa = False
        if job_args['sequencing_technology'] == 'single':
            new_job = VirusDiscoveryJob(user=job_args['email'], sample_name=job_args['sample_name'],
                                        genome=job_args['genome'], fasta=fa, paired=False,
                                        forward_file=job_args['single_file'], reverse_file=None)
        elif job_args['sequencing_technology'] == 'paired':
            new_job = VirusDiscoveryJob(user=job_args['email'], sample_name=job_args['sample_name'],
                                        genome=job_args['genome'], fasta=fa, paired=True,
                                        forward_file=job_args['forward_file'], reverse_file=job_args['reverse_file'])
        if new_job:
            db.session.add(new_job)
            db.session.commit()
            return True
        return False

    def is_scheduled(self):
        if self.stage == 0:
            if self.completed_analysis is not None:
                return True
        elif self.stage == 1:
            if self.completed_trimming is not None:
                return True
        else:
            if self.completed_discovery is not None:
                return True
        return False

    def schedule_for_trimming(self, job_args):
        if self.is_scheduled():
            self.stage = 1
            self.started_trimming = None
            self.completed_trimming = None
            self.adapter = job_args['adapter']
            self.window = job_args['window']
            self.min_len = job_args['min_len']
            db.session.commit()
            return True
        return False

    def schedule_for_discovery(self):
        if self.is_scheduled():
            self.stage = 2
            self.started_discovery = None
            self.completed_discovery = None
            db.session.commit()
            return True
        return False

    def verify_hash(self, job_hash, stage='analysis'):
        if stage == 'analysis':
            check_str = '{}%%{}%%{}'.format(self.id, self.user, self.completed_analysis).encode('utf-8')
        else:
            check_str = '{}%%{}%%{}'.format(self.id, self.user, self.completed_discovery).encode('utf-8')
        calc_hash = sha256(check_str).hexdigest()
        if calc_hash == job_hash:
            return True
        return False

    def get_analysis_reports(self):
        analysis_dir = os.path.join(config['app']['output_path'], str(self.id), 'fastqc_analysis')
        if os.path.exists(analysis_dir):
            if self.paired:
                forward_file_report = '{}_fastqc.html'.format(self.forward_file.split('.')[0])
                reverse_file_report = '{}_fastqc.html'.format(self.reverse_file.split('.')[0])
                forward_file_report_path = os.path.join(analysis_dir, forward_file_report)
                reverse_file_report_path = os.path.join(analysis_dir, reverse_file_report)
                return [get_text_output(forward_file_report_path), get_text_output(reverse_file_report_path)]
            else:
                forward_file_report = '{}_fastqc.html'.format(self.forward_file.split('.')[0])
                forward_file_report_path = os.path.join(analysis_dir, forward_file_report)
                return [get_text_output(forward_file_report_path)]
        return None

    def get_analysis_reports_zipped(self):
        if os.path.exists(config['app']['output_path']):
            tmp_file = os.path.join(config['app']['output_path'], str(uuid.uuid4()))
            zip_dir = os.path.join(config['app']['output_path'], str(self.id))
            zip_file = os.path.join(zip_dir, 'virus_discovery-{}-fastqc_analysis.zip'.format(str(self.id)))
            if os.path.exists(zip_file):
                os.remove(zip_file)
            shutil.make_archive(tmp_file, 'zip', root_dir=zip_dir, base_dir='fastqc_analysis')
            shutil.move('{}.zip'.format(tmp_file), zip_file)
            return zip_file
        return None

    def get_final_results(self):
        blast_file = os.path.join(config['app']['output_path'], str(self.id), 'discovery', 'output_blast.xml')
        if os.path.exists(blast_file):
            return parse_blast_xml(blast_file)
        return None

    def get_all_results_zipped(self):
        if os.path.exists(config['app']['output_path']):
            tmp_file = os.path.join(config['app']['output_path'], str(uuid.uuid4()))
            zip_dir = os.path.join(config['app']['output_path'], str(self.id))
            zip_file = os.path.join(zip_dir, 'virus_discovery-{}.zip'.format(str(self.id)))
            if os.path.exists(zip_file):
                os.remove(zip_file)
            shutil.make_archive(tmp_file, 'zip', root_dir=zip_dir)
            shutil.move('{}.zip'.format(tmp_file), zip_file)
            return zip_file
        return None
