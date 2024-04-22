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


class PocketomeJob(db.Model):
    __tablename__ = 'pocketome_jobs'

    id = db.Column('id', db.Integer, primary_key=True)
    submitted = db.Column('submitted', db.TIMESTAMP, nullable=False, default=datetime.now())
    started = db.Column('started', db.TIMESTAMP)
    completed = db.Column('completed', db.TIMESTAMP)
    user = db.Column('user', db.String(255), nullable=False)
    complex_pdb = db.Column('complex_pdb', db.String(255), nullable=False)
    chain_protein = db.Column('chain_protein', db.String(1), nullable=False)
    chain_ligand = db.Column('chain_ligand', db.String(1), nullable=False)
    complex_xtc = db.Column('complex_xtc', db.String(255))
    k_flag = db.Column('k_flag', db.Integer)
    dist = db.Column('dist', db.Numeric)
    sasa_threshold = db.Column('sasa_threshold', db.Numeric)
    dock_threshold = db.Column('dock_threshold', db.Numeric)
    extensive = db.Column('extensive', db.Integer, nullable=False, default=0)

    @staticmethod
    def get(job_id):
        return PocketomeJob.query.get(job_id)

    @staticmethod
    def create(job_args):
        complex_xtc = None
        if 'xtc-file' in job_args:
            complex_xtc = job_args['xtc-file']
        new_job = PocketomeJob(user=job_args['email'], complex_pdb=job_args['pdb-file'], complex_xtc=complex_xtc,
                               chain_protein=job_args['protein-chain'], chain_ligand=job_args['ligand-chain'],
                               k_flag=job_args['kflag'], dist=job_args['dist'],
                               sasa_threshold=job_args['sasa-threshold'], dock_threshold=job_args['dock-threshold'])
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
        root_path = config['app']['output_path'] + os.sep + str(self.id) + os.sep + '04-analysis' + os.sep
        visualisation_path = root_path + 'visualisation' + os.sep
        enrichment_path = root_path + 'enrichment_analysis' + os.sep
        centroids = []
        for f in sorted(os.listdir(visualisation_path)):
            if f.endswith('.pdb'):
                centroids.append(get_text_output(visualisation_path + f, False))
        return {
            'results': get_text_output(root_path + 'Results.txt'),
            'uni_prot_ids': get_text_output(root_path + 'UniProt_ids.txt'),
            'graphs': {
                'filtering': 'data:image/png;base64,'
                             + get_image_output(root_path + 'filtering.png'),
                'interactions': 'data:image/png;base64,'
                                + get_image_output(root_path + 'interactions.png'),
                'bp_bar': 'data:image/jpeg;base64,'
                          + get_image_output(enrichment_path + 'BP_barplot_third_level_plot.jpg'),
                'bp_goterm': 'data:image/jpeg;base64,'
                             + get_image_output(enrichment_path + 'GOTERM_BP_DIRECT_enrichment_plot.jpg'),
                'bp_gotree': 'data:image/jpeg;base64,'
                             + get_image_output(enrichment_path + 'BP_gotree_plot.jpg'),
                'cc_bar': 'data:image/jpeg;base64,'
                          + get_image_output(enrichment_path + 'CC_barplot_third_level_plot.jpg'),
                'cc_goterm': 'data:image/jpeg;base64,'
                             + get_image_output(enrichment_path + 'GOTERM_CC_DIRECT_enrichment_plot.jpg'),
                'cc_gotree': 'data:image/jpeg;base64,'
                             + get_image_output(enrichment_path + 'CC_gotree_plot.jpg'),
                'mf_bar': 'data:image/jpeg;base64,'
                             + get_image_output(enrichment_path + 'MF_barplot_third_level_plot.jpg'),
                'mf_goterm': 'data:image/jpeg;base64,'
                             + get_image_output(enrichment_path + 'GOTERM_MF_DIRECT_enrichment_plot.jpg'),
                'mf_gotree': 'data:image/jpeg;base64,'
                             + get_image_output(enrichment_path + 'MF_gotree_plot.jpg'),
                'pathway_enrichment': 'data:image/jpeg;base64,'
                             + get_image_output(enrichment_path + 'pathway_enrichment_plot.jpg'),
            },
            'centroids': centroids,
            'chains': {'protein': self.chain_protein, 'ligand': self.chain_ligand}
        }

    def get_zipped_results(self):
        root_path = config['app']['output_path'] + os.sep + str(self.id) + os.sep
        zip_path = root_path + 'pocketome-{}'.format(str(self.id))
        if not os.path.exists(zip_path + '.zip'):
            shutil.make_archive(zip_path, 'zip', root_path + '04-analysis')
        return zip_path + '.zip'
