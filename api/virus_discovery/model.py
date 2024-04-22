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
    genome = db.Column('genome', db.String(255), nullable=False)

    @staticmethod
    def get(job_id):
        return VirusDiscoveryJob.query.get(job_id)

    @staticmethod
    def create(job_args):
        new_job = VirusDiscoveryJob(user=job_args['email'], genome=job_args['genome'])
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
