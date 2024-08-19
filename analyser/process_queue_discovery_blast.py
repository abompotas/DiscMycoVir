import os
from json import load
from subprocess import Popen
from time import sleep

from helpers import notify_user
from virus_discovery_job import VirusDiscoveryJob

with open('config.json', encoding='utf-8') as fp:
    config = load(fp)


if __name__ == '__main__':
    discvir = VirusDiscoveryJob(config)
    while True:
        jobs_batch = discvir.get_blast_jobs(config['queue']['batch'])
        for job in jobs_batch:
            discvir.mark_as_started_blast(job['id'])
            cwd = os.path.dirname(os.path.abspath(__file__))
            discovery_exec = os.path.join(cwd, 'lib/discovery_blast.sh')
            discovery_dir = os.path.join(config['args']['output'], str(job['id']), 'discovery')
            Popen([discovery_exec, '-d', discovery_dir])

        jobs_batch = discvir.get_pending_blast_jobs(config['queue']['batch'])
        for job in jobs_batch:
            discovery_flag = os.path.join(config['args']['output'], str(job['id']), 'discovery', 'blast_finished')
            if os.path.exists(discovery_flag):
                timestamp = discvir.mark_as_completed_discovery(job['id'])
                if not notify_user(config, job['user'], job['id'], timestamp, 'discovery'):
                    print('[{}|{}]: Email (discovery) not sent to {}'.format(job['id'], job['sample_name'], job['user']))

        sleep(config['queue']['sleep'])
