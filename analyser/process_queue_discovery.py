import os
from json import load
from subprocess import run
from time import sleep

from helpers import build_args, check_sequencing_args, notify_user
from virus_discovery_job import VirusDiscoveryJob

with open('config.json', encoding='utf-8') as fp:
    config = load(fp)


def run_discovery(args):
    # TODO: .sh to python
    # step 1: fastq -> fasta NOT NEEDED
    # step 2: Trinity (paired-single)
    # step 3: BWA-MEM
    # step 4: SAMTOOLS
    # step 5: BLAST
    #
    # Usage (Single End): lib/discovery.sh [options] -g reference_genome -s file'
    # Usage (Paired End): lib/discovery.sh [options] -g reference_genome -f forward_file -r reverse_file
    if check_sequencing_args(args):
        cwd = os.path.dirname(os.path.abspath(__file__))
        discovery_exec = os.path.join(cwd, 'lib/discovery.sh')
        if args['single_paired'] == 'single':
            forward_file = os.path.join(args['output_job'], 'trimming', '{}.trimmed'.format(args['forward_file']))
            if not os.path.exists(forward_file):
                forward_file = os.path.join(config['args']['uploads'], args['forward_file'])
            res = run([discovery_exec,
                       '-n', args['sample_name'],
                       '-t', args['threads'],
                       '-m', args['max_memory'],
                       '-o', args['output_job'],
                       '-g', os.path.join(config['args']['uploads'], args['ref_genome']),
                       '-s', forward_file])
            if res.returncode == 0:
                print('Discovery executed.')
            else:
                print('Discovery failed.', res)

        elif args['single_paired'] == 'pair':
            forward_file = os.path.join(args['output_job'], 'trimming', '{}.trimmed'.format(args['forward_file']))
            reverse_file = os.path.join(args['output_job'], 'trimming', '{}.trimmed'.format(args['reverse_file']))
            if not os.path.exists(forward_file):
                forward_file = os.path.join(config['args']['uploads'], args['forward_file'])
                reverse_file = os.path.join(config['args']['uploads'], args['reverse_file'])
            res = run([discovery_exec,
                       '-n', args['sample_name'],
                       '-t', args['threads'],
                       '-m', args['max_memory'],
                       '-o', args['output_job'],
                       '-g', os.path.join(config['args']['uploads'], args['ref_genome']),
                       '-f', forward_file,
                       '-r', reverse_file])
            if res.returncode == 0:
                print('Discovery executed.')
            else:
                print('Discovery failed.', res)


if __name__ == '__main__':
    discvir = VirusDiscoveryJob(config)
    print('Processing queue...')
    while True:
        jobs_batch = discvir.get_discovery_jobs(config['queue']['batch'])
        if not jobs_batch:
            print('No discovery jobs found')
        for job in jobs_batch:
            discvir.mark_as_started_discovery(job['id'])
            job_args = build_args(config=config, job_id=job['id'], genome=job['genome'],
                                  paired=job['paired'], sample_name=job['sample_name'],
                                  forward_file=job['forward_file'], reverse_file=job['reverse_file'],
                                  adapter=job['adapter'], min_len=job['min_len'], window=job['window'])
            run_discovery(job_args)
            timestamp = discvir.mark_as_completed_discovery(job['id'])
            if not notify_user(config, job['user'], job['id'], timestamp, 'discovery'):
                print('Email (discovery) not sent to {}'.format(job['user']))
        sleep(config['queue']['sleep'])
