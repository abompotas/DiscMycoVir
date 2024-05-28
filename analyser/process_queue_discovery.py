import os
from json import load
from subprocess import run
from time import sleep

from helpers import build_args, check_sequencing_args, notify_user
from virus_discovery_job import VirusDiscoveryJob

with open('config.json', encoding='utf-8') as fp:
    config = load(fp)


def run_trimming(args, input_dir=''):
    if check_sequencing_args(args):
        # Usage (Single End): lib/trimming.sh [options] -s file'
        # Usage (Paired End): lib/trimming.sh [options] -f forward_file -r reverse_file
        cwd = os.path.dirname(os.path.abspath(__file__))
        trimming_exec = os.path.join(cwd, 'lib/trimming.sh')
        if args['single_paired'] == 'single':
            res = run([trimming_exec,
                       '-n', str(args['sample_name']),
                       '-t', str(args['threads']),
                       '-a', str(args['adapter']),
                       '-w', str(args['sliding_window']),
                       '-l', str(args['min_len']),
                       '-o', str(args['output_dir']),
                       '-s', os.path.join(input_dir, str(args['forward_file']))])
            if res.returncode == 0:
                print('Discovery executed.')
            else:
                print('Discovery failed.', res)

        elif args['single_paired'] == 'pair':
            res = run([trimming_exec,
                       '-n', str(args['sample_name']),
                       '-t', str(args['threads']),
                       '-a', str(args['adapter']),
                       '-w', str(args['sliding_window']),
                       '-l', str(args['min_len']),
                       '-o', str(args['output_dir']),
                       '-f', os.path.join(input_dir, str(args['forward_file'])),
                       '-r', os.path.join(input_dir, str(args['reverse_file']))])
            if res.returncode == 0:
                print('Discovery executed.')
            else:
                print('Discovery failed.', res)


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
        input_dir = os.path.join(config['args']['outputs'], 'trimming')
        if not os.path.exists(input_dir):
            input_dir = config['args']['uploads']
        cwd = os.path.dirname(os.path.abspath(__file__))
        discovery_exec = os.path.join(cwd, 'lib/discovery.sh')
        if args['single_paired'] == 'single':
            res = run([discovery_exec,
                       '-n', str(args['sample_name']),
                       '-t', str(args['threads']),
                       '-m', str(args['max_memory']),
                       '-o', str(args['output_dir']),
                       '-g', os.path.join(config['args']['uploads'], str(args['ref_genome'])),
                       '-s', os.path.join(input_dir, str(args['forward_file']))])
            if res.returncode == 0:
                print('Discovery executed.')
            else:
                print('Discovery failed.', res)

        elif args['single_paired'] == 'pair':
            res = run([discovery_exec,
                       '-n', str(args['sample_name']),
                       '-t', str(args['threads']),
                       '-m', str(args['max_memory']),
                       '-o', str(args['output_dir']),
                       '-g', os.path.join(config['args']['uploads'], str(args['ref_genome'])),
                       '-f', os.path.join(input_dir, str(args['forward_file'])),
                       '-r', os.path.join(input_dir, str(args['reverse_file']))])
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
