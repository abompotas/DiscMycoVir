import os
from json import load
from subprocess import run
from time import sleep

from helpers import build_args, check_sequencing_args, notify_user
from virus_discovery_job import VirusDiscoveryJob

with open('config.json', encoding='utf-8') as fp:
    config = load(fp)


def run_trimming(args):
    if check_sequencing_args(args):
        # Usage (Single End): lib/trimming.sh [options] -s file'
        # Usage (Paired End): lib/trimming.sh [options] -f forward_file -r reverse_file
        cwd = os.path.dirname(os.path.abspath(__file__))
        trimming_exec = os.path.join(cwd, 'lib/trimming.sh')
        if args['single_paired'] == 'single':
            res = run([trimming_exec,
                       '-n', args['sample_name'],
                       '-t', args['threads'],
                       '-a', os.path.join(config['args']['uploads'], args['adapter']),
                       '-w', args['sliding_window'],
                       '-l', args['min_len'],
                       '-o', args['output_job'],
                       '-i', config['args']['uploads'],
                       '-s', args['forward_file']])
            if res.returncode == 0:
                print('Trimming executed.')
            else:
                print('Trimming failed.', res)

        elif args['single_paired'] == 'pair':
            res = run([trimming_exec,
                       '-n', args['sample_name'],
                       '-t', args['threads'],
                       '-a', os.path.join(config['args']['uploads'], args['adapter']),
                       '-w', args['sliding_window'],
                       '-l', args['min_len'],
                       '-o', args['output_job'],
                       '-i', config['args']['uploads'],
                       '-f', args['forward_file'],
                       '-r', args['reverse_file']])
            if res.returncode == 0:
                print('Trimming executed.')
            else:
                print('Trimming failed.', res)


def run_analysis(args):
    if check_sequencing_args(args):
        # Usage (Single End): lib/analysis.sh [options] -s file'
        # Usage (Paired End): lib/analysis.sh [options] -f forward_file -r reverse_file
        cwd = os.path.dirname(os.path.abspath(__file__))
        analysis_exec = os.path.join(cwd, 'lib/analysis.sh')
        if args['single_paired'] == 'single':
            res = run([analysis_exec,
                       '-n', args['sample_name'],
                       '-t', args['threads'],
                       '-o', args['output_job'],
                       '-s', os.path.join(args['output_job'], 'trimming', args['forward_file'])])
            if res.returncode == 0:
                print('Analysis executed.')
            else:
                print('Analysis failed.', res)

        elif args['single_paired'] == 'pair':
            res = run([analysis_exec,
                       '-n', args['sample_name'],
                       '-t', args['threads'],
                       '-o', args['output_job'],
                       '-f', os.path.join(args['output_job'], 'trimming', args['forward_file']),
                       '-r', os.path.join(args['output_job'], 'trimming', args['reverse_file'])])
            if res.returncode == 0:
                print('Analysis executed.')
            else:
                print('Analysis failed.', res)


if __name__ == '__main__':
    discvir = VirusDiscoveryJob(config)
    print('Processing queue...')
    while True:
        jobs_batch = discvir.get_trimming_jobs(config['queue']['batch'])
        if not jobs_batch:
            print('No trimming jobs found')
        for job in jobs_batch:
            discvir.mark_as_started_trimming(job['id'])
            job_args = build_args(config=config, job_id=job['id'], genome=job['genome'],
                                  paired=job['paired'], sample_name=job['sample_name'],
                                  forward_file=job['forward_file'], reverse_file=job['reverse_file'],
                                  adapter=job['adapter'], min_len=job['min_len'], window=job['window'])
            run_trimming(job_args)
            discvir.mark_as_completed_trimming(job['id'])
            discvir.mark_as_started_analysis(job['id'])
            run_analysis(job_args)
            timestamp = discvir.mark_as_completed_analysis(job['id'])
            if not notify_user(config, job['user'], job['id'], timestamp, 'analysis'):
                print('Email (trimming) not sent to {}'.format(job['user']))
        sleep(config['queue']['sleep'])
