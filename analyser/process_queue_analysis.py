import os
from json import load
from subprocess import run
from time import sleep

from helpers import build_args, check_sequencing_args, notify_user
from virus_discovery_job import VirusDiscoveryJob

with open('config.json', encoding='utf-8') as fp:
    config = load(fp)


def run_analysis(args):
    if check_sequencing_args(args):
        # Usage (Single End): lib/analysis.sh [options] -s file'
        # Usage (Paired End): lib/analysis.sh [options] -f forward_file -r reverse_file
        cwd = os.path.dirname(os.path.abspath(__file__))
        analysis_exec = os.path.join(cwd, 'lib/analysis.sh')
        if args['single_paired'] == 'single':
            res = run([analysis_exec,
                       '-n', str(args['sample_name']),
                       '-t', str(args['threads']),
                       '-o', str(args['output_dir']),
                       '-s', os.path.join(config['args']['uploads'], str(args['forward_file']))])

            if res.returncode == 0:
                print('Analysis executed.')
            else:
                print('Analysis failed.', res)

        elif args['single_paired'] == 'pair':
            print(f"will run: {analysis_exec} '-n',"
                  f" {str(args['sample_name'])}, "
                  f"'-t', {str(args['threads'])}, "
                  f"'-o', {str(args['output_dir'])}, "
                  f"'-f', {os.path.join(config['args']['uploads'], str(args['forward_file']))}, "
                  f"'-r', {os.path.join(config['args']['uploads'], str(args['reverse_file']))}")
            res = run([analysis_exec,
                       '-n', str(args['sample_name']),
                       '-t', str(args['threads']),
                       '-o', str(args['output_dir']),
                       '-f', os.path.join(config['args']['uploads'], str(args['forward_file'])),
                       '-r', os.path.join(config['args']['uploads'], str(args['reverse_file']))])
            if res.returncode == 0:
                print('Analysis executed.')
            else:
                print('Analysis failed.', res)


if __name__ == '__main__':
    discvir = VirusDiscoveryJob(config)
    print('Processing queue...')
    while True:
        jobs_batch = discvir.get_analysis_jobs(config['queue']['batch'])
        if not jobs_batch:
            print('No analysis jobs found')
        for job in jobs_batch:
            discvir.mark_as_started_analysis(job['id'])
            job_args = build_args(config=config, job_id=job['id'], genome=job['genome'],
                                  paired=job['paired'], sample_name=job['sample_name'],
                                  forward_file=job['forward_file'], reverse_file=job['reverse_file'],
                                  adapter=job['adapter'], min_len=job['min_len'], window=job['window'])
            run_analysis(job_args)
            timestamp = discvir.mark_as_completed_analysis(job['id'])
            if not notify_user(config, job['user'], job['id'], timestamp, 'analysis'):
                print('Email (analysis) not sent to {}'.format(job['user']))
        sleep(config['queue']['sleep'])
