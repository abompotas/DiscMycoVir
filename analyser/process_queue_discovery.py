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
    #
    # Usage (Single End): lib/discovery.sh [options] -g reference_genome -s file'
    # Usage (Paired End): lib/discovery.sh [options] -g reference_genome -f forward_file -r reverse_file
    if check_sequencing_args(args):
        cwd = os.path.dirname(os.path.abspath(__file__))
        discovery_exec = os.path.join(cwd, 'lib/discovery.sh')
        if args['single_paired'] == 'single':
            input_format = 'fq'
            forward_file = os.path.join(args['output_job'], 'trimming', '{}.trimmed'.format(args['forward_file']))
            if not os.path.exists(forward_file):
                input_format = args['input_format']
                forward_file = os.path.join(config['args']['uploads'], args['forward_file'])
            res = run([discovery_exec,
                       '-n', args['sample_name'],
                       '-t', args['threads'],
                       '-m', args['max_memory'],
                       '-o', args['output_job'],
                       '-q', input_format,
                       '-g', os.path.join(config['args']['uploads'], args['ref_genome']),
                       '-s', forward_file])
            if res.returncode == 0:
                print('[{}|{}]: Discovery executed.'.format(job['id'], job['sample_name']))
            else:
                print('[{}|{}]: Discovery failed: '.format(job['id'], job['sample_name'], str(res)))

        elif args['single_paired'] == 'pair':
            input_format = 'fq'
            forward_file = os.path.join(args['output_job'], 'trimming', '{}.trimmed'.format(args['forward_file']))
            reverse_file = os.path.join(args['output_job'], 'trimming', '{}.trimmed'.format(args['reverse_file']))
            if not os.path.exists(forward_file):
                input_format = args['input_format']
                forward_file = os.path.join(config['args']['uploads'], args['forward_file'])
                reverse_file = os.path.join(config['args']['uploads'], args['reverse_file'])
            res = run([discovery_exec,
                       '-n', args['sample_name'],
                       '-t', args['threads'],
                       '-m', args['max_memory'],
                       '-o', args['output_job'],
                       '-q', input_format,
                       '-g', os.path.join(config['args']['uploads'], args['ref_genome']),
                       '-f', forward_file,
                       '-r', reverse_file])
            if res.returncode == 0:
                print('[{}|{}]: Discovery executed.'.format(job['id'], job['sample_name']))
            else:
                print('[{}|{}]: Discovery failed: '.format(job['id'], job['sample_name'], str(res)))


if __name__ == '__main__':
    discvir = VirusDiscoveryJob(config)
    while True:
        jobs_batch = discvir.get_discovery_jobs(config['queue']['batch'])
        for job in jobs_batch:
            discvir.mark_as_started_discovery(job['id'])
            job_args = build_args(config=config, job_id=job['id'], genome=job['genome'],
                                  fasta=job['fasta'], paired=job['paired'], sample_name=job['sample_name'],
                                  forward_file=job['forward_file'], reverse_file=job['reverse_file'],
                                  adapter=job['adapter'], min_len=job['min_len'], window=job['window'])
            run_discovery(job_args)
            timestamp = discvir.mark_as_completed_assembly(job['id'])
        sleep(config['queue']['sleep'])
