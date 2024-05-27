import os
from json import load
from subprocess import run
from time import sleep

from helpers import VirusDiscoveryJob, notify_user

with open('config.json', encoding='utf-8') as fp:
    config = load(fp)


def build_args(job_id, genome, paired=0, sample_name='', forward_file='', reverse_file='',
               adapter=None, min_len=None, window=None):
    if paired:
        single_paired = 'pair'
    else:
        single_paired = 'single'
    if adapter is None:
        adapter = config['args']['adapter']
    if min_len is None:
        min_len = config['args']['min_length']
    if window is None:
        window = config['args']['sliding_window']
    args = {
        # Number of threads for parallel execution
        'threads': config['args']['threads'],
        # Max memory to be used by Trinity
        'max_memory': config['args']['max_memory'],
        'single_paired': single_paired,
        # Adapter for Trimmomatic
        'adapter': adapter,
        # Sliding window for Trimmomatic
        'sliding_window': window,
        # Minimum length for Trimmomatic
        'min_len': min_len,
        # Sample name
        'sample_name': sample_name,
        # Reference genome
        'ref_genome': genome,
        # '-f	Paired forward input file' '-s	Single end input file'
        'forward_file': forward_file,
        # Paired reverse input file
        'reverse_file': reverse_file,
        # Output directory
        'output_dir': os.path.join(config['args']['output'], str(job_id))
    }
    return args


def check_sequencing_args(args):
    if args['single_paired'] == 'single':
        if not args['forward_file']:
            print('ERROR: MISSING single input file')
            return False
    elif args['single_paired'] == 'pair':
        if not (args['reverse_file'] and args['forward_file']):
            print('ERROR: MISSING Reverse and forward input files')
            return False
    else:
        print('ERROR: Unknown sequencing technology')
        return False
    return True


def run_analysis(args, input_dir=''):
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
                       '-s', os.path.join(input_dir, str(args['forward_file']))])

            if res.returncode == 0:
                print('Analysis executed.')
            else:
                print('Analysis failed.', res)

        elif args['single_paired'] == 'pair':
            res = run([analysis_exec,
                       '-n', str(args['sample_name']),
                       '-t', str(args['threads']),
                       '-o', str(args['output_dir']),
                       '-f', os.path.join(input_dir, str(args['forward_file'])),
                       '-r', os.path.join(input_dir, str(args['reverse_file']))])
            if res.returncode == 0:
                print('Analysis executed.')
            else:
                print('Analysis failed.', res)


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


def run_discovery(args, input_dir=''):
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
            res = run([discovery_exec,
                       '-n', str(args['sample_name']),
                       '-t', str(args['threads']),
                       '-m', str(args['max_memory']),
                       '-o', str(args['output_dir']),
                       '-g', os.path.join(input_dir, str(args['ref_genome'])),
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
                       '-g', os.path.join(input_dir, str(args['ref_genome'])),
                       '-f', os.path.join(input_dir, str(args['forward_file'])),
                       '-r', os.path.join(input_dir, str(args['reverse_file']))])
            if res.returncode == 0:
                print('Discovery executed.')
            else:
                print('Discovery failed.', res)


def process_analysis_queue(discvir):
    jobs_batch = discvir.get_analysis_jobs(config['queue']['batch'])
    if not jobs_batch:
        print('No analysis jobs found')
    for job in jobs_batch:
        discvir.mark_as_started_analysis(job['id'])
        job_args = build_args(job_id=job['id'], genome=job['genome'],
                              paired=job['paired'], sample_name=job['sample_name'],
                              forward_file=job['forward_file'], reverse_file=job['reverse_file'],
                              adapter=job['adapter'], min_len=job['min_len'], window=job['window'])
        run_analysis(job_args, config['args']['uploads'])
        timestamp = discvir.mark_as_completed_analysis(job['id'])
        if not notify_user(config, job['user'], job['id'], timestamp, 'analysis'):
            print('Email (analysis) not sent to {}'.format(job['user']))


def process_trimming_queue(discvir):
    jobs_batch = discvir.get_trimming_jobs(config['queue']['batch'])
    if not jobs_batch:
        print('No trimming jobs found')
    for job in jobs_batch:
        discvir.mark_as_started_trimming(job['id'])
        job_args = build_args(job_id=job['id'], genome=job['genome'],
                              paired=job['paired'], sample_name=job['sample_name'],
                              forward_file=job['forward_file'], reverse_file=job['reverse_file'],
                              adapter=job['adapter'], min_len=job['min_len'], window=job['window'])
        run_trimming(job_args, config['args']['uploads'])
        run_analysis(job_args, os.path.join(config['args']['outputs'], 'trimming'))
        timestamp = discvir.mark_as_completed_trimming(job['id'])
        if not notify_user(config, job['user'], job['id'], timestamp, 'analysis'):
            print('Email (trimming) not sent to {}'.format(job['user']))


def process_discovery_queue(discvir):
    jobs_batch = discvir.get_discovery_jobs(config['queue']['batch'])
    if not jobs_batch:
        print('No discovery jobs found')
    for job in jobs_batch:
        discvir.mark_as_started_discovery(job['id'])
        job_args = build_args(job_id=job['id'], genome=job['genome'],
                              paired=job['paired'], sample_name=job['sample_name'],
                              forward_file=job['forward_file'], reverse_file=job['reverse_file'],
                              adapter=job['adapter'], min_len=job['min_len'], window=job['window'])
        run_discovery(job_args, os.path.join(config['args']['outputs'], 'trimming'))
        timestamp = discvir.mark_as_completed_discovery(job['id'])
        if not notify_user(config, job['user'], job['id'], timestamp, 'discovery'):
            print('Email (discovery) not sent to {}'.format(job['user']))
        sleep(config['queue']['sleep'])


if __name__ == '__main__':
    jobs_helper = VirusDiscoveryJob(config)
    print('Processing queue...')
    while True:
        process_analysis_queue(jobs_helper)
        process_trimming_queue(jobs_helper)
        process_discovery_queue(jobs_helper)
