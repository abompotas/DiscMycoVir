import os

from json import load
from time import sleep
from helpers import VirusDiscoveryJob, notify_user


with open('config.json', encoding='utf-8') as fp:
    config = load(fp)


def build_args(job_id, genome, paired=0, sample_name="", forward_file="", reverse_file="",
               adapter=None, min_len=None, window=None):
    """Usage (Single End): ./virus_discovery_pipeline.sh [options] -c reference_genome -s file"
    "Usage (Paired End): ./virus_discovery_pipeline.sh [options] -c reference_genome -f forward_file -r reverse_file\
    """
    if paired:
        single_paired = "pair"
    else:
        single_paired = "single"
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
        # "-f	Paired forward input file" "-s	Single end input file"
        'forward_file': forward_file,
        # Paired reverse input file
        'reverse_file': reverse_file,
        # Output directory
        'output_dir': os.path.join(config['args']['output'], str(job_id))
    }
    return args


def run_pipeline(args):
    print(args)
    if args['single_paired'] == 'single':
        if not args['forward_file']:
            print('ERROR: MISSING single input file')
    elif args['single_paired'] == 'pair':
        if not (args['reverse_file'] and args['forward_file']):
            print('ERROR: MISSING Reverse and forward input files')
    else:
        print('ERROR: Unknown sequencing technology')
    #step 1: TrimmomaticSE (paired-single)
    #step 2: fastq -> fasta
    #step 3: Trinity (paired-single)
    #step 4: BWA-MEM
    #step 5: SAMTOOLS
    #step 6: BLAST


if __name__ == '__main__':
    discvir_job = VirusDiscoveryJob(config)
    print('Processing queue...')
    while True:
        jobs_batch = discvir_job.get_jobs(config['queue']['batch'])
        if not jobs_batch:
            print('No jobs found')
        for job in jobs_batch:
            discvir_job.mark_as_started(job['id'])
            job_args = build_args(job_id=job['id'], genome=job['genome'],
                                  paired=job['paired'], sample_name=job['sample_name'],
                                  forward_file=job['forward_file'], reverse_file=job['reverse_file'],
                                  adapter=job['adapter'], min_len=job['min_len'], window=job['window'])
            run_pipeline(job_args)
            discvir_job.mark_as_completed(job['id'])
            # if not notify_user(config, job['user'], job['id']):
            #     print('Email not sent to {}'.format(job['user']))
        sleep(config['queue']['sleep'])
