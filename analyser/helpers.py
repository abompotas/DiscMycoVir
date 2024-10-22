import os
import smtplib
import ssl
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from hashlib import sha256

import mysql.connector


def build_args(config, job_id, genome, fasta=0, paired=0, sample_name='', forward_file='', reverse_file='',
               adapter=None, min_len=None, window=None):
    if fasta:
        input_format = 'fa'
    else:
        input_format = 'fq'
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
        'threads': str(config['args']['threads']),
        # Max memory to be used by Trinity
        'max_memory': str(config['args']['max_memory']),
        'input_format': input_format,
        'single_paired': single_paired,
        # Adapter for Trimmomatic
        'adapter': str(adapter),
        # Sliding window for Trimmomatic
        'sliding_window': str(window),
        # Minimum length for Trimmomatic
        'min_len': str(min_len),
        # Sample name
        'sample_name': str(sample_name),
        # Reference genome
        'ref_genome': str(genome),
        # '-f	Paired forward input file' '-s	Single end input file'
        'forward_file': str(forward_file),
        # Paired reverse input file
        'reverse_file': str(reverse_file),
        # Output directory
        'output_job': os.path.join(config['args']['output'], str(job_id))
    }
    return args


def check_sequencing_args(args, stage='analysis'):
    if args['single_paired'] == 'single':
        if not args['forward_file']:
            print('[{}] ERROR: MISSING single input file'.format(stage))
            return False
    elif args['single_paired'] == 'pair':
        if not (args['reverse_file'] and args['forward_file']):
            print('[{}] ERROR: MISSING Reverse and forward input files'.format(stage))
            return False
    else:
        print('[{}] ERROR: Unknown sequencing technology'.format(stage))
        return False
    return True


def notify_user(config, email, job_id, timestamp, stage='analysis'):
    notified = False
    context = ssl.create_default_context()
    try:
        smtp = smtplib.SMTP_SSL(config['smtp']['server'], config['smtp']['port'], context=context)
    except Exception as e:
        print(str(e))
        return notified
    try:
        smtp.login(config['smtp']['username'], config['smtp']['password'])
        message = MIMEMultipart()
        message['From'] = config['smtp']['address']
        message['To'] = email
        if stage == 'analysis':
            message['Subject'] = 'Virus Discovery - FastQC analysis completed'
            results_hash = sha256('{}%%{}%%{}'.format(job_id, email, timestamp).encode('utf-8')).hexdigest()
            results_url = '{}/trimming/{}/{}'.format(config['webbapp'], job_id, results_hash)
            results_text = '''<html><body><p>Virus Discovery has analyzed your request</p>
                <p>View the results: <a href="{0}">{0}</a></p></body></html>'''.format(results_url)
        else:
            message['Subject'] = 'Virus Discovery - Search completed'
            results_hash = sha256('{}%%{}%%{}'.format(job_id, email, timestamp).encode('utf-8')).hexdigest()
            results_url = '{}/results/{}/{}'.format(config['webbapp'], job_id, results_hash)
            results_text = '''<html><body><p>Virus Discovery has analyzed your request</p>
                <p>View the results: <a href="{0}">{0}</a></p></body></html>'''.format(results_url)
        message.attach(MIMEText(results_text, 'html'))
        smtp.sendmail(config['smtp']['address'], email, message.as_string())
        notified = True
    except Exception as e:
        print(str(e))
    finally:
        smtp.quit()
    return notified
