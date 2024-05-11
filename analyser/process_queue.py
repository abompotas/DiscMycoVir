import os
import smtplib
import ssl
from hashlib import sha256
from json import load
from time import sleep
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

import mysql.connector

from virus_detection import run_pipeline

with open('config.json',encoding='utf-8') as fp:
    config = load(fp)
    db = mysql.connector.connect(
        host=config['db']['host'],
        user=config['db']['username'],
        password=config['db']['password'],
        database=config['db']['schema']
    )
def build_args(job_id, genome, paired = False, adapter = "NexteraPE-PE.fa:2:30:10",
               sample_name="", min_len="50", window = "5:20", forward_file="", reverse_file=""):
    """Usage (Single End): ./virus_discovery_pipeline.sh [options] -c reference_genome -s file"
    "Usage (Paired End): ./virus_discovery_pipeline.sh [options] -c reference_genome -f forward_file -r reverse_file\
    """
    if paired:
        single_paired = "pair"
    else:
        single_paired = "single"
    args = {
    # Number of threads for parallel execution
        'threads': config['args']['cpus'],
        'single_paired': single_paired,
    # Adapter for Trimmomatic
        'adapter' : adapter,
    # Sliding window for Trimmomatic
        'sliding_window' : window,
    # Minimum length for Trimmomatic
        'min_len' : min_len,
    # Max memory to be used by Trinity
        'max_memory': config['args']['max_mem'],
    # Sequence type to be used by Trinity
        'seq_type' : config['args']['seq_type'],
    # Sample name
        'sample_name' : sample_name,
    # Reference genome
        'ref_genome' : genome,
    #"-f	Paired forward input file" "-s	Single end input file"
        'forward_file' : forward_file,
    # Paired reverse input file
        'reverse_file' : reverse_file,
    # Output directory
        'output_dir' : config['args']['outfolder'] + os.sep + str(job_id) + os.sep}
    return args


def get_jobs(limit):
    cursor = db.cursor()
    cursor.execute(f"SELECT * FROM virus_discovery_jobs WHERE started IS NULL ORDER BY id ASC LIMIT 0,{limit}")
    jobs = cursor.fetchall()
    cursor.close()
    db.commit()
    return jobs


def mark_as_started(job_id):
    cursor = db.cursor()
    cursor.execute('UPDATE virus_discovery_jobs SET started=CURRENT_TIMESTAMP() WHERE id={}'.format(job_id))
    cursor.close()
    db.commit()


def mark_as_completed(job_id):
    cursor = db.cursor()
    cursor.execute('UPDATE virus_discovery_jobs SET completed=CURRENT_TIMESTAMP() WHERE id={}'.format(job_id))
    cursor.close()
    db.commit()


def notify_user(email, job_id):
    context = ssl.create_default_context()
    try:
        smtp = smtplib.SMTP_SSL(config['smtp']['server'], config['smtp']['port'], context=context)
        smtp.login(config['smtp']['username'], config['smtp']['password'])
        message = MIMEMultipart()
        message['From'] = config['smtp']['address']
        message['To'] = email
        message['Subject'] = 'Virtuous Pocketome - Search completed'
        results_hash = sha256('{}%%{}'.format(job_id, email).encode('utf-8')).hexdigest()
        results_url = '{}/{}/{}'.format(config['webbapp'], job_id, results_hash)
        results_text = '<html><body><p>VirtuousPocketome has analyzed your request</p>' \
                       '<p>View the results: <a href="{0}">{0}</a></p></body></html>'.format(results_url)
        message.attach(MIMEText(results_text, 'html'))
        smtp.sendmail(config['smtp']['address'], email, message.as_string())
    except Exception as e:
        print(str(e))
    finally:
        smtp.quit()


if __name__ == '__main__':
    print('Processing queue...')
    while True:
        jobs = get_jobs(config['queue']['batch'])
        if not jobs:
            print("No jobs found")
        for job in jobs:

            job_id, submitted, started, completed, user, genome, paired, adapter, \
                sample_name, min_len, window,forward_file, reverse_file = job
            mark_as_started(job_id)
            job_args = build_args(job_id,genome,paired,adapter,sample_name,
                                  min_len,window,forward_file,reverse_file)

            run_pipeline(job_args)
            mark_as_completed(job_id)
            # notify_user(job[4], job[0])
        sleep(config['queue']['sleep'])
