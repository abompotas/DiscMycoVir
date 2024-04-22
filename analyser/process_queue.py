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

with open('config.json') as fp:
    config = load(fp)
    db = mysql.connector.connect(
        host=config['db']['host'],
        user=config['db']['username'],
        password=config['db']['password'],
        database=config['db']['schema']
    )


def build_args(job_id, complex_pdb, chain_protein, chain_ligand,
               complex_xtc=None, k_flag=None, dist=None, sasa_threshold=None, dock_threshold=None, extensive=None):
    args = {
        'complexPDB': config['args']['uploads'] + os.sep + complex_pdb,
        'chainProtein': chain_protein,
        'chainLigand': chain_ligand,
        'db2bcompared': config['args']['db2bcompared'],
        'outfolder': config['args']['outfolder'] + os.sep + str(job_id) + os.sep,
        'cpus': config['args']['cpus']
    }
    if complex_xtc is not None:
        args['complexXTC'] = config['args']['uploads'] + os.sep + complex_xtc
    if k_flag is not None:
        args['kflag'] = k_flag
    else:
        args['kflag'] = config['args']['kflag']
    if dist is not None:
        args['dist'] = dist
    else:
        args['dist'] = config['args']['dist']
    if sasa_threshold is not None:
        args['sasaThreshold'] = sasa_threshold
    else:
        args['sasaThreshold'] = config['args']['sasaThreshold']
    if dock_threshold is not None:
        args['dockThreshold'] = dock_threshold
    else:
        args['dockThreshold'] = config['args']['dockThreshold']
    if extensive:
        args['extensive'] = extensive
    else:
        args['extensive'] = config['args']['extensive']
    return args


def get_jobs(limit):
    cursor = db.cursor()
    cursor.execute('SELECT * FROM pocketome_jobs WHERE started IS NULL ORDER BY id ASC LIMIT 0,{}'.format(limit))
    jobs = cursor.fetchall()
    cursor.close()
    db.commit()
    return jobs


def mark_as_started(job_id):
    cursor = db.cursor()
    cursor.execute('UPDATE pocketome_jobs SET started=CURRENT_TIMESTAMP() WHERE id={}'.format(job_id))
    cursor.close()
    db.commit()


def mark_as_completed(job_id):
    cursor = db.cursor()
    cursor.execute('UPDATE pocketome_jobs SET completed=CURRENT_TIMESTAMP() WHERE id={}'.format(job_id))
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
        for job in get_jobs(config['queue']['batch']):
            mark_as_started(job[0])
            job_args = build_args(job[0], job[5], job[6], job[7], job[8], job[9], job[10], job[11], job[12], job[13])
            run_pipeline(job_args)
            mark_as_completed(job[0])
            notify_user(job[4], job[0])
        sleep(config['queue']['sleep'])
