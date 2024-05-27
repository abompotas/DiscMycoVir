import smtplib
import ssl
import time
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from hashlib import sha256

import mysql.connector


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


class VirusDiscoveryJob:
    def __init__(self, config):
        self.config = config
        self.db = mysql.connector.connect(
            host=config['db']['host'],
            user=config['db']['username'],
            password=config['db']['password'],
            database=config['db']['schema']
        )

    def get_analysis_jobs(self, limit):
        cursor = self.db.cursor(dictionary=True)
        cursor.execute('''SELECT * FROM virus_discovery_jobs 
            WHERE started_analysis IS NULL ORDER BY id LIMIT 0,{}'''.format(limit))
        jobs = cursor.fetchall()
        cursor.close()
        self.db.commit()
        return jobs

    def mark_as_started_analysis(self, job_id):
        timestamp = time.strftime('%Y-%m-%d %H:%M:%S')
        cursor = self.db.cursor()
        cursor.execute('''UPDATE virus_discovery_jobs 
            SET started_analysis={} WHERE id={}'''.format(timestamp, job_id))
        cursor.close()
        self.db.commit()
        return timestamp

    def mark_as_completed_analysis(self, job_id):
        timestamp = time.strftime('%Y-%m-%d %H:%M:%S')
        cursor = self.db.cursor()
        cursor.execute('''UPDATE virus_discovery_jobs 
            SET completed_analysis={} WHERE id={}'''.format(timestamp, job_id))
        cursor.close()
        self.db.commit()
        return timestamp

    def get_trimming_jobs(self, limit):
        cursor = self.db.cursor(dictionary=True)
        cursor.execute('''SELECT * FROM virus_discovery_jobs 
            WHERE trimming_ready=0 AND started_trimming IS NULL ORDER BY id LIMIT 0,{}'''.format(limit))
        jobs = cursor.fetchall()
        cursor.close()
        self.db.commit()
        return jobs

    def mark_as_started_trimming(self, job_id):
        timestamp = time.strftime('%Y-%m-%d %H:%M:%S')
        cursor = self.db.cursor()
        cursor.execute('''UPDATE virus_discovery_jobs 
            SET started_trimming={} WHERE id={}'''.format(timestamp, job_id))
        cursor.close()
        self.db.commit()
        return timestamp

    def mark_as_completed_trimming(self, job_id):
        timestamp = time.strftime('%Y-%m-%d %H:%M:%S')
        cursor = self.db.cursor()
        cursor.execute('''UPDATE `virus_discovery_jobs` 
            SET completed_trimming={} WHERE id={}'''.format(timestamp, job_id))
        cursor.close()
        self.db.commit()
        return timestamp

    def get_discovery_jobs(self, limit):
        cursor = self.db.cursor(dictionary=True)
        cursor.execute('''SELECT * FROM virus_discovery_jobs 
            WHERE trimming_ready=1 AND started_discovery IS NULL ORDER BY id LIMIT 0,{}'''.format(limit))
        jobs = cursor.fetchall()
        cursor.close()
        self.db.commit()
        return jobs

    def mark_as_started_discovery(self, job_id):
        timestamp = time.strftime('%Y-%m-%d %H:%M:%S')
        cursor = self.db.cursor()
        cursor.execute('''UPDATE virus_discovery_jobs 
            SET started_discovery={} WHERE id={}'''.format(timestamp, job_id))
        cursor.close()
        self.db.commit()
        return timestamp

    def mark_as_completed_discovery(self, job_id):
        timestamp = time.strftime('%Y-%m-%d %H:%M:%S')
        cursor = self.db.cursor()
        cursor.execute('''UPDATE `virus_discovery_jobs` 
            SET completed_discovery={} WHERE id={}'''.format(timestamp, job_id))
        cursor.close()
        self.db.commit()
        return timestamp
