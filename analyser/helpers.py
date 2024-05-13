import smtplib
import ssl
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from hashlib import sha256
import mysql.connector


def notify_user(config, email, job_id):
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
        message['Subject'] = 'Virtuous Pocketome - Search completed'
        results_hash = sha256('{}%%{}'.format(job_id, email).encode('utf-8')).hexdigest()
        results_url = '{}/{}/{}'.format(config['webbapp'], job_id, results_hash)
        results_text = '<html><body><p>VirtuousPocketome has analyzed your request</p>' \
                       '<p>View the results: <a href="{0}">{0}</a></p></body></html>'.format(results_url)
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

    def get_jobs(self, limit):
        cursor = self.db.cursor(dictionary=True)
        cursor.execute("SELECT * FROM virus_discovery_jobs WHERE started IS NULL ORDER BY id ASC LIMIT 0,{}".format(limit))
        jobs = cursor.fetchall()
        cursor.close()
        self.db.commit()
        return jobs

    def mark_as_started(self, job_id):
        cursor = self.db.cursor()
        cursor.execute('UPDATE `virus_discovery_jobs` SET started=CURRENT_TIMESTAMP() WHERE id={}'.format(job_id))
        cursor.close()
        self.db.commit()

    def mark_as_completed(self, job_id):
        cursor = self.db.cursor()
        cursor.execute('UPDATE `virus_discovery_jobs` SET completed=CURRENT_TIMESTAMP() WHERE id={}'.format(job_id))
        cursor.close()
        self.db.commit()
