import time

import mysql.connector


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
            WHERE trimming_ready=1 AND started_trimming IS NULL ORDER BY id LIMIT 0,{}'''.format(limit))
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
            WHERE trimming_ready=2 AND started_discovery IS NULL ORDER BY id LIMIT 0,{}'''.format(limit))
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
