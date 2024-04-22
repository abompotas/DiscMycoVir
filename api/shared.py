from json import load
from flask_sqlalchemy import SQLAlchemy

with open('config.json') as fp:
    config = load(fp)

db = SQLAlchemy()
