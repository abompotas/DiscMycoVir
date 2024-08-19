import os

from flask import Flask
from flask_cors import CORS
from shared import config, db
from urllib.parse import quote_plus

from virus_discovery.routes import virus_discovery_blueprint


if not os.path.exists(config['app']['uploads_path']):
    os.makedirs(config['app']['uploads_path'])

app = Flask(__name__)
app.secret_key = 'PFtdbiPdkJevKACdS1eyVMsuKNCVRObt'
app.config['DEBUG'] = True
app.config['UPLOAD_FOLDER'] = config['app']['uploads_path']
app.config['MAX_CONTENT_LENGTH'] = config['app']['uploads_size'] * 1000 * 1000
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
app.config['SQLALCHEMY_DATABASE_URI'] = 'mysql+pymysql://{}:{}@{}/{}'.format(config['db']['username'],
                                                                             quote_plus(config['db']['password']),
                                                                             config['db']['host'],
                                                                             config['db']['schema'])

CORS(app)
db.init_app(app)

app.register_blueprint(virus_discovery_blueprint)


@app.route('/', methods=['GET'])
def home():
    return '<h1> DISCVIR: Virus Discovery API</h1>'


@app.after_request
def after_request(response):
    response.headers.add('Access-Control-Allow-Origin', '*')
    response.headers.add('Access-Control-Allow-Headers', 'Content-Type,Authorization')
    response.headers.add('Access-Control-Allow-Methods', 'GET,PUT,POST,DELETE')
    return response


if __name__ == '__main__':
    app.run(host=config['app']['host'], port=config['app']['port'])
