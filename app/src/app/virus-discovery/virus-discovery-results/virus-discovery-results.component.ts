import {Component, OnInit} from '@angular/core';
import {HttpClient} from '@angular/common/http';
import {ActivatedRoute, Router} from '@angular/router';
import {AlertController, LoadingController} from '@ionic/angular';
import {environment} from '../../../environments/environment';
import {VirusDiscoveryResponse} from '../../interfaces';


// noinspection DuplicatedCode
@Component({
  selector: 'app-virus-discovery-results',
  templateUrl: './virus-discovery-results.component.html',
  styleUrls: ['./virus-discovery-results.component.scss'],
})
export class VirusDiscoveryResultsComponent implements OnInit {

  jobId: number;
  hash: string;
  pocketomeURL: string;
  outputs: VirusDiscoveryResponse;
  resultsTable: Array<Array<string>>;

  constructor(private http: HttpClient, private route: ActivatedRoute, private router: Router,
              private alertController: AlertController, private loadingController: LoadingController) {
    this.pocketomeURL = '';
    this.outputs = null;
    this.resultsTable = [];
    this.jobId = 0;
    this.hash = '';
    this.route.params.subscribe(params => {
      if(params.hasOwnProperty('job')) {
        if(params.hasOwnProperty('hash')) {
          this.jobId = params.job;
          this.hash = params.hash;
        }
      }
    });
  }

  ngOnInit() {
    if((this.jobId) === 0 || (this.hash === '')) {
      this.router.navigate(['/']);
    }
    else {
      this.pocketomeURL = environment.discvirAPI + '/virus_discovery/' + this.jobId + '/' + this.hash;
      this.getResults();
    }
  }

  getResults() {
    this.loading().then(() => {
      this.http.get<VirusDiscoveryResponse>(this.pocketomeURL, {responseType: 'json'}).subscribe(
        x => this.parseResults(x),
        e => this.resultsError(e.error),
        () => this.loadingController.dismiss().then(null)
      );
    });
  }

  parseResults(resp) {
    this.resultsTable = [['PDBid', 'Matching results', 'Query residues', 'Hit residues', 'RMSD', 'SASA', 'Docking']];
    const colNames = ['PDBid', 'matchSize', 'query', 'hit', 'RMSD', 'SASA', 'Docking'];
    const cols = [];
    for(let r = 0; r < resp.results.results.length; r++) {
      let row = resp.results.results[r].replace('\n', '');
      row = row.split('\t');
      if(!r) {
        for(let c = 0; c < row.length; c++) {
          if(colNames.includes(row[c])) {
            cols.push(c);
          }
        }
      }
      else {
        const resultsRow = [];
        for(const c of cols) {
          resultsRow.push(row[c]);
        }
        this.resultsTable.push(resultsRow);
      }
    }
    this.outputs = resp.results;
  }

  async loading() {
    const loading = await this.loadingController.create({
      message: 'Please wait...',
    });
    await loading.present();
  }

  async resultsError(resp) {
    this.loadingController.dismiss().then(() => {
      let msg = 'Could not find any results matching your request.';
      if(resp.hasOwnProperty('error')) {
        msg = resp.error;
      }
      this.alertError(msg);
    });
  }

  async alertError(msg) {
    const alert = await this.alertController.create({
      header: 'Error!',
      message: msg,
      buttons: ['OK']
    });
    await alert.present();
  }

}
