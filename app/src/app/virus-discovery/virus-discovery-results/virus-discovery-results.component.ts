import {Component, OnInit} from '@angular/core';
import {HttpClient} from '@angular/common/http';
import {ActivatedRoute, Router} from '@angular/router';
import {AlertController, LoadingController} from '@ionic/angular';
import {environment} from '../../../environments/environment';
import {VirusDiscoveryResponse, VirusDiscoveryResults} from '../../interfaces';


// noinspection DuplicatedCode
@Component({
  selector: 'app-virus-discovery-results',
  templateUrl: './virus-discovery-results.component.html',
  styleUrls: ['./virus-discovery-results.component.scss'],
})
export class VirusDiscoveryResultsComponent implements OnInit {

  jobId: number;
  hash: string;
  resultsURL: string;
  downloadURL: string;
  results: Array<VirusDiscoveryResults>;

  constructor(private http: HttpClient, private route: ActivatedRoute, private router: Router,
              private alertController: AlertController, private loadingController: LoadingController) {
    this.resultsURL = '';
    this.downloadURL = '';
    this.jobId = 0;
    this.hash = '';
    this.results = [];
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
      this.resultsURL = environment.discvirAPI + '/results/' + this.jobId + '/' + this.hash;
      this.downloadURL = this.resultsURL + '/download';
      this.getResults();
    }
  }

  getResults() {
    this.loading().then(() => {
      this.http.get<VirusDiscoveryResponse>(this.resultsURL, {responseType: 'json'}).subscribe(
        x => this.parseResults(x),
        e => this.resultsError(e.error),
        () => this.loadingController.dismiss().then(null)
      );
    });
  }

  parseResults(resp) {
    for(const r of resp.results) {
      if(r.alignments.length > 0) {
        this.results.push(r);
      }
    }
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
