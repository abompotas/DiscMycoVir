import {Component, OnInit, Pipe, PipeTransform} from '@angular/core';
import {HttpClient} from '@angular/common/http';
import {ActivatedRoute, Router} from '@angular/router';
import {DomSanitizer} from '@angular/platform-browser';
import {AlertController, LoadingController} from '@ionic/angular';
import {VirusDiscoveryResponse} from '../../interfaces';
import {environment} from '../../../environments/environment';

@Pipe({name: 'safeIFrame'})
export class SafeHtmlPipe implements PipeTransform {
  constructor(private sanitizer: DomSanitizer) {
  }

  transform(inHtml) {
    let outHtml = inHtml.replaceAll(/ href="#M\d+/g, '');
    outHtml = outHtml.replaceAll('href=', 'target="_blank" href=');
    return this.sanitizer.bypassSecurityTrustHtml(outHtml);
  }
}

// noinspection DuplicatedCode
@Component({
  selector: 'app-virus-discovery-trimming',
  templateUrl: './virus-discovery-trimming.component.html',
  styleUrls: ['./virus-discovery-trimming.component.scss'],
})
export class VirusDiscoveryTrimmingComponent implements OnInit {

  jobId: number;
  hash: string;
  analysisReports: Array<string>;
  sequencingTechnology: string | null;
  slidingWindow: string | null;
  minLength: string | null;
  adapter: File | null;

  constructor(private http: HttpClient, private route: ActivatedRoute, private router: Router,
              private alertController: AlertController, private loadingController: LoadingController) {
    this.jobId = 0;
    this.hash = '';
    this.analysisReports = [];
    this.route.params.subscribe(params => {
      if(params.hasOwnProperty('job')) {
        if(params.hasOwnProperty('hash')) {
          this.jobId = params.job;
          this.hash = params.hash;
        }
      }
    });
    this.initForm();
  }

  ngOnInit() {
    if((this.jobId) === 0 || (this.hash === '')) {
      this.router.navigate(['/']).then(null);
    }
    else {
      this.fetchAnalysis();
    }
  }

  initForm() {
    this.sequencingTechnology = null;
    this.slidingWindow = null;
    this.minLength = null;
    this.adapter = null;
  }

  fetchAnalysis() {
    this.loading().then(() => {
      this.http.get<VirusDiscoveryResponse>(environment.discvirAPI + '/analysis/' + this.jobId + '/' + this.hash, {responseType: 'json'}).subscribe(
        x => {
          this.analysisReports = [...x.reports];
        },
        e => this.error(e.error),
        () => this.loadingController.dismiss().then(null)
      );
    });
  }

  onAdapterChange(event) {
    this.adapter = event.target.children['adapter'].files[0];
  }

  trim() {
    this.loading().then(() => {
      const formData = new FormData();
      if(this.adapter !== null) {
        formData.append('adapter', this.adapter, this.adapter.name);
      }
      if(this.slidingWindow !== null) {
        formData.append('sliding_window', this.slidingWindow);
      }
      if(this.minLength !== null) {
        formData.append('min_length', this.minLength);
      }
      this.http.put<VirusDiscoveryResponse>(environment.discvirAPI + '/trimming/' + this.jobId + '/' + this.hash, formData, {responseType: 'json'}).subscribe(
        x => this.response(x),
        e => this.error(e.error),
        () => {
          this.initForm();
          this.loadingController.dismiss().then(null);
        }
      );
    });
  }

  proceed() {
    this.loading().then(() => {
      this.http.put<VirusDiscoveryResponse>(environment.discvirAPI + '/discovery/' + this.jobId + '/' + this.hash, {}, {responseType: 'json'}).subscribe(
        x => this.response(x),
        e => this.error(e.error),
        () => {
          this.initForm();
          this.loadingController.dismiss().then(null);
        }
      );
    });
  }

  response(resp) {
    if(resp.status === 'success') {
      this.alertSuccess().then(null);
    }
    else {
      this.error(resp).then(null);
    }
  }

  async loading() {
    const loading = await this.loadingController.create({
      message: 'Please wait...',
    });
    await loading.present();
  }

  async error(resp) {
    this.loadingController.dismiss().then(() => {
      let msg = 'Oops! Something went wrong...';
      if(resp.hasOwnProperty('error')) {
        msg = resp.error;
      }
      this.alertError(msg);
    });
  }

  async alertSuccess() {
    const alert = await this.alertController.create({
      header: 'Success!',
      message: 'Your query has been submitted. Once the search is completed you will receive an email containing the results.',
      buttons: ['OK']
    });
    await alert.present();
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
