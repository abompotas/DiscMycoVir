import {Component, OnInit} from '@angular/core';
import {HttpClient} from '@angular/common/http';
import {Router} from '@angular/router';
import {AlertController, LoadingController} from '@ionic/angular';
import {environment} from '../../../environments/environment'
import {VirusDiscoveryResult} from '../../interfaces';


@Component({
  selector: 'app-virus-discovery-form',
  templateUrl: './virus-discovery-form.component.html',
  styleUrls: ['./virus-discovery-form.component.scss']
})
export class VirusDiscoveryFormComponent implements OnInit {

  email: string | null;
  sampleName: string | null;
  sequencingTechnology: string | null;
  adapter: string | null;
  slidingWindow: string | null;
  minLength: string | null;
  private singleFile: File | null;
  private forwardFile: File | null;
  private reverseFile: File | null;
  private referenceGenome: File | null;

  constructor(private http: HttpClient, private router: Router,
              private alertController: AlertController, private loadingController: LoadingController) {
    this.initForm();
  }

  initForm() {
    this.email = null;
    this.sampleName = null;
    this.sequencingTechnology = null;
    this.adapter = null;
    this.slidingWindow = null;
    this.minLength = null;
    this.singleFile = null;
    this.forwardFile = null;
    this.reverseFile = null;
    this.referenceGenome = null;
  }

  ngOnInit() {
  }

  onSingleFileChange(event) {
    this.singleFile = event.target.children['single_file'].files[0];
  }

  onForwardFileChange(event) {
    this.forwardFile = event.target.children['forward_file'].files[0];
  }

  onReverseFileChange(event) {
    this.reverseFile = event.target.children['reverse_file'].files[0];
  }

  onGenomeFileChange(event) {
    this.referenceGenome = event.target.children['reference_genome'].files[0];
  }

  search() {
    this.loading().then(() => {
      if(this.validateForm()) {
        const formData = new FormData();
        formData.append('email', this.email);
        formData.append('sample_name', this.sampleName);
        formData.append('sequencing_technology', this.sequencingTechnology);
        if(this.sequencingTechnology === 'single') {
          formData.append('single_file', this.singleFile, this.singleFile.name);
        }
        else if(this.sequencingTechnology === 'paired') {
          formData.append('forward_file', this.forwardFile, this.forwardFile.name);
          formData.append('reverse_file', this.reverseFile, this.reverseFile.name);
        }
        formData.append('reference_genome', this.referenceGenome, this.referenceGenome.name);
        this.http.post<VirusDiscoveryResult>(environment.discvirAPI + '/virus-discovery',
          formData, {responseType: 'json'}).subscribe(
          x => this.searchResponse(x),
          e => this.searchError(e.error),
          () => {
            this.initForm();
            this.loadingController.dismiss().then(null);
          }
        );
      }
    });
  }

  validateForm() {
    if(this.sequencingTechnology === null) {
      this.searchError({error: 'Please select the sequencing technology type.'}).then(null);
      return false;
    }
    else {
      if(this.sequencingTechnology === 'single') {
        if(this.singleFile === null) {
          this.searchError({error: 'Please select the input file to search in.'}).then(null);
          return false;
        }
      }
      else if(this.sequencingTechnology === 'paired') {
        if(this.forwardFile === null) {
          this.searchError({error: 'Please select the forward read input file to search in.'}).then(null);
          return false;
        }
        if(this.reverseFile === null) {
          this.searchError({error: 'Please select the reverse read input file to search in.'}).then(null);
          return false;
        }
      }
      else {
        this.searchError({error: 'Unknown sequencing technology type.'}).then(null);
        return false;
      }
    }
    if(this.referenceGenome === null) {
      this.searchError({error: 'Please select a file for the reference genome.'}).then(null);
      return false;
    }
    return true;
  }

  searchResponse(response) {
    if(response.status === 'success') {
      this.alertSuccess().then(null);
    }
    else {
      this.searchError(response).then(null);
    }
  }

  async loading() {
    const loading = await this.loadingController.create({
      message: 'Please wait...',
    });
    await loading.present();
  }

  async searchError(resp) {
    this.loadingController.dismiss().then(() => {
      let msg = 'Check your input for missing values.';
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
