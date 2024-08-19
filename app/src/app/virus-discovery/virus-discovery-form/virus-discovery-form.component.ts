import {Component, OnInit} from '@angular/core';
import {HttpClient} from '@angular/common/http';
import {Router} from '@angular/router';
import {AlertController, LoadingController} from '@ionic/angular';
import {environment} from '../../../environments/environment'
import {VirusDiscoveryResponse} from '../../interfaces';


@Component({
  selector: 'app-virus-discovery-form',
  templateUrl: './virus-discovery-form.component.html',
  styleUrls: ['./virus-discovery-form.component.scss']
})
export class VirusDiscoveryFormComponent implements OnInit {

  email: string | null;
  sampleName: string | null;
  inputFormat: string | null;
  sequencingTechnology: string | null;
  singleFile: File | null;
  forwardFile: File | null;
  reverseFile: File | null;
  referenceGenome: File | null;

  constructor(private http: HttpClient, private router: Router,
              private alertController: AlertController, private loadingController: LoadingController) {
    this.initForm();
  }

  initForm() {
    this.email = null;
    this.sampleName = null;
    this.inputFormat = null;
    this.sequencingTechnology = null;
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
        formData.append('input_format', this.inputFormat);
        formData.append('sequencing_technology', this.sequencingTechnology);
        if(this.sequencingTechnology === 'single') {
          formData.append('single_file', this.singleFile, this.singleFile.name);
        }
        else if(this.sequencingTechnology === 'paired') {
          formData.append('forward_file', this.forwardFile, this.forwardFile.name);
          formData.append('reverse_file', this.reverseFile, this.reverseFile.name);
        }
        formData.append('reference_genome', this.referenceGenome, this.referenceGenome.name);
        this.http.post<VirusDiscoveryResponse>(environment.discvirAPI + '/job', formData, {responseType: 'json'}).subscribe(
          x => this.response(x),
          e => this.error(e.error),
          () => {
            this.initForm();
            this.loadingController.dismiss().then(null);
          }
        );
      }
    });
  }

  validateForm() {
    let validEmail = true;
    if(this.email === null) {
      validEmail = false;
    }
    else {
      // noinspection RegExpRedundantEscape
      const matches = this.email.match(
        /^(([^<>()[\]\\.,;:\s@\"]+(\.[^<>()[\]\\.,;:\s@\"]+)*)|(\".+\"))@((\[[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}\])|(([a-zA-Z\-0-9]+\.)+[a-zA-Z]{2,}))$/
      );
      if(!matches) {
        validEmail = false;
      }
    }
    if(!validEmail) {
      this.error({error: 'Please fill in a valid email address.'}).then(null);
      return false;
    }
    let validName = true;
    if(this.sampleName === null) {
      validName = false;
    }
    else {
      const matches = this.sampleName.match(/^\w+$/);
      if(!matches) {
        validName = false;
      }
    }
    if(!validName) {
      this.error({error: 'Please fill in a name for your experiment (only alphanumeric characters are permitted).'}).then(null);
      return false;
    }
    return this.validateInputFiles();
  }

  validateInputFiles() {
    if(this.inputFormat === null) {
      this.error({error: 'Please select format of the input files.'}).then(null);
      return false;
    }
    if(this.sequencingTechnology === null) {
      this.error({error: 'Please select the sequencing technology type.'}).then(null);
      return false;
    }
    else {
      if(this.sequencingTechnology === 'single') {
        if(this.singleFile === null) {
          this.error({error: 'Please select the input file to search in.'}).then(null);
          return false;
        }
      }
      else if(this.sequencingTechnology === 'paired') {
        if(this.forwardFile === null) {
          this.error({error: 'Please select the forward read input file to search in.'}).then(null);
          return false;
        }
        if(this.reverseFile === null) {
          this.error({error: 'Please select the reverse read input file to search in.'}).then(null);
          return false;
        }
      }
      else {
        this.error({error: 'Unknown sequencing technology type.'}).then(null);
        return false;
      }
    }
    if(this.referenceGenome === null) {
      this.error({error: 'Please select a file for the reference genome.'}).then(null);
      return false;
    }
    return true;
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
