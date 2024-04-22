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

  private genome: File | null;
  private email: string | null;

  constructor(private http: HttpClient, private router: Router,
              private alertController: AlertController, private loadingController: LoadingController) {
    this.initForm();
  }

  ngOnInit() {
  }

  onGenomeFileChange(event) {
    this.genome = event.target.children['reference-genome'].files[0];
  }

  search() {
    this.loading().then(() => {
      if(this.genome === null) {
        this.searchError({error: 'Please select a file for the reference genome.'}).then(null);
      }
      else {
        const formData = new FormData();
        formData.append('pdb-file', this.genome, this.genome.name);
        formData.append('email', this.email);
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

  initForm() {
    this.genome = null;
    this.email = null;
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
