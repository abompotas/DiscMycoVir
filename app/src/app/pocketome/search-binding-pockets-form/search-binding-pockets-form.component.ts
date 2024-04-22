import {Component, OnInit} from '@angular/core';
import {HttpClient} from '@angular/common/http';
import {Router} from '@angular/router';
import {AlertController, LoadingController} from '@ionic/angular';
import {environment} from '../../../environments/environment'
import {CompoundAnalysisResult} from '../../interfaces';


@Component({
  selector: 'app-search-binding-pockets-form',
  templateUrl: './search-binding-pockets-form.component.html',
  styleUrls: ['./search-binding-pockets-form.component.scss']
})
export class SearchBindingPocketsFormComponent implements OnInit {

  private pdbFile: File | null;
  private xtcFile: File | null;
  private email: string | null;
  private proteinChain: string | null;
  private ligandChain: string | null;
  private kFlag: number | null;
  private dist: number | null;
  private sasaThreshold: number | null;
  private dockThreshold: number | null;

  constructor(private http: HttpClient, private router: Router,
              private alertController: AlertController, private loadingController: LoadingController) {
    this.initForm();
  }

  ngOnInit() {
  }

  onPDBFileChange(event) {
    this.pdbFile = event.target.children['pdb-file'].files[0];
  }

  onXTCFileChange(event) {
    this.xtcFile = event.target.children['xtc-file'].files[0];
  }

  search() {
    this.loading().then(() => {
      if(this.pdbFile === null) {
        this.searchError({error: 'Please select a PDB file.'}).then(null);
      }
      else {
        const formData = new FormData();
        formData.append('pdb-file', this.pdbFile, this.pdbFile.name);
        if(this.xtcFile !== null) {
          formData.append('xtc-file', this.xtcFile, this.xtcFile.name);
        }
        formData.append('email', this.email);
        formData.append('protein-chain', this.proteinChain);
        formData.append('ligand-chain', this.ligandChain);
        formData.append('kflag', this.kFlag.toString());
        formData.append('dist', this.dist.toString());
        formData.append('sasa-threshold', this.sasaThreshold.toString());
        formData.append('dock-threshold', this.dockThreshold.toString());
        this.http.post<CompoundAnalysisResult>(environment.virtuousPocketome + '/pocketome',
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
    this.pdbFile = null;
    this.xtcFile = null;
    this.email = null;
    this.proteinChain = null;
    this.ligandChain = null;
    this.kFlag = 0;
    this.dist = 10;
    this.sasaThreshold = 0.75;
    this.dockThreshold = 0.1;
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
