import {Component, EventEmitter, OnInit, Output} from '@angular/core';
import {HttpClient} from '@angular/common/http';
import {Router} from '@angular/router';
import {AlertController, LoadingController} from '@ionic/angular';
import {environment} from '../../../environments/environment'
import {CompoundAnalysis, CompoundAnalysisResult} from '../../interfaces';


@Component({
  selector: 'app-analyse-umami-compound-form',
  templateUrl: './analyse-umami-compound-form.component.html',
  styleUrls: ['./analyse-umami-compound-form.component.scss']
})
export class AnalyseUmamiCompoundFormComponent implements OnInit {

  @Output() analysisResults = new EventEmitter<Array<CompoundAnalysis>>();

  defaultQuery: string;
  smiles: string;
  compoundsType: string | null;
  compoundsTypeFile: string | null;
  private smilesFile: File | null;
  private compounds: Array<CompoundAnalysis>;

  constructor(private http: HttpClient, private router: Router,
              private alertController: AlertController, private loadingController: LoadingController) {
    this.smiles = '';
    this.smilesFile = null;
    this.compoundsType = null;
    this.compoundsTypeFile = null;
    this.defaultQuery = 'N=C(N)NCCC[C@H](N)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCC(=O)O)C(=O)O\n' + 'GNPDIEHP'
  }

  ngOnInit() {
  }

  analyse() {
    this.loading().then(() => {
      if(this.smiles === '') {
        this.smiles = this.defaultQuery;
      }
      if(this.compoundsType === 'auto') {
        this.compoundsType = null;
      }
      this.http.post<CompoundAnalysisResult>(environment.virtuousAPI + '/analysis/umami/taste',
        {data: this.smiles, type: this.compoundsType}, {responseType: 'json'}).subscribe(
        x => this.analysisResponse(x),
        e => this.analysisError(e.error),
        () => this.loadingController.dismiss().then(() => this.analysisResults.emit(this.compounds))
      );
    });
  }

  onFileChange(event) {
    this.smilesFile = event.target.children['smiles-file'].files[0];
  }

  analyseFile() {
    this.loading().then(() => {
      if(this.smilesFile === null) {
        this.analysisError({error: 'Please select a file.'}).then(null);
      }
      else {
        if(this.compoundsTypeFile === 'auto') {
          this.compoundsTypeFile = null;
        }
        const formData = new FormData();
        formData.append('type', this.compoundsTypeFile);
        formData.append('smiles-file', this.smilesFile, this.smilesFile.name);
        this.http.post<CompoundAnalysisResult>(environment.virtuousAPI + '/analysis/umami/taste/file',
          formData, {responseType: 'json'}).subscribe(
          x => this.analysisResponse(x),
          e => this.analysisError(e.error),
          () => this.loadingController.dismiss().then(() => this.analysisResults.emit(this.compounds))
        );
      }
    });
  }

  analysisResponse(response) {
    if(response.status === 'success') {
      this.compounds = response.result;
      if(response.hasOwnProperty('error')) {
        this.alertMsg('Warning!', response.error).then(null);
      }
    }
    else {
      this.analysisError(response).then(null);
    }
  }

  async loading() {
    const loading = await this.loadingController.create({
      message: 'Please wait...',
    });
    await loading.present();
  }

  async analysisError(resp) {
    this.loadingController.dismiss().then(() => {
      let msg = 'Check your molecules! You may also need to specify their type. Note that allowed query types are SMILES, FASTA, Inchi, Sequence, Smarts or pubchem name.';
      if(resp.hasOwnProperty('error')) {
        msg = resp.error;
      }
      this.alertMsg('Error!', msg);
    });
  }

  async alertMsg(header, msg) {
    const alert = await this.alertController.create({
      header: header,
      message: msg,
      buttons: ['OK']
    });
    await alert.present();
  }

}
