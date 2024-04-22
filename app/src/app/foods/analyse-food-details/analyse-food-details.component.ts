import {Component, Input, OnInit} from '@angular/core';
import {ActivatedRoute, Router} from '@angular/router';
import {HttpClient} from '@angular/common/http';
import {AlertController, LoadingController} from '@ionic/angular';
import {environment} from '../../../environments/environment';
import {CompoundAnalysis, CompoundAnalysisResult, FoodAnalysisResult, FooDBContent, FooDBContentsCountResult, FooDBContentsResult, FooDBDetailed, FooDBDetailsResult, Taste} from '../../interfaces';
import {PaginationService} from '../../pagination.service';


@Component({
  selector: 'app-analyse-food-details',
  templateUrl: './analyse-food-details.component.html',
  styleUrls: ['./analyse-food-details.component.scss'],
  providers: [PaginationService]
})
export class AnalyseFoodDetailsComponent implements OnInit {

  @Input() foodId: number;

  foodDetails: FooDBDetailed;
  foodContents: Array<FooDBContent>;
  private foodTaste: Taste;
  private compoundAnalysis: Array<CompoundAnalysis>;
  private readonly pagesLimit: number;

  constructor(private route: ActivatedRoute, private router: Router, private http: HttpClient, public pagination: PaginationService,
              private alertController: AlertController, private loadingController: LoadingController) {
    this.foodId = 0;
    this.foodDetails = null;
    this.pagesLimit = 10;
    this.foodContents = [];
    this.compoundAnalysis = [];
  }

  ngOnInit() {
    this.pagination.initialize();
    this.getFoodDetails();
    this.getFoodCompoundsCount(true);
  }

  getFoodDetails() {
    this.http.get<FooDBDetailsResult>(environment.virtuousAPI + '/foods/' + this.foodId, {responseType: 'json'}).subscribe(
      x => {
        if(x.status === 'success') {
          this.foodDetails = x.result;
        }
        else {
          console.log(x.status)
        }
      },
      err => console.error(err)
    );
  }

  getFoodCompounds() {
    this.http.get<FooDBContentsResult>(environment.virtuousAPI + '/compounds',
      {params: {fid: this.foodId, page: this.pagination.newPage, limit: this.pagesLimit}, responseType: 'json'}).subscribe(
      x => {
        if(x.status === 'success') {
          this.foodContents = x.result;
        }
        else {
          console.log(x.status)
        }
      },
      err => console.error(err),
      () => this.pagination.update()
    );
  }

  getFoodCompoundsCount(both) {
    this.http.get<FooDBContentsCountResult>(environment.virtuousAPI + '/compounds/count',
      {params: {fid: this.foodId}, responseType: 'json'}).subscribe(
      x => {
        if(x.status === 'success') {
          this.pagination.count(x.result, this.pagesLimit);
          if(both) {
            this.getFoodCompounds();
          }
        }
        else {
          console.log(x.status);
        }
      },
      err => console.error(err),
      () => this.pagination.update()
    );
  }

  downloadCompounds() {
    return environment.virtuousAPI + '/compounds/download/' + this.foodId;
  }

  analyseFood() {
    this.loading().then(() => {
      const url = environment.virtuousAPI + '/analysis/multitaste/taste/food/' + this.foodId;
      this.http.get<FoodAnalysisResult>(url, {responseType: 'json'}).subscribe(x => {
          if(x.status === 'success') {
            this.foodTaste = x.result;
          }
          else {
            this.analysisError(x).then(null);
          }
        },
        e => this.analysisError(e.error),
        () => {
          this.loadingController.dismiss().then(() => {
            this.router.navigate(['taste'], {
              relativeTo: this.route,
              state: {taste: this.foodTaste}
            }).then(null);
          });
        }
      );
    });
  }

  analyseCompound(smiles, mode) {
    this.loading().then(() => {
      let url = environment.virtuousAPI;
      if(mode === 'umami') {
        url += '/analysis/umami/taste';
      }
      else if(mode === 'sweetbitter') {
        url += '/analysis/sweetbitter/taste';
      }
      else if(mode === 'multitaste') {
        url += '/analysis/multitaste/taste';
      }
      this.http.post<CompoundAnalysisResult>(url, {data: smiles}, {responseType: 'json'}).subscribe(x => {
          if(x.status === 'success') {
            this.compoundAnalysis = x.result;
          }
          else {
            this.analysisError(x).then(null);
          }
        },
        e => this.analysisError(e.error),
        () => {
          this.loadingController.dismiss().then(() => {
            this.router.navigate(['results'], {
              relativeTo: this.route,
              state: {
                view: mode,
                compounds: this.compoundAnalysis
              }
            }).then(null);
          });
        }
      );
    });
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

  toFixed(num) {
    return Number.parseFloat(num).toFixed(2);
  }

  goToPage(ev: Event, page: number | string) {
    ev.preventDefault();
    this.pagination.select(page);
    this.getFoodCompoundsCount(true);
  }

}
