import {Component, Input, OnInit} from '@angular/core';
import {ActivatedRoute, Router} from '@angular/router';
import {HttpClient} from '@angular/common/http';
import {AlertController, LoadingController} from '@ionic/angular';
import {environment} from '../../../environments/environment';
import {FoodTrialResult, FoodTrialField, FooDBDetailsResult} from '../../interfaces';

@Component({
  selector: 'app-organoleptic-trial-taste',
  templateUrl: './organoleptic-trial-taste.component.html',
  styleUrls: ['./organoleptic-trial-taste.component.scss'],
})
export class OrganolepticTrialTasteComponent implements OnInit {

  @Input() foodId: number;
  @Input() origin: string;

  foodName: string;
  values: Array<number>;
  fields: Array<FoodTrialField>;

  constructor(private router: Router, private route: ActivatedRoute, private http: HttpClient,
              private loadingController: LoadingController, private alertController: AlertController) {
    this.foodName = '';
    this.values = [3, 3, 3, 3];
    this.fields = [
      {index: 0, name: 'sweetbitter', tasteLabel: 'Bitter/Sweet', minLabel: 'Bitter', maxLabel: 'Sweet'},
      {index: 1, name: 'umami', tasteLabel: 'Umami', minLabel: 'Not at all', maxLabel: 'Very'},
      {index: 2, name: 'salty', tasteLabel: 'Salty', minLabel: 'Not at all', maxLabel: 'Very'},
      {index: 3, name: 'sour', tasteLabel: 'Sour', minLabel: 'Not at all', maxLabel: 'Very'}
    ];
  }

  ngOnInit() {
    this.getFoodDetails();
  }

  getFoodDetails() {
    this.http.get<FooDBDetailsResult>(environment.virtuousAPI + '/foods/' + this.foodId, {responseType: 'json'}).subscribe(
      x => {
        if(x.status === 'success') {
          this.foodName = x.result.food_name;
        }
        else {
          console.log(x.status)
        }
      },
      err => console.error(err)
    );
  }

  submit(ev) {
    ev.preventDefault();
    this.loading().then(() => {
      const trialValues = {};
      for(const f of this.fields) {
        trialValues[f.name] = this.values[f.index];
      }
      this.http.post<FoodTrialResult>(environment.virtuousAPI + '/trials/' + this.foodId,
        trialValues, {responseType: 'json'}).subscribe(
        x => {
          if(x.status === 'success') {
            this.loadingController.dismiss().then(() => this.alertSuccess());
          }
          else {
            this.loadingController.dismiss().then(() => this.alertError(x.error));
          }
        },
        err => this.loadingController.dismiss().then(() => this.alertError(err.statusText))
      );
    });
  }

  async loading() {
    const loading = await this.loadingController.create({
      message: 'Please wait...',
    });
    await loading.present();
  }

  async alertError(error) {
    const alert = await this.alertController.create({
      header: 'Error!',
      message: error,
      buttons: ['OK']
    });
    await alert.present();
  }

  async alertSuccess() {
    const alert = await this.alertController.create({
      header: 'Thank you!',
      message: 'Your feedback is valuable to us.',
      buttons: [{
        text: 'OK',
        role: 'confirm',
        handler: () => this.router.navigate(['..'], {relativeTo: this.route})
      }]
    });
    await alert.present();
  }

}
