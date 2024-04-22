import {ApplicationRef, Component, EventEmitter, OnInit, Output} from '@angular/core';
import {HttpClient} from '@angular/common/http';
import {LoadingController} from '@ionic/angular';
import {FooDBListResult, FooDBShort} from '../../interfaces';
import {environment} from '../../../environments/environment';


@Component({
  selector: 'app-organoleptic-trial-select',
  templateUrl: './organoleptic-trial-select.component.html',
  styleUrls: ['./organoleptic-trial-select.component.scss'],
})
export class OrganolepticTrialSelectComponent implements OnInit {

  @Output() selectedFood = new EventEmitter<number>();

  private foodId: number;
  private foods: Array<FooDBShort>

  constructor(private ref: ApplicationRef, private http: HttpClient, private loadingController: LoadingController) {
    this.foodId = 0;
    this.foods = [];
  }

  ngOnInit() {
    this.getFoods()
  }

  getFoods() {
    this.loading().then(() => {
      this.http.get<FooDBListResult>(environment.virtuousAPI + '/foods?ab=1', {responseType: 'json'}).subscribe(
        x => {
          if(x.status === 'success') {
            this.foods = x.result;
          }
        },
        err => console.error(err),
        () => {
          this.ref.tick();
          this.loadingController.dismiss().then(null);
        }
      );
    })
  }

  async loading() {
    const loading = await this.loadingController.create({
      message: 'Please wait...',
    });
    await loading.present();
  }

  changeSelection(ev) {
    this.foodId = ev.target.value;
  }

  emitSelection() {
    if(this.foodId) {
      this.selectedFood.emit(this.foodId);
    }
  }

}
