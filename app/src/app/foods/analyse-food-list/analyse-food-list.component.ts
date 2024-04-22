import {Component, EventEmitter, OnInit, Output} from '@angular/core';
import {FormBuilder} from '@angular/forms';
import {HttpClient} from '@angular/common/http';
import {environment} from '../../../environments/environment';
import {FooDBCountResult, FooDBListResult, FooDBShort} from '../../interfaces';
import {PaginationService} from '../../pagination.service';


@Component({
  selector: 'app-analyse-food-list',
  templateUrl: './analyse-food-list.component.html',
  styleUrls: ['./analyse-food-list.component.scss'],
  providers: [PaginationService]
})
export class AnalyseFoodListComponent implements OnInit {

  @Output() selectFood = new EventEmitter<number>();

  searchForm = this.formBuilder.group({search: ''});
  private foods: Array<FooDBShort>;
  private search: boolean;
  private readonly pagesLimit: number;

  constructor(private formBuilder: FormBuilder, private http: HttpClient, private pagination: PaginationService) {
    this.search = false;
    this.pagesLimit = 20;
    this.foods = [];
  }

  ngOnInit() {
    this.foodsInit();
  }

  foodsInit() {
    this.search = false;
    this.pagination.initialize();
    this.getFoodsCount(true);
  }

  getFoods() {
    this.http.get<FooDBListResult>(environment.virtuousAPI + '/foods',
      {params: {page: this.pagination.newPage, limit: this.pagesLimit}, responseType: 'json'}).subscribe(
      x => {
        if(x.status === 'success') {
          this.foods = x.result;
        }
        else {
          console.log(x.status);
        }
      },
      err => console.error(err),
      () => this.pagination.update()
    );
  }

  getFoodsCount(both) {
    this.http.get<FooDBCountResult>(environment.virtuousAPI + '/foods/count', {responseType: 'json'}).subscribe(
      x => {
        if(x.status === 'success') {
          this.pagination.count(x.result, this.pagesLimit);
          if(both) {
            this.getFoods();
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

  searchFoods() {
    this.search = true;
    this.http.get<FooDBListResult>(environment.virtuousAPI + '/foods/search/' + this.searchForm.value.search,
      {params: {page: this.pagination.newPage, limit: this.pagesLimit}, responseType: 'json'}).subscribe(
      x => {
        if(x.status === 'success') {
          this.foods = x.result;
        }
        else {
          console.log(x.status);
        }
      },
      err => console.error(err),
      () => this.pagination.update()
    );
  }

  searchFoodsCount(both) {
    this.http.get<FooDBCountResult>(environment.virtuousAPI + '/foods/search/' + this.searchForm.value.search + '/count', {responseType: 'json'}).subscribe(
      x => {
        if(x.status === 'success') {
          this.pagination.count(x.result, this.pagesLimit);
          if(both) {
            this.searchFoods();
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

  searchFormSubmit() {
    this.pagination.initialize();
    this.searchFoodsCount(true);
  }

  goToPage(ev: Event, page: number | string) {
    ev.preventDefault();
    this.pagination.select(page);
    if(this.search) {
      this.searchFoodsCount(true);
    }
    else {
      this.getFoodsCount(true);
    }
  }

  showFoodDetails(ev: Event, fid: number) {
    ev.preventDefault();
    this.selectFood.emit(fid);
  }

}
