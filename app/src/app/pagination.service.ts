import {ApplicationRef, Injectable} from '@angular/core';
import {PaginationItem} from './interfaces';


@Injectable({
  providedIn: 'root'
})
export class PaginationService {

  public page: number;
  public newPage: number;
  public pagesCount: number;
  public pagesList: Array<PaginationItem>;

  constructor(private ref: ApplicationRef) {
    this.initialize();
  }

  initialize() {
    this.page = 0;
    this.newPage = 0;
    this.pagesCount = 0;
    this.pagesList = [];
  }

  count(items, limit) {
    this.pagesCount = Math.ceil(items / limit);
  }

  select(page) {
    if(isNaN(page)) {
      if(page === 'prev') {
        if(this.page > 0) {
          this.newPage--;
        }
      }
      else if(page === 'next') {
        if(this.page < (this.pagesCount - 1)) {
          this.newPage++;
        }
      }
      else if(page === 'first') {
        this.newPage = 0;
      }
      else if(page === 'last') {
        this.newPage = this.pagesCount - 1;
      }
    }
    else {
      if(page > -1) {
        this.newPage = page;
      }
    }
  }

  update() {
    this.pagesList = [];
    if(this.pagesCount) {
      let skipped = false;
      for(let i = 0; i < this.pagesCount; i++) {
        if(i < 2) {
          this.pagesList.push({title: (i + 1).toString(), page: i, class: 'ion-hide-md-down'});
        }
        else if((i > (this.newPage - 3)) && (i < (this.newPage + 3))) {
          skipped = false;
          if(i === this.newPage) {
            this.pagesList.push({title: (i + 1).toString(), page: i});
          }
          else {
            this.pagesList.push({title: (i + 1).toString(), page: i, class: 'ion-hide-md-down'});
          }
        }
        else if(i > (this.pagesCount - 3)) {
          this.pagesList.push({title: (i + 1).toString(), page: i, class: 'ion-hide-md-down'});
        }
        else {
          if(!skipped) {
            skipped = true;
            this.pagesList.push({title: '...', page: -1, class: 'ion-hide-md-down'});
          }
        }
      }
      this.page = this.newPage;
      this.ref.tick();
    }
  }
}
