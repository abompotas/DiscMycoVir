import {Injectable} from '@angular/core';
import {ActivatedRoute, NavigationEnd, Router} from '@angular/router';
import {filter} from 'rxjs/operators';
import {CompoundAnalysis, Taste} from './interfaces';

@Injectable({
  providedIn: 'root'
})
export class ViewSelectorService {

  public view: string;
  public taste: Taste;
  public compounds: Array<CompoundAnalysis>;

  constructor(private route: ActivatedRoute, private router: Router) {
    this.view = '';
    this.compounds = [];
  }

  selectView(defaultView) {
    this.router.events.pipe(filter(e => e instanceof NavigationEnd)).subscribe(() => {
      const segment = this.router.url.split('/').pop();
      if(segment === 'taste') {
        this.selectTasteOrResults('taste');
      }
      else if(segment === 'results') {
        this.selectTasteOrResults('results');
      }
      else if(segment === 'organoleptic-trial') {
        this.taste = null;
        this.compounds = [];
        this.view = 'organoleptic-trial';
      }
      else {
        this.taste = null;
        this.compounds = [];
        this.view = defaultView;
      }
    });
  }

  selectTasteOrResults(mode) {
    let flag = false;
    const nav = this.router.getCurrentNavigation();
    if(nav) {
      if(nav.extras) {
        if(nav.extras.hasOwnProperty('state')) {
          if(mode === 'taste') {
            if(nav.extras.state.hasOwnProperty('taste')) {
              this.taste = nav.extras.state.taste;
            }
            this.view = 'taste';
          }
          else {
            if(nav.extras.state.hasOwnProperty('compounds')) {
              this.compounds = nav.extras.state.compounds;
            }
            if(nav.extras.state.hasOwnProperty('view')) {
              this.view = nav.extras.state.view;
            }
            else {
              this.view = 'results';
            }
          }
          flag = true;
        }
      }
    }
    if(!flag) {
      this.router.navigate(['..'], {relativeTo: this.route});
    }
  }

  showResults(compounds: Array<CompoundAnalysis>, view: string) {
    this.router.navigate(['results'], {
      relativeTo: this.route,
      state: {
        compounds: compounds,
        view: view
      }
    });
  }

}
