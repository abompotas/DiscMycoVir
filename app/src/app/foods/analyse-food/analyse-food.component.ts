import {Component, OnInit} from '@angular/core';
import {ActivatedRoute, Router} from '@angular/router';
import {ViewSelectorService} from '../../view-selector.service';
import {Taste} from '../../interfaces';

@Component({
  selector: 'app-analyse-food',
  templateUrl: './analyse-food.component.html',
  styleUrls: ['./analyse-food.component.scss'],
  providers: [ViewSelectorService]
})
export class AnalyseFoodComponent implements OnInit {

  private foodId: number;

  constructor(private route: ActivatedRoute, private router: Router, private selector: ViewSelectorService) {
    this.foodId = 0;
  }

  ngOnInit() {
    this.route.params.subscribe(params => {
      if(params.hasOwnProperty('food')) {
        this.foodId = params.food;
        this.selector.selectView('details');
      }
      else {
        this.showList();
      }
    })
  }

  showList() {
    this.foodId = 0;
    this.selector.view = 'list';
  }

  showDetails(fid) {
    this.router.navigate([fid], {relativeTo: this.route});
  }

}
