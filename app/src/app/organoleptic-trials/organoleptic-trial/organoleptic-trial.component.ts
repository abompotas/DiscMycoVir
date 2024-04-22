import {Component, Input, OnInit} from '@angular/core';
import {ActivatedRoute, Router} from '@angular/router';


@Component({
  selector: 'app-organoleptic-trial',
  templateUrl: './organoleptic-trial.component.html',
  styleUrls: ['./organoleptic-trial.component.scss'],
})
export class OrganolepticTrialComponent implements OnInit {

  @Input() foodId: number;

  constructor(private route: ActivatedRoute, private router: Router) {
    this.route.params.subscribe(params => {
      if(params.hasOwnProperty('food')) {
        this.foodId = params.food;
      }
      else {
        this.foodId = 0;
      }
    });
  }

  ngOnInit() {}

  selectFood(food) {
    this.router.navigate([food], {relativeTo: this.route});
  }

}
