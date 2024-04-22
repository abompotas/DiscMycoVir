import {Component, OnInit} from '@angular/core';
import {ActivatedRoute, Router} from '@angular/router';

@Component({
  selector: 'app-search-binding-pockets',
  templateUrl: './search-binding-pockets.component.html',
  styleUrls: ['./search-binding-pockets.component.scss']
})
export class SearchBindingPocketsComponent implements OnInit {

  private jobId: number;
  private hash: string;

  constructor(private route: ActivatedRoute, private router: Router) {
    this.jobId = 0;
    this.hash = '';
    this.route.params.subscribe(params => {
      if(params.hasOwnProperty('job')) {
        if(params.hasOwnProperty('hash')) {
          this.jobId = params.job;
          this.hash = params.hash;
        }
        else {
          this.router.navigate(['..'], {relativeTo: this.route});
        }
      }
      else {
        this.jobId = 0;
        this.hash = '';
      }
    });
  }

  ngOnInit() {}

}
