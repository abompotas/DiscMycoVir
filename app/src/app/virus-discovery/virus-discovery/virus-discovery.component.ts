import {Component, OnInit} from '@angular/core';
import {ActivatedRoute, Router} from '@angular/router';

@Component({
  selector: 'app-virus-discovery',
  templateUrl: './virus-discovery.component.html',
  styleUrls: ['./virus-discovery.component.scss']
})
export class VirusDiscoveryComponent implements OnInit {

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
