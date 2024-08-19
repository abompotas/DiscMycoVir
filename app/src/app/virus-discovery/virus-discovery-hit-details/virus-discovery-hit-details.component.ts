import {Component, Input, OnInit} from '@angular/core';


@Component({
  selector: 'app-virus-discovery-hit-details',
  templateUrl: './virus-discovery-hit-details.component.html',
  styleUrls: ['./virus-discovery-hit-details.component.scss'],
})
export class VirusDiscoveryHitDetailsComponent implements OnInit {

  @Input() hsp;

  constructor() {
  }

  ngOnInit() {
  }

}
