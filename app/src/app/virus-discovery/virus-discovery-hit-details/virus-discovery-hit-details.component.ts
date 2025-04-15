import {Component, Input, OnChanges, OnInit, SimpleChanges} from '@angular/core';


@Component({
  selector: 'app-virus-discovery-hit-details',
  templateUrl: './virus-discovery-hit-details.component.html',
  styleUrls: ['./virus-discovery-hit-details.component.scss'],
})
export class VirusDiscoveryHitDetailsComponent implements OnInit, OnChanges {

  @Input() hsp;

  constructor() {
  }

  ngOnInit() {
  }

  ngOnChanges(changes: SimpleChanges) {
  }

}
