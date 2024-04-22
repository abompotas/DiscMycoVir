import {Component, Input, OnInit} from '@angular/core';

@Component({
  selector: 'app-virtuous-topbar',
  templateUrl: './virtuous-topbar.component.html',
  styleUrls: ['./virtuous-topbar.component.scss'],
})
export class VirtuousTopbarComponent implements OnInit {

  @Input() buttons: boolean = true;

  constructor() {
  }

  ngOnInit() {
  }

}
