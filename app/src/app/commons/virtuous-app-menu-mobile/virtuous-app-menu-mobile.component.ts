import {Component, Input, OnInit} from '@angular/core';

@Component({
  selector: 'app-virtuous-app-menu-mobile',
  templateUrl: './virtuous-app-menu-mobile.component.html',
  styleUrls: ['./virtuous-app-menu-mobile.component.scss'],
})
export class VirtuousAppMenuMobileComponent implements OnInit {

  @Input() active;

  constructor() {
  }

  ngOnInit() {
  }

  isActive(segment) {
    if(segment == this.active) {
      return 'active';
    }
    return '';
  }

}
