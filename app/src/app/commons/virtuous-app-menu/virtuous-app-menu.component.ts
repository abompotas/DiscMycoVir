import {Component, Input, OnInit} from '@angular/core';

@Component({
  selector: 'app-virtuous-app-menu',
  templateUrl: './virtuous-app-menu.component.html',
  styleUrls: ['./virtuous-app-menu.component.scss'],
})
export class VirtuousAppMenuComponent implements OnInit {

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
