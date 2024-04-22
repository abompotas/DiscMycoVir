import {Component, Input, OnInit} from '@angular/core';

@Component({
  selector: 'app-virtuous-footer',
  templateUrl: './virtuous-footer.component.html',
  styleUrls: ['./virtuous-footer.component.scss'],
})
export class VirtuousFooterComponent implements OnInit {

  @Input() spacer = false;

  constructor() {
  }

  ngOnInit() {
  }

}
