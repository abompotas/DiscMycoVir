import {Component, OnInit} from '@angular/core';
import {ViewSelectorService} from '../../view-selector.service';

@Component({
  selector: 'app-analyse-sweetbitter-compound',
  templateUrl: './analyse-sweetbitter-compound.component.html',
  styleUrls: ['./analyse-sweetbitter-compound.component.scss'],
  providers: [ViewSelectorService]
})
export class AnalyseSweetBitterCompoundComponent implements OnInit {

  constructor(private selector: ViewSelectorService) {
  }

  ngOnInit() {
    this.selector.selectView('form');
  }

}
