import {Component, OnInit} from '@angular/core';
import {ViewSelectorService} from '../../view-selector.service';

@Component({
  selector: 'app-analyse-umami-compound',
  templateUrl: './analyse-umami-compound.component.html',
  styleUrls: ['./analyse-umami-compound.component.scss'],
  providers: [ViewSelectorService]
})
export class AnalyseUmamiCompoundComponent implements OnInit {

  constructor(private selector: ViewSelectorService) {
  }

  ngOnInit() {
    this.selector.selectView('form');
  }

}
