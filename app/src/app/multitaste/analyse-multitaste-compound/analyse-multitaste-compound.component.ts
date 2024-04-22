import {Component, OnInit} from '@angular/core';
import {ViewSelectorService} from '../../view-selector.service';

@Component({
  selector: 'app-analyse-multitaste-compound',
  templateUrl: './analyse-multitaste-compound.component.html',
  styleUrls: ['./analyse-multitaste-compound.component.scss'],
  providers: [ViewSelectorService]
})
export class AnalyseMultitasteCompoundComponent implements OnInit {

  constructor(private selector: ViewSelectorService) {
  }

  ngOnInit() {
    this.selector.selectView('form');
  }

}
