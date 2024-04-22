import {Component, Input, OnInit} from '@angular/core';
import {CompoundAnalysis} from '../../interfaces';
import {environment} from '../../../environments/environment';

@Component({
  selector: 'app-analyse-sweetbitter-compound-results',
  templateUrl: './analyse-sweetbitter-compound-results.component.html',
  styleUrls: ['./analyse-sweetbitter-compound-results.component.scss'],
})
export class AnalyseSweetBitterCompoundResultsComponent implements OnInit {

  @Input() origin: String;
  @Input() compounds: Array<CompoundAnalysis>;

  constructor() {
    this.compounds = [];
  }

  ngOnInit() {
  }

  goBack(ev) {
    ev.preventDefault();
    history.back();
  }

  formatADResults(checkad) {
    let str = 'False';
    if(checkad) {
      str = 'True';
    }
    return str;
  }

  formatProbabilityResults(prob) {
    let str = 'Bitter';
    if(prob > 0) {
      str = 'Sweet';
    }
    return str;
  }

  downloadDescriptors(smiles, best) {
    if(best) {
      return environment.virtuousAPI + '/analysis/sweetbitter/descriptors/' + smiles + '/best';
    }
    return environment.virtuousAPI + '/analysis/sweetbitter/descriptors/' + smiles;
  }

}
