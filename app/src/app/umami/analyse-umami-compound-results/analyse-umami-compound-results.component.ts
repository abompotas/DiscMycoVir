import {Component, Input, OnInit} from '@angular/core';
import {CompoundAnalysis} from '../../interfaces';
import {environment} from '../../../environments/environment';

@Component({
  selector: 'app-analyse-umami-compound-results',
  templateUrl: './analyse-umami-compound-results.component.html',
  styleUrls: ['./analyse-umami-compound-results.component.scss'],
})
export class AnalyseUmamiCompoundResultsComponent implements OnInit {

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
    let str = 'No';
    if(prob > 0) {
      str = 'Yes';
    }
    return str;
  }

  downloadDescriptors(smiles, best) {
    if(best) {
      return environment.virtuousAPI + '/analysis/umami/descriptors/' + smiles + '/best';
    }
    return environment.virtuousAPI + '/analysis/umami/descriptors/' + smiles;
  }

}
