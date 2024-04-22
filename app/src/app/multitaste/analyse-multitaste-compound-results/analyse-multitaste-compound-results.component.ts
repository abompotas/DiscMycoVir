import {Component, Input, OnInit} from '@angular/core';
import {CompoundAnalysis, Taste} from '../../interfaces';
import {environment} from '../../../environments/environment';

@Component({
  selector: 'app-analyse-multitaste-compound-results',
  templateUrl: './analyse-multitaste-compound-results.component.html',
  styleUrls: ['./analyse-multitaste-compound-results.component.scss']
})
export class AnalyseMultitasteCompoundResultsComponent implements OnInit {

  @Input() origin: String;
  @Input() compounds: Array<CompoundAnalysis>;

  private dominantTaste: Array<Array<string>>;

  constructor() {
    this.dominantTaste = [];
  }

  ngOnInit() {
    for(let c in this.compounds) {
      let max = 0;
      this.dominantTaste[c] = [];
      for(let t in this.compounds[c].taste) {
        if(max < this.compounds[c].taste[t]) {
          max = this.compounds[c].taste[t];
          this.dominantTaste[c] = [t];
        }
        else if(max === this.compounds[c].taste[t]) {
          this.dominantTaste[c].push(t);
        }
      }
    }
  }

  generateId(part1, part2) {
    const secs = Math.round(Date.now() / 1000)
    return part1 + '-' + part2 + '-' + secs;
  }

  goBack(ev) {
    ev.preventDefault();
    history.back();
  }

  formatProbabilityResults(prob) {
    return prob + '%';
  }

  downloadDescriptors(smiles, best) {
    if(best) {
      return environment.virtuousAPI + '/analysis/multitaste/descriptors/' + smiles + '/best';
    }
    return environment.virtuousAPI + '/analysis/multitaste/descriptors/' + smiles;
  }

  isDominantTaste(compound, taste) {
    if(this.dominantTaste[compound].includes(taste)) {
      return 'dominant-taste';
    }
    return '';
  }

}
