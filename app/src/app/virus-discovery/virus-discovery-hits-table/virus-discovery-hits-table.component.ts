import {AfterViewInit, Component, Input, OnInit} from '@angular/core';
import {VirusDiscoveryResultHSP, VirusDiscoveryResults} from '../../interfaces';
import DataTable from 'datatables.net-dt';


@Component({
  selector: 'app-virus-discovery-hits-table',
  templateUrl: './virus-discovery-hits-table.component.html',
  styleUrls: ['./virus-discovery-hits-table.component.scss'],
})
export class VirusDiscoveryHitsTableComponent implements OnInit, AfterViewInit {

  @Input() qid: number;
  @Input() qData: VirusDiscoveryResults;
  alignmentsData: Array<VirusDiscoveryResultHSP>;
  matchesData: Array<any>;
  selectedHSP: any;

  constructor() {
    this.alignmentsData = [];
    this.matchesData = [];
    this.selectedHSP = null
  }

  ngOnInit() {
    for(let a of this.qData.alignments) {
      for(let h of a.hsps) {
        h.title = a.hitId + a.hitDef;
        h.length = a.length;
        this.alignmentsData.push(h);
        this.matchesData.push(this.parseHSPMatch(h, 50));
      }
    }
  }

  ngAfterViewInit() {
    new DataTable('#alignments-' + this.qid);
  }

  parseHSPMatch(hsp, limit) {
    const parsedData = {
      'queryStart': [], 'queryEnd': [], 'query': [],
      'sbjctStart': [], 'sbjctEnd': [], 'sbjct': [],
      'match': [], 'title': hsp.title
    }
    let step = 1;
    if(hsp.strand[1] === 'Minus') {
      step = -1;
    }
    let query = '';
    let sbjct = '';
    let match = '';
    let queryStart = hsp.queryStart - 1;
    let queryEnd = queryStart;
    let sbjctStart = hsp.sbjctStart - step;
    let sbjctEnd = sbjctStart;
    for(let i = 0; i < hsp.query.length; i++) {
      query += hsp.query[i];
      sbjct += hsp.sbjct[i];
      match += hsp.match[i];
      if(hsp.query[i] !== '-') {
        queryEnd++;
      }
      if(hsp.sbjct[i] !== '-') {
        sbjctEnd += step;
      }
      if((i % limit) === 0) {
        if(hsp.query[i] !== '-') {
          queryStart++;
        }
        if(hsp.sbjct[i] !== '-') {
          sbjctStart += step;
        }
      }
      if((i % limit) === (limit - 1)) {
        parsedData.queryStart.push(queryStart);
        parsedData.queryEnd.push(queryEnd);
        parsedData.query.push(query);
        parsedData.sbjctStart.push(sbjctStart);
        parsedData.sbjctEnd.push(sbjctEnd);
        parsedData.sbjct.push(sbjct);
        parsedData.match.push(match);
        query = '';
        sbjct = '';
        match = '';
        queryStart = queryEnd;
        sbjctStart = sbjctEnd;
      }
    }
    if(query !== '') {
      parsedData.queryStart.push(queryStart);
      parsedData.queryEnd.push(hsp.queryEnd);
      parsedData.query.push(query);
      parsedData.sbjctStart.push(sbjctStart);
      parsedData.sbjctEnd.push(hsp.sbjctEnd);
      parsedData.sbjct.push(sbjct);
      parsedData.match.push(match);
    }
    return parsedData
  }

  showHSPDetails(i: number) {
    this.selectedHSP = this.matchesData[i];
    return false;
  }

}
