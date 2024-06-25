import {AfterViewInit, Component, Input, OnInit} from '@angular/core';
import {Chart} from 'chart.js/auto';
import {VirusDiscoveryResults} from '../../interfaces';


@Component({
  selector: 'app-virus-discovery-hits-graph',
  templateUrl: './virus-discovery-hits-graph.component.html',
  styleUrls: ['./virus-discovery-hits-graph.component.scss'],
})
export class VirusDiscoveryHitsGraphComponent implements OnInit, AfterViewInit {

  @Input() qid: number;
  @Input() qData: VirusDiscoveryResults;
  canvas: HTMLCanvasElement;
  rootStyle: CSSStyleDeclaration;
  graphData: any;
  datasets: any;

  constructor() {
    this.canvas = null;
    this.rootStyle = null;
    this.graphData = null;
    this.datasets = {'200': [], '80-200': [], '50-80': [], '40-50': [], '0-40': []};
  }

  ngOnInit() {
    this.rootStyle = getComputedStyle(document.body);
    this.createDatasets();
    this.graphData = {
      datasets: [{
        label: '>=200',
        data: this.datasets['200'],
        barThickness: 15,
        borderWidth: 2,
        borderColor: this.rootStyle.getPropertyValue('--ion-color-danger-shade'),
        backgroundColor: this.rootStyle.getPropertyValue('--ion-color-danger-tint'),
        stack: 'stack-0',
      }, {
        label: '80-200',
        data: this.datasets['80-200'],
        barThickness: 15,
        borderWidth: 2,
        borderColor: this.rootStyle.getPropertyValue('--ion-color-warning-shade'),
        backgroundColor: this.rootStyle.getPropertyValue('--ion-color-warning-tint'),
        stack: 'stack-0',
      }, {
        label: '50-80',
        data: this.datasets['50-80'],
        barThickness: 15,
        borderWidth: 2,
        borderColor: this.rootStyle.getPropertyValue('--ion-color-secondary-shade'),
        backgroundColor: this.rootStyle.getPropertyValue('--ion-color-secondary-tint'),
        stack: 'stack-0',
      }, {
        label: '40-50',
        data: this.datasets['40-50'],
        barThickness: 15,
        borderWidth: 2,
        borderColor: this.rootStyle.getPropertyValue('--ion-color-primary-shade'),
        backgroundColor: this.rootStyle.getPropertyValue('--ion-color-primary-tint'),
        stack: 'stack-0',
      }, {
        label: '<40',
        data: this.datasets['0-40'],
        barThickness: 15,
        borderWidth: 2,
        borderColor: this.rootStyle.getPropertyValue('--ion-color-dark-shade'),
        backgroundColor: this.rootStyle.getPropertyValue('--ion-color-dark'),
        stack: 'stack-0',
      }]
    };
  }

  ngAfterViewInit() {
    this.canvas = <HTMLCanvasElement>document.getElementById('hits-graph-' + this.qid);
    this.canvas.height = this.qData.alignments.length * 6 + 15;
    new Chart(this.canvas.getContext('2d'), {
      type: 'bar',
      data: this.graphData,
      options: {
        indexAxis: 'y',
        scales: {
          y: {display: false},
          x: {
            title: {
              text: 'Query',
              display: true,
              padding: {top: 10, bottom: 5},
              color: '#ffffff'
            },
            ticks: {
              color: '#ffffff',
              align: 'inner',
              padding: 0
            },
            position: 'top',
            backgroundColor: this.rootStyle.getPropertyValue('--ion-color-tertiary-tint'),
            max: this.qData.queryLetters
          }
        },
        plugins: {
          legend: {position: 'top'}
        }
      }
    });
  }

  createDatasets() {
    for(let a of this.qData.alignments) {
      for(let h of a.hsps) {
        const dataPoint = {
          y: a.hitId + ' Score: ' + h.score,
          x: h.alignLength
        }
        if(h.score >= 200) {
          this.datasets['200'].push(dataPoint);
        }
        else if(h.score < 200 && h.score >= 80) {
          this.datasets['80-200'].push(dataPoint);
        }
        else if(h.score < 80 && h.score >= 50) {
          this.datasets['50-80'].push(dataPoint);
        }
        else if(h.score < 50 && h.score >= 40) {
          this.datasets['40-50'].push(dataPoint);
        }
        else {
          this.datasets['0-40'].push(dataPoint);
        }
      }
    }
  }

}
