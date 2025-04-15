import {AfterViewInit, Component, Input, OnChanges, OnInit, SimpleChanges} from '@angular/core';
import {Chart} from 'chart.js/auto';
import {VirusDiscoveryResultHSP} from '../../interfaces';


@Component({
  selector: 'app-virus-discovery-hit-coverage',
  templateUrl: './virus-discovery-hit-coverage.component.html',
  styleUrls: ['./virus-discovery-hit-coverage.component.scss'],
})
export class VirusDiscoveryHitCoverageComponent implements OnInit, OnChanges, AfterViewInit {

  @Input() qid: number;
  @Input() hsp: VirusDiscoveryResultHSP;
  coverageGraph: any;
  canvas: HTMLCanvasElement;
  rootStyle: CSSStyleDeclaration;
  graphData: any;
  datasets: any;

  constructor() {
    this.coverageGraph = null;
    this.canvas = null;
    this.rootStyle = null;
    this.graphData = null;
    this.datasets = {'200': [], '80-200': [], '50-80': [], '40-50': [], '0-40': []};
  }

  ngOnInit() {
    this.rootStyle = getComputedStyle(document.body);
    this.updateData();
  }

  ngOnChanges(changes: SimpleChanges) {
    if(this.rootStyle !== null) {
      this.updateData();
      this.drawGraph();
    }
  }

  ngAfterViewInit() {
    this.drawGraph();
  }

  updateData() {
    this.createDatasets();
    this.graphData = {
      datasets: [{
        label: 'Score >=200',
        data: this.datasets['200'],
        barThickness: 15,
        borderWidth: 2,
        borderColor: this.rootStyle.getPropertyValue('--ion-color-danger-shade'),
        backgroundColor: this.rootStyle.getPropertyValue('--ion-color-danger-tint'),
        stack: 'stack-1'
      }, {
        label: '80<= Score <200',
        data: this.datasets['80-200'],
        barThickness: 15,
        borderWidth: 2,
        borderColor: this.rootStyle.getPropertyValue('--ion-color-warning-shade'),
        backgroundColor: this.rootStyle.getPropertyValue('--ion-color-warning-tint'),
        stack: 'stack-1'
      }, {
        label: '50<= Score <80',
        data: this.datasets['50-80'],
        barThickness: 15,
        borderWidth: 2,
        borderColor: this.rootStyle.getPropertyValue('--ion-color-secondary-shade'),
        backgroundColor: this.rootStyle.getPropertyValue('--ion-color-secondary-tint'),
        stack: 'stack-1'
      }, {
        label: '40<= Score <50',
        data: this.datasets['40-50'],
        barThickness: 15,
        borderWidth: 2,
        borderColor: this.rootStyle.getPropertyValue('--ion-color-primary-shade'),
        backgroundColor: this.rootStyle.getPropertyValue('--ion-color-primary-tint'),
        stack: 'stack-1'
      }, {
        label: 'Score <40',
        data: this.datasets['0-40'],
        barThickness: 15,
        borderWidth: 2,
        borderColor: this.rootStyle.getPropertyValue('--ion-color-dark-shade'),
        backgroundColor: this.rootStyle.getPropertyValue('--ion-color-dark'),
        stack: 'stack-1'
      }]
    };
  }

  drawGraph() {
    if(this.coverageGraph !== null) {
      this.coverageGraph.destroy();
    }
    this.canvas = <HTMLCanvasElement>document.getElementById('hit-coverage-' + this.qid);
    this.canvas.height = 21;
    this.coverageGraph = new Chart(this.canvas.getContext('2d'), {
      type: 'bar',
      data: this.graphData,
      options: {
        indexAxis: 'y',
        scales: {
          y: {display: false},
          x: {
            title: {
              text: 'Sequence',
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
            backgroundColor: this.rootStyle.getPropertyValue('--ion-color-medium'),
            max: this.hsp.length
          }
        },
        plugins: {
          legend: {position: 'top'}
        }
      }
    });
  }

  createDatasets() {
    const per = Math.round(100 * this.hsp.alignLength / this.hsp.length).toFixed(2);
    const dataPoint = {
      y: 'Score: ' + this.hsp.score + ', Alignment length: ' + this.hsp.alignLength + ' (' + per + '%), E-value: ' + this.hsp.expect,
      x: [this.hsp.sbjctStart, this.hsp.sbjctEnd]
    }
    this.datasets = {'200': [], '80-200': [], '50-80': [], '40-50': [], '0-40': []};
    if(this.hsp.score >= 200) {
      this.datasets['200'].push(dataPoint);
    }
    else if(this.hsp.score < 200 && this.hsp.score >= 80) {
      this.datasets['80-200'].push(dataPoint);
    }
    else if(this.hsp.score < 80 && this.hsp.score >= 50) {
      this.datasets['50-80'].push(dataPoint);
    }
    else if(this.hsp.score < 50 && this.hsp.score >= 40) {
      this.datasets['40-50'].push(dataPoint);
    }
    else {
      this.datasets['0-40'].push(dataPoint);
    }
  }

}
