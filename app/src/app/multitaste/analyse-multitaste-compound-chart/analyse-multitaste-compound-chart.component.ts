import {AfterViewChecked, Component, Input, OnInit} from '@angular/core';
import {Chart, ChartConfiguration} from 'chart.js/auto';
import {Taste} from '../../interfaces';

@Component({
  selector: 'app-analyse-multitaste-compound-chart',
  templateUrl: './analyse-multitaste-compound-chart.component.html',
  styleUrls: ['./analyse-multitaste-compound-chart.component.scss'],
})
export class AnalyseMultitasteCompoundChartComponent implements OnInit, AfterViewChecked {

  @Input() canvasId: string;
  @Input() chartData: Taste;

  private initFlag: boolean;

  constructor() {
    this.initFlag = false;
  }

  ngOnInit() {
  }

  async ngAfterViewChecked() {
    if(!this.initFlag) {
      if(document.getElementById(this.canvasId)) {
        this.initFlag = true;
        try {
          await new Chart(this.canvasId, this.chartConfig());
        }
        catch(e) {
        }
      }
    }
  }

  chartConfig(): ChartConfiguration {
    return {
      type: 'radar',
      data: {
        labels: ['Bitter', 'Other', 'Sweet', 'Umami'],
        datasets: [{
          data: [this.chartData.bitter, this.chartData.other, this.chartData.sweet, this.chartData.umami],
          fill: true,
          backgroundColor: 'rgba(100, 58, 144, 0.2)',
          borderColor: 'rgb(100, 58, 144)',
          pointBackgroundColor: 'rgb(100, 58, 144)',
          pointBorderColor: '#ffffff',
          pointHoverBackgroundColor: '#ffffff',
          pointHoverBorderColor: 'rgb(100, 58, 144)'
        }]
      },
      options: {
        elements: {line: {borderWidth: 2}},
        scales: {r: {suggestedMin: 0, suggestedMax: 100}},
        plugins: {legend: {display: false}}
      }
    };
  }

}
