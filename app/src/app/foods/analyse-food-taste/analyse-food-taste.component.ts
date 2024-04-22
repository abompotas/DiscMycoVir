import {AfterViewChecked, Component, Input, OnInit} from '@angular/core';
import {HttpClient} from '@angular/common/http';
import {FooDBDetailsResult, Taste} from '../../interfaces';
import {Chart, ChartConfiguration} from 'chart.js/auto';
import {environment} from '../../../environments/environment';

@Component({
  selector: 'app-analyse-food-taste',
  templateUrl: './analyse-food-taste.component.html',
  styleUrls: ['./analyse-food-taste.component.scss']
})
export class AnalyseFoodTasteComponent implements OnInit, AfterViewChecked {

  @Input() foodId: number;
  @Input() taste: Taste;

  private initFlag: boolean;
  private canvasId: string;
  private foodName: string;
  private dominantTaste: Array<string>;

  constructor(private http: HttpClient) {
    this.initFlag = false;
    this.canvasId = '';
    this.foodName = '';
    this.dominantTaste = [];
  }

  ngOnInit() {
    this.canvasId = this.generateId('bar', 0);
    this.getFoodDetails();
    this.calculateDominantTaste();
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
      type: 'bar',
      data: {
        labels: ['Bitter', 'Sweet', 'Other', 'Umami'],
        datasets: [{
          data: [this.taste.bitter, this.taste.sweet, this.taste.other, this.taste.umami]
        }]
      },
      options: {
        elements: {line: {borderWidth: 2}},
        scales: {y: {beginAtZero: true}},
        plugins: {legend: {display: false}}
      }
    };
  }

  calculateDominantTaste() {
    let max = 0;
    this.dominantTaste = [];
    for(let t in this.taste) {
      if(max < this.taste[t]) {
        max = this.taste[t];
        this.dominantTaste = [t];
      }
      else if(max === this.taste[t]) {
        this.dominantTaste.push(t);
      }
    }
  }

  getFoodDetails() {
    this.http.get<FooDBDetailsResult>(environment.virtuousAPI + '/foods/' + this.foodId, {responseType: 'json'}).subscribe(
      x => {
        if(x.status === 'success') {
          this.foodName = x.result.food_name;
        }
        else {
          console.log(x.status)
        }
      },
      err => console.error(err)
    );
  }

  generateId(part1, part2) {
    const secs = Math.round(Date.now() / 1000)
    return part1 + '-' + part2 + '-' + secs;
  }

  goBack(ev) {
    ev.preventDefault();
    history.back();
  }

  isDominantTaste(taste) {
    if(this.dominantTaste.includes(taste)) {
      return 'dominant-taste';
    }
    return '';
  }

}
