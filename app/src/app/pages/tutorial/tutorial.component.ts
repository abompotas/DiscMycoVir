import {Component, OnInit, ViewChild} from '@angular/core';
import {IonAccordionGroup} from '@ionic/angular';

@Component({
  selector: 'app-tutorial',
  templateUrl: './tutorial.component.html',
  styleUrls: ['./tutorial.component.scss'],
})
export class TutorialPageComponent implements OnInit {

  @ViewChild('step1', { static: true }) accordionGroup!: IonAccordionGroup;

  constructor() {
  }

  ngOnInit() {
    const nativeEl = this.accordionGroup;
    nativeEl.value = 'step1';
  }

}
