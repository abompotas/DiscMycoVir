import {ComponentFixture, TestBed, waitForAsync} from '@angular/core/testing';
import {IonicModule} from '@ionic/angular';

import {AnalyseMultitasteCompoundChartComponent} from './analyse-multitaste-compound-chart.component';

describe('AnalyseMultitasteCompoundChartComponent', () => {
  let component: AnalyseMultitasteCompoundChartComponent;
  let fixture: ComponentFixture<AnalyseMultitasteCompoundChartComponent>;

  beforeEach(waitForAsync(() => {
    TestBed.configureTestingModule({
      declarations: [AnalyseMultitasteCompoundChartComponent],
      imports: [IonicModule.forRoot()]
    }).compileComponents();

    fixture = TestBed.createComponent(AnalyseMultitasteCompoundChartComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  }));

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
