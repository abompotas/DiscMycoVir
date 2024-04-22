import {ComponentFixture, TestBed, waitForAsync} from '@angular/core/testing';
import {IonicModule} from '@ionic/angular';

import {AnalyseMultitasteCompoundResultsComponent} from './analyse-multitaste-compound-results.component';

describe('AnalyseUmamiCompoundResultsComponent', () => {
  let component: AnalyseMultitasteCompoundResultsComponent;
  let fixture: ComponentFixture<AnalyseMultitasteCompoundResultsComponent>;

  beforeEach(waitForAsync(() => {
    TestBed.configureTestingModule({
      declarations: [AnalyseMultitasteCompoundResultsComponent],
      imports: [IonicModule.forRoot()]
    }).compileComponents();

    fixture = TestBed.createComponent(AnalyseMultitasteCompoundResultsComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  }));

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
