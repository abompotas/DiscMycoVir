import {ComponentFixture, TestBed, waitForAsync} from '@angular/core/testing';
import {IonicModule} from '@ionic/angular';

import {AnalyseSweetBitterCompoundResultsComponent} from './analyse-sweetbitter-compound-results.component';

describe('AnalyseSweetBitterCompoundResultsComponent', () => {
  let component: AnalyseSweetBitterCompoundResultsComponent;
  let fixture: ComponentFixture<AnalyseSweetBitterCompoundResultsComponent>;

  beforeEach(waitForAsync(() => {
    TestBed.configureTestingModule({
      declarations: [AnalyseSweetBitterCompoundResultsComponent],
      imports: [IonicModule.forRoot()]
    }).compileComponents();

    fixture = TestBed.createComponent(AnalyseSweetBitterCompoundResultsComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  }));

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
