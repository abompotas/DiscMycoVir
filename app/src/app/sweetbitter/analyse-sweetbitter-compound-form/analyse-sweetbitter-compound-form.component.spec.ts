import {ComponentFixture, TestBed, waitForAsync} from '@angular/core/testing';
import {IonicModule} from '@ionic/angular';

import {AnalyseSweetBitterCompoundFormComponent} from './analyse-sweetbitter-compound-form.component';

describe('AnalyseSweetBitterCompoundFormComponent', () => {
  let component: AnalyseSweetBitterCompoundFormComponent;
  let fixture: ComponentFixture<AnalyseSweetBitterCompoundFormComponent>;

  beforeEach(waitForAsync(() => {
    TestBed.configureTestingModule({
      declarations: [AnalyseSweetBitterCompoundFormComponent],
      imports: [IonicModule.forRoot()]
    }).compileComponents();

    fixture = TestBed.createComponent(AnalyseSweetBitterCompoundFormComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  }));

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
