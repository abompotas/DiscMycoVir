import {ComponentFixture, TestBed, waitForAsync} from '@angular/core/testing';
import {IonicModule} from '@ionic/angular';

import {AnalyseSweetBitterCompoundComponent} from './analyse-sweetbitter-compound.component';

describe('AnalyseSweetBitterCompoundComponent', () => {
  let component: AnalyseSweetBitterCompoundComponent;
  let fixture: ComponentFixture<AnalyseSweetBitterCompoundComponent>;

  beforeEach(waitForAsync(() => {
    TestBed.configureTestingModule({
      declarations: [AnalyseSweetBitterCompoundComponent],
      imports: [IonicModule.forRoot()]
    }).compileComponents();

    fixture = TestBed.createComponent(AnalyseSweetBitterCompoundComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  }));

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
