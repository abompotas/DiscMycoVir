import {ComponentFixture, TestBed, waitForAsync} from '@angular/core/testing';
import {IonicModule} from '@ionic/angular';

import {AnalyseUmamiCompoundResultsComponent} from './analyse-umami-compound-results.component';

describe('AnalyseUmamiCompoundResultsComponent', () => {
  let component: AnalyseUmamiCompoundResultsComponent;
  let fixture: ComponentFixture<AnalyseUmamiCompoundResultsComponent>;

  beforeEach(waitForAsync(() => {
    TestBed.configureTestingModule({
      declarations: [AnalyseUmamiCompoundResultsComponent],
      imports: [IonicModule.forRoot()]
    }).compileComponents();

    fixture = TestBed.createComponent(AnalyseUmamiCompoundResultsComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  }));

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
