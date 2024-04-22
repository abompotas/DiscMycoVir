import {ComponentFixture, TestBed, waitForAsync} from '@angular/core/testing';
import {IonicModule} from '@ionic/angular';

import {AnalyseFoodTasteComponent} from './analyse-food-taste.component';

describe('AnalyseUmamiCompoundResultsComponent', () => {
  let component: AnalyseFoodTasteComponent;
  let fixture: ComponentFixture<AnalyseFoodTasteComponent>;

  beforeEach(waitForAsync(() => {
    TestBed.configureTestingModule({
      declarations: [AnalyseFoodTasteComponent],
      imports: [IonicModule.forRoot()]
    }).compileComponents();

    fixture = TestBed.createComponent(AnalyseFoodTasteComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  }));

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
