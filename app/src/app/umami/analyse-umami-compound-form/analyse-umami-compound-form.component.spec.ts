import {ComponentFixture, TestBed, waitForAsync} from '@angular/core/testing';
import {IonicModule} from '@ionic/angular';

import {AnalyseUmamiCompoundFormComponent} from './analyse-umami-compound-form.component';

describe('AnalyseUmamiCompoundFormComponent', () => {
  let component: AnalyseUmamiCompoundFormComponent;
  let fixture: ComponentFixture<AnalyseUmamiCompoundFormComponent>;

  beforeEach(waitForAsync(() => {
    TestBed.configureTestingModule({
      declarations: [AnalyseUmamiCompoundFormComponent],
      imports: [IonicModule.forRoot()]
    }).compileComponents();

    fixture = TestBed.createComponent(AnalyseUmamiCompoundFormComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  }));

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
