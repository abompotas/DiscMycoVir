import {ComponentFixture, TestBed, waitForAsync} from '@angular/core/testing';
import {IonicModule} from '@ionic/angular';

import {OrganolepticTrialComponent} from './organoleptic-trial.component';

describe('SamplingReaserchFormComponent', () => {
  let component: OrganolepticTrialComponent;
  let fixture: ComponentFixture<OrganolepticTrialComponent>;

  beforeEach(waitForAsync(() => {
    TestBed.configureTestingModule({
      declarations: [OrganolepticTrialComponent],
      imports: [IonicModule.forRoot()]
    }).compileComponents();

    fixture = TestBed.createComponent(OrganolepticTrialComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  }));

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
