import {ComponentFixture, TestBed, waitForAsync} from '@angular/core/testing';
import {IonicModule} from '@ionic/angular';

import {VirtuousOrganolepticTrialComponent} from './virtuous-organoleptic-trial.component';

describe('VirtuousOrganolepticTrialComponent', () => {
  let component: VirtuousOrganolepticTrialComponent;
  let fixture: ComponentFixture<VirtuousOrganolepticTrialComponent>;

  beforeEach(waitForAsync(() => {
    TestBed.configureTestingModule({
      declarations: [VirtuousOrganolepticTrialComponent],
      imports: [IonicModule.forRoot()]
    }).compileComponents();

    fixture = TestBed.createComponent(VirtuousOrganolepticTrialComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  }));

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
