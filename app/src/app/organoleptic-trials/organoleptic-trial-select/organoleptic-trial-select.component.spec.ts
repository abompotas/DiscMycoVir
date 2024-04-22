import {ComponentFixture, TestBed, waitForAsync} from '@angular/core/testing';
import {IonicModule} from '@ionic/angular';

import {OrganolepticTrialSelectComponent} from './organoleptic-trial-select.component';

describe('OrganolepticTrialSelectComponent', () => {
  let component: OrganolepticTrialSelectComponent;
  let fixture: ComponentFixture<OrganolepticTrialSelectComponent>;

  beforeEach(waitForAsync(() => {
    TestBed.configureTestingModule({
      declarations: [OrganolepticTrialSelectComponent],
      imports: [IonicModule.forRoot()]
    }).compileComponents();

    fixture = TestBed.createComponent(OrganolepticTrialSelectComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  }));

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
