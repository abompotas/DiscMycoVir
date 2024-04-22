import {ComponentFixture, TestBed, waitForAsync} from '@angular/core/testing';
import {IonicModule} from '@ionic/angular';

import {OrganolepticTrialTasteComponent} from './organoleptic-trial-taste.component';

describe('OrganolepticTrialTasteComponent', () => {
  let component: OrganolepticTrialTasteComponent;
  let fixture: ComponentFixture<OrganolepticTrialTasteComponent>;

  beforeEach(waitForAsync(() => {
    TestBed.configureTestingModule({
      declarations: [OrganolepticTrialTasteComponent],
      imports: [IonicModule.forRoot()]
    }).compileComponents();

    fixture = TestBed.createComponent(OrganolepticTrialTasteComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  }));

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
