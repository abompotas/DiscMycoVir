import {ComponentFixture, TestBed, waitForAsync} from '@angular/core/testing';
import {IonicModule} from '@ionic/angular';

import {AnalyseMultitasteCompoundFormComponent} from './analyse-multitaste-compound-form.component';

describe('AnalyseUmamiCompoundFormComponent', () => {
  let component: AnalyseMultitasteCompoundFormComponent;
  let fixture: ComponentFixture<AnalyseMultitasteCompoundFormComponent>;

  beforeEach(waitForAsync(() => {
    TestBed.configureTestingModule({
      declarations: [AnalyseMultitasteCompoundFormComponent],
      imports: [IonicModule.forRoot()]
    }).compileComponents();

    fixture = TestBed.createComponent(AnalyseMultitasteCompoundFormComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  }));

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
