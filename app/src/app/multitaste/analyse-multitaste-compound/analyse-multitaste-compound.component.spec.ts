import {ComponentFixture, TestBed, waitForAsync} from '@angular/core/testing';
import {IonicModule} from '@ionic/angular';

import {AnalyseMultitasteCompoundComponent} from './analyse-multitaste-compound.component';

describe('AnalyseUmamiCompoundComponent', () => {
  let component: AnalyseMultitasteCompoundComponent;
  let fixture: ComponentFixture<AnalyseMultitasteCompoundComponent>;

  beforeEach(waitForAsync(() => {
    TestBed.configureTestingModule({
      declarations: [AnalyseMultitasteCompoundComponent],
      imports: [IonicModule.forRoot()]
    }).compileComponents();

    fixture = TestBed.createComponent(AnalyseMultitasteCompoundComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  }));

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
