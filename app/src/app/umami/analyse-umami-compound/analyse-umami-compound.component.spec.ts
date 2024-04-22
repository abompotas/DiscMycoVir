import {ComponentFixture, TestBed, waitForAsync} from '@angular/core/testing';
import {IonicModule} from '@ionic/angular';

import {AnalyseUmamiCompoundComponent} from './analyse-umami-compound.component';

describe('AnalyseUmamiCompoundComponent', () => {
  let component: AnalyseUmamiCompoundComponent;
  let fixture: ComponentFixture<AnalyseUmamiCompoundComponent>;

  beforeEach(waitForAsync(() => {
    TestBed.configureTestingModule({
      declarations: [AnalyseUmamiCompoundComponent],
      imports: [IonicModule.forRoot()]
    }).compileComponents();

    fixture = TestBed.createComponent(AnalyseUmamiCompoundComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  }));

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
