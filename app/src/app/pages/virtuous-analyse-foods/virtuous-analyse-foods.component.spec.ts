import {ComponentFixture, TestBed, waitForAsync} from '@angular/core/testing';
import {IonicModule} from '@ionic/angular';

import {VirtuousAnalyseFoodsComponent} from './virtuous-analyse-foods.component';

describe('VirtuousAnalyseFoodsComponent', () => {
  let component: VirtuousAnalyseFoodsComponent;
  let fixture: ComponentFixture<VirtuousAnalyseFoodsComponent>;

  beforeEach(waitForAsync(() => {
    TestBed.configureTestingModule({
      declarations: [VirtuousAnalyseFoodsComponent],
      imports: [IonicModule.forRoot()]
    }).compileComponents();

    fixture = TestBed.createComponent(VirtuousAnalyseFoodsComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  }));

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
