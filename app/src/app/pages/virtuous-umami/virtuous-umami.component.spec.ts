import {ComponentFixture, TestBed, waitForAsync} from '@angular/core/testing';
import {IonicModule} from '@ionic/angular';

import {VirtuousUmamiComponent} from './virtuous-umami.component';

describe('VirtuousUmamiComponent', () => {
  let component: VirtuousUmamiComponent;
  let fixture: ComponentFixture<VirtuousUmamiComponent>;

  beforeEach(waitForAsync(() => {
    TestBed.configureTestingModule({
      declarations: [VirtuousUmamiComponent],
      imports: [IonicModule.forRoot()]
    }).compileComponents();

    fixture = TestBed.createComponent(VirtuousUmamiComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  }));

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
