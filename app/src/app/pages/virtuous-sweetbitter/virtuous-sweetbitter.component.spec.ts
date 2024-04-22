import {ComponentFixture, TestBed, waitForAsync} from '@angular/core/testing';
import {IonicModule} from '@ionic/angular';

import {VirtuousSweetBitterComponent} from './virtuous-sweetbitter.component';

describe('VirtuousSweetBitterComponent', () => {
  let component: VirtuousSweetBitterComponent;
  let fixture: ComponentFixture<VirtuousSweetBitterComponent>;

  beforeEach(waitForAsync(() => {
    TestBed.configureTestingModule({
      declarations: [VirtuousSweetBitterComponent],
      imports: [IonicModule.forRoot()]
    }).compileComponents();

    fixture = TestBed.createComponent(VirtuousSweetBitterComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  }));

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
