import {ComponentFixture, TestBed, waitForAsync} from '@angular/core/testing';
import {IonicModule} from '@ionic/angular';

import {VirtuousAppMenuMobileComponent} from './virtuous-app-menu-mobile.component';

describe('VirtuousAppMenuMobileComponent', () => {
  let component: VirtuousAppMenuMobileComponent;
  let fixture: ComponentFixture<VirtuousAppMenuMobileComponent>;

  beforeEach(waitForAsync(() => {
    TestBed.configureTestingModule({
      declarations: [VirtuousAppMenuMobileComponent],
      imports: [IonicModule.forRoot()]
    }).compileComponents();

    fixture = TestBed.createComponent(VirtuousAppMenuMobileComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  }));

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
