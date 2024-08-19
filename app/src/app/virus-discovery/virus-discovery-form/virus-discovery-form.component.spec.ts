import {ComponentFixture, TestBed, waitForAsync} from '@angular/core/testing';
import {IonicModule} from '@ionic/angular';

import {VirusDiscoveryFormComponent} from './virus-discovery-form.component';

describe('VirusDiscoveryFormComponent', () => {
  let component: VirusDiscoveryFormComponent;
  let fixture: ComponentFixture<VirusDiscoveryFormComponent>;

  beforeEach(waitForAsync(() => {
    TestBed.configureTestingModule({
      declarations: [VirusDiscoveryFormComponent],
      imports: [IonicModule.forRoot()]
    }).compileComponents();

    fixture = TestBed.createComponent(VirusDiscoveryFormComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  }));

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
