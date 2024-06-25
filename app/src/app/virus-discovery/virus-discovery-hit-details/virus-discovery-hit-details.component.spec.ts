import { ComponentFixture, TestBed, waitForAsync } from '@angular/core/testing';
import { IonicModule } from '@ionic/angular';

import { VirusDiscoveryHitDetailsComponent } from './virus-discovery-hit-details.component';

describe('VirusDiscoveryHitDetailsComponent', () => {
  let component: VirusDiscoveryHitDetailsComponent;
  let fixture: ComponentFixture<VirusDiscoveryHitDetailsComponent>;

  beforeEach(waitForAsync(() => {
    TestBed.configureTestingModule({
      declarations: [ VirusDiscoveryHitDetailsComponent ],
      imports: [IonicModule.forRoot()]
    }).compileComponents();

    fixture = TestBed.createComponent(VirusDiscoveryHitDetailsComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  }));

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
