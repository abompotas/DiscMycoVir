import { ComponentFixture, TestBed, waitForAsync } from '@angular/core/testing';
import { IonicModule } from '@ionic/angular';

import { VirusDiscoveryHitCoverageComponent } from './virus-discovery-hit-coverage.component';

describe('VirusDiscoveryHitCoverageComponent', () => {
  let component: VirusDiscoveryHitCoverageComponent;
  let fixture: ComponentFixture<VirusDiscoveryHitCoverageComponent>;

  beforeEach(waitForAsync(() => {
    TestBed.configureTestingModule({
      declarations: [ VirusDiscoveryHitCoverageComponent ],
      imports: [IonicModule.forRoot()]
    }).compileComponents();

    fixture = TestBed.createComponent(VirusDiscoveryHitCoverageComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  }));

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
