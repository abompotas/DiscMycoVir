import { ComponentFixture, TestBed, waitForAsync } from '@angular/core/testing';
import { IonicModule } from '@ionic/angular';

import { VirusDiscoveryHitsGraphComponent } from './virus-discovery-hits-graph.component';

describe('VirusDiscoveryHitsGraphComponent', () => {
  let component: VirusDiscoveryHitsGraphComponent;
  let fixture: ComponentFixture<VirusDiscoveryHitsGraphComponent>;

  beforeEach(waitForAsync(() => {
    TestBed.configureTestingModule({
      declarations: [ VirusDiscoveryHitsGraphComponent ],
      imports: [IonicModule.forRoot()]
    }).compileComponents();

    fixture = TestBed.createComponent(VirusDiscoveryHitsGraphComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  }));

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
