import { ComponentFixture, TestBed, waitForAsync } from '@angular/core/testing';
import { IonicModule } from '@ionic/angular';

import { VirusDiscoveryTrimmingComponent } from './virus-discovery-trimming.component';

describe('VirusDiscoveryTrimmingComponent', () => {
  let component: VirusDiscoveryTrimmingComponent;
  let fixture: ComponentFixture<VirusDiscoveryTrimmingComponent>;

  beforeEach(waitForAsync(() => {
    TestBed.configureTestingModule({
      declarations: [ VirusDiscoveryTrimmingComponent ],
      imports: [IonicModule.forRoot()]
    }).compileComponents();

    fixture = TestBed.createComponent(VirusDiscoveryTrimmingComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  }));

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
