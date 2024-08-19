import { ComponentFixture, TestBed, waitForAsync } from '@angular/core/testing';
import { IonicModule } from '@ionic/angular';

import { VirusDiscoveryHitsTableComponent } from './virus-discovery-hits-table.component';

describe('VirusDiscoveryHitsTableComponent', () => {
  let component: VirusDiscoveryHitsTableComponent;
  let fixture: ComponentFixture<VirusDiscoveryHitsTableComponent>;

  beforeEach(waitForAsync(() => {
    TestBed.configureTestingModule({
      declarations: [ VirusDiscoveryHitsTableComponent ],
      imports: [IonicModule.forRoot()]
    }).compileComponents();

    fixture = TestBed.createComponent(VirusDiscoveryHitsTableComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  }));

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
