import {ComponentFixture, TestBed, waitForAsync} from '@angular/core/testing';
import {IonicModule} from '@ionic/angular';

import {VirusDiscoveryResultsComponent} from './search-binding-pockets-results.component';

describe('SearchBindingPocketsResultsComponent', () => {
  let component: VirusDiscoveryResultsComponent;
  let fixture: ComponentFixture<VirusDiscoveryResultsComponent>;

  beforeEach(waitForAsync(() => {
    TestBed.configureTestingModule({
      declarations: [VirusDiscoveryResultsComponent],
      imports: [IonicModule.forRoot()]
    }).compileComponents();

    fixture = TestBed.createComponent(VirusDiscoveryResultsComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  }));

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
