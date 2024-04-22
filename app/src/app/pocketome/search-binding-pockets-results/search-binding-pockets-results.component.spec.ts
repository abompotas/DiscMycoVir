import {ComponentFixture, TestBed, waitForAsync} from '@angular/core/testing';
import {IonicModule} from '@ionic/angular';

import {SearchBindingPocketsResultsComponent} from './search-binding-pockets-results.component';

describe('SearchBindingPocketsResultsComponent', () => {
  let component: SearchBindingPocketsResultsComponent;
  let fixture: ComponentFixture<SearchBindingPocketsResultsComponent>;

  beforeEach(waitForAsync(() => {
    TestBed.configureTestingModule({
      declarations: [SearchBindingPocketsResultsComponent],
      imports: [IonicModule.forRoot()]
    }).compileComponents();

    fixture = TestBed.createComponent(SearchBindingPocketsResultsComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  }));

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
