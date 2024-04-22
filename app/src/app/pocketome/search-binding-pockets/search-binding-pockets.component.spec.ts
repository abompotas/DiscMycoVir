import {ComponentFixture, TestBed, waitForAsync} from '@angular/core/testing';
import {IonicModule} from '@ionic/angular';

import {SearchBindingPocketsComponent} from './search-binding-pockets.component';

describe('SearchBindingPocketsComponent', () => {
  let component: SearchBindingPocketsComponent;
  let fixture: ComponentFixture<SearchBindingPocketsComponent>;

  beforeEach(waitForAsync(() => {
    TestBed.configureTestingModule({
      declarations: [SearchBindingPocketsComponent],
      imports: [IonicModule.forRoot()]
    }).compileComponents();

    fixture = TestBed.createComponent(SearchBindingPocketsComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  }));

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
