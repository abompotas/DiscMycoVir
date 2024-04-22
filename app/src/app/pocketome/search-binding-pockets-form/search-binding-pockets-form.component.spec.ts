import {ComponentFixture, TestBed, waitForAsync} from '@angular/core/testing';
import {IonicModule} from '@ionic/angular';

import {SearchBindingPocketsFormComponent} from './search-binding-pockets-form.component';

describe('SearchBindingPocketsFormComponent', () => {
  let component: SearchBindingPocketsFormComponent;
  let fixture: ComponentFixture<SearchBindingPocketsFormComponent>;

  beforeEach(waitForAsync(() => {
    TestBed.configureTestingModule({
      declarations: [SearchBindingPocketsFormComponent],
      imports: [IonicModule.forRoot()]
    }).compileComponents();

    fixture = TestBed.createComponent(SearchBindingPocketsFormComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  }));

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
