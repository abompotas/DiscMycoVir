import { ComponentFixture, TestBed, waitForAsync } from '@angular/core/testing';
import { IonicModule } from '@ionic/angular';

import { SearchBindingPocketsVisualizationComponent } from './search-binding-pockets-visualization.component';

describe('SearchBindingPocketsVisualizationComponent', () => {
  let component: SearchBindingPocketsVisualizationComponent;
  let fixture: ComponentFixture<SearchBindingPocketsVisualizationComponent>;

  beforeEach(waitForAsync(() => {
    TestBed.configureTestingModule({
      declarations: [ SearchBindingPocketsVisualizationComponent ],
      imports: [IonicModule.forRoot()]
    }).compileComponents();

    fixture = TestBed.createComponent(SearchBindingPocketsVisualizationComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  }));

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
