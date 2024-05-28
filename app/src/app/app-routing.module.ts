import {NgModule} from '@angular/core';
import {RouterModule, Routes} from '@angular/router';
import {HomePageComponent} from './pages/home/home-page.component';
import {TrimmingPageComponent} from './pages/trimming/trimming-page.component';
import {ResultsPageComponent} from './pages/results/results-page.component';

const routes: Routes = [
  {
    path: '',
    component: HomePageComponent
  },
  {
    path: 'trimming',
    component: TrimmingPageComponent
  },
  {
    path: 'trimming/:job',
    component: TrimmingPageComponent
  },
  {
    path: 'trimming/:job/:hash',
    component: TrimmingPageComponent
  },
  {
    path: 'results',
    component: ResultsPageComponent
  },
  {
    path: 'results/:job',
    component: ResultsPageComponent
  },
  {
    path: 'results/:job/:hash',
    component: ResultsPageComponent
  }
];

@NgModule({
  imports: [
    RouterModule.forRoot(routes, {useHash: false})
  ],
  exports: [RouterModule]
})
export class AppRoutingModule {
}
