import {NgModule} from '@angular/core';
import {RouterModule, Routes} from '@angular/router';
import {FoodsComponent} from './pages/foods/foods.component';
import {UmamiComponent} from './pages/umami/umami.component';
import {MultitasteComponent} from './pages/multitaste/multitaste.component';
import {SweetBitterComponent} from './pages/sweetbitter/sweetbitter.component';
import {VirtuousMultitasteComponent} from './pages/virtuous-multitaste/virtuous-multitaste.component';
import {VirtuousSweetBitterComponent} from './pages/virtuous-sweetbitter/virtuous-sweetbitter.component';
import {VirtuousUmamiComponent} from './pages/virtuous-umami/virtuous-umami.component';
import {VirtuousAnalyseFoodsComponent} from './pages/virtuous-analyse-foods/virtuous-analyse-foods.component';
import {VirtuousOrganolepticTrialComponent} from './pages/virtuous-organoleptic-trial/virtuous-organoleptic-trial.component';
import {VirtuousPocketomeComponent} from './pages/virtuous-pocketome/virtuous-pocketome.component';

const routes: Routes = [
  {
    path: '',
    redirectTo: 'foods',
    pathMatch: 'full'
  },
  {
    path: 'foods',
    component: FoodsComponent,
  },
  {
    path: 'foods/:food',
    component: FoodsComponent,
  },
  {
    path: 'foods/:food/taste',
    component: FoodsComponent,
  },
  {
    path: 'foods/:food/results',
    component: FoodsComponent,
  },
  {
    path: 'foods/:food/organoleptic-trial',
    component: FoodsComponent,
  },
  {
    path: 'umami',
    component: UmamiComponent,
  },
  {
    path: 'umami/results',
    component: UmamiComponent,
  },
  {
    path: 'sweetbitter',
    component: SweetBitterComponent,
  },
  {
    path: 'sweetbitter/results',
    component: SweetBitterComponent,
  },
  {
    path: 'multitaste',
    component: MultitasteComponent,
  },
  {
    path: 'multitaste/results',
    component: MultitasteComponent,
  },
  {
    path: 'virtuous-analyse-foods',
    component: VirtuousAnalyseFoodsComponent,
  },
  {
    path: 'virtuous-umami',
    component: VirtuousUmamiComponent,
  },
  {
    path: 'virtuous-umami/results',
    component: VirtuousUmamiComponent,
  },
  {
    path: 'virtuous-sweetbitter',
    component: VirtuousSweetBitterComponent,
  },
  {
    path: 'virtuous-sweetbitter/results',
    component: VirtuousSweetBitterComponent,
  },
  {
    path: 'virtuous-multitaste',
    component: VirtuousMultitasteComponent,
  },
  {
    path: 'virtuous-multitaste/results',
    component: VirtuousMultitasteComponent,
  },
  {
    path: 'virtuous-organoleptic-trial',
    component: VirtuousOrganolepticTrialComponent,
  },
  {
    path: 'virtuous-organoleptic-trial/:food',
    component: VirtuousOrganolepticTrialComponent
  },
  {
    path: 'virtuous-pocketome',
    component: VirtuousPocketomeComponent
  },
  {
    path: 'virtuous-pocketome/:job',
    component: VirtuousPocketomeComponent
  },
  {
    path: 'virtuous-pocketome/:job/:hash',
    component: VirtuousPocketomeComponent
  }
];

@NgModule({
  imports: [
    RouterModule.forRoot(routes, {useHash: true})
  ],
  exports: [RouterModule]
})
export class AppRoutingModule {
}
