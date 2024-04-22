import {NgModule} from '@angular/core';
import {BrowserModule} from '@angular/platform-browser';
import {RouteReuseStrategy} from '@angular/router';
import {HttpClientModule} from '@angular/common/http';
import {FormsModule, ReactiveFormsModule} from '@angular/forms';

import {IonicModule, IonicRouteStrategy} from '@ionic/angular';

import {AppComponent} from './app.component';
import {AppRoutingModule} from './app-routing.module';

import {AnalyseFoodListComponent} from './foods/analyse-food-list/analyse-food-list.component';
import {FoodsComponent} from './pages/foods/foods.component';
import {UmamiComponent} from './pages/umami/umami.component';
import {MultitasteComponent} from './pages/multitaste/multitaste.component';
import {SweetBitterComponent} from './pages/sweetbitter/sweetbitter.component';
import {OrganolepticTrialComponent} from './organoleptic-trials/organoleptic-trial/organoleptic-trial.component';
import {AnalyseFoodComponent} from './foods/analyse-food/analyse-food.component';
import {AnalyseFoodDetailsComponent} from './foods/analyse-food-details/analyse-food-details.component';
import {AnalyseFoodTasteComponent} from './foods/analyse-food-taste/analyse-food-taste.component';
import {AnalyseUmamiCompoundComponent} from './umami/analyse-umami-compound/analyse-umami-compound.component';
import {AnalyseUmamiCompoundFormComponent} from './umami/analyse-umami-compound-form/analyse-umami-compound-form.component';
import {AnalyseUmamiCompoundResultsComponent} from './umami/analyse-umami-compound-results/analyse-umami-compound-results.component';
import {AnalyseSweetBitterCompoundComponent} from './sweetbitter/analyse-sweetbitter-compound/analyse-sweetbitter-compound.component';
import {AnalyseSweetBitterCompoundFormComponent} from './sweetbitter/analyse-sweetbitter-compound-form/analyse-sweetbitter-compound-form.component';
import {AnalyseSweetBitterCompoundResultsComponent} from './sweetbitter/analyse-sweetbitter-compound-results/analyse-sweetbitter-compound-results.component';
import {AnalyseMultitasteCompoundComponent} from './multitaste/analyse-multitaste-compound/analyse-multitaste-compound.component';
import {AnalyseMultitasteCompoundFormComponent} from './multitaste/analyse-multitaste-compound-form/analyse-multitaste-compound-form.component';
import {AnalyseMultitasteCompoundResultsComponent} from './multitaste/analyse-multitaste-compound-results/analyse-multitaste-compound-results.component';
import {AnalyseMultitasteCompoundChartComponent} from './multitaste/analyse-multitaste-compound-chart/analyse-multitaste-compound-chart.component';
import {SearchBindingPocketsComponent} from './pocketome/search-binding-pockets/search-binding-pockets.component';
import {SearchBindingPocketsFormComponent} from './pocketome/search-binding-pockets-form/search-binding-pockets-form.component';
import {SearchBindingPocketsResultsComponent} from './pocketome/search-binding-pockets-results/search-binding-pockets-results.component';
import {SearchBindingPocketsVisualizationComponent} from './pocketome/search-binding-pockets-visualization/search-binding-pockets-visualization.component';
import {OrganolepticTrialSelectComponent} from './organoleptic-trials/organoleptic-trial-select/organoleptic-trial-select.component';
import {OrganolepticTrialTasteComponent} from './organoleptic-trials/organoleptic-trial-taste/organoleptic-trial-taste.component';
import {VirtuousTopbarComponent} from './commons/virtuous-topbar/virtuous-topbar.component';
import {VirtuousFooterComponent} from './commons/virtuous-footer/virtuous-footer.component';
import {VirtuousAppMenuComponent} from './commons/virtuous-app-menu/virtuous-app-menu.component';
import {VirtuousAppMenuMobileComponent} from './commons/virtuous-app-menu-mobile/virtuous-app-menu-mobile.component';
import {VirtuousUmamiComponent} from './pages/virtuous-umami/virtuous-umami.component';
import {VirtuousSweetBitterComponent} from './pages/virtuous-sweetbitter/virtuous-sweetbitter.component';
import {VirtuousMultitasteComponent} from './pages/virtuous-multitaste/virtuous-multitaste.component';
import {VirtuousAnalyseFoodsComponent} from './pages/virtuous-analyse-foods/virtuous-analyse-foods.component';
import {VirtuousOrganolepticTrialComponent} from './pages/virtuous-organoleptic-trial/virtuous-organoleptic-trial.component';
import {VirtuousPocketomeComponent} from './pages/virtuous-pocketome/virtuous-pocketome.component';


@NgModule({
  declarations: [
    AppComponent,
    FoodsComponent,
    UmamiComponent,
    MultitasteComponent,
    SweetBitterComponent,
    AnalyseUmamiCompoundComponent,
    AnalyseUmamiCompoundFormComponent,
    AnalyseUmamiCompoundResultsComponent,
    AnalyseSweetBitterCompoundComponent,
    AnalyseSweetBitterCompoundFormComponent,
    AnalyseSweetBitterCompoundResultsComponent,
    AnalyseMultitasteCompoundComponent,
    AnalyseMultitasteCompoundFormComponent,
    AnalyseMultitasteCompoundResultsComponent,
    AnalyseMultitasteCompoundChartComponent,
    SearchBindingPocketsComponent,
    SearchBindingPocketsFormComponent,
    SearchBindingPocketsResultsComponent,
    SearchBindingPocketsVisualizationComponent,
    AnalyseFoodComponent,
    AnalyseFoodListComponent,
    AnalyseFoodDetailsComponent,
    AnalyseFoodTasteComponent,
    OrganolepticTrialComponent,
    OrganolepticTrialSelectComponent,
    OrganolepticTrialTasteComponent,
    VirtuousTopbarComponent,
    VirtuousFooterComponent,
    VirtuousAppMenuComponent,
    VirtuousAppMenuMobileComponent,
    VirtuousUmamiComponent,
    VirtuousSweetBitterComponent,
    VirtuousMultitasteComponent,
    VirtuousAnalyseFoodsComponent,
    VirtuousOrganolepticTrialComponent,
    VirtuousPocketomeComponent
  ],
  entryComponents: [],
  imports: [
    BrowserModule,
    IonicModule.forRoot(),
    AppRoutingModule,
    HttpClientModule,
    FormsModule,
    ReactiveFormsModule
  ],
  providers: [{provide: RouteReuseStrategy, useClass: IonicRouteStrategy}],
  bootstrap: [AppComponent]
})
export class AppModule {
}
