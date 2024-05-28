import {NgModule} from '@angular/core';
import {BrowserModule} from '@angular/platform-browser';
import {RouteReuseStrategy} from '@angular/router';
import {HttpClientModule} from '@angular/common/http';
import {FormsModule, ReactiveFormsModule} from '@angular/forms';

import {IonicModule, IonicRouteStrategy} from '@ionic/angular';

import {AppComponent} from './app.component';
import {AppRoutingModule} from './app-routing.module';

import {HomePageComponent} from './pages/home/home-page.component';
import {TrimmingPageComponent} from './pages/trimming/trimming-page.component';
import {ResultsPageComponent} from './pages/results/results-page.component';
import {FooterComponent} from './commons/footer/footer.component';
import {TopbarComponent} from './commons/topbar/topbar.component';
import {VirusDiscoveryFormComponent} from './virus-discovery/virus-discovery-form/virus-discovery-form.component';
import {SafeHtmlPipe, VirusDiscoveryTrimmingComponent} from './virus-discovery/virus-discovery-trimming/virus-discovery-trimming.component';
import {VirusDiscoveryResultsComponent} from './virus-discovery/virus-discovery-results/virus-discovery-results.component';


@NgModule({
  declarations: [
    SafeHtmlPipe,
    AppComponent,
    HomePageComponent,
    TrimmingPageComponent,
    ResultsPageComponent,
    FooterComponent,
    TopbarComponent,
    VirusDiscoveryFormComponent,
    VirusDiscoveryTrimmingComponent,
    VirusDiscoveryResultsComponent
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
