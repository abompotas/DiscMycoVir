import {NgModule} from '@angular/core';
import {BrowserModule} from '@angular/platform-browser';
import {RouteReuseStrategy} from '@angular/router';
import {HttpClientModule} from '@angular/common/http';
import {FormsModule, ReactiveFormsModule} from '@angular/forms';

import {IonicModule, IonicRouteStrategy} from '@ionic/angular';

import {AppComponent} from './app.component';
import {AppRoutingModule} from './app-routing.module';

import {HomePageComponent} from './pages/home/home-page.component';
import {FooterComponent} from './commons/footer/footer.component';
import {TopbarComponent} from './commons/topbar/topbar.component';
import {VirusDiscoveryComponent} from './virus-discovery/virus-discovery/virus-discovery.component';
import {VirusDiscoveryFormComponent} from './virus-discovery/virus-discovery-form/virus-discovery-form.component';
import {VirusDiscoveryResultsComponent} from './virus-discovery/virus-discovery-results/virus-discovery-results.component';


@NgModule({
  declarations: [
    AppComponent,
    HomePageComponent,
    FooterComponent,
    TopbarComponent,
    VirusDiscoveryComponent,
    VirusDiscoveryFormComponent,
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
