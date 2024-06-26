(self["webpackChunkVirusDiscovery"] = self["webpackChunkVirusDiscovery"] || []).push([["main"],{

/***/ 8255:
/*!*******************************************************!*\
  !*** ./$_lazy_route_resources/ lazy namespace object ***!
  \*******************************************************/
/***/ ((module) => {

function webpackEmptyAsyncContext(req) {
	// Here Promise.resolve().then() is used instead of new Promise() to prevent
	// uncaught exception popping up in devtools
	return Promise.resolve().then(() => {
		var e = new Error("Cannot find module '" + req + "'");
		e.code = 'MODULE_NOT_FOUND';
		throw e;
	});
}
webpackEmptyAsyncContext.keys = () => ([]);
webpackEmptyAsyncContext.resolve = webpackEmptyAsyncContext;
webpackEmptyAsyncContext.id = 8255;
module.exports = webpackEmptyAsyncContext;

/***/ }),

/***/ 158:
/*!***************************************!*\
  !*** ./src/app/app-routing.module.ts ***!
  \***************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "AppRoutingModule": () => (/* binding */ AppRoutingModule)
/* harmony export */ });
/* harmony import */ var tslib__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! tslib */ 5353);
/* harmony import */ var _angular_core__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! @angular/core */ 7716);
/* harmony import */ var _angular_router__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! @angular/router */ 9895);
/* harmony import */ var _pages_home_home_page_component__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./pages/home/home-page.component */ 2249);
/* harmony import */ var _pages_trimming_trimming_page_component__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./pages/trimming/trimming-page.component */ 102);
/* harmony import */ var _pages_results_results_page_component__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./pages/results/results-page.component */ 3209);






const routes = [
    {
        path: '',
        component: _pages_home_home_page_component__WEBPACK_IMPORTED_MODULE_0__.HomePageComponent
    },
    {
        path: 'trimming',
        component: _pages_trimming_trimming_page_component__WEBPACK_IMPORTED_MODULE_1__.TrimmingPageComponent
    },
    {
        path: 'trimming/:job',
        component: _pages_trimming_trimming_page_component__WEBPACK_IMPORTED_MODULE_1__.TrimmingPageComponent
    },
    {
        path: 'trimming/:job/:hash',
        component: _pages_trimming_trimming_page_component__WEBPACK_IMPORTED_MODULE_1__.TrimmingPageComponent
    },
    {
        path: 'results',
        component: _pages_results_results_page_component__WEBPACK_IMPORTED_MODULE_2__.ResultsPageComponent
    },
    {
        path: 'results/:job',
        component: _pages_results_results_page_component__WEBPACK_IMPORTED_MODULE_2__.ResultsPageComponent
    },
    {
        path: 'results/:job/:hash',
        component: _pages_results_results_page_component__WEBPACK_IMPORTED_MODULE_2__.ResultsPageComponent
    }
];
let AppRoutingModule = class AppRoutingModule {
};
AppRoutingModule = (0,tslib__WEBPACK_IMPORTED_MODULE_3__.__decorate)([
    (0,_angular_core__WEBPACK_IMPORTED_MODULE_4__.NgModule)({
        imports: [
            _angular_router__WEBPACK_IMPORTED_MODULE_5__.RouterModule.forRoot(routes, { useHash: true })
        ],
        exports: [_angular_router__WEBPACK_IMPORTED_MODULE_5__.RouterModule]
    })
], AppRoutingModule);



/***/ }),

/***/ 5041:
/*!**********************************!*\
  !*** ./src/app/app.component.ts ***!
  \**********************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "AppComponent": () => (/* binding */ AppComponent)
/* harmony export */ });
/* harmony import */ var tslib__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! tslib */ 5353);
/* harmony import */ var _raw_loader_app_component_html__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! !raw-loader!./app.component.html */ 1106);
/* harmony import */ var _app_component_scss__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./app.component.scss */ 3069);
/* harmony import */ var _angular_core__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! @angular/core */ 7716);
/* harmony import */ var _ionic_angular__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! @ionic/angular */ 9122);





let AppComponent = class AppComponent {
    constructor(platform, alertController) {
        this.platform = platform;
        this.alertController = alertController;
        this.platform.backButton.subscribe(() => (0,tslib__WEBPACK_IMPORTED_MODULE_2__.__awaiter)(this, void 0, void 0, function* () {
            const alert = yield this.alertController.create({
                header: 'Exit the Mycovirus Discovery App?',
                message: 'You are going to exit the Mycovirus Discovery App',
                buttons: [
                    { text: 'Cancel', role: 'cancel' },
                    {
                        text: 'OK',
                        role: 'confirm',
                        handler: () => {
                            if (navigator.hasOwnProperty('app')) {
                                // @ts-ignore
                                navigator.app.exitApp();
                            }
                            else if (navigator.hasOwnProperty('device')) {
                                // @ts-ignore
                                navigator.device.exitApp();
                            }
                        },
                    }
                ]
            });
            yield alert.present();
        }));
    }
};
AppComponent.ctorParameters = () => [
    { type: _ionic_angular__WEBPACK_IMPORTED_MODULE_3__.Platform },
    { type: _ionic_angular__WEBPACK_IMPORTED_MODULE_3__.AlertController }
];
AppComponent = (0,tslib__WEBPACK_IMPORTED_MODULE_2__.__decorate)([
    (0,_angular_core__WEBPACK_IMPORTED_MODULE_4__.Component)({
        selector: 'app-root',
        template: _raw_loader_app_component_html__WEBPACK_IMPORTED_MODULE_0__.default,
        styles: [_app_component_scss__WEBPACK_IMPORTED_MODULE_1__.default]
    })
], AppComponent);



/***/ }),

/***/ 6747:
/*!*******************************!*\
  !*** ./src/app/app.module.ts ***!
  \*******************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "AppModule": () => (/* binding */ AppModule)
/* harmony export */ });
/* harmony import */ var tslib__WEBPACK_IMPORTED_MODULE_13__ = __webpack_require__(/*! tslib */ 5353);
/* harmony import */ var _angular_core__WEBPACK_IMPORTED_MODULE_14__ = __webpack_require__(/*! @angular/core */ 7716);
/* harmony import */ var _angular_platform_browser__WEBPACK_IMPORTED_MODULE_15__ = __webpack_require__(/*! @angular/platform-browser */ 9075);
/* harmony import */ var _angular_router__WEBPACK_IMPORTED_MODULE_19__ = __webpack_require__(/*! @angular/router */ 9895);
/* harmony import */ var _angular_common_http__WEBPACK_IMPORTED_MODULE_17__ = __webpack_require__(/*! @angular/common/http */ 1841);
/* harmony import */ var _angular_forms__WEBPACK_IMPORTED_MODULE_18__ = __webpack_require__(/*! @angular/forms */ 3679);
/* harmony import */ var _ionic_angular__WEBPACK_IMPORTED_MODULE_16__ = __webpack_require__(/*! @ionic/angular */ 9122);
/* harmony import */ var _app_component__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./app.component */ 5041);
/* harmony import */ var _app_routing_module__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./app-routing.module */ 158);
/* harmony import */ var _pages_home_home_page_component__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./pages/home/home-page.component */ 2249);
/* harmony import */ var _pages_trimming_trimming_page_component__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./pages/trimming/trimming-page.component */ 102);
/* harmony import */ var _pages_results_results_page_component__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ./pages/results/results-page.component */ 3209);
/* harmony import */ var _commons_footer_footer_component__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! ./commons/footer/footer.component */ 3758);
/* harmony import */ var _commons_topbar_topbar_component__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(/*! ./commons/topbar/topbar.component */ 2608);
/* harmony import */ var _virus_discovery_virus_discovery_form_virus_discovery_form_component__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(/*! ./virus-discovery/virus-discovery-form/virus-discovery-form.component */ 3743);
/* harmony import */ var _virus_discovery_virus_discovery_trimming_virus_discovery_trimming_component__WEBPACK_IMPORTED_MODULE_8__ = __webpack_require__(/*! ./virus-discovery/virus-discovery-trimming/virus-discovery-trimming.component */ 7575);
/* harmony import */ var _virus_discovery_virus_discovery_results_virus_discovery_results_component__WEBPACK_IMPORTED_MODULE_9__ = __webpack_require__(/*! ./virus-discovery/virus-discovery-results/virus-discovery-results.component */ 1126);
/* harmony import */ var _virus_discovery_virus_discovery_hits_graph_virus_discovery_hits_graph_component__WEBPACK_IMPORTED_MODULE_10__ = __webpack_require__(/*! ./virus-discovery/virus-discovery-hits-graph/virus-discovery-hits-graph.component */ 5017);
/* harmony import */ var _virus_discovery_virus_discovery_hits_table_virus_discovery_hits_table_component__WEBPACK_IMPORTED_MODULE_11__ = __webpack_require__(/*! ./virus-discovery/virus-discovery-hits-table/virus-discovery-hits-table.component */ 125);
/* harmony import */ var _virus_discovery_virus_discovery_hit_details_virus_discovery_hit_details_component__WEBPACK_IMPORTED_MODULE_12__ = __webpack_require__(/*! ./virus-discovery/virus-discovery-hit-details/virus-discovery-hit-details.component */ 2365);




















let AppModule = class AppModule {
};
AppModule = (0,tslib__WEBPACK_IMPORTED_MODULE_13__.__decorate)([
    (0,_angular_core__WEBPACK_IMPORTED_MODULE_14__.NgModule)({
        declarations: [
            _virus_discovery_virus_discovery_trimming_virus_discovery_trimming_component__WEBPACK_IMPORTED_MODULE_8__.SafeHtmlPipe,
            _app_component__WEBPACK_IMPORTED_MODULE_0__.AppComponent,
            _pages_home_home_page_component__WEBPACK_IMPORTED_MODULE_2__.HomePageComponent,
            _pages_trimming_trimming_page_component__WEBPACK_IMPORTED_MODULE_3__.TrimmingPageComponent,
            _pages_results_results_page_component__WEBPACK_IMPORTED_MODULE_4__.ResultsPageComponent,
            _commons_footer_footer_component__WEBPACK_IMPORTED_MODULE_5__.FooterComponent,
            _commons_topbar_topbar_component__WEBPACK_IMPORTED_MODULE_6__.TopbarComponent,
            _virus_discovery_virus_discovery_form_virus_discovery_form_component__WEBPACK_IMPORTED_MODULE_7__.VirusDiscoveryFormComponent,
            _virus_discovery_virus_discovery_trimming_virus_discovery_trimming_component__WEBPACK_IMPORTED_MODULE_8__.VirusDiscoveryTrimmingComponent,
            _virus_discovery_virus_discovery_results_virus_discovery_results_component__WEBPACK_IMPORTED_MODULE_9__.VirusDiscoveryResultsComponent,
            _virus_discovery_virus_discovery_hits_table_virus_discovery_hits_table_component__WEBPACK_IMPORTED_MODULE_11__.VirusDiscoveryHitsTableComponent,
            _virus_discovery_virus_discovery_hit_details_virus_discovery_hit_details_component__WEBPACK_IMPORTED_MODULE_12__.VirusDiscoveryHitDetailsComponent,
            _virus_discovery_virus_discovery_hits_graph_virus_discovery_hits_graph_component__WEBPACK_IMPORTED_MODULE_10__.VirusDiscoveryHitsGraphComponent
        ],
        entryComponents: [],
        imports: [
            _angular_platform_browser__WEBPACK_IMPORTED_MODULE_15__.BrowserModule,
            _ionic_angular__WEBPACK_IMPORTED_MODULE_16__.IonicModule.forRoot(),
            _app_routing_module__WEBPACK_IMPORTED_MODULE_1__.AppRoutingModule,
            _angular_common_http__WEBPACK_IMPORTED_MODULE_17__.HttpClientModule,
            _angular_forms__WEBPACK_IMPORTED_MODULE_18__.FormsModule,
            _angular_forms__WEBPACK_IMPORTED_MODULE_18__.ReactiveFormsModule
        ],
        providers: [{ provide: _angular_router__WEBPACK_IMPORTED_MODULE_19__.RouteReuseStrategy, useClass: _ionic_angular__WEBPACK_IMPORTED_MODULE_16__.IonicRouteStrategy }],
        bootstrap: [_app_component__WEBPACK_IMPORTED_MODULE_0__.AppComponent]
    })
], AppModule);



/***/ }),

/***/ 3758:
/*!****************************************************!*\
  !*** ./src/app/commons/footer/footer.component.ts ***!
  \****************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "FooterComponent": () => (/* binding */ FooterComponent)
/* harmony export */ });
/* harmony import */ var tslib__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! tslib */ 5353);
/* harmony import */ var _raw_loader_footer_component_html__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! !raw-loader!./footer.component.html */ 5302);
/* harmony import */ var _footer_component_scss__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./footer.component.scss */ 2696);
/* harmony import */ var _angular_core__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! @angular/core */ 7716);




let FooterComponent = class FooterComponent {
    constructor() {
        this.spacer = false;
    }
    ngOnInit() {
    }
};
FooterComponent.ctorParameters = () => [];
FooterComponent.propDecorators = {
    spacer: [{ type: _angular_core__WEBPACK_IMPORTED_MODULE_2__.Input }]
};
FooterComponent = (0,tslib__WEBPACK_IMPORTED_MODULE_3__.__decorate)([
    (0,_angular_core__WEBPACK_IMPORTED_MODULE_2__.Component)({
        selector: 'app-footer',
        template: _raw_loader_footer_component_html__WEBPACK_IMPORTED_MODULE_0__.default,
        styles: [_footer_component_scss__WEBPACK_IMPORTED_MODULE_1__.default]
    })
], FooterComponent);



/***/ }),

/***/ 2608:
/*!****************************************************!*\
  !*** ./src/app/commons/topbar/topbar.component.ts ***!
  \****************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "TopbarComponent": () => (/* binding */ TopbarComponent)
/* harmony export */ });
/* harmony import */ var tslib__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! tslib */ 5353);
/* harmony import */ var _raw_loader_topbar_component_html__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! !raw-loader!./topbar.component.html */ 4503);
/* harmony import */ var _topbar_component_scss__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./topbar.component.scss */ 772);
/* harmony import */ var _angular_core__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! @angular/core */ 7716);




let TopbarComponent = class TopbarComponent {
    constructor() {
        this.buttons = true;
    }
    ngOnInit() {
    }
};
TopbarComponent.ctorParameters = () => [];
TopbarComponent.propDecorators = {
    buttons: [{ type: _angular_core__WEBPACK_IMPORTED_MODULE_2__.Input }]
};
TopbarComponent = (0,tslib__WEBPACK_IMPORTED_MODULE_3__.__decorate)([
    (0,_angular_core__WEBPACK_IMPORTED_MODULE_2__.Component)({
        selector: 'app-topbar',
        template: _raw_loader_topbar_component_html__WEBPACK_IMPORTED_MODULE_0__.default,
        styles: [_topbar_component_scss__WEBPACK_IMPORTED_MODULE_1__.default]
    })
], TopbarComponent);



/***/ }),

/***/ 2249:
/*!***************************************************!*\
  !*** ./src/app/pages/home/home-page.component.ts ***!
  \***************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "HomePageComponent": () => (/* binding */ HomePageComponent)
/* harmony export */ });
/* harmony import */ var tslib__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! tslib */ 5353);
/* harmony import */ var _raw_loader_home_page_component_html__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! !raw-loader!./home-page.component.html */ 4217);
/* harmony import */ var _home_page_component_scss__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./home-page.component.scss */ 4728);
/* harmony import */ var _angular_core__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! @angular/core */ 7716);




let HomePageComponent = class HomePageComponent {
    constructor() {
    }
    ngOnInit() {
    }
};
HomePageComponent.ctorParameters = () => [];
HomePageComponent = (0,tslib__WEBPACK_IMPORTED_MODULE_2__.__decorate)([
    (0,_angular_core__WEBPACK_IMPORTED_MODULE_3__.Component)({
        selector: 'app-home-page',
        template: _raw_loader_home_page_component_html__WEBPACK_IMPORTED_MODULE_0__.default,
        styles: [_home_page_component_scss__WEBPACK_IMPORTED_MODULE_1__.default]
    })
], HomePageComponent);



/***/ }),

/***/ 3209:
/*!*********************************************************!*\
  !*** ./src/app/pages/results/results-page.component.ts ***!
  \*********************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "ResultsPageComponent": () => (/* binding */ ResultsPageComponent)
/* harmony export */ });
/* harmony import */ var tslib__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! tslib */ 5353);
/* harmony import */ var _raw_loader_results_page_component_html__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! !raw-loader!./results-page.component.html */ 8069);
/* harmony import */ var _results_page_component_scss__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./results-page.component.scss */ 1308);
/* harmony import */ var _angular_core__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! @angular/core */ 7716);




let ResultsPageComponent = class ResultsPageComponent {
    constructor() { }
    ngOnInit() { }
};
ResultsPageComponent.ctorParameters = () => [];
ResultsPageComponent = (0,tslib__WEBPACK_IMPORTED_MODULE_2__.__decorate)([
    (0,_angular_core__WEBPACK_IMPORTED_MODULE_3__.Component)({
        selector: 'app-results-page',
        template: _raw_loader_results_page_component_html__WEBPACK_IMPORTED_MODULE_0__.default,
        styles: [_results_page_component_scss__WEBPACK_IMPORTED_MODULE_1__.default]
    })
], ResultsPageComponent);



/***/ }),

/***/ 102:
/*!***********************************************************!*\
  !*** ./src/app/pages/trimming/trimming-page.component.ts ***!
  \***********************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "TrimmingPageComponent": () => (/* binding */ TrimmingPageComponent)
/* harmony export */ });
/* harmony import */ var tslib__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! tslib */ 5353);
/* harmony import */ var _raw_loader_trimming_page_component_html__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! !raw-loader!./trimming-page.component.html */ 4210);
/* harmony import */ var _trimming_page_component_scss__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./trimming-page.component.scss */ 4422);
/* harmony import */ var _angular_core__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! @angular/core */ 7716);




let TrimmingPageComponent = class TrimmingPageComponent {
    constructor() { }
    ngOnInit() { }
};
TrimmingPageComponent.ctorParameters = () => [];
TrimmingPageComponent = (0,tslib__WEBPACK_IMPORTED_MODULE_2__.__decorate)([
    (0,_angular_core__WEBPACK_IMPORTED_MODULE_3__.Component)({
        selector: 'app-trimming-page',
        template: _raw_loader_trimming_page_component_html__WEBPACK_IMPORTED_MODULE_0__.default,
        styles: [_trimming_page_component_scss__WEBPACK_IMPORTED_MODULE_1__.default]
    })
], TrimmingPageComponent);



/***/ }),

/***/ 3743:
/*!****************************************************************************************!*\
  !*** ./src/app/virus-discovery/virus-discovery-form/virus-discovery-form.component.ts ***!
  \****************************************************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "VirusDiscoveryFormComponent": () => (/* binding */ VirusDiscoveryFormComponent)
/* harmony export */ });
/* harmony import */ var tslib__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! tslib */ 5353);
/* harmony import */ var _raw_loader_virus_discovery_form_component_html__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! !raw-loader!./virus-discovery-form.component.html */ 8790);
/* harmony import */ var _virus_discovery_form_component_scss__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./virus-discovery-form.component.scss */ 8827);
/* harmony import */ var _angular_core__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(/*! @angular/core */ 7716);
/* harmony import */ var _angular_common_http__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! @angular/common/http */ 1841);
/* harmony import */ var _angular_router__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! @angular/router */ 9895);
/* harmony import */ var _ionic_angular__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(/*! @ionic/angular */ 9122);
/* harmony import */ var _environments_environment__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ../../../environments/environment */ 2340);








let VirusDiscoveryFormComponent = class VirusDiscoveryFormComponent {
    constructor(http, router, alertController, loadingController) {
        this.http = http;
        this.router = router;
        this.alertController = alertController;
        this.loadingController = loadingController;
        this.initForm();
    }
    initForm() {
        this.email = null;
        this.sampleName = null;
        this.sequencingTechnology = null;
        this.singleFile = null;
        this.forwardFile = null;
        this.reverseFile = null;
        this.referenceGenome = null;
    }
    ngOnInit() {
    }
    onSingleFileChange(event) {
        this.singleFile = event.target.children['single_file'].files[0];
    }
    onForwardFileChange(event) {
        this.forwardFile = event.target.children['forward_file'].files[0];
    }
    onReverseFileChange(event) {
        this.reverseFile = event.target.children['reverse_file'].files[0];
    }
    onGenomeFileChange(event) {
        this.referenceGenome = event.target.children['reference_genome'].files[0];
    }
    search() {
        this.loading().then(() => {
            if (this.validateForm()) {
                const formData = new FormData();
                formData.append('email', this.email);
                formData.append('sample_name', this.sampleName);
                formData.append('sequencing_technology', this.sequencingTechnology);
                if (this.sequencingTechnology === 'single') {
                    formData.append('single_file', this.singleFile, this.singleFile.name);
                }
                else if (this.sequencingTechnology === 'paired') {
                    formData.append('forward_file', this.forwardFile, this.forwardFile.name);
                    formData.append('reverse_file', this.reverseFile, this.reverseFile.name);
                }
                formData.append('reference_genome', this.referenceGenome, this.referenceGenome.name);
                this.http.post(_environments_environment__WEBPACK_IMPORTED_MODULE_2__.environment.discvirAPI + '/job', formData, { responseType: 'json' }).subscribe(x => this.response(x), e => this.error(e.error), () => {
                    this.initForm();
                    this.loadingController.dismiss().then(null);
                });
            }
        });
    }
    validateForm() {
        let validEmail = true;
        if (this.email === null) {
            validEmail = false;
        }
        else {
            // noinspection RegExpRedundantEscape
            const matches = this.email.match(/^(([^<>()[\]\\.,;:\s@\"]+(\.[^<>()[\]\\.,;:\s@\"]+)*)|(\".+\"))@((\[[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}\])|(([a-zA-Z\-0-9]+\.)+[a-zA-Z]{2,}))$/);
            if (!matches) {
                validEmail = false;
            }
        }
        if (!validEmail) {
            this.error({ error: 'Please fill in a valid email address.' }).then(null);
            return false;
        }
        let validName = true;
        if (this.sampleName === null) {
            validName = false;
        }
        else {
            const matches = this.sampleName.match(/^\w+$/);
            if (!matches) {
                validName = false;
            }
        }
        if (!validName) {
            this.error({ error: 'Please fill in a name for your experiment (only alphanumeric characters are permitted).' }).then(null);
            return false;
        }
        return this.validateInputFiles();
    }
    validateInputFiles() {
        if (this.sequencingTechnology === null) {
            this.error({ error: 'Please select the sequencing technology type.' }).then(null);
            return false;
        }
        else {
            if (this.sequencingTechnology === 'single') {
                if (this.singleFile === null) {
                    this.error({ error: 'Please select the input file to search in.' }).then(null);
                    return false;
                }
            }
            else if (this.sequencingTechnology === 'paired') {
                if (this.forwardFile === null) {
                    this.error({ error: 'Please select the forward read input file to search in.' }).then(null);
                    return false;
                }
                if (this.reverseFile === null) {
                    this.error({ error: 'Please select the reverse read input file to search in.' }).then(null);
                    return false;
                }
            }
            else {
                this.error({ error: 'Unknown sequencing technology type.' }).then(null);
                return false;
            }
        }
        if (this.referenceGenome === null) {
            this.error({ error: 'Please select a file for the reference genome.' }).then(null);
            return false;
        }
        return true;
    }
    response(resp) {
        if (resp.status === 'success') {
            this.alertSuccess().then(null);
        }
        else {
            this.error(resp).then(null);
        }
    }
    loading() {
        return (0,tslib__WEBPACK_IMPORTED_MODULE_3__.__awaiter)(this, void 0, void 0, function* () {
            const loading = yield this.loadingController.create({
                message: 'Please wait...',
            });
            yield loading.present();
        });
    }
    error(resp) {
        return (0,tslib__WEBPACK_IMPORTED_MODULE_3__.__awaiter)(this, void 0, void 0, function* () {
            this.loadingController.dismiss().then(() => {
                let msg = 'Check your input for missing values.';
                if (resp.hasOwnProperty('error')) {
                    msg = resp.error;
                }
                this.alertError(msg);
            });
        });
    }
    alertSuccess() {
        return (0,tslib__WEBPACK_IMPORTED_MODULE_3__.__awaiter)(this, void 0, void 0, function* () {
            const alert = yield this.alertController.create({
                header: 'Success!',
                message: 'Your query has been submitted. Once the search is completed you will receive an email containing the results.',
                buttons: ['OK']
            });
            yield alert.present();
        });
    }
    alertError(msg) {
        return (0,tslib__WEBPACK_IMPORTED_MODULE_3__.__awaiter)(this, void 0, void 0, function* () {
            const alert = yield this.alertController.create({
                header: 'Error!',
                message: msg,
                buttons: ['OK']
            });
            yield alert.present();
        });
    }
};
VirusDiscoveryFormComponent.ctorParameters = () => [
    { type: _angular_common_http__WEBPACK_IMPORTED_MODULE_4__.HttpClient },
    { type: _angular_router__WEBPACK_IMPORTED_MODULE_5__.Router },
    { type: _ionic_angular__WEBPACK_IMPORTED_MODULE_6__.AlertController },
    { type: _ionic_angular__WEBPACK_IMPORTED_MODULE_6__.LoadingController }
];
VirusDiscoveryFormComponent = (0,tslib__WEBPACK_IMPORTED_MODULE_3__.__decorate)([
    (0,_angular_core__WEBPACK_IMPORTED_MODULE_7__.Component)({
        selector: 'app-virus-discovery-form',
        template: _raw_loader_virus_discovery_form_component_html__WEBPACK_IMPORTED_MODULE_0__.default,
        styles: [_virus_discovery_form_component_scss__WEBPACK_IMPORTED_MODULE_1__.default]
    })
], VirusDiscoveryFormComponent);



/***/ }),

/***/ 2365:
/*!******************************************************************************************************!*\
  !*** ./src/app/virus-discovery/virus-discovery-hit-details/virus-discovery-hit-details.component.ts ***!
  \******************************************************************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "VirusDiscoveryHitDetailsComponent": () => (/* binding */ VirusDiscoveryHitDetailsComponent)
/* harmony export */ });
/* harmony import */ var tslib__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! tslib */ 5353);
/* harmony import */ var _raw_loader_virus_discovery_hit_details_component_html__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! !raw-loader!./virus-discovery-hit-details.component.html */ 8546);
/* harmony import */ var _virus_discovery_hit_details_component_scss__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./virus-discovery-hit-details.component.scss */ 3407);
/* harmony import */ var _angular_core__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! @angular/core */ 7716);




let VirusDiscoveryHitDetailsComponent = class VirusDiscoveryHitDetailsComponent {
    constructor() {
    }
    ngOnInit() {
    }
};
VirusDiscoveryHitDetailsComponent.ctorParameters = () => [];
VirusDiscoveryHitDetailsComponent.propDecorators = {
    hsp: [{ type: _angular_core__WEBPACK_IMPORTED_MODULE_2__.Input }]
};
VirusDiscoveryHitDetailsComponent = (0,tslib__WEBPACK_IMPORTED_MODULE_3__.__decorate)([
    (0,_angular_core__WEBPACK_IMPORTED_MODULE_2__.Component)({
        selector: 'app-virus-discovery-hit-details',
        template: _raw_loader_virus_discovery_hit_details_component_html__WEBPACK_IMPORTED_MODULE_0__.default,
        styles: [_virus_discovery_hit_details_component_scss__WEBPACK_IMPORTED_MODULE_1__.default]
    })
], VirusDiscoveryHitDetailsComponent);



/***/ }),

/***/ 5017:
/*!****************************************************************************************************!*\
  !*** ./src/app/virus-discovery/virus-discovery-hits-graph/virus-discovery-hits-graph.component.ts ***!
  \****************************************************************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "VirusDiscoveryHitsGraphComponent": () => (/* binding */ VirusDiscoveryHitsGraphComponent)
/* harmony export */ });
/* harmony import */ var tslib__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! tslib */ 5353);
/* harmony import */ var _raw_loader_virus_discovery_hits_graph_component_html__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! !raw-loader!./virus-discovery-hits-graph.component.html */ 7566);
/* harmony import */ var _virus_discovery_hits_graph_component_scss__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./virus-discovery-hits-graph.component.scss */ 1129);
/* harmony import */ var _angular_core__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! @angular/core */ 7716);
/* harmony import */ var chart_js_auto__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! chart.js/auto */ 5649);





let VirusDiscoveryHitsGraphComponent = class VirusDiscoveryHitsGraphComponent {
    constructor() {
        this.canvas = null;
        this.rootStyle = null;
        this.graphData = null;
        this.datasets = { '200': [], '80-200': [], '50-80': [], '40-50': [], '0-40': [] };
    }
    ngOnInit() {
        this.rootStyle = getComputedStyle(document.body);
        this.createDatasets();
        this.graphData = {
            datasets: [{
                    label: '>=200',
                    data: this.datasets['200'],
                    barThickness: 15,
                    borderWidth: 2,
                    borderColor: this.rootStyle.getPropertyValue('--ion-color-danger-shade'),
                    backgroundColor: this.rootStyle.getPropertyValue('--ion-color-danger-tint'),
                    stack: 'stack-0',
                }, {
                    label: '80-200',
                    data: this.datasets['80-200'],
                    barThickness: 15,
                    borderWidth: 2,
                    borderColor: this.rootStyle.getPropertyValue('--ion-color-warning-shade'),
                    backgroundColor: this.rootStyle.getPropertyValue('--ion-color-warning-tint'),
                    stack: 'stack-0',
                }, {
                    label: '50-80',
                    data: this.datasets['50-80'],
                    barThickness: 15,
                    borderWidth: 2,
                    borderColor: this.rootStyle.getPropertyValue('--ion-color-secondary-shade'),
                    backgroundColor: this.rootStyle.getPropertyValue('--ion-color-secondary-tint'),
                    stack: 'stack-0',
                }, {
                    label: '40-50',
                    data: this.datasets['40-50'],
                    barThickness: 15,
                    borderWidth: 2,
                    borderColor: this.rootStyle.getPropertyValue('--ion-color-primary-shade'),
                    backgroundColor: this.rootStyle.getPropertyValue('--ion-color-primary-tint'),
                    stack: 'stack-0',
                }, {
                    label: '<40',
                    data: this.datasets['0-40'],
                    barThickness: 15,
                    borderWidth: 2,
                    borderColor: this.rootStyle.getPropertyValue('--ion-color-dark-shade'),
                    backgroundColor: this.rootStyle.getPropertyValue('--ion-color-dark'),
                    stack: 'stack-0',
                }]
        };
    }
    ngAfterViewInit() {
        this.canvas = document.getElementById('hits-graph-' + this.qid);
        this.canvas.height = this.qData.alignments.length * 6 + 15;
        new chart_js_auto__WEBPACK_IMPORTED_MODULE_2__.Chart(this.canvas.getContext('2d'), {
            type: 'bar',
            data: this.graphData,
            options: {
                indexAxis: 'y',
                scales: {
                    y: { display: false },
                    x: {
                        title: {
                            text: 'Query',
                            display: true,
                            padding: { top: 10, bottom: 5 },
                            color: '#ffffff'
                        },
                        ticks: {
                            color: '#ffffff',
                            align: 'inner',
                            padding: 0
                        },
                        position: 'top',
                        backgroundColor: this.rootStyle.getPropertyValue('--ion-color-tertiary-tint'),
                        max: this.qData.queryLetters
                    }
                },
                plugins: {
                    legend: { position: 'top' }
                }
            }
        });
    }
    createDatasets() {
        for (let a of this.qData.alignments) {
            for (let h of a.hsps) {
                const dataPoint = {
                    y: 'Score: ' + h.score + ', Alignment length: ' + h.alignLength + ', E-value: ' + h.expect,
                    x: [h.queryStart, h.queryEnd]
                };
                if (h.score >= 200) {
                    this.datasets['200'].push(dataPoint);
                }
                else if (h.score < 200 && h.score >= 80) {
                    this.datasets['80-200'].push(dataPoint);
                }
                else if (h.score < 80 && h.score >= 50) {
                    this.datasets['50-80'].push(dataPoint);
                }
                else if (h.score < 50 && h.score >= 40) {
                    this.datasets['40-50'].push(dataPoint);
                }
                else {
                    this.datasets['0-40'].push(dataPoint);
                }
            }
        }
    }
};
VirusDiscoveryHitsGraphComponent.ctorParameters = () => [];
VirusDiscoveryHitsGraphComponent.propDecorators = {
    qid: [{ type: _angular_core__WEBPACK_IMPORTED_MODULE_3__.Input }],
    qData: [{ type: _angular_core__WEBPACK_IMPORTED_MODULE_3__.Input }]
};
VirusDiscoveryHitsGraphComponent = (0,tslib__WEBPACK_IMPORTED_MODULE_4__.__decorate)([
    (0,_angular_core__WEBPACK_IMPORTED_MODULE_3__.Component)({
        selector: 'app-virus-discovery-hits-graph',
        template: _raw_loader_virus_discovery_hits_graph_component_html__WEBPACK_IMPORTED_MODULE_0__.default,
        styles: [_virus_discovery_hits_graph_component_scss__WEBPACK_IMPORTED_MODULE_1__.default]
    })
], VirusDiscoveryHitsGraphComponent);



/***/ }),

/***/ 125:
/*!****************************************************************************************************!*\
  !*** ./src/app/virus-discovery/virus-discovery-hits-table/virus-discovery-hits-table.component.ts ***!
  \****************************************************************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "VirusDiscoveryHitsTableComponent": () => (/* binding */ VirusDiscoveryHitsTableComponent)
/* harmony export */ });
/* harmony import */ var tslib__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! tslib */ 5353);
/* harmony import */ var _raw_loader_virus_discovery_hits_table_component_html__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! !raw-loader!./virus-discovery-hits-table.component.html */ 8723);
/* harmony import */ var _virus_discovery_hits_table_component_scss__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./virus-discovery-hits-table.component.scss */ 3697);
/* harmony import */ var _angular_core__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! @angular/core */ 7716);
/* harmony import */ var datatables_net_dt__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! datatables.net-dt */ 2661);





let VirusDiscoveryHitsTableComponent = class VirusDiscoveryHitsTableComponent {
    constructor() {
        this.alignmentsData = [];
        this.matchesData = [];
        this.selectedHSP = null;
    }
    ngOnInit() {
        for (let a of this.qData.alignments) {
            for (let h of a.hsps) {
                h.title = a.hitId + a.hitDef;
                h.length = a.length;
                this.alignmentsData.push(h);
                this.matchesData.push(this.parseHSPMatch(h, 50));
            }
        }
    }
    ngAfterViewInit() {
        new datatables_net_dt__WEBPACK_IMPORTED_MODULE_2__.default('#alignments-' + this.qid, {
            autoWidth: false,
            order: [[2, 'desc']],
            columns: [
                { width: '30%' },
                { width: '7%' },
                { width: '7%' },
                { width: '7%' },
                { width: '7%' },
                { width: '7%' },
                { width: '7%' },
                { width: '7%' },
                { width: '7%' },
                { width: '7%' },
                { width: '7%' }
            ]
        });
    }
    parseHSPMatch(hsp, limit) {
        const parsedData = {
            'queryStart': [], 'queryEnd': [], 'query': [],
            'sbjctStart': [], 'sbjctEnd': [], 'sbjct': [],
            'match': [], 'title': hsp.title
        };
        let step = 1;
        if (hsp.strand[1] === 'Minus') {
            step = -1;
        }
        let query = '';
        let sbjct = '';
        let match = '';
        let queryStart = hsp.queryStart - 1;
        let queryEnd = queryStart;
        let sbjctStart = hsp.sbjctStart - step;
        let sbjctEnd = sbjctStart;
        for (let i = 0; i < hsp.query.length; i++) {
            query += hsp.query[i];
            sbjct += hsp.sbjct[i];
            match += hsp.match[i];
            if (hsp.query[i] !== '-') {
                queryEnd++;
            }
            if (hsp.sbjct[i] !== '-') {
                sbjctEnd += step;
            }
            if ((i % limit) === 0) {
                if (hsp.query[i] !== '-') {
                    queryStart++;
                }
                if (hsp.sbjct[i] !== '-') {
                    sbjctStart += step;
                }
            }
            if ((i % limit) === (limit - 1)) {
                parsedData.queryStart.push(queryStart);
                parsedData.queryEnd.push(queryEnd);
                parsedData.query.push(query);
                parsedData.sbjctStart.push(sbjctStart);
                parsedData.sbjctEnd.push(sbjctEnd);
                parsedData.sbjct.push(sbjct);
                parsedData.match.push(match);
                query = '';
                sbjct = '';
                match = '';
                queryStart = queryEnd;
                sbjctStart = sbjctEnd;
            }
        }
        if (query !== '') {
            parsedData.queryStart.push(queryStart);
            parsedData.queryEnd.push(hsp.queryEnd);
            parsedData.query.push(query);
            parsedData.sbjctStart.push(sbjctStart);
            parsedData.sbjctEnd.push(hsp.sbjctEnd);
            parsedData.sbjct.push(sbjct);
            parsedData.match.push(match);
        }
        return parsedData;
    }
    showHSPDetails(i) {
        this.selectedHSP = this.matchesData[i];
        return false;
    }
};
VirusDiscoveryHitsTableComponent.ctorParameters = () => [];
VirusDiscoveryHitsTableComponent.propDecorators = {
    qid: [{ type: _angular_core__WEBPACK_IMPORTED_MODULE_3__.Input }],
    qData: [{ type: _angular_core__WEBPACK_IMPORTED_MODULE_3__.Input }]
};
VirusDiscoveryHitsTableComponent = (0,tslib__WEBPACK_IMPORTED_MODULE_4__.__decorate)([
    (0,_angular_core__WEBPACK_IMPORTED_MODULE_3__.Component)({
        selector: 'app-virus-discovery-hits-table',
        template: _raw_loader_virus_discovery_hits_table_component_html__WEBPACK_IMPORTED_MODULE_0__.default,
        styles: [_virus_discovery_hits_table_component_scss__WEBPACK_IMPORTED_MODULE_1__.default]
    })
], VirusDiscoveryHitsTableComponent);



/***/ }),

/***/ 1126:
/*!**********************************************************************************************!*\
  !*** ./src/app/virus-discovery/virus-discovery-results/virus-discovery-results.component.ts ***!
  \**********************************************************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "VirusDiscoveryResultsComponent": () => (/* binding */ VirusDiscoveryResultsComponent)
/* harmony export */ });
/* harmony import */ var tslib__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! tslib */ 5353);
/* harmony import */ var _raw_loader_virus_discovery_results_component_html__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! !raw-loader!./virus-discovery-results.component.html */ 9064);
/* harmony import */ var _virus_discovery_results_component_scss__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./virus-discovery-results.component.scss */ 8081);
/* harmony import */ var _angular_core__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(/*! @angular/core */ 7716);
/* harmony import */ var _angular_common_http__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! @angular/common/http */ 1841);
/* harmony import */ var _angular_router__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! @angular/router */ 9895);
/* harmony import */ var _ionic_angular__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(/*! @ionic/angular */ 9122);
/* harmony import */ var _environments_environment__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ../../../environments/environment */ 2340);








// noinspection DuplicatedCode
let VirusDiscoveryResultsComponent = class VirusDiscoveryResultsComponent {
    constructor(http, route, router, alertController, loadingController) {
        this.http = http;
        this.route = route;
        this.router = router;
        this.alertController = alertController;
        this.loadingController = loadingController;
        this.resultsURL = '';
        this.downloadURL = '';
        this.jobId = 0;
        this.hash = '';
        this.results = null;
        this.route.params.subscribe(params => {
            if (params.hasOwnProperty('job')) {
                if (params.hasOwnProperty('hash')) {
                    this.jobId = params.job;
                    this.hash = params.hash;
                }
            }
        });
    }
    ngOnInit() {
        if ((this.jobId) === 0 || (this.hash === '')) {
            this.router.navigate(['/']);
        }
        else {
            this.resultsURL = _environments_environment__WEBPACK_IMPORTED_MODULE_2__.environment.discvirAPI + '/results/' + this.jobId + '/' + this.hash;
            this.downloadURL = this.resultsURL + '/download';
            this.getResults();
        }
    }
    getResults() {
        this.loading().then(() => {
            this.http.get(this.resultsURL, { responseType: 'json' }).subscribe(x => this.parseResults(x), e => this.resultsError(e.error), () => this.loadingController.dismiss().then(null));
        });
    }
    parseResults(resp) {
        this.results = [...resp.results];
    }
    loading() {
        return (0,tslib__WEBPACK_IMPORTED_MODULE_3__.__awaiter)(this, void 0, void 0, function* () {
            const loading = yield this.loadingController.create({
                message: 'Please wait...',
            });
            yield loading.present();
        });
    }
    resultsError(resp) {
        return (0,tslib__WEBPACK_IMPORTED_MODULE_3__.__awaiter)(this, void 0, void 0, function* () {
            this.loadingController.dismiss().then(() => {
                let msg = 'Could not find any results matching your request.';
                if (resp.hasOwnProperty('error')) {
                    msg = resp.error;
                }
                this.alertError(msg);
            });
        });
    }
    alertError(msg) {
        return (0,tslib__WEBPACK_IMPORTED_MODULE_3__.__awaiter)(this, void 0, void 0, function* () {
            const alert = yield this.alertController.create({
                header: 'Error!',
                message: msg,
                buttons: ['OK']
            });
            yield alert.present();
        });
    }
};
VirusDiscoveryResultsComponent.ctorParameters = () => [
    { type: _angular_common_http__WEBPACK_IMPORTED_MODULE_4__.HttpClient },
    { type: _angular_router__WEBPACK_IMPORTED_MODULE_5__.ActivatedRoute },
    { type: _angular_router__WEBPACK_IMPORTED_MODULE_5__.Router },
    { type: _ionic_angular__WEBPACK_IMPORTED_MODULE_6__.AlertController },
    { type: _ionic_angular__WEBPACK_IMPORTED_MODULE_6__.LoadingController }
];
VirusDiscoveryResultsComponent = (0,tslib__WEBPACK_IMPORTED_MODULE_3__.__decorate)([
    (0,_angular_core__WEBPACK_IMPORTED_MODULE_7__.Component)({
        selector: 'app-virus-discovery-results',
        template: _raw_loader_virus_discovery_results_component_html__WEBPACK_IMPORTED_MODULE_0__.default,
        styles: [_virus_discovery_results_component_scss__WEBPACK_IMPORTED_MODULE_1__.default]
    })
], VirusDiscoveryResultsComponent);



/***/ }),

/***/ 7575:
/*!************************************************************************************************!*\
  !*** ./src/app/virus-discovery/virus-discovery-trimming/virus-discovery-trimming.component.ts ***!
  \************************************************************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "SafeHtmlPipe": () => (/* binding */ SafeHtmlPipe),
/* harmony export */   "VirusDiscoveryTrimmingComponent": () => (/* binding */ VirusDiscoveryTrimmingComponent)
/* harmony export */ });
/* harmony import */ var tslib__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! tslib */ 5353);
/* harmony import */ var _raw_loader_virus_discovery_trimming_component_html__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! !raw-loader!./virus-discovery-trimming.component.html */ 3352);
/* harmony import */ var _virus_discovery_trimming_component_scss__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./virus-discovery-trimming.component.scss */ 1400);
/* harmony import */ var _angular_core__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! @angular/core */ 7716);
/* harmony import */ var _angular_common_http__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(/*! @angular/common/http */ 1841);
/* harmony import */ var _angular_router__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(/*! @angular/router */ 9895);
/* harmony import */ var _angular_platform_browser__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! @angular/platform-browser */ 9075);
/* harmony import */ var _ionic_angular__WEBPACK_IMPORTED_MODULE_8__ = __webpack_require__(/*! @ionic/angular */ 9122);
/* harmony import */ var _environments_environment__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ../../../environments/environment */ 2340);









let SafeHtmlPipe = class SafeHtmlPipe {
    constructor(sanitizer) {
        this.sanitizer = sanitizer;
    }
    transform(inHtml) {
        let outHtml = inHtml.replaceAll(/ href="#M\d+/g, '');
        outHtml = outHtml.replaceAll('href=', 'target="_blank" href=');
        return this.sanitizer.bypassSecurityTrustHtml(outHtml);
    }
};
SafeHtmlPipe.ctorParameters = () => [
    { type: _angular_platform_browser__WEBPACK_IMPORTED_MODULE_3__.DomSanitizer }
];
SafeHtmlPipe = (0,tslib__WEBPACK_IMPORTED_MODULE_4__.__decorate)([
    (0,_angular_core__WEBPACK_IMPORTED_MODULE_5__.Pipe)({ name: 'safeIFrame' })
], SafeHtmlPipe);

// noinspection DuplicatedCode
let VirusDiscoveryTrimmingComponent = class VirusDiscoveryTrimmingComponent {
    constructor(http, route, router, alertController, loadingController) {
        this.http = http;
        this.route = route;
        this.router = router;
        this.alertController = alertController;
        this.loadingController = loadingController;
        this.jobId = 0;
        this.hash = '';
        this.analysisReports = [];
        this.route.params.subscribe(params => {
            if (params.hasOwnProperty('job')) {
                if (params.hasOwnProperty('hash')) {
                    this.jobId = params.job;
                    this.hash = params.hash;
                }
            }
        });
        this.initForm();
    }
    ngOnInit() {
        if ((this.jobId) === 0 || (this.hash === '')) {
            this.router.navigate(['/']).then(null);
        }
        else {
            this.fetchAnalysis();
        }
    }
    initForm() {
        this.sequencingTechnology = null;
        this.slidingWindow = null;
        this.minLength = null;
        this.adapter = null;
    }
    fetchAnalysis() {
        this.loading().then(() => {
            this.http.get(_environments_environment__WEBPACK_IMPORTED_MODULE_2__.environment.discvirAPI + '/analysis/' + this.jobId + '/' + this.hash, { responseType: 'json' }).subscribe(x => {
                this.analysisReports = [...x.reports];
            }, e => this.error(e.error), () => this.loadingController.dismiss().then(null));
        });
    }
    onAdapterChange(event) {
        this.adapter = event.target.children['adapter'].files[0];
    }
    trim() {
        this.loading().then(() => {
            const formData = new FormData();
            if (this.adapter !== null) {
                formData.append('adapter', this.adapter, this.adapter.name);
            }
            if (this.slidingWindow !== null) {
                formData.append('sliding_window', this.slidingWindow);
            }
            if (this.minLength !== null) {
                formData.append('min_length', this.minLength);
            }
            this.http.put(_environments_environment__WEBPACK_IMPORTED_MODULE_2__.environment.discvirAPI + '/trimming/' + this.jobId + '/' + this.hash, formData, { responseType: 'json' }).subscribe(x => this.response(x), e => this.error(e.error), () => {
                this.initForm();
                this.loadingController.dismiss().then(null);
            });
        });
    }
    proceed() {
        this.loading().then(() => {
            this.http.put(_environments_environment__WEBPACK_IMPORTED_MODULE_2__.environment.discvirAPI + '/discovery/' + this.jobId + '/' + this.hash, {}, { responseType: 'json' }).subscribe(x => this.response(x), e => this.error(e.error), () => {
                this.initForm();
                this.loadingController.dismiss().then(null);
            });
        });
    }
    response(resp) {
        if (resp.status === 'success') {
            this.alertSuccess().then(null);
        }
        else {
            this.error(resp).then(null);
        }
    }
    loading() {
        return (0,tslib__WEBPACK_IMPORTED_MODULE_4__.__awaiter)(this, void 0, void 0, function* () {
            const loading = yield this.loadingController.create({
                message: 'Please wait...',
            });
            yield loading.present();
        });
    }
    error(resp) {
        return (0,tslib__WEBPACK_IMPORTED_MODULE_4__.__awaiter)(this, void 0, void 0, function* () {
            this.loadingController.dismiss().then(() => {
                let msg = 'Oops! Something went wrong...';
                if (resp.hasOwnProperty('error')) {
                    msg = resp.error;
                }
                this.alertError(msg);
            });
        });
    }
    alertSuccess() {
        return (0,tslib__WEBPACK_IMPORTED_MODULE_4__.__awaiter)(this, void 0, void 0, function* () {
            const alert = yield this.alertController.create({
                header: 'Success!',
                message: 'Your query has been submitted. Once the search is completed you will receive an email containing the results.',
                buttons: ['OK']
            });
            yield alert.present();
        });
    }
    alertError(msg) {
        return (0,tslib__WEBPACK_IMPORTED_MODULE_4__.__awaiter)(this, void 0, void 0, function* () {
            const alert = yield this.alertController.create({
                header: 'Error!',
                message: msg,
                buttons: ['OK']
            });
            yield alert.present();
        });
    }
};
VirusDiscoveryTrimmingComponent.ctorParameters = () => [
    { type: _angular_common_http__WEBPACK_IMPORTED_MODULE_6__.HttpClient },
    { type: _angular_router__WEBPACK_IMPORTED_MODULE_7__.ActivatedRoute },
    { type: _angular_router__WEBPACK_IMPORTED_MODULE_7__.Router },
    { type: _ionic_angular__WEBPACK_IMPORTED_MODULE_8__.AlertController },
    { type: _ionic_angular__WEBPACK_IMPORTED_MODULE_8__.LoadingController }
];
VirusDiscoveryTrimmingComponent = (0,tslib__WEBPACK_IMPORTED_MODULE_4__.__decorate)([
    (0,_angular_core__WEBPACK_IMPORTED_MODULE_5__.Component)({
        selector: 'app-virus-discovery-trimming',
        template: _raw_loader_virus_discovery_trimming_component_html__WEBPACK_IMPORTED_MODULE_0__.default,
        styles: [_virus_discovery_trimming_component_scss__WEBPACK_IMPORTED_MODULE_1__.default]
    })
], VirusDiscoveryTrimmingComponent);



/***/ }),

/***/ 2340:
/*!*****************************************!*\
  !*** ./src/environments/environment.ts ***!
  \*****************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "environment": () => (/* binding */ environment)
/* harmony export */ });
// This file can be replaced during build by using the `fileReplacements` array.
// `ng build --prod` replaces `environment.ts` with `environment.prod.ts`.
// The list of file replacements can be found in `angular.json`.
const environment = {
    production: false,
    discvirAPI: 'http://localhost:5000',
    discvir: 'http://localhost:8100'
};
/*
 * For easier debugging in development mode, you can import the following file
 * to ignore zone related error stack frames such as `zone.run`, `zoneDelegate.invokeTask`.
 *
 * This import should be commented out in production mode because it will have a negative impact
 * on performance if an error is thrown.
 */
// import 'zone.js/dist/zone-error';  // Included with Angular CLI.


/***/ }),

/***/ 4431:
/*!*********************!*\
  !*** ./src/main.ts ***!
  \*********************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony import */ var _angular_core__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! @angular/core */ 7716);
/* harmony import */ var _angular_platform_browser_dynamic__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! @angular/platform-browser-dynamic */ 4608);
/* harmony import */ var _app_app_module__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./app/app.module */ 6747);
/* harmony import */ var _environments_environment__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./environments/environment */ 2340);




if (_environments_environment__WEBPACK_IMPORTED_MODULE_1__.environment.production) {
    (0,_angular_core__WEBPACK_IMPORTED_MODULE_2__.enableProdMode)();
}
(0,_angular_platform_browser_dynamic__WEBPACK_IMPORTED_MODULE_3__.platformBrowserDynamic)().bootstrapModule(_app_app_module__WEBPACK_IMPORTED_MODULE_0__.AppModule)
    .catch(err => console.log(err));


/***/ }),

/***/ 863:
/*!******************************************************************************************************************************************!*\
  !*** ./node_modules/@ionic/core/dist/esm/ lazy ^\.\/.*\.entry\.js$ include: \.entry\.js$ exclude: \.system\.entry\.js$ namespace object ***!
  \******************************************************************************************************************************************/
/***/ ((module, __unused_webpack_exports, __webpack_require__) => {

var map = {
	"./ion-accordion_2.entry.js": [
		8359,
		"common",
		"node_modules_ionic_core_dist_esm_ion-accordion_2_entry_js"
	],
	"./ion-action-sheet.entry.js": [
		7321,
		"common",
		"node_modules_ionic_core_dist_esm_ion-action-sheet_entry_js"
	],
	"./ion-alert.entry.js": [
		6108,
		"common",
		"node_modules_ionic_core_dist_esm_ion-alert_entry_js"
	],
	"./ion-app_8.entry.js": [
		1489,
		"common",
		"node_modules_ionic_core_dist_esm_ion-app_8_entry_js"
	],
	"./ion-avatar_3.entry.js": [
		305,
		"node_modules_ionic_core_dist_esm_ion-avatar_3_entry_js"
	],
	"./ion-back-button.entry.js": [
		5830,
		"common",
		"node_modules_ionic_core_dist_esm_ion-back-button_entry_js"
	],
	"./ion-backdrop.entry.js": [
		7757,
		"node_modules_ionic_core_dist_esm_ion-backdrop_entry_js"
	],
	"./ion-breadcrumb_2.entry.js": [
		4355,
		"common",
		"node_modules_ionic_core_dist_esm_ion-breadcrumb_2_entry_js"
	],
	"./ion-button_2.entry.js": [
		392,
		"node_modules_ionic_core_dist_esm_ion-button_2_entry_js"
	],
	"./ion-card_5.entry.js": [
		6911,
		"node_modules_ionic_core_dist_esm_ion-card_5_entry_js"
	],
	"./ion-checkbox.entry.js": [
		937,
		"node_modules_ionic_core_dist_esm_ion-checkbox_entry_js"
	],
	"./ion-chip.entry.js": [
		8695,
		"node_modules_ionic_core_dist_esm_ion-chip_entry_js"
	],
	"./ion-col_3.entry.js": [
		6034,
		"node_modules_ionic_core_dist_esm_ion-col_3_entry_js"
	],
	"./ion-datetime-button.entry.js": [
		1135,
		"default-node_modules_ionic_core_dist_esm_data-caf38df0_js-node_modules_ionic_core_dist_esm_th-3ec348",
		"node_modules_ionic_core_dist_esm_ion-datetime-button_entry_js"
	],
	"./ion-datetime_3.entry.js": [
		8837,
		"default-node_modules_ionic_core_dist_esm_data-caf38df0_js-node_modules_ionic_core_dist_esm_th-3ec348",
		"common",
		"node_modules_ionic_core_dist_esm_ion-datetime_3_entry_js"
	],
	"./ion-fab_3.entry.js": [
		4195,
		"common",
		"node_modules_ionic_core_dist_esm_ion-fab_3_entry_js"
	],
	"./ion-img.entry.js": [
		1709,
		"node_modules_ionic_core_dist_esm_ion-img_entry_js"
	],
	"./ion-infinite-scroll_2.entry.js": [
		3087,
		"common",
		"node_modules_ionic_core_dist_esm_ion-infinite-scroll_2_entry_js"
	],
	"./ion-input.entry.js": [
		4513,
		"common",
		"node_modules_ionic_core_dist_esm_ion-input_entry_js"
	],
	"./ion-item-option_3.entry.js": [
		8056,
		"common",
		"node_modules_ionic_core_dist_esm_ion-item-option_3_entry_js"
	],
	"./ion-item_8.entry.js": [
		862,
		"common",
		"node_modules_ionic_core_dist_esm_ion-item_8_entry_js"
	],
	"./ion-loading.entry.js": [
		7509,
		"node_modules_ionic_core_dist_esm_ion-loading_entry_js"
	],
	"./ion-menu_3.entry.js": [
		6272,
		"common",
		"node_modules_ionic_core_dist_esm_ion-menu_3_entry_js"
	],
	"./ion-modal.entry.js": [
		1855,
		"common",
		"node_modules_ionic_core_dist_esm_ion-modal_entry_js"
	],
	"./ion-nav_2.entry.js": [
		8708,
		"common",
		"node_modules_ionic_core_dist_esm_ion-nav_2_entry_js"
	],
	"./ion-picker-column-internal.entry.js": [
		1349,
		"common",
		"node_modules_ionic_core_dist_esm_ion-picker-column-internal_entry_js"
	],
	"./ion-picker-internal.entry.js": [
		7915,
		"node_modules_ionic_core_dist_esm_ion-picker-internal_entry_js"
	],
	"./ion-popover.entry.js": [
		3527,
		"common",
		"node_modules_ionic_core_dist_esm_ion-popover_entry_js"
	],
	"./ion-progress-bar.entry.js": [
		4694,
		"node_modules_ionic_core_dist_esm_ion-progress-bar_entry_js"
	],
	"./ion-radio_2.entry.js": [
		9222,
		"node_modules_ionic_core_dist_esm_ion-radio_2_entry_js"
	],
	"./ion-range.entry.js": [
		5277,
		"common",
		"node_modules_ionic_core_dist_esm_ion-range_entry_js"
	],
	"./ion-refresher_2.entry.js": [
		9921,
		"common",
		"node_modules_ionic_core_dist_esm_ion-refresher_2_entry_js"
	],
	"./ion-reorder_2.entry.js": [
		3122,
		"common",
		"node_modules_ionic_core_dist_esm_ion-reorder_2_entry_js"
	],
	"./ion-ripple-effect.entry.js": [
		1602,
		"node_modules_ionic_core_dist_esm_ion-ripple-effect_entry_js"
	],
	"./ion-route_4.entry.js": [
		5174,
		"node_modules_ionic_core_dist_esm_ion-route_4_entry_js"
	],
	"./ion-searchbar.entry.js": [
		7895,
		"common",
		"node_modules_ionic_core_dist_esm_ion-searchbar_entry_js"
	],
	"./ion-segment_2.entry.js": [
		6164,
		"common",
		"node_modules_ionic_core_dist_esm_ion-segment_2_entry_js"
	],
	"./ion-select_3.entry.js": [
		592,
		"node_modules_ionic_core_dist_esm_ion-select_3_entry_js"
	],
	"./ion-slide_2.entry.js": [
		7162,
		"node_modules_ionic_core_dist_esm_ion-slide_2_entry_js"
	],
	"./ion-spinner.entry.js": [
		1374,
		"common",
		"node_modules_ionic_core_dist_esm_ion-spinner_entry_js"
	],
	"./ion-split-pane.entry.js": [
		7896,
		"node_modules_ionic_core_dist_esm_ion-split-pane_entry_js"
	],
	"./ion-tab-bar_2.entry.js": [
		5043,
		"common",
		"node_modules_ionic_core_dist_esm_ion-tab-bar_2_entry_js"
	],
	"./ion-tab_2.entry.js": [
		7802,
		"common",
		"node_modules_ionic_core_dist_esm_ion-tab_2_entry_js"
	],
	"./ion-text.entry.js": [
		9072,
		"node_modules_ionic_core_dist_esm_ion-text_entry_js"
	],
	"./ion-textarea.entry.js": [
		2191,
		"node_modules_ionic_core_dist_esm_ion-textarea_entry_js"
	],
	"./ion-toast.entry.js": [
		801,
		"node_modules_ionic_core_dist_esm_ion-toast_entry_js"
	],
	"./ion-toggle.entry.js": [
		7110,
		"common",
		"node_modules_ionic_core_dist_esm_ion-toggle_entry_js"
	],
	"./ion-virtual-scroll.entry.js": [
		431,
		"node_modules_ionic_core_dist_esm_ion-virtual-scroll_entry_js"
	]
};
function webpackAsyncContext(req) {
	if(!__webpack_require__.o(map, req)) {
		return Promise.resolve().then(() => {
			var e = new Error("Cannot find module '" + req + "'");
			e.code = 'MODULE_NOT_FOUND';
			throw e;
		});
	}

	var ids = map[req], id = ids[0];
	return Promise.all(ids.slice(1).map(__webpack_require__.e)).then(() => {
		return __webpack_require__(id);
	});
}
webpackAsyncContext.keys = () => (Object.keys(map));
webpackAsyncContext.id = 863;
module.exports = webpackAsyncContext;

/***/ }),

/***/ 3069:
/*!************************************!*\
  !*** ./src/app/app.component.scss ***!
  \************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "default": () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = ("#main-container {\n  background: #262626;\n}\n/*# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbImFwcC5jb21wb25lbnQuc2NzcyJdLCJuYW1lcyI6W10sIm1hcHBpbmdzIjoiQUFBQTtFQUNFLG1CQUFBO0FBQ0YiLCJmaWxlIjoiYXBwLmNvbXBvbmVudC5zY3NzIiwic291cmNlc0NvbnRlbnQiOlsiI21haW4tY29udGFpbmVyIHtcbiAgYmFja2dyb3VuZDogIzI2MjYyNjtcbn1cbiJdfQ== */");

/***/ }),

/***/ 2696:
/*!******************************************************!*\
  !*** ./src/app/commons/footer/footer.component.scss ***!
  \******************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "default": () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (":host {\n  font-size: 0.75em;\n}\n\n.footer {\n  background: #171717;\n}\n\n.spacer {\n  padding-top: 5.5em;\n  background: #171717;\n}\n\na ion-icon {\n  margin: 0;\n  padding: 0 0.25em 0;\n  font-size: 2em;\n  color: var(--ion-color-light);\n}\n\n.disclaimer {\n  display: flex;\n  justify-content: center;\n  align-items: center;\n  color: var(--ion-color-light);\n}\n/*# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbImZvb3Rlci5jb21wb25lbnQuc2NzcyJdLCJuYW1lcyI6W10sIm1hcHBpbmdzIjoiQUFBQTtFQUNFLGlCQUFBO0FBQ0Y7O0FBRUE7RUFDRSxtQkFBQTtBQUNGOztBQUVBO0VBQ0Usa0JBQUE7RUFDQSxtQkFBQTtBQUNGOztBQUVBO0VBQ0UsU0FBQTtFQUNBLG1CQUFBO0VBQ0EsY0FBQTtFQUNBLDZCQUFBO0FBQ0Y7O0FBRUE7RUFDRSxhQUFBO0VBQ0EsdUJBQUE7RUFDQSxtQkFBQTtFQUNBLDZCQUFBO0FBQ0YiLCJmaWxlIjoiZm9vdGVyLmNvbXBvbmVudC5zY3NzIiwic291cmNlc0NvbnRlbnQiOlsiOmhvc3Qge1xuICBmb250LXNpemU6IDAuNzVlbTtcbn1cblxuLmZvb3RlciB7XG4gIGJhY2tncm91bmQ6ICMxNzE3MTc7XG59XG5cbi5zcGFjZXIge1xuICBwYWRkaW5nLXRvcDogNS41ZW07XG4gIGJhY2tncm91bmQ6ICMxNzE3MTc7XG59XG5cbmEgaW9uLWljb24ge1xuICBtYXJnaW46IDA7XG4gIHBhZGRpbmc6IDAgMC4yNWVtIDA7XG4gIGZvbnQtc2l6ZTogMmVtO1xuICBjb2xvcjogdmFyKC0taW9uLWNvbG9yLWxpZ2h0KTtcbn1cblxuLmRpc2NsYWltZXIge1xuICBkaXNwbGF5OiBmbGV4O1xuICBqdXN0aWZ5LWNvbnRlbnQ6IGNlbnRlcjtcbiAgYWxpZ24taXRlbXM6IGNlbnRlcjtcbiAgY29sb3I6IHZhcigtLWlvbi1jb2xvci1saWdodCk7XG59XG4iXX0= */");

/***/ }),

/***/ 772:
/*!******************************************************!*\
  !*** ./src/app/commons/topbar/topbar.component.scss ***!
  \******************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "default": () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = ("ion-toolbar {\n  --padding-top: 5px;\n  --padding-bottom: 5px;\n  --background: var(--ion-color-primary);\n  --min-height: 80px;\n}\n\nion-buttons {\n  margin: 0 5px;\n  padding: 0;\n}\n\nion-buttons a {\n  margin: 0;\n  padding: 0;\n}\n\nion-buttons a ion-icon {\n  margin-top: 5px;\n  font-size: 42px;\n  color: #eeeeee;\n}\n\nion-img.large {\n  height: 75px;\n  margin: 2.5px;\n}\n\nion-img.small {\n  height: 60px;\n}\n\nion-chip {\n  background-color: #eeeeee;\n}\n\nion-chip ion-icon {\n  color: #333333;\n  font-size: 25px;\n}\n/*# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbInRvcGJhci5jb21wb25lbnQuc2NzcyJdLCJuYW1lcyI6W10sIm1hcHBpbmdzIjoiQUFBQTtFQUNFLGtCQUFBO0VBQ0EscUJBQUE7RUFDQSxzQ0FBQTtFQUNBLGtCQUFBO0FBQ0Y7O0FBRUE7RUFDRSxhQUFBO0VBQ0EsVUFBQTtBQUNGOztBQUVBO0VBQ0UsU0FBQTtFQUNBLFVBQUE7QUFDRjs7QUFFQTtFQUNFLGVBQUE7RUFDQSxlQUFBO0VBQ0EsY0FBQTtBQUNGOztBQUVBO0VBQ0UsWUFBQTtFQUNBLGFBQUE7QUFDRjs7QUFFQTtFQUNFLFlBQUE7QUFDRjs7QUFFQTtFQUNFLHlCQUFBO0FBQ0Y7O0FBRUE7RUFDRSxjQUFBO0VBQ0EsZUFBQTtBQUNGIiwiZmlsZSI6InRvcGJhci5jb21wb25lbnQuc2NzcyIsInNvdXJjZXNDb250ZW50IjpbImlvbi10b29sYmFyIHtcbiAgLS1wYWRkaW5nLXRvcDogNXB4O1xuICAtLXBhZGRpbmctYm90dG9tOiA1cHg7XG4gIC0tYmFja2dyb3VuZDogdmFyKC0taW9uLWNvbG9yLXByaW1hcnkpO1xuICAtLW1pbi1oZWlnaHQ6IDgwcHg7XG59XG5cbmlvbi1idXR0b25zIHtcbiAgbWFyZ2luOiAwIDVweDtcbiAgcGFkZGluZzogMDtcbn1cblxuaW9uLWJ1dHRvbnMgYSB7XG4gIG1hcmdpbjogMDtcbiAgcGFkZGluZzogMDtcbn1cblxuaW9uLWJ1dHRvbnMgYSBpb24taWNvbiB7XG4gIG1hcmdpbi10b3A6IDVweDtcbiAgZm9udC1zaXplOiA0MnB4O1xuICBjb2xvcjogI2VlZWVlZTtcbn1cblxuaW9uLWltZy5sYXJnZSB7XG4gIGhlaWdodDogNzVweDtcbiAgbWFyZ2luOiAyLjVweDtcbn1cblxuaW9uLWltZy5zbWFsbCB7XG4gIGhlaWdodDogNjBweDtcbn1cblxuaW9uLWNoaXAge1xuICBiYWNrZ3JvdW5kLWNvbG9yOiAjZWVlZWVlO1xufVxuXG5pb24tY2hpcCBpb24taWNvbiB7XG4gIGNvbG9yOiAjMzMzMzMzO1xuICBmb250LXNpemU6IDI1cHg7XG59XG4iXX0= */");

/***/ }),

/***/ 4728:
/*!*****************************************************!*\
  !*** ./src/app/pages/home/home-page.component.scss ***!
  \*****************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "default": () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = ("\n/*# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbXSwibmFtZXMiOltdLCJtYXBwaW5ncyI6IiIsImZpbGUiOiJob21lLXBhZ2UuY29tcG9uZW50LnNjc3MifQ== */");

/***/ }),

/***/ 1308:
/*!***********************************************************!*\
  !*** ./src/app/pages/results/results-page.component.scss ***!
  \***********************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "default": () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = ("\n/*# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbXSwibmFtZXMiOltdLCJtYXBwaW5ncyI6IiIsImZpbGUiOiJyZXN1bHRzLXBhZ2UuY29tcG9uZW50LnNjc3MifQ== */");

/***/ }),

/***/ 4422:
/*!*************************************************************!*\
  !*** ./src/app/pages/trimming/trimming-page.component.scss ***!
  \*************************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "default": () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = ("\n/*# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbXSwibmFtZXMiOltdLCJtYXBwaW5ncyI6IiIsImZpbGUiOiJ0cmltbWluZy1wYWdlLmNvbXBvbmVudC5zY3NzIn0= */");

/***/ }),

/***/ 8827:
/*!******************************************************************************************!*\
  !*** ./src/app/virus-discovery/virus-discovery-form/virus-discovery-form.component.scss ***!
  \******************************************************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "default": () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = ("ion-col.vcenter {\n  display: flex;\n  align-items: center;\n}\n\nion-col.section-head {\n  padding: 20px 0 10px;\n}\n\nion-input, ion-select {\n  background-color: #fdfdfd;\n  border: 1px solid var(--ion-color-medium);\n  border-radius: 0.75em;\n}\n\nion-button {\n  margin: 0;\n}\n/*# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbInZpcnVzLWRpc2NvdmVyeS1mb3JtLmNvbXBvbmVudC5zY3NzIl0sIm5hbWVzIjpbXSwibWFwcGluZ3MiOiJBQUFBO0VBQ0UsYUFBQTtFQUNBLG1CQUFBO0FBQ0Y7O0FBRUE7RUFDRSxvQkFBQTtBQUNGOztBQUVBO0VBQ0UseUJBQUE7RUFDQSx5Q0FBQTtFQUNBLHFCQUFBO0FBQ0Y7O0FBRUE7RUFDRSxTQUFBO0FBQ0YiLCJmaWxlIjoidmlydXMtZGlzY292ZXJ5LWZvcm0uY29tcG9uZW50LnNjc3MiLCJzb3VyY2VzQ29udGVudCI6WyJpb24tY29sLnZjZW50ZXIge1xuICBkaXNwbGF5OiBmbGV4O1xuICBhbGlnbi1pdGVtczogY2VudGVyO1xufVxuXG5pb24tY29sLnNlY3Rpb24taGVhZCB7XG4gIHBhZGRpbmc6IDIwcHggMCAxMHB4O1xufVxuXG5pb24taW5wdXQsIGlvbi1zZWxlY3QgIHtcbiAgYmFja2dyb3VuZC1jb2xvcjogI2ZkZmRmZDtcbiAgYm9yZGVyOiAxcHggc29saWQgdmFyKC0taW9uLWNvbG9yLW1lZGl1bSk7XG4gIGJvcmRlci1yYWRpdXM6IDAuNzVlbTtcbn1cblxuaW9uLWJ1dHRvbiB7XG4gIG1hcmdpbjogMDtcbn1cbiJdfQ== */");

/***/ }),

/***/ 3407:
/*!********************************************************************************************************!*\
  !*** ./src/app/virus-discovery/virus-discovery-hit-details/virus-discovery-hit-details.component.scss ***!
  \********************************************************************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "default": () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (".hsp-data > ion-col > div {\n  padding: 0.5em 0;\n}\n\n.hsp-data .alignments {\n  display: flex;\n  justify-content: center;\n}\n\n.hsp-data .bio-base {\n  display: block;\n  width: 1.5em;\n}\n\n.hsp-data .bio-base.base-A {\n  background-color: var(--ion-color-tertiary-tint);\n}\n\n.hsp-data .bio-base.base-T {\n  background-color: var(--ion-color-secondary-tint);\n}\n\n.hsp-data .bio-base.base-G {\n  background-color: var(--ion-color-warning-tint);\n}\n\n.hsp-data .bio-base.base-C {\n  background-color: var(--ion-color-danger-tint);\n}\n/*# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbInZpcnVzLWRpc2NvdmVyeS1oaXQtZGV0YWlscy5jb21wb25lbnQuc2NzcyJdLCJuYW1lcyI6W10sIm1hcHBpbmdzIjoiQUFBQTtFQUNFLGdCQUFBO0FBQ0Y7O0FBRUE7RUFDRSxhQUFBO0VBQ0EsdUJBQUE7QUFDRjs7QUFDQTtFQUNFLGNBQUE7RUFDQSxZQUFBO0FBRUY7O0FBQ0E7RUFDRSxnREFBQTtBQUVGOztBQUNBO0VBQ0UsaURBQUE7QUFFRjs7QUFDQTtFQUNFLCtDQUFBO0FBRUY7O0FBQ0E7RUFDRSw4Q0FBQTtBQUVGIiwiZmlsZSI6InZpcnVzLWRpc2NvdmVyeS1oaXQtZGV0YWlscy5jb21wb25lbnQuc2NzcyIsInNvdXJjZXNDb250ZW50IjpbIi5oc3AtZGF0YSA+IGlvbi1jb2wgPiBkaXYge1xuICBwYWRkaW5nOiAwLjVlbSAwO1xufVxuXG4uaHNwLWRhdGEgLmFsaWdubWVudHMge1xuICBkaXNwbGF5OiBmbGV4O1xuICBqdXN0aWZ5LWNvbnRlbnQ6IGNlbnRlcjtcbn1cbi5oc3AtZGF0YSAuYmlvLWJhc2Uge1xuICBkaXNwbGF5OiBibG9jaztcbiAgd2lkdGg6IDEuNWVtO1xufVxuXG4uaHNwLWRhdGEgLmJpby1iYXNlLmJhc2UtQSB7XG4gIGJhY2tncm91bmQtY29sb3I6IHZhcigtLWlvbi1jb2xvci10ZXJ0aWFyeS10aW50KTtcbn1cblxuLmhzcC1kYXRhIC5iaW8tYmFzZS5iYXNlLVQge1xuICBiYWNrZ3JvdW5kLWNvbG9yOiB2YXIoLS1pb24tY29sb3Itc2Vjb25kYXJ5LXRpbnQpO1xufVxuXG4uaHNwLWRhdGEgLmJpby1iYXNlLmJhc2UtRyB7XG4gIGJhY2tncm91bmQtY29sb3I6IHZhcigtLWlvbi1jb2xvci13YXJuaW5nLXRpbnQpO1xufVxuXG4uaHNwLWRhdGEgLmJpby1iYXNlLmJhc2UtQyB7XG4gIGJhY2tncm91bmQtY29sb3I6IHZhcigtLWlvbi1jb2xvci1kYW5nZXItdGludCk7XG59XG4iXX0= */");

/***/ }),

/***/ 1129:
/*!******************************************************************************************************!*\
  !*** ./src/app/virus-discovery/virus-discovery-hits-graph/virus-discovery-hits-graph.component.scss ***!
  \******************************************************************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "default": () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = ("h3 {\n  font-size: 1.5em;\n  font-weight: bold;\n}\n/*# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbInZpcnVzLWRpc2NvdmVyeS1oaXRzLWdyYXBoLmNvbXBvbmVudC5zY3NzIl0sIm5hbWVzIjpbXSwibWFwcGluZ3MiOiJBQUFBO0VBQ0UsZ0JBQUE7RUFDQSxpQkFBQTtBQUNGIiwiZmlsZSI6InZpcnVzLWRpc2NvdmVyeS1oaXRzLWdyYXBoLmNvbXBvbmVudC5zY3NzIiwic291cmNlc0NvbnRlbnQiOlsiaDMge1xuICBmb250LXNpemU6IDEuNWVtO1xuICBmb250LXdlaWdodDogYm9sZDtcbn1cbiJdfQ== */");

/***/ }),

/***/ 3697:
/*!******************************************************************************************************!*\
  !*** ./src/app/virus-discovery/virus-discovery-hits-table/virus-discovery-hits-table.component.scss ***!
  \******************************************************************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "default": () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = ("h3 {\n  font-size: 1.5em;\n  font-weight: bold;\n}\n/*# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbInZpcnVzLWRpc2NvdmVyeS1oaXRzLXRhYmxlLmNvbXBvbmVudC5zY3NzIl0sIm5hbWVzIjpbXSwibWFwcGluZ3MiOiJBQUFBO0VBQ0UsZ0JBQUE7RUFDQSxpQkFBQTtBQUNGIiwiZmlsZSI6InZpcnVzLWRpc2NvdmVyeS1oaXRzLXRhYmxlLmNvbXBvbmVudC5zY3NzIiwic291cmNlc0NvbnRlbnQiOlsiaDMge1xuICBmb250LXNpemU6IDEuNWVtO1xuICBmb250LXdlaWdodDogYm9sZDtcbn1cbiJdfQ== */");

/***/ }),

/***/ 8081:
/*!************************************************************************************************!*\
  !*** ./src/app/virus-discovery/virus-discovery-results/virus-discovery-results.component.scss ***!
  \************************************************************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "default": () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = ("ion-item {\n  width: 100%;\n}\n\nion-text {\n  size: 5em;\n}\n\nion-button {\n  min-width: 11em;\n}\n\nh2 {\n  margin: 1em 0 0.5em;\n  font-size: 24px;\n  font-weight: 500;\n}\n\nh3 {\n  font-size: 1.5em;\n  font-weight: bold;\n}\n\np {\n  margin: 1em 0;\n  font-size: 1em;\n}\n/*# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbInZpcnVzLWRpc2NvdmVyeS1yZXN1bHRzLmNvbXBvbmVudC5zY3NzIl0sIm5hbWVzIjpbXSwibWFwcGluZ3MiOiJBQUFBO0VBQ0UsV0FBQTtBQUNGOztBQUVBO0VBQ0UsU0FBQTtBQUNGOztBQUVBO0VBQ0UsZUFBQTtBQUNGOztBQUVBO0VBQ0UsbUJBQUE7RUFDQSxlQUFBO0VBQ0EsZ0JBQUE7QUFDRjs7QUFFQTtFQUNFLGdCQUFBO0VBQ0EsaUJBQUE7QUFDRjs7QUFFQTtFQUNFLGFBQUE7RUFDQSxjQUFBO0FBQ0YiLCJmaWxlIjoidmlydXMtZGlzY292ZXJ5LXJlc3VsdHMuY29tcG9uZW50LnNjc3MiLCJzb3VyY2VzQ29udGVudCI6WyJpb24taXRlbSB7XG4gIHdpZHRoOiAxMDAlO1xufVxuXG5pb24tdGV4dCB7XG4gIHNpemU6IDVlbTtcbn1cblxuaW9uLWJ1dHRvbiB7XG4gIG1pbi13aWR0aDogMTFlbTtcbn1cblxuaDIge1xuICBtYXJnaW46IDFlbSAwIDAuNWVtO1xuICBmb250LXNpemU6IDI0cHg7XG4gIGZvbnQtd2VpZ2h0OiA1MDA7XG59XG5cbmgzIHtcbiAgZm9udC1zaXplOiAxLjVlbTtcbiAgZm9udC13ZWlnaHQ6IGJvbGQ7XG59XG5cbnAge1xuICBtYXJnaW46IDFlbSAwO1xuICBmb250LXNpemU6IDFlbTtcbn1cbiJdfQ== */");

/***/ }),

/***/ 1400:
/*!**************************************************************************************************!*\
  !*** ./src/app/virus-discovery/virus-discovery-trimming/virus-discovery-trimming.component.scss ***!
  \**************************************************************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "default": () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = ("ion-col.vcenter {\n  display: flex;\n  align-items: center;\n}\n\nion-col.section-head {\n  padding: 20px 0 10px;\n}\n\nion-input, ion-select {\n  background-color: #fdfdfd;\n  border: 1px solid var(--ion-color-medium);\n  border-radius: 0.75em;\n}\n\nion-button {\n  margin: 0;\n}\n\niframe {\n  width: 100%;\n  min-height: calc(100vh - 100px);\n  border: 0;\n}\n/*# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbInZpcnVzLWRpc2NvdmVyeS10cmltbWluZy5jb21wb25lbnQuc2NzcyJdLCJuYW1lcyI6W10sIm1hcHBpbmdzIjoiQUFBQTtFQUNFLGFBQUE7RUFDQSxtQkFBQTtBQUNGOztBQUVBO0VBQ0Usb0JBQUE7QUFDRjs7QUFFQTtFQUNFLHlCQUFBO0VBQ0EseUNBQUE7RUFDQSxxQkFBQTtBQUNGOztBQUVBO0VBQ0UsU0FBQTtBQUNGOztBQUVBO0VBQ0UsV0FBQTtFQUNBLCtCQUFBO0VBQ0EsU0FBQTtBQUNGIiwiZmlsZSI6InZpcnVzLWRpc2NvdmVyeS10cmltbWluZy5jb21wb25lbnQuc2NzcyIsInNvdXJjZXNDb250ZW50IjpbImlvbi1jb2wudmNlbnRlciB7XG4gIGRpc3BsYXk6IGZsZXg7XG4gIGFsaWduLWl0ZW1zOiBjZW50ZXI7XG59XG5cbmlvbi1jb2wuc2VjdGlvbi1oZWFkIHtcbiAgcGFkZGluZzogMjBweCAwIDEwcHg7XG59XG5cbmlvbi1pbnB1dCwgaW9uLXNlbGVjdCAge1xuICBiYWNrZ3JvdW5kLWNvbG9yOiAjZmRmZGZkO1xuICBib3JkZXI6IDFweCBzb2xpZCB2YXIoLS1pb24tY29sb3ItbWVkaXVtKTtcbiAgYm9yZGVyLXJhZGl1czogMC43NWVtO1xufVxuXG5pb24tYnV0dG9uIHtcbiAgbWFyZ2luOiAwO1xufVxuXG5pZnJhbWUge1xuICB3aWR0aDogMTAwJTtcbiAgbWluLWhlaWdodDogY2FsYygxMDB2aCAtIDEwMHB4KTtcbiAgYm9yZGVyOiAwO1xufVxuIl19 */");

/***/ }),

/***/ 1106:
/*!**************************************************************************!*\
  !*** ./node_modules/raw-loader/dist/cjs.js!./src/app/app.component.html ***!
  \**************************************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "default": () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = ("<ion-app>\n  <ion-router-outlet id=\"main-container\" animated=\"false\"></ion-router-outlet>\n</ion-app>\n");

/***/ }),

/***/ 5302:
/*!********************************************************************************************!*\
  !*** ./node_modules/raw-loader/dist/cjs.js!./src/app/commons/footer/footer.component.html ***!
  \********************************************************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "default": () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = ("<div class=\"footer\">\n  <ion-grid>\n    <ion-row class=\"ion-align-items-center ion-hide-md-down\">\n      <ion-col>\n        <div class=\"disclaimer\">\n          <p class=\"ion-text-center\">&copy; Agorakis Bompotas, Nikitas Kalogeropoulos, Maria Giachali, Ioly Kotta-Loizou, Christos Makris</p>\n        </div>\n      </ion-col>\n    </ion-row>\n  </ion-grid>\n</div>\n<div *ngIf=\"spacer\" class=\"spacer ion-hide-lg-up\"></div>\n");

/***/ }),

/***/ 4503:
/*!********************************************************************************************!*\
  !*** ./node_modules/raw-loader/dist/cjs.js!./src/app/commons/topbar/topbar.component.html ***!
  \********************************************************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "default": () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = ("<ion-header [translucent]=\"true\">\n  <ion-toolbar>\n    <ion-buttons *ngIf=\"buttons === true\" slot=\"start\">\n      <a href=\"https://www.imslab.gr/\" target=\"_blank\">\n        <ion-icon name=\"information-circle\"></ion-icon>\n      </a>\n    </ion-buttons>\n    <ion-title>\n      <a routerLink=\"/\">\n        <ion-img class=\"large ion-hide-lg-down\" src=\"/assets/discvir.png\"></ion-img>\n      </a>\n      <a routerLink=\"/\">\n        <ion-img class=\"small ion-hide-lg-up\" src=\"/assets/discvir.png\"></ion-img>\n      </a>\n    </ion-title>\n  </ion-toolbar>\n</ion-header>\n");

/***/ }),

/***/ 4217:
/*!*******************************************************************************************!*\
  !*** ./node_modules/raw-loader/dist/cjs.js!./src/app/pages/home/home-page.component.html ***!
  \*******************************************************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "default": () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = ("<app-topbar buttons=\"false\"></app-topbar>\n\n<ion-content [fullscreen]=\"true\">\n  <div class=\"page-container\">\n    <ion-grid>\n      <ion-row>\n        <ion-col size=\"0\" sizeXl=\"1\"></ion-col>\n        <ion-col size=\"12\" sizeXl=\"10\">\n          <app-virus-discovery-form></app-virus-discovery-form>\n        </ion-col>\n        <ion-col size=\"0\" sizeXl=\"1\"></ion-col>\n      </ion-row>\n    </ion-grid>\n\n    <app-footer></app-footer>\n  </div>\n</ion-content>\n");

/***/ }),

/***/ 8069:
/*!*************************************************************************************************!*\
  !*** ./node_modules/raw-loader/dist/cjs.js!./src/app/pages/results/results-page.component.html ***!
  \*************************************************************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "default": () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = ("<app-topbar buttons=\"false\"></app-topbar>\n\n<ion-content [fullscreen]=\"true\">\n  <div class=\"page-container\">\n    <ion-grid>\n      <ion-row>\n        <ion-col size=\"0\" sizeXl=\"1\"></ion-col>\n        <ion-col size=\"12\" sizeXl=\"10\">\n          <app-virus-discovery-results></app-virus-discovery-results>\n        </ion-col>\n        <ion-col size=\"0\" sizeXl=\"1\"></ion-col>\n      </ion-row>\n    </ion-grid>\n\n    <app-footer></app-footer>\n  </div>\n</ion-content>\n");

/***/ }),

/***/ 4210:
/*!***************************************************************************************************!*\
  !*** ./node_modules/raw-loader/dist/cjs.js!./src/app/pages/trimming/trimming-page.component.html ***!
  \***************************************************************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "default": () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = ("<app-topbar buttons=\"false\"></app-topbar>\n\n<ion-content [fullscreen]=\"true\">\n  <div class=\"page-container\">\n    <ion-grid>\n      <ion-row>\n        <ion-col size=\"0\" sizeXl=\"1\"></ion-col>\n        <ion-col size=\"12\" sizeXl=\"10\">\n          <app-virus-discovery-trimming></app-virus-discovery-trimming>\n        </ion-col>\n        <ion-col size=\"0\" sizeXl=\"1\"></ion-col>\n      </ion-row>\n    </ion-grid>\n\n    <app-footer></app-footer>\n  </div>\n</ion-content>\n");

/***/ }),

/***/ 8790:
/*!********************************************************************************************************************************!*\
  !*** ./node_modules/raw-loader/dist/cjs.js!./src/app/virus-discovery/virus-discovery-form/virus-discovery-form.component.html ***!
  \********************************************************************************************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "default": () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = ("<ion-card>\n  <ion-card-header class=\"ion-margin-top ion-text-center\">\n    <ion-card-title>Mycovirus Discovery</ion-card-title>\n  </ion-card-header>\n  <ion-card-content class=\"ion-text-center\">\n    <ion-label>\n      Search genomes for discovering mycoviruses\n    </ion-label>\n    <br/><br/>\n    <ion-grid>\n      <ion-row>\n        <ion-col></ion-col>\n        <ion-col size=\"12\" sizeLg=\"10\" sizeXl=\"8\">\n          <!--suppress AngularUndefinedBinding -->\n          <form (ngSubmit)=\"search()\">\n            <ion-row>\n              <ion-col class=\"section-head\">\n                <ion-label>\n                  <strong>Provide the genome:</strong>\n                </ion-label>\n              </ion-col>\n            </ion-row>\n            <ion-row class=\"ion-text-start\">\n              <ion-col sizeSm=\"3\" sizeXs=\"4\" class=\"vcenter\">\n                <ion-label>E-mail:</ion-label>\n              </ion-col>\n              <ion-col sizeSm=\"9\" sizeXs=\"8\">\n                <ion-input [(ngModel)]=\"email\" name=\"email\" type=\"email\" placeholder=\"johndoe@email.com\" required=\"true\"></ion-input>\n              </ion-col>\n            </ion-row>\n            <ion-row class=\"ion-text-start\">\n              <ion-col sizeSm=\"3\" sizeXs=\"4\" class=\"vcenter\">\n                <ion-label>Sample name:</ion-label>\n              </ion-col>\n              <ion-col sizeSm=\"9\" sizeXs=\"8\">\n                <ion-input [(ngModel)]=\"sampleName\" name=\"sample_name\" type=\"text\" placeholder=\"Sample 1\" required=\"true\"></ion-input>\n              </ion-col>\n            </ion-row>\n            <ion-row class=\"ion-text-start\">\n              <ion-col sizeSm=\"3\" sizeXs=\"4\" class=\"vcenter\">\n                <ion-label>Sequencing technology:</ion-label>\n              </ion-col>\n              <ion-col sizeSm=\"9\" sizeXs=\"8\">\n                <ion-select [(ngModel)]=\"sequencingTechnology\" name=\"sequencing_technology\"  placeholder=\"Single or Pair end\" required=\"true\">\n                  <ion-select-option value=\"single\">Single end</ion-select-option>\n                  <ion-select-option value=\"paired\">Paired end</ion-select-option>\n                </ion-select>\n              </ion-col>\n            </ion-row>\n            <ion-row class=\"ion-text-start\" *ngIf=\"sequencingTechnology === 'single'\">\n              <ion-col sizeSm=\"3\" sizeXs=\"4\" class=\"vcenter\">\n                <ion-label>Input file:</ion-label>\n              </ion-col>\n              <ion-col sizeSm=\"9\" sizeXs=\"8\">\n                <ion-input name=\"single_file\" type=\"file\" (ionChange)=\"onSingleFileChange($event)\" required=\"true\"></ion-input>\n              </ion-col>\n            </ion-row>\n            <ion-row class=\"ion-text-start\" *ngIf=\"sequencingTechnology === 'paired'\">\n              <ion-col sizeSm=\"3\" sizeXs=\"4\" class=\"vcenter\">\n                <ion-label>Input file (forward):</ion-label>\n              </ion-col>\n              <ion-col sizeSm=\"9\" sizeXs=\"8\">\n                <ion-input name=\"forward_file\" type=\"file\" (ionChange)=\"onForwardFileChange($event)\" required=\"true\"></ion-input>\n              </ion-col>\n            </ion-row>\n            <ion-row class=\"ion-text-start\" *ngIf=\"sequencingTechnology === 'paired'\">\n              <ion-col sizeSm=\"3\" sizeXs=\"4\" class=\"vcenter\">\n                <ion-label>Input file (reverse):</ion-label>\n              </ion-col>\n              <ion-col sizeSm=\"9\" sizeXs=\"8\">\n                <ion-input name=\"reverse_file\" type=\"file\" (ionChange)=\"onReverseFileChange($event)\" required=\"true\"></ion-input>\n              </ion-col>\n            </ion-row>\n            <ion-row class=\"ion-text-start\">\n              <ion-col sizeSm=\"3\" sizeXs=\"4\" class=\"vcenter\">\n                <ion-label>Reference Genome file:</ion-label>\n              </ion-col>\n              <ion-col sizeSm=\"9\" sizeXs=\"8\">\n                <ion-input name=\"reference_genome\" type=\"file\" (ionChange)=\"onGenomeFileChange($event)\" required=\"true\"></ion-input>\n              </ion-col>\n            </ion-row>\n            <br/><br/>\n            <ion-row>\n              <ion-col></ion-col>\n              <ion-col size=\"10\" sizeMd=\"4\" class=\"ion-text-center\">\n                <ion-button type=\"submit\" expand=\"block\" color=\"secondary\">Submit</ion-button>\n              </ion-col>\n              <ion-col></ion-col>\n            </ion-row>\n          </form>\n        </ion-col>\n        <ion-col></ion-col>\n      </ion-row>\n    </ion-grid>\n    <br/><br/>\n  </ion-card-content>\n</ion-card>\n");

/***/ }),

/***/ 8546:
/*!**********************************************************************************************************************************************!*\
  !*** ./node_modules/raw-loader/dist/cjs.js!./src/app/virus-discovery/virus-discovery-hit-details/virus-discovery-hit-details.component.html ***!
  \**********************************************************************************************************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "default": () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = ("<ion-row class=\"hsp-data\">\n  <ion-col size=\"1\" class=\"ion-text-end\">\n    <div *ngFor=\"let h of hsp.query; index as i\">\n      <div><strong>Query Start</strong></div>\n      <div>{{hsp.queryStart[i]}}</div>\n      <div>&nbsp;</div>\n      <div>{{hsp.sbjctStart[i]}}</div>\n      <div><strong>Subject Start</strong></div>\n    </div>\n  </ion-col>\n  <ion-col size=\"10\" class=\"ion-text-center\">\n    <div *ngFor=\"let h of hsp.query; index as i\">\n      <div><strong>Query</strong></div>\n      <div class=\"alignments\">\n        <span *ngFor=\"let b of hsp.query[i].split('');\" class=\"bio-base base-{{b}}\">{{b}}</span>\n      </div>\n      <div class=\"alignments\">\n        <span *ngFor=\"let b of hsp.match[i].split('');\" class=\"bio-base\"><strong>{{b}}</strong></span>\n      </div>\n      <div class=\"alignments\">\n        <span *ngFor=\"let b of hsp.sbjct[i].split('');\" class=\"bio-base base-{{b}}\">{{b}}</span>\n      </div>\n      <div><strong>Subject</strong></div>\n    </div>\n  </ion-col>\n  <ion-col size=\"1\" class=\"ion-text-start\">\n    <div *ngFor=\"let h of hsp.query; index as i\">\n      <div><strong>Query End</strong></div>\n      <div>{{hsp.queryEnd[i]}}</div>\n      <div>&nbsp;</div>\n      <div>{{hsp.sbjctEnd[i]}}</div>\n      <div><strong>Subject End</strong></div>\n    </div>\n  </ion-col>\n</ion-row>\n");

/***/ }),

/***/ 7566:
/*!********************************************************************************************************************************************!*\
  !*** ./node_modules/raw-loader/dist/cjs.js!./src/app/virus-discovery/virus-discovery-hits-graph/virus-discovery-hits-graph.component.html ***!
  \********************************************************************************************************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "default": () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = ("<ion-grid>\n  <ion-row>\n    <ion-col class=\"ion-text-center\">\n      <h3>Graphic Summary</h3>\n    </ion-col>\n  </ion-row>\n  <ion-row>\n    <ion-col>\n      <canvas id=\"hits-graph-{{qid}}\"></canvas>\n    </ion-col>\n  </ion-row>\n</ion-grid>\n");

/***/ }),

/***/ 8723:
/*!********************************************************************************************************************************************!*\
  !*** ./node_modules/raw-loader/dist/cjs.js!./src/app/virus-discovery/virus-discovery-hits-table/virus-discovery-hits-table.component.html ***!
  \********************************************************************************************************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "default": () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = ("<ion-grid>\n  <ion-row>\n    <ion-col class=\"ion-text-center\">\n      <ion-label>\n        <h3>Hits</h3>\n      </ion-label>\n    </ion-col>\n  </ion-row>\n  <ion-row>\n    <ion-col>\n      <table id=\"alignments-{{qid}}\" class=\"display\" style=\"width:100%\">\n        <thead>\n        <tr>\n          <th>Title</th>\n          <th>Length</th>\n          <th>Score</th>\n          <th>Bits</th>\n          <th>Expect</th>\n          <th>Alignment Length</th>\n          <th>Gaps\n          <th>Identities</th>\n          <th>Positives</th>\n          <th>Strand</th>\n          <th>Details</th>\n        </tr>\n        </thead>\n        <tbody>\n        <tr *ngFor=\"let a of alignmentsData; index as i\">\n          <td>{{a.title}}</td>\n          <td>{{a.length}}</td>\n          <td>{{a.score}}</td>\n          <td>{{a.bits}}</td>\n          <td>{{a.expect}}</td>\n          <td>{{a.alignLength}}</td>\n          <td>{{a.gaps}}</td>\n          <td>{{a.identities}}</td>\n          <td>{{a.positives}}</td>\n          <td>{{a.strand.join('/')}}</td>\n          <td>\n            <ion-button (click)=\"showHSPDetails(i)\" color=\"secondary\">Show</ion-button>\n          </td>\n        </tr>\n        </tbody>\n      </table>\n    </ion-col>\n  </ion-row>\n  <ion-row *ngIf=\"selectedHSP !== null\">\n    <ion-col class=\"ion-text-center\">\n      <ion-label>\n        <h3>High Scoring Pairs</h3>\n        <p><strong>{{selectedHSP.title}}</strong></p>\n      </ion-label>\n    </ion-col>\n  </ion-row>\n  <app-virus-discovery-hit-details *ngIf=\"selectedHSP !== null\" [hsp]=\"selectedHSP\"></app-virus-discovery-hit-details>\n</ion-grid>\n");

/***/ }),

/***/ 9064:
/*!**************************************************************************************************************************************!*\
  !*** ./node_modules/raw-loader/dist/cjs.js!./src/app/virus-discovery/virus-discovery-results/virus-discovery-results.component.html ***!
  \**************************************************************************************************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "default": () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = ("<ion-card>\n  <ion-card-header class=\"ion-margin-top ion-text-center\">\n    <ion-card-title>Analysis Results</ion-card-title>\n  </ion-card-header>\n  <ion-card-content class=\"ion-text-center\">\n    <br/>\n    <ion-button color=\"secondary\" [href]=\"downloadURL\">Download all</ion-button>\n    <br/><br/>\n    <ion-grid>\n      <ion-row *ngFor=\"let row of results; index as i\">\n        <ion-col>\n          <ion-accordion-group>\n            <ion-accordion value=\"first\">\n              <ion-item slot=\"header\" color=\"primary\">\n                <ion-label>Query: {{row.queryName}}</ion-label>\n              </ion-item>\n              <div class=\"ion-padding ion-text-start\" slot=\"content\">\n                <ion-grid>\n                  <ion-row>\n                    <ion-col class=\"ion-text-center\">\n                      <ion-label>\n                        <h3>Query Info</h3>\n                      </ion-label>\n                      <ion-label>\n                        <strong>Name:</strong> {{row.queryName}}\n                        &nbsp;&nbsp;\n                        <strong>Letters:</strong> {{row.queryLetters}}\n                      </ion-label>\n                    </ion-col>\n                  </ion-row>\n                </ion-grid>\n                <br/>\n                <app-virus-discovery-hits-graph [qid]=\"i\" [qData]=\"row\"></app-virus-discovery-hits-graph>\n                <br/>\n                <app-virus-discovery-hits-table [qid]=\"i\" [qData]=\"row\"></app-virus-discovery-hits-table>\n              </div>\n            </ion-accordion>\n          </ion-accordion-group>\n        </ion-col>\n      </ion-row>\n    </ion-grid>\n    <br/><br/>\n  </ion-card-content>\n</ion-card>\n");

/***/ }),

/***/ 3352:
/*!****************************************************************************************************************************************!*\
  !*** ./node_modules/raw-loader/dist/cjs.js!./src/app/virus-discovery/virus-discovery-trimming/virus-discovery-trimming.component.html ***!
  \****************************************************************************************************************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "default": () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = ("<ion-card>\n  <ion-card-header class=\"ion-margin-top ion-text-center\">\n    <ion-card-title>Mycovirus Discovery</ion-card-title>\n  </ion-card-header>\n  <ion-card-content class=\"ion-text-center\">\n    <ion-label>\n      Search genomes for discovering mycoviruses\n    </ion-label>\n    <br/><br/>\n    <!--suppress AngularUndefinedBinding -->\n    <form (ngSubmit)=\"trim()\">\n      <ion-grid>\n        <ion-row *ngFor=\"let report of analysisReports; index as i\">\n          <ion-col>\n            <ion-accordion-group>\n              <ion-accordion value=\"first\">\n                <ion-item slot=\"header\" color=\"primary\">\n                  <ion-label>Report file #{{i + 1}}</ion-label>\n                </ion-item>\n                <div class=\"ion-padding ion-text-start\" slot=\"content\">\n                  <iframe [srcdoc]=\"report | safeIFrame\"></iframe>\n                </div>\n              </ion-accordion>\n            </ion-accordion-group>\n          </ion-col>\n        </ion-row>\n        <br/><br/>\n        <ion-row>\n          <ion-col></ion-col>\n          <ion-col size=\"12\" sizeLg=\"10\" sizeXl=\"8\">\n            <ion-row class=\"ion-text-start\">\n              <ion-col sizeSm=\"3\" sizeXs=\"4\" class=\"vcenter\">\n                <ion-label>Trimmomatic adapter:</ion-label>\n              </ion-col>\n              <ion-col sizeSm=\"9\" sizeXs=\"8\">\n                <ion-input name=\"adapter\" type=\"file\" (ionChange)=\"onAdapterChange($event)\" required=\"true\"></ion-input>\n              </ion-col>\n            </ion-row>\n            <ion-row class=\"ion-text-start\">\n              <ion-col sizeSm=\"3\" sizeXs=\"4\" class=\"vcenter\">\n                <ion-label>Trimmomatic sliding window:</ion-label>\n              </ion-col>\n              <ion-col sizeSm=\"9\" sizeXs=\"8\">\n                <ion-input [(ngModel)]=\"slidingWindow\" name=\"sliding_window\" type=\"text\" required=\"true\"></ion-input>\n              </ion-col>\n            </ion-row>\n            <ion-row class=\"ion-text-start\">\n              <ion-col sizeSm=\"3\" sizeXs=\"4\" class=\"vcenter\">\n                <ion-label>Trimmomatic minimum length:</ion-label>\n              </ion-col>\n              <ion-col sizeSm=\"9\" sizeXs=\"8\">\n                <ion-input [(ngModel)]=\"minLength\" name=\"min_length\" type=\"text\" required=\"true\"></ion-input>\n              </ion-col>\n            </ion-row>\n            <br/><br/>\n            <ion-row>\n              <ion-col></ion-col>\n              <ion-col size=\"10\" sizeMd=\"4\" class=\"ion-text-center\">\n                <ion-button type=\"submit\" expand=\"block\" color=\"primary\">Trim</ion-button>\n              </ion-col>\n              <ion-col size=\"1\"></ion-col>\n              <ion-col size=\"1\"></ion-col>\n              <ion-col size=\"10\" sizeMd=\"4\" class=\"ion-text-center\">\n                <ion-button (click)=\"proceed()\" expand=\"block\" color=\"secondary\">Proceed</ion-button>\n              </ion-col>\n              <ion-col></ion-col>\n            </ion-row>\n          </ion-col>\n          <ion-col></ion-col>\n        </ion-row>\n      </ion-grid>\n    </form>\n    <br/><br/>\n  </ion-card-content>\n</ion-card>\n");

/***/ })

},
/******/ __webpack_require__ => { // webpackRuntimeModules
/******/ "use strict";
/******/ 
/******/ var __webpack_exec__ = (moduleId) => (__webpack_require__(__webpack_require__.s = moduleId))
/******/ __webpack_require__.O(0, ["vendor"], () => (__webpack_exec__(4431)));
/******/ var __webpack_exports__ = __webpack_require__.O();
/******/ }
]);
//# sourceMappingURL=main.js.map