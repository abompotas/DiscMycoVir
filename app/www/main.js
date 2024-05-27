(self["webpackChunkVirus_Discovery"] = self["webpackChunkVirus_Discovery"] || []).push([["main"],{

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
/* harmony import */ var tslib__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! tslib */ 4762);
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
/* harmony import */ var tslib__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! tslib */ 4762);
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
/* harmony import */ var tslib__WEBPACK_IMPORTED_MODULE_10__ = __webpack_require__(/*! tslib */ 4762);
/* harmony import */ var _angular_core__WEBPACK_IMPORTED_MODULE_11__ = __webpack_require__(/*! @angular/core */ 7716);
/* harmony import */ var _angular_platform_browser__WEBPACK_IMPORTED_MODULE_12__ = __webpack_require__(/*! @angular/platform-browser */ 9075);
/* harmony import */ var _angular_router__WEBPACK_IMPORTED_MODULE_16__ = __webpack_require__(/*! @angular/router */ 9895);
/* harmony import */ var _angular_common_http__WEBPACK_IMPORTED_MODULE_14__ = __webpack_require__(/*! @angular/common/http */ 1841);
/* harmony import */ var _angular_forms__WEBPACK_IMPORTED_MODULE_15__ = __webpack_require__(/*! @angular/forms */ 3679);
/* harmony import */ var _ionic_angular__WEBPACK_IMPORTED_MODULE_13__ = __webpack_require__(/*! @ionic/angular */ 9122);
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

















let AppModule = class AppModule {
};
AppModule = (0,tslib__WEBPACK_IMPORTED_MODULE_10__.__decorate)([
    (0,_angular_core__WEBPACK_IMPORTED_MODULE_11__.NgModule)({
        declarations: [
            _app_component__WEBPACK_IMPORTED_MODULE_0__.AppComponent,
            _pages_home_home_page_component__WEBPACK_IMPORTED_MODULE_2__.HomePageComponent,
            _pages_trimming_trimming_page_component__WEBPACK_IMPORTED_MODULE_3__.TrimmingPageComponent,
            _pages_results_results_page_component__WEBPACK_IMPORTED_MODULE_4__.ResultsPageComponent,
            _commons_footer_footer_component__WEBPACK_IMPORTED_MODULE_5__.FooterComponent,
            _commons_topbar_topbar_component__WEBPACK_IMPORTED_MODULE_6__.TopbarComponent,
            _virus_discovery_virus_discovery_form_virus_discovery_form_component__WEBPACK_IMPORTED_MODULE_7__.VirusDiscoveryFormComponent,
            _virus_discovery_virus_discovery_trimming_virus_discovery_trimming_component__WEBPACK_IMPORTED_MODULE_8__.VirusDiscoveryTrimmingComponent,
            _virus_discovery_virus_discovery_results_virus_discovery_results_component__WEBPACK_IMPORTED_MODULE_9__.VirusDiscoveryResultsComponent
        ],
        entryComponents: [],
        imports: [
            _angular_platform_browser__WEBPACK_IMPORTED_MODULE_12__.BrowserModule,
            _ionic_angular__WEBPACK_IMPORTED_MODULE_13__.IonicModule.forRoot(),
            _app_routing_module__WEBPACK_IMPORTED_MODULE_1__.AppRoutingModule,
            _angular_common_http__WEBPACK_IMPORTED_MODULE_14__.HttpClientModule,
            _angular_forms__WEBPACK_IMPORTED_MODULE_15__.FormsModule,
            _angular_forms__WEBPACK_IMPORTED_MODULE_15__.ReactiveFormsModule
        ],
        providers: [{ provide: _angular_router__WEBPACK_IMPORTED_MODULE_16__.RouteReuseStrategy, useClass: _ionic_angular__WEBPACK_IMPORTED_MODULE_13__.IonicRouteStrategy }],
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
/* harmony import */ var tslib__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! tslib */ 4762);
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
/* harmony import */ var tslib__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! tslib */ 4762);
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
/* harmony import */ var tslib__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! tslib */ 4762);
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
/* harmony import */ var tslib__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! tslib */ 4762);
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
/* harmony import */ var tslib__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! tslib */ 4762);
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
/* harmony import */ var tslib__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! tslib */ 4762);
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
                this.http.post(_environments_environment__WEBPACK_IMPORTED_MODULE_2__.environment.discvirAPI + '/virus-discovery', formData, { responseType: 'json' }).subscribe(x => this.searchResponse(x), e => this.searchError(e.error), () => {
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
            const matches = this.email.match(/^(([^<>()[\]\\.,;:\s@\"]+(\.[^<>()[\]\\.,;:\s@\"]+)*)|(\".+\"))@((\[[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}\])|(([a-zA-Z\-0-9]+\.)+[a-zA-Z]{2,}))$/);
            if (!matches) {
                validEmail = false;
            }
        }
        if (!validEmail) {
            this.searchError({ error: 'Please fill in a valid email address.' }).then(null);
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
            this.searchError({ error: 'Please fill in a name for your experiment (only alphanumeric characters are permitted).' }).then(null);
            return false;
        }
        return this.validateInputFiles();
    }
    validateInputFiles() {
        if (this.sequencingTechnology === null) {
            this.searchError({ error: 'Please select the sequencing technology type.' }).then(null);
            return false;
        }
        else {
            if (this.sequencingTechnology === 'single') {
                if (this.singleFile === null) {
                    this.searchError({ error: 'Please select the input file to search in.' }).then(null);
                    return false;
                }
            }
            else if (this.sequencingTechnology === 'paired') {
                if (this.forwardFile === null) {
                    this.searchError({ error: 'Please select the forward read input file to search in.' }).then(null);
                    return false;
                }
                if (this.reverseFile === null) {
                    this.searchError({ error: 'Please select the reverse read input file to search in.' }).then(null);
                    return false;
                }
            }
            else {
                this.searchError({ error: 'Unknown sequencing technology type.' }).then(null);
                return false;
            }
        }
        if (this.referenceGenome === null) {
            this.searchError({ error: 'Please select a file for the reference genome.' }).then(null);
            return false;
        }
        return true;
    }
    searchResponse(response) {
        if (response.status === 'success') {
            this.alertSuccess().then(null);
        }
        else {
            this.searchError(response).then(null);
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
    searchError(resp) {
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
/* harmony import */ var tslib__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! tslib */ 4762);
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
        this.pocketomeURL = '';
        this.outputs = null;
        this.resultsTable = [];
        this.jobId = 0;
        this.hash = '';
        this.route.params.subscribe(params => {
            if (params.hasOwnProperty('job')) {
                if (params.hasOwnProperty('hash')) {
                    this.jobId = params.job;
                    this.hash = params.hash;
                }
            }
        }, () => router.navigate(['/']), () => {
            if ((this.jobId) === 0 || (this.hash === '')) {
                router.navigate(['/']);
            }
        });
    }
    ngOnInit() {
        this.pocketomeURL = _environments_environment__WEBPACK_IMPORTED_MODULE_2__.environment.discvirAPI + '/virus_discovery/' + this.jobId + '/' + this.hash;
        this.getResults();
    }
    getResults() {
        this.loading().then(() => {
            this.http.get(this.pocketomeURL, { responseType: 'json' }).subscribe(x => this.parseResults(x), e => this.resultsError(e.error), () => this.loadingController.dismiss().then(null));
        });
    }
    parseResults(resp) {
        this.resultsTable = [['PDBid', 'Matching results', 'Query residues', 'Hit residues', 'RMSD', 'SASA', 'Docking']];
        const colNames = ['PDBid', 'matchSize', 'query', 'hit', 'RMSD', 'SASA', 'Docking'];
        const cols = [];
        for (let r = 0; r < resp.results.results.length; r++) {
            let row = resp.results.results[r].replace('\n', '');
            row = row.split('\t');
            if (!r) {
                for (let c = 0; c < row.length; c++) {
                    if (colNames.includes(row[c])) {
                        cols.push(c);
                    }
                }
            }
            else {
                const resultsRow = [];
                for (const c of cols) {
                    resultsRow.push(row[c]);
                }
                this.resultsTable.push(resultsRow);
            }
        }
        this.outputs = resp.results;
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
/* harmony export */   "VirusDiscoveryTrimmingComponent": () => (/* binding */ VirusDiscoveryTrimmingComponent)
/* harmony export */ });
/* harmony import */ var tslib__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! tslib */ 4762);
/* harmony import */ var _raw_loader_virus_discovery_trimming_component_html__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! !raw-loader!./virus-discovery-trimming.component.html */ 3352);
/* harmony import */ var _virus_discovery_trimming_component_scss__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./virus-discovery-trimming.component.scss */ 1400);
/* harmony import */ var _angular_core__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(/*! @angular/core */ 7716);
/* harmony import */ var _angular_common_http__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! @angular/common/http */ 1841);
/* harmony import */ var _angular_router__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! @angular/router */ 9895);
/* harmony import */ var _ionic_angular__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(/*! @ionic/angular */ 9122);
/* harmony import */ var _environments_environment__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ../../../environments/environment */ 2340);








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
        this.analysisResults = [];
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
            this.router.navigate(['/']);
        }
        else {
            this.loading().then(() => {
                this.http.get(_environments_environment__WEBPACK_IMPORTED_MODULE_2__.environment.discvirAPI + '/virus-discovery/analysis/' + this.jobId + '/' + this.hash, { responseType: 'json' }).subscribe(x => {
                    this.analysisResults = [...x.results];
                }, e => this.resultsError(e.error), () => this.loadingController.dismiss().then(null));
            });
        }
    }
    initForm() {
        this.sequencingTechnology = null;
        this.slidingWindow = null;
        this.minLength = null;
        this.adapter = null;
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
            this.http.post(_environments_environment__WEBPACK_IMPORTED_MODULE_2__.environment.discvirAPI + '/virus-discovery/trimming', formData, { responseType: 'json' }).subscribe(x => this.trimResponse(x), e => this.resultsError(e.error), () => {
                this.initForm();
                this.loadingController.dismiss().then(null);
            });
        });
    }
    trimResponse(response) {
        if (response.status === 'success') {
            this.alertSuccess().then(null);
        }
        else {
            this.resultsError(response).then(null);
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
    resultsError(resp) {
        return (0,tslib__WEBPACK_IMPORTED_MODULE_3__.__awaiter)(this, void 0, void 0, function* () {
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
VirusDiscoveryTrimmingComponent.ctorParameters = () => [
    { type: _angular_common_http__WEBPACK_IMPORTED_MODULE_4__.HttpClient },
    { type: _angular_router__WEBPACK_IMPORTED_MODULE_5__.ActivatedRoute },
    { type: _angular_router__WEBPACK_IMPORTED_MODULE_5__.Router },
    { type: _ionic_angular__WEBPACK_IMPORTED_MODULE_6__.AlertController },
    { type: _ionic_angular__WEBPACK_IMPORTED_MODULE_6__.LoadingController }
];
VirusDiscoveryTrimmingComponent = (0,tslib__WEBPACK_IMPORTED_MODULE_3__.__decorate)([
    (0,_angular_core__WEBPACK_IMPORTED_MODULE_7__.Component)({
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
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = ("ion-item {\n  width: 100%;\n}\n\nion-item .error-ad {\n  color: var(--ion-color-danger);\n}\n\nion-text {\n  size: 5em;\n}\n\nion-icon {\n  margin-right: 0.5em;\n}\n\nion-grid.seven-cols {\n  --ion-grid-columns: 7;\n}\n\nion-col {\n  text-align: center;\n}\n\nion-col.bg-white {\n  background-color: #ffffff;\n}\n\nion-button {\n  min-width: 11em;\n}\n\nh2 {\n  margin: 1em 0 0.5em;\n  font-size: 24px;\n  font-weight: 500;\n}\n\np {\n  margin: 1em 0;\n  font-size: 1em;\n}\n\nimg {\n  position: relative;\n  top: 50%;\n  transform: translateY(-50%);\n}\n/*# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbInZpcnVzLWRpc2NvdmVyeS1yZXN1bHRzLmNvbXBvbmVudC5zY3NzIl0sIm5hbWVzIjpbXSwibWFwcGluZ3MiOiJBQUFBO0VBQ0UsV0FBQTtBQUNGOztBQUVBO0VBQ0UsOEJBQUE7QUFDRjs7QUFFQTtFQUNFLFNBQUE7QUFDRjs7QUFFQTtFQUNFLG1CQUFBO0FBQ0Y7O0FBRUE7RUFDRSxxQkFBQTtBQUNGOztBQUVBO0VBQ0Usa0JBQUE7QUFDRjs7QUFFQTtFQUNFLHlCQUFBO0FBQ0Y7O0FBRUE7RUFDRSxlQUFBO0FBQ0Y7O0FBRUE7RUFDRSxtQkFBQTtFQUNBLGVBQUE7RUFDQSxnQkFBQTtBQUNGOztBQUVBO0VBQ0UsYUFBQTtFQUNBLGNBQUE7QUFDRjs7QUFFQTtFQUNFLGtCQUFBO0VBQ0EsUUFBQTtFQUNBLDJCQUFBO0FBQ0YiLCJmaWxlIjoidmlydXMtZGlzY292ZXJ5LXJlc3VsdHMuY29tcG9uZW50LnNjc3MiLCJzb3VyY2VzQ29udGVudCI6WyJpb24taXRlbSB7XG4gIHdpZHRoOiAxMDAlO1xufVxuXG5pb24taXRlbSAuZXJyb3ItYWQge1xuICBjb2xvcjogdmFyKC0taW9uLWNvbG9yLWRhbmdlcik7XG59XG5cbmlvbi10ZXh0IHtcbiAgc2l6ZTogNWVtO1xufVxuXG5pb24taWNvbiB7XG4gIG1hcmdpbi1yaWdodDogMC41ZW07XG59XG5cbmlvbi1ncmlkLnNldmVuLWNvbHMge1xuICAtLWlvbi1ncmlkLWNvbHVtbnM6IDc7XG59XG5cbmlvbi1jb2wge1xuICB0ZXh0LWFsaWduOiBjZW50ZXI7XG59XG5cbmlvbi1jb2wuYmctd2hpdGUge1xuICBiYWNrZ3JvdW5kLWNvbG9yOiAjZmZmZmZmO1xufVxuXG5pb24tYnV0dG9uIHtcbiAgbWluLXdpZHRoOiAxMWVtO1xufVxuXG5oMiB7XG4gIG1hcmdpbjogMWVtIDAgMC41ZW07XG4gIGZvbnQtc2l6ZTogMjRweDtcbiAgZm9udC13ZWlnaHQ6IDUwMDtcbn1cblxucCB7XG4gIG1hcmdpbjogMWVtIDA7XG4gIGZvbnQtc2l6ZTogMWVtO1xufVxuXG5pbWcge1xuICBwb3NpdGlvbjogcmVsYXRpdmU7XG4gIHRvcDogNTAlO1xuICB0cmFuc2Zvcm06IHRyYW5zbGF0ZVkoLTUwJSk7XG59XG4iXX0= */");

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
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = ("ion-col.vcenter {\n  display: flex;\n  align-items: center;\n}\n\nion-col.section-head {\n  padding: 20px 0 10px;\n}\n\nion-input, ion-select {\n  background-color: #fdfdfd;\n  border: 1px solid var(--ion-color-medium);\n  border-radius: 0.75em;\n}\n\nion-button {\n  margin: 0;\n}\n/*# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbInZpcnVzLWRpc2NvdmVyeS10cmltbWluZy5jb21wb25lbnQuc2NzcyJdLCJuYW1lcyI6W10sIm1hcHBpbmdzIjoiQUFBQTtFQUNFLGFBQUE7RUFDQSxtQkFBQTtBQUNGOztBQUVBO0VBQ0Usb0JBQUE7QUFDRjs7QUFFQTtFQUNFLHlCQUFBO0VBQ0EseUNBQUE7RUFDQSxxQkFBQTtBQUNGOztBQUVBO0VBQ0UsU0FBQTtBQUNGIiwiZmlsZSI6InZpcnVzLWRpc2NvdmVyeS10cmltbWluZy5jb21wb25lbnQuc2NzcyIsInNvdXJjZXNDb250ZW50IjpbImlvbi1jb2wudmNlbnRlciB7XG4gIGRpc3BsYXk6IGZsZXg7XG4gIGFsaWduLWl0ZW1zOiBjZW50ZXI7XG59XG5cbmlvbi1jb2wuc2VjdGlvbi1oZWFkIHtcbiAgcGFkZGluZzogMjBweCAwIDEwcHg7XG59XG5cbmlvbi1pbnB1dCwgaW9uLXNlbGVjdCAge1xuICBiYWNrZ3JvdW5kLWNvbG9yOiAjZmRmZGZkO1xuICBib3JkZXI6IDFweCBzb2xpZCB2YXIoLS1pb24tY29sb3ItbWVkaXVtKTtcbiAgYm9yZGVyLXJhZGl1czogMC43NWVtO1xufVxuXG5pb24tYnV0dG9uIHtcbiAgbWFyZ2luOiAwO1xufVxuIl19 */");

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
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = ("<ion-card>\n  <ion-card-header class=\"ion-margin-top ion-text-center\">\n    <ion-card-title>Pocketome Analysis</ion-card-title>\n  </ion-card-header>\n  <ion-card-content class=\"ion-text-center\">\n    <ion-list>\n      <ion-grid>\n        <ion-row>\n          <ion-item class=\"ion-no-padding\">\n            <ion-col><h2>Results</h2></ion-col>\n          </ion-item>\n        </ion-row>\n\n        <div *ngIf=\"outputs !== null\">\n        </div>\n\n        <br/>\n\n      </ion-grid>\n    </ion-list>\n    <br/>\n  </ion-card-content>\n</ion-card>\n");

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
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = ("<ion-card>\n  <ion-card-header class=\"ion-margin-top ion-text-center\">\n    <ion-card-title>Mycovirus Discovery</ion-card-title>\n  </ion-card-header>\n  <ion-card-content class=\"ion-text-center\">\n    <ion-label>\n      Search genomes for discovering mycoviruses\n    </ion-label>\n    <br/><br/>\n    <ion-grid>\n      <ion-row>\n        <ion-col></ion-col>\n        <ion-col size=\"12\" sizeLg=\"10\" sizeXl=\"8\">\n          <!--suppress AngularUndefinedBinding -->\n          <form (ngSubmit)=\"trim()\">\n            <ion-accordion-group *ngFor=\"let report of analysisResults; index as i\">\n              <ion-accordion value=\"first\">\n                <ion-item slot=\"header\" color=\"primary\">\n                  <ion-label>Report file #{{i+1}}</ion-label>\n                </ion-item>\n                <div class=\"ion-padding ion-text-start\" slot=\"content\">\n                  <iframe srcdoc=\"{{report}}\"></iframe>\n                </div>\n              </ion-accordion>\n            </ion-accordion-group>\n            <br/><br/>\n            <ion-row class=\"ion-text-start\">\n              <ion-col sizeSm=\"3\" sizeXs=\"4\" class=\"vcenter\">\n                <ion-label>Trimmomatic adapter:</ion-label>\n              </ion-col>\n              <ion-col sizeSm=\"9\" sizeXs=\"8\">\n                <ion-input name=\"adapter\" type=\"file\" (ionChange)=\"onAdapterChange($event)\" required=\"true\"></ion-input>\n              </ion-col>\n            </ion-row>\n            <ion-row class=\"ion-text-start\">\n              <ion-col sizeSm=\"3\" sizeXs=\"4\" class=\"vcenter\">\n                <ion-label>Trimmomatic sliding window:</ion-label>\n              </ion-col>\n              <ion-col sizeSm=\"9\" sizeXs=\"8\">\n                <ion-input [(ngModel)]=\"slidingWindow\" name=\"sliding_window\" type=\"text\" required=\"true\"></ion-input>\n              </ion-col>\n            </ion-row>\n            <ion-row class=\"ion-text-start\">\n              <ion-col sizeSm=\"3\" sizeXs=\"4\" class=\"vcenter\">\n                <ion-label>Trimmomatic minimum length:</ion-label>\n              </ion-col>\n              <ion-col sizeSm=\"9\" sizeXs=\"8\">\n                <ion-input [(ngModel)]=\"minLength\" name=\"min_length\" type=\"text\" required=\"true\"></ion-input>\n              </ion-col>\n            </ion-row>\n            <br/><br/>\n            <ion-row>\n              <ion-col></ion-col>\n              <ion-col size=\"10\" sizeMd=\"4\" class=\"ion-text-center\">\n                <ion-button type=\"submit\" expand=\"block\" color=\"primary\">Trim</ion-button>\n              </ion-col>\n              <ion-col size=\"1\"></ion-col>\n              <ion-col size=\"1\"></ion-col>\n              <ion-col size=\"10\" sizeMd=\"4\" class=\"ion-text-center\">\n                <ion-button type=\"submit\" expand=\"block\" color=\"secondary\">Submit</ion-button>\n              </ion-col>\n              <ion-col></ion-col>\n            </ion-row>\n          </form>\n        </ion-col>\n        <ion-col></ion-col>\n      </ion-row>\n    </ion-grid>\n    <br/><br/>\n  </ion-card-content>\n</ion-card>\n");

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