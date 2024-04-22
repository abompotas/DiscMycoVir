import {AfterViewInit, Component, ElementRef, Input, OnInit, ViewChild} from '@angular/core';
import * as $3Dmol from '../../../assets/js/3Dmol-min.js';

@Component({
  selector: 'app-search-binding-pockets-visualization',
  templateUrl: './search-binding-pockets-visualization.component.html',
  styleUrls: ['./search-binding-pockets-visualization.component.scss'],
})
export class SearchBindingPocketsVisualizationComponent implements OnInit, AfterViewInit {

  @ViewChild('containerCentroid') container: ElementRef;
  @Input() id;
  @Input() pdb;
  @Input() proteinChain;
  @Input() ligandChain;

  constructor() {
  }

  ngOnInit() {
  }

  ngAfterViewInit() {
    let viewer = $3Dmol.createViewer(this.container.nativeElement);
    viewer.addModel(this.pdb, 'pdb');
    viewer.setStyle({chain: this.proteinChain}, {cartoon: {color: '#cccccc'}});
    viewer.setStyle({chain: this.ligandChain}, {stick: {color: '#6495ed', radius: 0.33}});
    viewer.zoomTo();
    viewer.render();
    viewer.zoom(1.2, 1000);
  }
}
