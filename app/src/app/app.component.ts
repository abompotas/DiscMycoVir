import {Component} from '@angular/core';
import {AlertController, Platform} from '@ionic/angular';

@Component({
  selector: 'app-root',
  templateUrl: 'app.component.html',
  styleUrls: ['app.component.scss'],
})
export class AppComponent {
  constructor(private platform: Platform, private alertController: AlertController) {
    this.platform.backButton.subscribe(async () => {
      const alert = await this.alertController.create({
        header: 'Exit the Mycovirus Discovery App?',
        message: 'You are going to exit the Mycovirus Discovery App',
        buttons: [
          {text: 'Cancel', role: 'cancel'},
          {
            text: 'OK',
            role: 'confirm',
            handler: () => {
              if(navigator.hasOwnProperty('app')) {
                // @ts-ignore
                navigator.app.exitApp();
              }
              else if(navigator.hasOwnProperty('device')) {
                // @ts-ignore
                navigator.device.exitApp();
              }
            },
          }]
      });
      await alert.present();
    });
  }
}
