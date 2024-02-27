clc
close all
clear all

%% Ucitavanje signala - vremenski i frekvencijski domen

load eegdata.mat

% treći EEG kanal prvog subjekta tokom zadatka relaksacije - P3 kanal
eeg = data{1}{4}(3,:);
fs = 250;
N = length(eeg);
t = 0:1/fs:(N-1)/fs;

figure
    plot(t, eeg, 'Color', [0 0 0]);
      title('Treći EEG kanal prvog subjekta tokom zadatka relaksacije - P3 kanal');
      %xlim([0, 5]);
      xlabel('Vreme [s]');
      ylabel('Amplituda [\muV]' ); 
      grid on;
      
eegfreq = fftshift(abs(fft(eeg)));
f = (fs/N)*((-N/2+1):(N/2));      

figure
     plot(f, eegfreq,'k');
       %xlim([0, 125]);
       title('Furijeova transformacija signala');
       xlabel('Frekvencija [Hz]');
       ylabel('Magnituda [a.u.]');
       grid on;

       
%% Eleminisanje suma napajanja na 60Hz

[b, a] = butter(3, [59.5/(fs/2) 60.5/(fs/2)], 'stop');

eeg_fil = filter(b,a,eeg);
eegfreq_fil = fftshift(abs(fft(eeg_fil)));
figure
    plot(t, eeg_fil, 'Color', [0 0 0]);
      title('Signal posle odstranjivanja šuma napajanja');
      %xlim([0, 5]);
      xlabel('Vreme [s]');
      ylabel('Amplituda [\muV]' ); 
      grid on;
figure
     plot(f, eegfreq_fil,'k');
       %xlim([0, 125]);
       title('Furijeova transformacija signala posle eleminisanja šuma napajanja');
       xlabel('Frekvencija [Hz]');
       ylabel('Magnituda [a.u.]');
       grid on;

%% Dodavanje artifakta pokreta

f0 = 0.5;
artifact =5* sin(2*pi*f0*t);

eeg_noise =artifact + eeg_fil;

figure
  subplot(2,1,1)
    plot(t, eeg_fil,'Color', [0 0 0]);
      %xlim([0, 10]);
      %ylim([-30 30]);
      xlabel('Vreme [s]');
      ylabel('Amplituda [\muV]' ); 
      title('Vremenski prikaz signala pre dodavanja suma pokreta');
      grid on;
  subplot(2,1,2)
    plot(t, eeg_noise,'Color', [0 0 0]);
      %xlim([0, 10]);
      %ylim([-30 30]);
      xlabel('Vreme [s]');
      ylabel('Amplituda [\muV]' ); 
      title('Vremenski prikaz signala posle dodavanja suma pokreta');
      grid on;
      
      
eegfreq_noise = fftshift(abs(fft(eeg_noise)));

figure
  subplot(2,1,1)
    plot(f, eegfreq_fil,'Color', [0 0 0]);
      xlim([f(1), f(end)]);
      %ylim([0, 4500]);
      xlabel('Frekvencija [Hz]');
      ylabel('Magnituda [a.u.]');
      title('Spektar signala pre dodavanja suma pokreta');
      grid on;
  subplot(2,1,2)
    plot(f, eegfreq_noise,'Color', [0 0 0]);
      xlim([f(1), f(end)]);
      %ylim([0, 4500]);
      xlabel('Frekvencija [Hz]');
      ylabel('Magnituda [a.u.]');
      title('Spektar signala posle dodavanja suma pokreta');
      grid on;

%% Polinomijalno fitovanje

indices = floor(linspace(1, N, 20));
t_dec = t(indices);
eeg_dec = eeg_noise(indices);

eeg_dec = [0 eeg_dec 0];


pp = spline(t_dec, eeg_dec);


yy = ppval(pp, t);

yy_freq = fftshift(abs(fft(yy)));

figure
 
    plot(f, yy_freq,'Color', [0 0 0]);
      xlim([f(1), f(end)]);
      %ylim([0, 4500]);
      xlabel('Frekvencija [Hz]');
      ylabel('Magnituda [a.u.]');
      title('Spektar signala dobijen polinomijalnm fitovanjem');
      grid on;

figure
    plot(t,eeg_noise, 'Color', [0 0 0])
    hold on;
    plot(t, yy, 'r--')
      xlabel('Vreme [s]');
      ylabel('Amplituda [\muV]' ); 
      title('Polinomijalno fitovanje');
      grid on;

eeg_filt = eeg_noise-yy;
% do the baseline removal:subtract the spline estimate from the raw data
figure
    plot(t, eeg_filt, 'Color', [0 0 0])
 
    
      xlabel('Vreme [s]');
      ylabel('Amplituda [\muV]' ); 
      title('Filtriran signal');
      grid on;

eegfreq_filt = fftshift(abs(fft(eeg_filt)));

figure
  subplot(2,1,1)
    plot(t, eeg_noise,'Color', [0 0 0]);
      %xlim([0, 10]);
      %ylim([-30 30]);
      xlabel('Vreme [s]');
      ylabel('Amplituda [\muV]' ); 
      title('Vremenski prikaz signala pre filtriranja');
      grid on;
  subplot(2,1,2)
    plot(t, eeg_filt,'Color', [0 0 0]);
      %xlim([0, 10]);
      %ylim([-30 30]);
      xlabel('Vreme [s]');
      ylabel('Amplituda [\muV]' ); 
      title('Vremenski prikaz signala posle filtriranja');
      grid on;
figure
  subplot(2,1,1)
    plot(f, eegfreq_noise,'Color', [0 0 0]);
      xlim([f(1), f(end)]);
      %ylim([0, 4500]);
      xlabel('Frekvencija [Hz]');
      ylabel('Magnituda [a.u.]');
      title('Spektar signala pre filtriranja');
      grid on;
  subplot(2,1,2)
    plot(f, eegfreq_filt,'Color', [0 0 0]);
      xlim([f(1), f(end)]);
      %ylim([0, 4500]);
      xlabel('Frekvencija [Hz]');
      ylabel('Magnituda [a.u.]');
      title('Spektar signala posle filtriranja');
      grid on;
      
figure
  subplot(2,1,1)
    plot(t, eeg_fil,'Color', [0 0 0]);
      %xlim([0, 10]);
      %ylim([-30 30]);
      xlabel('Vreme [s]');
      ylabel('Amplituda [\muV]' ); 
      title('Vremenski prikaz signala pre dodavanja suma pokreta');
      grid on;
  subplot(2,1,2)
    plot(t, eeg_filt,'Color', [0 0 0]);
      %xlim([0, 10]);
      %ylim([-30 30]);
      xlabel('Vreme [s]');
      ylabel('Amplituda [\muV]' ); 
      title('Vremenski prikaz signala posle filtriranja');
      grid on;
      

%% Procenti EEG talasa u ukupnoj spektralnoj snazi pre dodavanja suma, posle i nakon filtriranja

power_data_1 = eegfreq_fil.^2;
power_data_2 = eegfreq_noise.^2;
power_data_3 = eegfreq_filt.^2;

delta = 1:20;
theta = 20:70;
alpha = 70:130;
beta = 130:640;

center = 1250;

power_total_1 = sum(power_data_1(610:1890)); % suma odbiraka
power_delta_1 = (sum(power_data_1(center-delta))+sum(power_data_1(center+delta))) / power_total_1;
power_theta_1 = (sum(power_data_1(center-theta))+sum(power_data_1(center+theta))) / power_total_1;
power_alpha_1 = (sum(power_data_1(center-alpha))+sum(power_data_1(center+alpha))) / power_total_1;
power_beta_1 = (sum(power_data_1(center-beta))+sum(power_data_1(center+beta))) / power_total_1;

power_1 = power_delta_1 + power_theta_1 + power_alpha_1 + power_beta_1; % preklapanje podataka

x_1 = {'\alpha'; '\beta'; '\delta'; '\theta'};
y_1 = [power_alpha_1/power_1*100; power_beta_1/power_1*100; power_delta_1/power_1*100; power_theta_1/power_1*100];

power_total_2 = sum(power_data_2(610:1890)); % suma odbiraka
power_delta_2 = (sum(power_data_2(center-delta))+sum(power_data_2(center+delta))) / power_total_2;
power_theta_2 = (sum(power_data_2(center-theta))+sum(power_data_2(center+theta))) / power_total_2;
power_alpha_2 = (sum(power_data_2(center-alpha))+sum(power_data_2(center+alpha))) / power_total_2;
power_beta_2= (sum(power_data_2(center-beta))+sum(power_data_2(center+beta))) / power_total_2;

power_2 = power_delta_2 + power_theta_2 + power_alpha_2 + power_beta_2; % preklapanje podataka

x_2 = {'\alpha'; '\beta'; '\delta'; '\theta'};
y_2 = [power_alpha_2/power_2*100; power_beta_2/power_2*100; power_delta_2/power_2*100; power_theta_2/power_2*100];

power_total_3 = sum(power_data_3(610:1890)); % suma odbiraka
power_delta_3 = (sum(power_data_3(center-delta))+sum(power_data_3(center+delta))) / power_total_3;
power_theta_3 = (sum(power_data_3(center-theta))+sum(power_data_3(center+theta))) / power_total_3;
power_alpha_3 = (sum(power_data_3(center-alpha))+sum(power_data_3(center+alpha))) / power_total_3;
power_beta_3= (sum(power_data_3(center-beta))+sum(power_data_3(center+beta))) / power_total_3;

power_3 = power_delta_3 + power_theta_3 + power_alpha_3 + power_beta_3; % preklapanje podataka

x_3 = {'\alpha'; '\beta'; '\delta'; '\theta'};
y_3 = [power_alpha_3/power_3*100; power_beta_3/power_3*100; power_delta_3/power_3*100; power_theta_3/power_3*100];

figure
 subplot(3,1,1);
   h = barh (y_1);
   set (h, "facecolor", "k");
   title ("Udeo snaga pojedinačnih talasa za originalni signal");
   set(gca, 'YTickLabel', x_1, 'YTick', 1:numel(x_1))
   xlabel('Udeo snage [%]');
   xlim([0 80])
 subplot(3,1,2);
   h = barh (y_2);
   set (h, "facecolor", "k");
   title ("Udeo snaga pojedinačnih talasa za signal sa artefaktom šuma");
   set(gca, 'YTickLabel', x_2, 'YTick', 1:numel(x_2))
   xlabel('Udeo snage [%]');
   xlim([0 80])
 subplot(3,1,3);
   h = barh (y_3);
   set (h, "facecolor", "k");
   title ("Udeo snaga pojedinačnih talasa za filtrirani signal");
   set(gca, 'YTickLabel', x_3, 'YTick', 1:numel(x_3))
   xlabel('Udeo snage [%]');  
   xlim([0 80])
   
   
y_sum_relativno = abs(y_1-y_2)./y_1*100;
y_filtar_relativno = abs(y_1-y_3)./y_1*100;

figure
 subplot(2,1,1);
   h = barh (y_sum_relativno);
   set (h, "facecolor", "k");
   title ("Relativna razlika udela u spektru originalnog signala i onog posle dodatog šuma");
   set(gca, 'YTickLabel', x_1, 'YTick', 1:numel(x_1))
   xlabel('Udeo snage [%]');
  
 subplot(2,1,2);
   h = barh (y_filtar_relativno);
   set (h, "facecolor", "k");
   title ("Relativna razlika udela u spektru originalnog signala i posle filtriranja");
   set(gca, 'YTickLabel', x_2, 'YTick', 1:numel(x_2))
   xlabel('Udeo snage [%]');
  










