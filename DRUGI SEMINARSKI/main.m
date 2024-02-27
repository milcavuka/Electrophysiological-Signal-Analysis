clear;
close all; 
clc; 

Fs = 256;                          % frekvencija odabiranja EEG signala
T = 30;                            % trajanje signala u sekundama
time = 0:1/Fs:T-1/Fs;     

% ucitavanje EEG signala

EEG = importdata('EEG_clean.mat');
EEG = EEG(:, 1:T*256); 
N = length(EEG);

for i=1:20
    
    eegfreqshift = fftshift(abs(fft(EEG(i,:))));
    f = (Fs/N)*((-N/2+1):(N/2));
    
    figure(i)
      subplot(2,1,1)
        plot(time, EEG(i,:),'Color', [0 0 0]);
          %xlim([0, 10]);
          xlabel('Vreme [s]');
          ylabel('Amplituda [\muV]' ); 
          title(['Vremenski prikaz signala pre filtriranja kanal ', num2str(i)]);
          grid on;
      subplot(2,1,2)
        plot(f, eegfreqshift,'Color', [0 0 0]);
          xlim([-50, 50]);  
          xlabel('Frekvencija [Hz]');
          ylabel('Magnituda [a.u.]');
          title(['Spektar signala pre fitriranja kanal ', num2str(i)]);
          grid on;
end 

% FILTRIRANJE NIJE POTREBNO
% filtriran_eeg = zeros(20,N);
% for i=1:20
%     
%    
%       
%      [b,a] = butter(0.5,1/Fs/2,'high');
%      filtriran_eeg(i,:) = filter(b,a,EEG(i,:));     
%      eegfreqshift_filt = fftshift(abs(fft(filtriran_eeg(i,:))));     
%           
%      figure()
%      
%       subplot(2,1,1)
%         plot(time, filtriran_eeg(i,:),'Color', [0 0 0]);
%           xlim([0, 10]);
%           ylim([-30 30]);
%           xlabel('Vreme [s]');
%           ylabel('Amplituda [\muV]' ); 
%           title(['Vremenski prikaz signala posle filtriranja kanal ', num2str(i)]);
%           grid on;
%       subplot(2,1,2)
%         plot(f, eegfreqshift_filt,'Color', [0 0 0]);
%           xlim([-50, 50]);  
%           xlabel('Frekvencija [Hz]');
%           ylabel('Magnituda [a.u.]');
%           title(['Spektar signala posle fitriranja kanal ' , num2str(i)]);
%           grid on;
%      
% end 
%%

for i = 1:20
    
    
    window = 512;
    noverlap= 512/2;
    % 
    [S,F,T,P]=spectrogram(EEG(i,:),window,noverlap,[],Fs);
    normalized_spectrogram = 10*log10(abs(S) / max(abs(S(:))));
    figure(40+i)
    subplot(2,1,1)
    colormap bone
    imagesc(T,F,normalized_spectrogram); 
    axis xy
    ylim([0 25]); 
    colorbar('vert')
    xlabel('Vreme[s]'); ylabel('Frekvencija[Hz]');
    title(['Spektrogram za sirinu prozora 512 kanal ', num2str(i)])
    
    window = 1024;
    noverlap= 1024/2;
    [S,F,T,P]=spectrogram(EEG(i,:),window,noverlap,[],Fs, 'yaxis');
    normalized_spectrogram = 10*log10(abs(S) / max(abs(S(:))));
    subplot(2,1,2)
    colormap bone
    imagesc(T,F,normalized_spectrogram);
    axis xy
    ylim([0 25]); 
    colorbar('vert')
    xlabel('Vreme[s]'); ylabel('Frekvencija[Hz]');
    title(['Spektrogram za sirinu prozora 1024 kanal ' , num2str(i)])
    
    figure(60+i)
    [cfs,f] = cwt(EEG(i,:),Fs,'FrequencyLimits',[0 30]);
   
   
     normalized_cwt = abs(cfs) / max(abs(cfs(:)));
    colormap bone
    surface(time,f,10*log10(normalized_cwt));
    shading flat
    
    colorbar('vert')
    xlabel('Vreme[s]'); ylabel('Frekvencija[Hz]');
    %set(gca,"yscale","log")
    title(['Skalogram za kanal ' ,num2str(i)])
    
end
%%
bpm = 90;
duration = 30;
fs = 256;
amp = 100; 

s = ECGwaveGen(bpm,duration,fs,amp)./10;

figure(81)

    subplot(2,1,1)
    plot(time, s,'Color', [0 0 0])
    xlim([0, 10]);
    ylim([-10 60]);
    xlabel('Vreme [s]');
    ylabel('Amplituda [\muV]' ); 
    title(['Vremenski prikaz suma EKG-a ']);
    grid on;
    
    subplot(2,1,2)
    ecgfreqshift = fftshift(abs(fft(s)));
    f = (Fs/N)*((-N/2+1):(N/2));
    plot(f, ecgfreqshift,'Color', [0 0 0]);
          xlabel('Frekvencija [Hz]');
          ylabel('Magnituda [a.u.]');
          title(['Spektar signala EKG-a suma']);
          grid on;
          xlim([-50,50])

eeg_example = EEG(1,:) + s;

window = 512;
noverlap= 512/2;
%
[S,F,T,P]=spectrogram(eeg_example,window,noverlap,[],Fs);
normalized_spectrogram = 10*log10(abs(S) / max(abs(S(:))));
figure(83)
subplot(2,1,1)
colormap bone
imagesc(T,F,normalized_spectrogram); 
ylim([0 25]); 
colorbar('vert')
xlabel('Vreme[s]'); ylabel('Frekvencija[Hz]');
title(['Spektrogram za sirinu prozora 512 zasumljeni kanal 1'])

window = 1024;
noverlap= 1024/2;
[S,F,T,P]=spectrogram(eeg_example,window,noverlap,[],Fs);
normalized_spectrogram = 10*log10(abs(S) / max(abs(S(:))));
subplot(2,1,2)
colormap bone
imagesc(T,F,normalized_spectrogram);
ylim([0 25]); 
colorbar('vert')
axis xy
xlabel('Vreme[s]'); ylabel('Frekvencija[Hz]');
title(['Spektrogram za sirinu prozora 1024 zasumljeni kanal 1'])

figure(84)
[cfs,f] = cwt(eeg_example,Fs,'FrequencyLimits',[0 30]);
normalized_cwt = abs(cfs) / max(abs(cfs(:)));

colormap bone
surface(time,f,(10*log10(normalized_cwt)));
shading flat

colorbar('vert')
xlabel('Vreme[s]'); ylabel('Frekvencija[Hz]');
%set(gca,"yscale","log")
title(['Skalogram za zasumljeni kanal 1' ])



% 
%     [cfs,f] = cwt(filtriran_eeg(i,:),Fs);
%     colormap bone
%     imagesc(time,f,10*log10(abs(cfs))); 
%     ylim([0 25]); 
%     colorbar('vert')
%     xlabel('Vreme[s]'); ylabel('Frekvencija[Hz]');