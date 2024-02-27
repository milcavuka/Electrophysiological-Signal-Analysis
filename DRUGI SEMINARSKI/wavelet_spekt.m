function wavelet_spekt(x, fs, time)
    res = [];
    sigm = 0.7;
    
    fmin = input('Unesite f_min: ');
    fmax = input('Unesite f_max: ');
    
    ns = 0;
    fo = 1;
    pas = 0.5; 
    pass = 1;
    number = fix(length(x)/2)*2;
    
    data = x';
    fsamp = fs;
    nfft = number;
    nsaut = ns;

% function res = wt_mor1(data,fsamp,nfft,nsaut,fo,fmin,fmax,pas,pass,sigm);
% time frequency analysis by wavelet analysis of Morlet
% x-signal vector
% Fsamp-sample frequency
% nfft-number of points of the signal
% nsaut-number of points not used in the begining of signal (=0)
% fo- fo of the wavelet (=1)
% fmin-fmin of the representation
% fmax-fmax
% pas-step of the scale between fmin and fmax
% pass-step of the time of the result (step of the resolution)
% sigm- time frequency resolution general case =.7;
% res-result (it is better to plot abs(res))

    y = fft(data(1+nsaut:nfft+nsaut));
    brojac = 1;
    
    WT_menu = menu('Wavelet', 'Morlet', 'Window functions');
    
    if WT_menu == 2
        WTmenu2 = menu('menu', 'Block', 'Hamming', 'Haning', ...
            'Normalized Hanning');
    end
     
    if WT_menu == 2
        eta = 1e-5;
        h = log10(eta+y); y = y/max(y);
    end
    
        
 for fff = fmin:pas:fmax 
    
    if WT_menu == 1
        
    a = 1/fff;
    wavelet= exp(-2/sigm*(pi)^2*(([0:nfft/2-1,-nfft/2:-1]*(fsamp*a/nfft))-(fo)).^2);      
    size(wavelet.'); 
    g = wavelet.'.*y;
    t = ifft(g).';
    res(brojac,:) = (a.^0.5).*t(1:pass:nfft);
   
    % WT real time
    figure(55)
        subplot(3,1,1);
            plot(linspace(0,fs,number), abs(wavelet), 'r', ...
                linspace(0,fs,number), abs(y)/max(abs(y)), 'b'); 
                title('Wavelet u domenu FT');
                xlabel('frekvencija [Hz]');
                axis tight; xlim([fmin fmax]);
        subplot(3,1,2);
            % crtanje proizvoda prozorske f-je i signala
            plot(linspace(0,fs,number), abs(g));
                title('Proizvod FT wavelet fje i FT signala');
                xlabel('frekvencija [Hz]');
                axis tight; xlim([fmin fmax]);
        subplot(3,1,3);
            plot(linspace(time(1), time(end), length(res(1,:))), ...
                real(res(brojac,:))); 
                title('BP signal'); xlabel('vreme [s]'); axis tight;     
    end
    
    if WT_menu == 2
        a = 1/fff; %scale
       
        if WTmenu2 == 1 
            sigm2 = 0.3;
            % h = abs((([0:number/2-1,-number/2:-1]*(fs*a/number))-(fo)))<a;%block window fiksne duzine
            h = 1/sigm2*abs((([0:number/2-1,-number/2:-1]*(fs*a/number))-(fo)))<sigm2;%block window koji se siri i pomera
        end

        if WTmenu2 == 2 
            sigm2 = 0.2;
            h = cos((([0:number/2-1,-number/2:-1]*(fs*a/number))-(fo))*1/sigm2*pi/2).*(abs((([0:number/2-1,-number/2:-1]*(fs*a/number))-(fo)))<sigm2);%Hamming
        end

        if WTmenu2 == 3
            sigm2 = 0.2;
            h = (cos((([0:number/2-1,-number/2:-1]*(fs*a/number))-(fo))*1/sigm2*pi)+1)/2 .* (abs((([0:number/2-1,-number/2:-1]*(fs*a/number))-(fo)))<sigm2);%Haning
        end

        if WTmenu2 == 4 
            sigm2 = 0.2;
            h = sqrt(2)/2*(1+cos((([0:number/2-1,-number/2:-1]*(fs*a/number))-(fo))*1/sigm2*pi))./sqrt(1+cos((([0:number/2-1,-number/2:-1]*(fs*a/number))-(fo))*1/sigm2*pi).^2 ).*(abs((([0:number/2-1,-number/2:-1]*(fs*a/number))-(fo)))<sigm2);
        end  

        g = h.'.*y;
        t = ifft(g).';
        res(brojac,:) = (a.^0.5).*t(1:pass:number);
    
         % WT real time
         figure(66)
            subplot(3,1,1);
                plot(linspace(0,fs,number), (h), 'r', ...
                    linspace(0,fs,number), abs(y)/max(abs(y)), 'b'); 
                    title('Wavelet u domenu FT');
                    xlabel('frekvencija [Hz]');
                    axis tight; xlim([fmin fmax]);
            subplot(3,1,2);
                % crtanje proizvoda prozorske f-je i signala
                plot(linspace(0,fs,number), abs(g));
                    title('Proizvod FT wavelet f-je i FT signala');
                    xlabel('frekvencija [Hz]');
                    axis tight; xlim([fmin fmax]);
            subplot(3,1,3);
                plot(linspace(time(1),time(end),length(res(1,:))), ...
                    real(res(brojac,:))); 
                    title('BP signal'); xlabel('vreme [s]'); axis tight;
        %
    end

    brojac = brojac + 1;      
 end
  
%% KRAJ WAVELET FUNKCIJA    

[T,F] = size(res');
f_novo = linspace(fmin, fmax, F);
time_novo = linspace(time(1), time(end), T);

if WT_menu == 1
    figure
        subplot(2,1,1)
            plot(time, x);
                xlabel('vreme [s]'); ylabel('amplituda'); axis tight; 
                title('signal'); xlabel('vreme [s]'); axis('tight'); 
        subplot(2,1,2)
            imagesc(time_novo,f_novo,abs(res)); colorbar('vert');
                ylabel('frekvencija [Hz]'); xlabel('vreme [s]');
                ylim([0 30]); title('Wavelet t-f prezentacija');
end


if WT_menu == 2
    figure
        subplot(2,1,1)
            plot(time, x);
                xlabel('vreme [s]'); ylabel('amplituda'); axis tight; 
                title('signal'); xlabel('vreme [s]'); axis('tight'); 
        subplot(2,1,2)
            imagesc(time_novo, f_novo, abs(res));
                ylabel('frekvencija [Hz]'); xlabel('vreme [s]');
                ylim([0 30]); title('Wavelet f-t prezentacija');
end