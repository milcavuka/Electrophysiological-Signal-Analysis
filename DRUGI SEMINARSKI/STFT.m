function STFT(x, fs, time) 

N = length(x);
SPEKT = [];
izbor_win = menu('PROZORSKA FUNKCIJA', 'Block Window', ...
    'Hamming', 'Haning', 'Normalized Hanning');
   
sirina = input('\nUneti sirinu prozorske f-je u odbircima: ');
korak_vreme = input('\nUneti korak pomeranja prozora u odbircima: ');  

if (korak_vreme > sirina)
    korak_vreme = sirina; 
end

brojac = 1;
t = linspace(-1, 1, sirina);

% PROZORSKE FUNKCIJE    
% eta = 1e-5;

if izbor_win == 1
% block window has a sharp transition and thus a poor frequency
% localization.
    h = (abs(t) < 1);
    length(h);
    % hf = fftshift(abs(fft(h)));
    % hf = log10(eta + hf);
    % hf = hf / max(hf);
end
%------------------------------------
if izbor_win == 2
% Hamming window is smoother.
    h = cos(t*pi/2) .* (abs(t)<1);
    % hf = fftshift(abs(fft(h)));
    % hf = log10(eta + hf);
    % hf = hf / max(hf);
end
%--------------------------------------
if izbor_win == 3
% Haning window has continuous derivatives.
    h = (cos(t*pi)+1)/2 .* (abs(t)<1);
    % hf = fftshift(abs(fft(h)));
    % hf = log10(eta + hf);
    % hf = hf / max(hf);
end
%--------------------------------------
if izbor_win == 4
% normalized Haning window has a sharper transition.
% It has the advantage of generating a tight frame STFT,
% and is used in the following. 
    h = sqrt(2)/2 * (1+cos(t*pi)) ./ sqrt( 1+cos(t*pi).^2 ) .* (abs(t)<1);
    % hf = fftshift(abs(fft(h)));
    % hf = log10(eta+hf);
    % hf = hf / max(hf);
end
% KRAJ PROZORSKIH FUNKCIJA

dopuna_left = flip(x(1:sirina/2),2); 
dopuna_right = flip(x(N-sirina/2:N-1),2); 

data1 = [dopuna_left x dopuna_right]; 

time_osa = (-sirina/2:length(x)-1+sirina/2)/fs; 

N1 = length(data1); % nova duzina vektora podataka
%ff = (fs/N1)*((-N1/2+1):(N1/2)); % nova frekevencijska osa

%% PETLJA

for tt = 1:korak_vreme:length(time_osa)-(sirina-1)
     
    % ovo sluzi za iscrtavanje
    ttt = zeros(1, N1); 
    ttt(tt:tt+(sirina-1)) = h';
    %
    data_win = data1(tt:tt+(sirina-1));           
    g = data_win.*h; % proizvod prozorske f-je i originalnog signala 
    yy = abs(fft(g)); 
    spekt_i = yy(1:round(length(yy)/2));
    spekt_i = spekt_i.^2;  
    %
    SPEKT(brojac, :) = spekt_i; 
    % crtanje 
    figure(33)
        subplot(3,1,1);
            plot(time_osa, data1/(max(data1)), 'b', ...
                time_osa, ttt, 'r', 'LineWidth', 1); 
                title('Prozorska F-ja'); xlabel('vreme [s]');
                axis('tight');
        subplot(3,1,2);
            % crtanje proizvoda prozorske f-je i signala
            plot(time_osa(tt:tt+(sirina-1)), g);
                title('Proizvod prozora i signala');
                xlabel('vreme [s]'); axis('tight');
        subplot(3,1,3);
            plot(linspace(0,fs/2,length(spekt_i)), spekt_i); 
                title('Spektar proizvoda'); xlabel('frekvencija [Hz]');
                xlim([0 30]);
    %     
    brojac = brojac + 1; % brojac prolazaka kroz petlju
end
     
[T,F] = size(SPEKT);
% vremenska osa za spektrogram
time_novo = linspace(time(1), time(end),T);
% frekvencijska osa za spektrogram
f_novo = linspace(0, fs/2, F);

figure
    subplot(2,1,1)
        plot(time,x);
            xlabel('vreme [s]'); ylabel('amplituda'); axis tight; 
    subplot(2,1,2)
        imagesc(time_novo,f_novo,SPEKT'); colorbar('vert')
            ylabel('frekvencija [Hz]'); xlabel('vreme [s]');
            ylim([0 30]);
            title(['Spektrogram, sirina prozora: ' num2str(sirina) ...
                ', korak prozora: ' num2str(korak_vreme) ]);

end