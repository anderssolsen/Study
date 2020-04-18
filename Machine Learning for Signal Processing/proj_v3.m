% Proj v3

%% Creating and recreating sines
fs = 1000;
secs = 10;
phase = pi/2;
t = linspace(0,secs,secs*fs+1);
sin1 = 2*sin(2*pi*5 *t + phase);
sin2 = 1*sin(2*pi*10*t + phase+1);
figure;plot(t,sin1);hold on;plot(t,sin2)

fourier1 = fft(sin1);fourier2 = fft(sin2);
f = linspace(-fs/2,fs/2,secs*fs+1);
figure;plot(f,abs(fftshift(fourier1)));
hold on;plot(f,abs(fftshift(fourier2)));

% reconstructing
sin1re = ifft(fourier1);
sin2re = ifft(fourier2);
figure;plot(t,sin1re);hold on;plot(t,sin2re)
figure;plot(t,sin1re-sin1);hold on;plot(t,sin2re-sin2)

% Conclusion: should be possible for periodic signals.

%% Creating and recreating aperiodic signals
clear;
fs = 1000;
secs = 10;



%% Miniproject on sparsity aware learning
% Anders Olsen and Lena Mahler Nilsson

clear
close all

[data,Fs] = audioread('Guitar.m4a');    % Loading created audio signal
data = data(1:340001,1);                % Cutting out the desired time
% signal = downsample(data,20);           % Downsampling the signal
signal = data;
% signal = [signal;0];
clearvars data
% Fs = Fs/20;                             % Adjusting the sampling frequency

time = linspace(0,length(signal)/Fs,length(signal)); % Creating time axis


%% Visualization of signal
figure
plot(time,signal)
title('Sound signal','FontSize',16)
xlabel('Time [s]','FontSize',16)
ylabel('Amplitude','FontSize',16)
print -dpng signal


%% Power spectrum
N = length(signal)  ;
% w = hamming(N);
ft = fft(signal);                       % Fourier transform
ftshift = fftshift(ft);                      % zero-shifted (zero in center)
psd = abs(ftshift).^2/N;                     % Power spectral density (ABS)
freqaxis = (-N/2:N/2-1)*(Fs/N);         % zero-centered frequency range


figure;
plot(freqaxis,psd,'LineWidth',1.5);
title('Two-sided PSD','FontSize',16)
xlabel('Frequency [Hz]','FontSize',16)
ylabel('Power [energy/frequency]','FontSize',16)
xlim([-1500 1500])
 print -dpng PSD
clearvars ftshift


n = 0:(N-1);
k = n';
% M = exp((-1i*2*pi*k*n)/N); %only compute this once!  Takes disk space
% y = M*signal;
% figure;
% plot(freqaxis,abs(fftshift(y)).^2/N,'LineWidth',1.5);
% title('Two-sided PSD','FontSize',16)
% xlabel('Frequency [Hz]','FontSize',16)
% ylabel('Power [energy/frequency]','FontSize',16)

%% Lasso with sparse model
tic

X = sparse(N,N);                        % Setup model matrix
v = ones(N,1);
X(1:N+1:end) = v;                       % Ones in diagonal
% sol_lasso = lasso(X,psd, 'Lambda', lambda); % LASSO solution
toc

lambda = 10^(3);                     % (-3.5) selects 2 frequencies

% LSre = lsqminnorm(X,real(ft));
% LSim = lsqminnorm(X,imag(ft));
LS = lsqminnorm(X,abs(ft));toc
% sol_lassore = sign(LSre).*max(LSre-lambda/2,0);
% sol_lassoim = sign(LSim).*max(LSim-lambda/2,0);
% sol_lasso = complex(sol_lassore,sol_lassoim);
sol_lasso = sign(LS).*max(LS-lambda/2,0);
% figure;plot(abs(sol_lasso))
% hold on; plot(abs(sol_lasso))
% toc
notzero = find(sol_lasso>0);
for p = 1:length(notzero)
    if sol_lasso(notzero(p))<sol_lasso(notzero(p)+1)
        sol_lasso(notzero(p)) =0;
    elseif sol_lasso(notzero(p))<sol_lasso(notzero(p)-1)
        sol_lasso(notzero(p)) =0;
    end
end;toc
figure;
plot(freqaxis,fftshift(sol_lasso).^2/N,'LineWidth',1.5)
xlim([-1500 1500])
title('Sparse frequency spectrum', 'FontSize',16)
xlabel('Frequency [Hz]')
ylabel('Power [energy/frequency]','FontSize',16)
print -dpng sol_lasso
%% Reconstructing signal

figure;                                 % Plot "selected" frequencies
plot(freqaxis,fftshift(sol_lasso).^2/N,'LineWidth',1.5)
title('Reconstructed PSD','FontSize',16)
xlabel('Frequency [Hz]','FontSize',16)
ylabel('Power [energy/frequency]','FontSize',16)
xlim([-1500 1500])
 print -dpng sol_lasso_fftshift

% sol_lasso_shift = ifftshift(sol_lasso);  % Inverse zero-shift
recons = ifft(sol_lasso);          % Sound signal reconstruction
recons = abs(recons);                     % Absolute value
figure                                   % Plot signal reconstruction
plot(time,recons)
title('Reconstructed sound signal','FontSize',16)
xlabel('Time [s]','FontSize',16)
ylabel('Amplitude','FontSize',16)
 print -dpng Reconsound
soundsc(recons,Fs);                     % Listen to reconstruction

%%%%% Amplitudes in both reconstructed PSD and reconstructed signal are not
%%%%% the same - but the sound signal is good anyways


%% Spectrogram of signal
% figure
spectrogram(signal,100,[],[],Fs,'yaxis'); % Computing and visualizing the spectrogram of the signal
% % [s,f,t] = spectrogram(signal,600,550,900,Fs); % Computing and visualizing the spectrogram of the signal
% title('Spectrogram of signal','FontSize',16)
% print -dpng spec
% BlockSize = Fs/2;
% % window = hamming(BlockSize)
% % overlap = Blocksize
% 
% 
% % s = stft(signal,
% s = STFT_mns(signal,Fs/5,Fs/5);
%%
clearvars -except signal Fs time recons
B = round(Fs/10)+1; %Block size
Bs = Fs/20; %Block skip
nB = floor((length(signal)-B)/Bs);
y = zeros(B, nB);

ndft = 1:B;
k = ndft';
M = exp((-1i*2*pi*k*ndft)/B); %only compute this once!  Takes disk space


lambda = 10^(1.8);
N = size(y,1);
X = sparse(N,N);                        % Setup model matrix
v = ones(N,1);
X(1:N+1:end) = v;
tic
for b = 1:nB
    
    x = signal(Bs*(b-1)+(1:B));
    y(:,b) = M*x;
    
    LS = lsqminnorm(X,abs(y(:,b)));
    sol_lasso(:,b) = sign(LS).*max(LS-lambda/2,0);
    
    notzero = find(sol_lasso(:,b)>0);
    for p = 1:length(notzero) try
        if sol_lasso(:,notzero(p))<sol_lasso(:,notzero(p)+1)
            sol_lasso(:,notzero(p)) =0;
        elseif sol_lasso(:,notzero(p))<sol_lasso(:,notzero(p)-1)
            sol_lasso(:,notzero(p)) =0;
        end
        catch end
    end
    recon(:,b) = abs(ifft(sol_lasso(:,b)));
    recon2(:,b) = recon(1:Bs,b);
    toc
end
figure;
plot(sol_lasso)
title('Sparse frequency spectrum of STFT','FontSize',16)
xlabel('Frequency [Hz]','FontSize',16)
ylabel('Power [energy/frequency]','FontSize',16)
reconscon = reshape(recon2,[1,size(recon2,1)*size(recon2,2)]);
figure;
plot(time(1:length(reconscon)),reconscon)
title('Reconstructed sound signal from spectrogram','FontSize',16)
xlabel('Time [s]','FontSize',16)
ylabel('Amplitude','FontSize',16)
print -dpng ReconsoundSTFT
soundsc(reconscon,Fs)
%% incorporate prior knowledge - only one frequency at a time
clearvars -except signal Fs time recons reconscon
B = round(Fs/10)+1; %Block size
Bs = Fs/20; %Block skip
nB = floor((length(signal)-B)/Bs);
y = zeros(B, nB);

ndft = 1:B;
k = ndft';
F = exp((-1i*2*pi*k*ndft)/B); %only compute this once!  Takes disk space
sol_lasso2 = zeros(B,nB);

lambda = 10^(0.9);
N = size(y,1);
X = sparse(N,N);                        % Setup model matrix
v = ones(N,1);
X(1:N+1:end) = v;
tic
for b = 1:nB
    
    x = signal(Bs*(b-1)+(1:B));
    y(:,b) = F*x;
    
    LS = lsqminnorm(X,abs(y(:,b)));
    sol_lasso(:,b) = sign(LS).*max(LS-lambda/2,0);
    
    [sorted,index] = sort(sol_lasso(:,b),'descend');
    sol_lasso2(index(1),b) = sorted(1);
    sol_lasso2(index(2),b) = sorted(2);
    
    recon(:,b) = abs(ifft(sol_lasso2(:,b)));
    recon2(:,b) = recon(1:Bs,b);
    toc
end
figure;
plot(sol_lasso2)
title('Sparse frequency spectrum of STFT','FontSize',16)
xlabel('Frequency [Hz]','FontSize',16)
ylabel('Power [energy/frequency]','FontSize',16)
print -dpng sol_lasso_prior
reconscon2 = reshape(recon2,[1,size(recon2,1)*size(recon2,2)]);
figure
plot(time(1:length(reconscon2)),reconscon2)
title('Reconstructed sound signal from spectrogram','FontSize',16)
xlabel('Time [s]','FontSize',16)
ylabel('Amplitude','FontSize',16)
print -dpng ReconsoundSTFT_prior
soundsc(reconscon2,Fs)
%%


soundsc(signal,Fs) %slide 3
soundsc(recons,Fs) % slide 4
soundsc(reconscon,Fs) %slide 5
soundsc(reconscon2,Fs) %slide 6














