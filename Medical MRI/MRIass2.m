%%% Initialize
clear;close all
directory = 'C:\Users\ander\Documents\10. Semester\Medical MRI\ExerciseData\ass2';
datadir = strcat(directory,'\Dataold\');
addpath(genpath(directory))
cd(directory)

% Info from SPAR file:
Nsamples = 1024;
synthFreq = 127759137;  %Hz (?, must be)
fs = 2000; %Hz
TE = 0.144;               %s
TR = 3;              %s
numAvs = 48;            %48 averages for both scans
DwellTime = 0.0005000000237; %secs

% Constants:
ppmNAA   = 2.01; % Frequency shift in ppm
ppmCho   = 3.20;
ppmCr    = 3.04;
ppmWater = 4.70;
ppmLipid = 1.30;
Chemshift = [ppmWater,ppmNAA,ppmCho,ppmCr,ppmLipid]; % Chemical  shift
freqshift = (ppmWater - Chemshift)*synthFreq/10^6;   % Frequency shift

Nprot = [2,3,9,5]';
T1 = [4,1.5,1.2,1.6]';
T2     = [0.1,0.450,0.360,0.210,0.085]';
% Nmetabs = length(Chemshifts);



% Load data
FID_act = readSDAT('X35819_WIP_SVS_wref+wsup_integrated_6_2_raw_act.SDAT',0);
FID_ref = readSDAT('X35819_WIP_SVS_wref+wsup_integrated_6_2_raw_ref.SDAT',0);


%% Some initial plots
% Freqaxis from -980Hz to 1000Hz
freqaxis = (1:Nsamples)/(Nsamples*DwellTime)-1/(DwellTime*2);
% Freqaxis in ppm (slide 11) scaled to water
ppmaxis = freqaxis/(10^(-6)*synthFreq)+ppmWater;


act_fft = fftshift(fft(FID_act));
ref_fft = fftshift(fft(FID_ref));

%%%%% As PPM
figure;subplot(2,1,1);
plot(ppmaxis,abs(act_fft),'k');
set(gca, 'XDir','reverse')
xlabel('Frequency (ppm)')
ylabel('Spectrum')
title('Water-suppressed scan')
axis([-2 12 -inf inf])

subplot(2,1,2);
plot(ppmaxis,abs(ref_fft),'k');
set(gca, 'XDir','reverse')
xlabel('Frequency (ppm)')
ylabel('Spectrum')
title('Water reference scan')
axis([-2 12 -inf inf])


%%%%% As frequency
figure;subplot(2,1,1);
plot(freqaxis,abs(act_fft));
xlabel('Frequency (Hz)')
ylabel('Spectrum')

subplot(2,1,2);
plot(freqaxis,abs(ref_fft));
xlabel('Frequency (Hz)')
ylabel('Spectrum')

%% Building the design matrix

% Elements: 
% 1) Scaled unsuppressed water signals
t = ((0:Nsamples-1)*DwellTime)';

FIDs = FID_ref.*exp(-1i*2*pi*freqshift.*t);     % Frequency shift in time
FIDs = FIDs./exp(-t/T2(1)).*exp(-t./T2');       % Scale by T2 of water

% We need to normalize to avoid blow-up????
% We use plots to see when the constant c is right
c = 5;
FIDs = FIDs.*exp(-c*t);
% figure;
% subplot(3,1,1);plot(t,abs(FIDs(:,2)));
% subplot(3,1,2);plot(t,abs(FIDs(:,3)));
% subplot(3,1,3);plot(t,abs(FIDs(:,4)));

% FID_act = FID_act.*exp(c*t);
% FID_ref = FID_ref.*exp(c*t);

% Convert to spectra
X = fftshift(fft(FIDs),1);
% X = abs(X); %This is right??
X = X./sum(X);

% 2) Derivatives using symmetric difference quotient
der1 = symdifquot(X);
der2 = symdifquot(der1);
X = [X,der1,der2];

% 3) Baseline
X(:,size(X,2)+1) = ones(1,Nsamples);
X(:,size(X,2)+1) = linspace(0,1,Nsamples); %Is this right??
X(:,size(X,2)+1) = X(:,size(X,2)).^2;

% X = X./sum(X)
% act_fft=act_fft*exp(-c*t)
w = (X'*X)^(-1)*X'*(act_fft);
% w1 = (abs(X)'*abs(X))^(-1)*abs(X)'*abs(act_fft);
% Plot of regressors
figure;plot(ppmaxis,abs(X))
set(gca, 'XDir','reverse')
legend('Water','Naa','Cho','Cr','Lipid','dWater','dNaa','dCho','dCr','dLipid',...
'd2Water','d2Naa','d2Cho','d2Cr','d2Lipid','BL1','BL2','BL3')

figure;plot(ppmaxis,abs(X(:,1:5)))
set(gca, 'XDir','reverse')
legend('Water','Naa','Cho','Cr','Lipid')

%% Some intermediary plots
act_fft = fftshift(fft(FID_act));

figure;subplot(3,1,1);
plot(ppmaxis,abs(act_fft),'k');
set(gca, 'XDir','reverse')
% xlabel('Frequency (ppm)')
ylabel('Spectrum')
axis([0 8 -0.1 0.8])
title('Original signal')

subplot(3,1,2);
plot(ppmaxis,abs(X*w),'k');
set(gca, 'XDir','reverse')
% xlabel('Frequency (ppm)')
ylabel('Spectrum')
axis([0 8 -0.1 0.8])
title('Reconstructed signal after GLM')

subplot(3,1,3);
plot(ppmaxis,abs(act_fft-X*w),'k');
set(gca, 'XDir','reverse')
xlabel('Frequency (ppm)')
ylabel('Spectrum')
axis([0 8 -0.1 0.1])
title('Residual y-X\beta')

figure;
plot(freqaxis,abs(act_fft)+0.2);hold on
set(gca, 'XDir','reverse')
plot(freqaxis,abs(X*w));
set(gca, 'XDir','reverse')
plot(freqaxis,abs(act_fft-X*w)-0.2);
legend('signal','Model','Residual')
set(gca, 'XDir','reverse')
xlabel('Frequency (ppm)')
ylabel('Spectrum')

figure;
plot(ppmaxis,abs(X*w),'k');hold on
plot(ppmaxis,abs(act_fft),'r');
set(gca, 'XDir','reverse')
%  'Color', uint8([50 50 50]))
%%
figure;hold on
names = {'Water','NAA','Cho','CR','Lipids',...
         'dWater','dNAA','dCho','dCR','dLipids',...
         'd2Water','d2NAA','d2Cho','d2CR','d2Lipids'...
         'Constant','Linear','2nd-order'};
     for k = 1:18
         
         if k<16
         plot(ppmaxis,k+real(X(:,k))/max(abs(real(X(:,k))))/2);
         plot(ppmaxis,k+imag(X(:,k))/max(abs(imag(X(:,k))))/2);
         else
             plot(ppmaxis,k+X(:,k));
         end
     end  
set(gca, 'XDir','reverse')
xlabel('Frequency (ppm)')

%% Concentration estimation


% Sm(1) = abs(sum(X(:,[2])*w([2]))); %NAA
% Sm(2) = abs(sum(X(:,[3])*w([3]))); %Cho
% Sm(3) = abs(sum(X(:,[4])*w([4]))); %Cr
Sm(1) = abs(sum(X(:,[2,7,12])*w([2,6,11])));
Sm(2) = abs(sum(X(:,[3,8,13])*w([3,7,12])));
Sm(3) = abs(sum(X(:,[4,9,14])*w([4,8,13])));
% Sm(1) = abs(w(2)); %NAA
% Sm(2) = abs(w(3)); %Cho
% Sm(3) = abs(w(4)); %Cr


FID_wat = FID_ref.*exp(-c*t);
Sw = abs(sum(fft(FID_wat)));

Cw = 45;
Ww = exp(-TE/T2(1));
Wm = exp(-TE./T2(2:4)).*(1-exp(-TR./T1(2:4)));

C = Ww*Nprot(1)*Sm'*Cw./(Wm.*Nprot(2:4)*Sw)

%% with Ww2
Ww = exp(-TE/T2(1)).*(1-exp(-TR./T1(1)));
C = Ww*Nprot(1)*Sm'*Cw./(Wm.*Nprot(2:4)*Sw)
%% other way to test
Sw2 = abs(FID_ref(1));
Smn = ifft(ifftshift(X(:,2)*w(2)));
SmC = ifft(ifftshift(X(:,3)*w(3)));
Smc = ifft(ifftshift(X(:,4)*w(4)));
Sm2(1) = abs(Smn(1));
Sm2(2) = abs(SmC(2));
Sm2(3) = abs(Smc(3));

C = Ww*Nprot(1)*Sm2'*Cw./(Wm.*Nprot(2:4)*Sw2)


%% t-test
k = 18; %num columns
df = Nsamples-k;
D = (act_fft-X*w)'*(act_fft-X*w); %distance
var_beta = D/df; %variance estimate

tmp = (X'*X)^(-1);
t2 = abs(w(2)/(var_beta*sqrt(tmp(2,2))))
t3 = abs(w(3)/(var_beta*sqrt(tmp(3,3))))
t4 = abs(w(4)/(var_beta*sqrt(tmp(4,4))))

p2 = 1-tcdf(t2,df)
p3 = 1-tcdf(t3,df)
p4 = 1-tcdf(t4,df)









