%%%% Write a function, that given the following parameters returns a
% corresponding FID: A frequency shift in ppm, a relaxation time
% (e.g. 100 ms), the number of samples (e.g. 1024), the time between
% samples (called the dwell time, e.g. 1 ms).

%%% Function: S(t) = S(0)*exp(i(-2pi*f_rot*t+phi)-t/T2*)

% Chemical shift of water is 4.7
% Frequency shift of a metabolite is the ppm eq of the distance to water
% f_rot in the equation is the Hz equiv
% Relaxation time is equal to T2*

function [S,t] = FID(freqshift, T2, Nsamples,DwellTime,weights,S0,larmorfreq,phi)

if ~exist('S0')
    S0 = 1;end
if ~exist('larmorfreq')
    larmorfreq = 100;end
if ~exist('phi')
    phi = 0;end
if ~exist('weights')
    weights = 1;end


% t is modeled as including 0, otherwise use the one below
t = linspace(0,Nsamples*DwellTime,Nsamples);
% t = linspace(DwellTime,Nsamples*DwellTime,Nsamples);


f_rot = (4.7-freqshift)*larmorfreq;
S = 0; %for while-loop
while ~isempty(f_rot)
S = S + weights(1)*S0*exp(1i*(-2*pi*f_rot(1)*t+phi)-t/T2(1));
f_rot(1) = [];
weights(1) = [];
T2(1) = [];
end
S = S.';
S(1) = S(1)/2;
end



