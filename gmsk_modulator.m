close all; clear all; clc;
%% GMSK

%% Step 1: Calling the Gaussian LPF Function and Plotting Its frequency response
N=10000;%Number of symbols to transmit
Fc=800;%carrier frequency in Hertz
L =16; %oversampling factor,use L= Fs/Fc, where Fs >> 2xFc
Fs = L*Fc; Ts=1/Fs; Tb = L*Ts;
k=1;%truncation length for Gaussian LPF

[h1, t1] = gaussianLPF(0.3,Tb,L,k);
[h2, t2] = gaussianLPF(0.28,Tb,L,k);
[h3, t3] = gaussianLPF(0.25,Tb,L,k);
[h4, t4] = gaussianLPF(100,Tb,L,k);

figure()
hold on
p1 = plot(t1,h1); p2 = plot(t2,h2); p3 = plot(t3,h3); p4 = plot(t4,h4);
ylabel('h(t)'); xlabel('time(s)'); title('Impulse response for different BT values');
M1 = "BT = 0.3";
M2 = "BT = 0.28";
M3 = "BT = 0.25";
M4 = "BT = 100";
legend([p1; p2; p3; p4], [M1; M2; M3; M4]);
%% Step 2: Implementing the GSMK Modulator Function

% PSD of GMSK signals with various BT products
clearvars; clc;
N=10000;%Number of symbols to transmit
Fc=800;%carrier frequency in Hertz
L =16; %oversampling factor,use L= Fs/Fc, where Fs >> 2xFc
Fs = L*Fc;

%Random Number Generator
a = rand(N,1)>0.5; %random symbols input to modulator

s1 = gmsk_mod(a,Fc,L,0.3); %BT_b=0.3
s2 = gmsk_mod(a,Fc,L,0.5); %BT_b=0.5
s3 = gmsk_mod(a,Fc,L,0.7); %BT_b=0.7
s4 = gmsk_mod(a,Fc,L,10000); %BT_b=10000 (very large value-MSK)
%see section - 'Power Spectral Density plots' for definition
figure; plotWelchPSD(s1,Fs,Fc,'r'); hold on;
plotWelchPSD(s2,Fs,Fc,'b'); plotWelchPSD(s3,Fs,Fc,'m');
plotWelchPSD(s4,Fs,Fc,'k');
xlabel('f-f_c'); ylabel('PSD (dB/Hz)');
legend('BT_b=0.3','BT_b=0.5','BT_b=1','BT_b=\infty');














