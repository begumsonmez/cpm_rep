%% PSD of QPSK, MSK and GMSK signals 

clearvars; clc;
N=10000;%Number of symbols to transmit
Fc=800;%carrier frequency in Hertz
L =16; %oversampling factor,use L= Fs/Fc, where Fs >> 2xFc

Fs = L*Fc;
a = rand(N,1)>0.5; %random symbols input to modulator

s1= gmsk_mod(a,Fc,L,0.3); %BT_b=0.3
s2 = qpsk_mod(a,Fc,L);
s3 = gmsk_mod(a,Fc,L,10000); %BT_b=10000 (very large value-MSK)
s4 = bpsk_mod(a,L);

%see section - 'Power Spectral Density plots' for definition
figure; plotWelchPSD(s1,Fs,Fc,'r'); hold on;
plotWelchPSD(s2,Fs,Fc,'b'); 
plotWelchPSD(s3,Fs,Fc,'k');
plotWelchPSD(s4,Fs,Fc,'c');

legend('GMSK with BT_b=0.3','QPSK','MSK','BPSK');
