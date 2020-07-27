% Demonstration of Eb/N0 Vs BER for QPSK vs GMSK
close all; clear all; clc;

%% Start of Code
N=100000;%Number of symbols to transmit
EbN0dB = 0:16; % Eb/N0 range in dB for simulation
Fc = 800; % Carrier frequency in Hz (must be < fs/2 and > fg)
BT = 0.3; %Gaussian LPF's BT product
L = 16; %oversampling factor

EbN0lin = 10.^(EbN0dB/10); %converting dB values to linear scale
BER_gmsk = zeros(length(EbN0dB),1); %BER values for each Eb/N0 %GMSK
BER_qpsk = zeros(length(EbN0dB),1); %BER values for each Eb/N0 %QPSK
BER_msk = zeros(length(EbN0dB),1);  %BER values for each Eb/N0 %MSK
BER_bpsk = zeros(length(EbN0dB),1); %BER values for each Eb/N0 %BPSK

%-----------------Transmitter---------------------
%% GMSK
a = rand(N,1)>0.5; %random symbols for modulation
[s_gmsk,t_gmsk,s_complex] = gmsk_mod(a,Fc,L,BT);%GMSK modulation
%% MSK (BT=10000)
[s_msk,t_msk,s_complex2] = gmsk_mod(a,Fc,L,10000);%MSK modulation
%% QPSK
[s_qpsk,t_qpsk] = qpsk_mod(a,Fc,L);%QPSK modulation
%% BPSK
[s_bpsk,t_bpsk] = bpsk_mod(a,L);%QPSK modulation

%GMSK BER calc
for i=1:length(EbN0dB)
    
Eb = sum(abs(s_complex).^2)/(length(s_complex)); %compute Energy
N0 = Eb/EbN0lin(i); %required noise spectral density from Eb/N0
n = sqrt(N0/2)*(randn(size(s_complex))+1i*randn(size(s_complex))); %randn: gaussian distributed rand numbers
r = s_complex + n ; %noise added baseband GMSK signal

%---------------Receiver--------------------
a_cap_gmsk = gmsk_demod(r,L);%Baseband GMSK demodulation
BER_gmsk(i) = sum(a~=a_cap_gmsk)/N;%Bit Error Rate Computation
end

%QPSK BER calc
for j=1:length(EbN0dB)

Eb = L*sum(abs(s_qpsk).^2)/(length(s_qpsk)); %compute energy per bit
N0 = Eb/EbN0lin(j); %required noise spectral density from Eb/N0
n = sqrt(N0/2)*(randn(1,length(s_qpsk)));%computed noise
r_qpsk = s_qpsk + n;%add noise
a_cap_qpsk = qpsk_demod(r_qpsk,Fc,L); %QPSK demodulation
BER_qpsk(j) = sum(a~=a_cap_qpsk.')/N;%Bit Error Rate Computation
end

% MSK BER Calculation
for k=1:length(EbN0dB)
    
Eb = sum(abs(s_complex2).^2)/(length(s_complex2)); %compute Energy
N0 = Eb/EbN0lin(k); %required noise spectral density from Eb/N0
n = sqrt(N0/2)*(randn(size(s_complex2))+1i*randn(size(s_complex2))); %randn: gaussian distributed rand numbers
r_msk = s_complex2 + n ; %noise added baseband GMSK signal

%---------------Receiver--------------------
a_cap_msk = gmsk_demod(r_msk,L);%Baseband GMSK demodulation
BER_msk(k) = sum(a~=a_cap_msk)/N;%Bit Error Rate Computation
end

for g=1:length(EbN0dB)

Eb = L*sum(abs(s_bpsk).^2)/(length(s_bpsk)); %compute energy per bit
N0 = Eb/EbN0lin(g); %required noise spectral density from Eb/N0
n = sqrt(N0/2)*(randn(1,length(s_bpsk)));%computed noise
r_bpsk = s_bpsk + n;%add noise

%---------------Receiver--------------------
a_cap_bpsk = bpsk_demod(r_bpsk,L); %BPSK demodulation
BER_bpsk(g) = sum(a~=a_cap_bpsk)/N;%Bit Error Rate Computation
end


% figure;%Plot performance curves
% semilogy(EbN0dB,BER_gmsk,'g*-','LineWidth',1.5); hold on;
% title('Probability of Bit Error for GMSK, MSK, QPSK and BPSK');
% xlabel('E_b/N_0 (dB)');ylabel('Probability of Bit Error - P_b');
% hold on;
% semilogy(EbN0dB,BER_msk,'co--','LineWidth',1.5); 
% semilogy(EbN0dB,BER_qpsk,'bx-','LineWidth',1.5); 
% semilogy(EbN0dB,BER_bpsk,'m+-.','LineWidth',1.5); 
% legend('GMSK','MSK','QPSK','BPSK');grid on;

M=2;
EbNo = (0:16)';
berQ = zeros(length(EbNo));
berMSK = zeros(length(EbNo));
% Theroretical BER curves
berQ = berawgn(EbNo,'psk',M,'nondiff');
berMSK = berawgn(EbNo,'fsk',M,'coherent');

Fs = L*Fc; Ts=1/Fs; Tb = L*Ts;
B = BT/Tb;
alpha = 0.5887/B;
berGMSK = zeros(length(EbNo),1);
for len=1:length(EbNo)
% berGMSK(len) = qfunc(sqrt(2*EbN0lin(len)*alpha));
berGMSK(len) = 0.5*erfc(sqrt(EbN0lin(len)));
end


close all;
%Plot performance curves
figure(1)
semilogy(EbN0dB,BER_qpsk,'b*','LineWidth',1.5); hold on; semilogy(EbNo,berQ,'r-','LineWidth',1.5);
title('BER for QPSK'); xlabel('E_b/N_0 (dB)');ylabel('Probability of Bit Error - P_b');
legend('simulation','theoretical');

figure(2)
semilogy(EbN0dB,BER_msk,'b*','LineWidth',1.5); hold on; semilogy(EbNo,berMSK,'c-','LineWidth',1.5);
title('BER for MSK'); xlabel('E_b/N_0 (dB)');ylabel('Probability of Bit Error - P_b');
legend('simulation','theoretical');

figure(3)
semilogy(EbN0dB,BER_bpsk,'b*','LineWidth',1.5); hold on; semilogy(EbNo,berQ,'g-','LineWidth',1.5);
title('BER for BPSK'); xlabel('E_b/N_0 (dB)');ylabel('Probability of Bit Error - P_b');
legend('simulation','theoretical');

figure(4)
semilogy(EbN0dB,BER_gmsk,'b*','LineWidth',1.5); hold on; semilogy(EbNo,berGMSK,'m-','LineWidth',1.5);
title('BER for GMSK'); xlabel('E_b/N_0 (dB)');ylabel('Probability of Bit Error - P_b');
legend('simulation','theoretical');

% hold on;
% semilogy(EbN0dB,BER_msk,'co--','LineWidth',1.5); 
% semilogy(EbN0dB,BER_qpsk,'bx-','LineWidth',1.5); 
% semilogy(EbN0dB,BER_bpsk,'m+-.','LineWidth',1.5); 
% legend('GMSK','MSK','QPSK','BPSK');grid on;


