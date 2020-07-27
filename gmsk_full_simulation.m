% Demonstration of Eb/N0 Vs SER for baseband GMSK modulation scheme
close all; clear all; clc;

%% Steps for the GMSK Modulation:
% MODULATOR
% 1. Create bitstream: a
% 2. Turn the bistream into a series of rectangular pulses (of amplitude 0 and 1V)
% 3. Convert this from unipolar to binary NRZ format ( -1V for "0" and +1V for "1"): c(t)
% 4. Apply Gaussian filter: b(t)
% 5. Integrate the waveform: phi(t)
% 6. I/Q Modulation: Obtain the in-phase and quadrature components of the signal: I(t) and Q(t)
% 7. Multiply I and Q with the carrier (cos and a (pi/2 shift + cos = ) sin) and add together: s(t)
% *****TRANSMIT THROUGH CHANNEL****** : AWGN or Fading
% DEMODULATOR
% 8. 



%% Start of Code
N=100000;%Number of symbols to transmit
EbN0dB = 0:2:18; % Eb/N0 range in dB for simulation
Fc = 800; % Carrier frequency in Hz (must be < fs/2 and > fg)
BT= 0.3; %Gaussian LPF's BT product
L = 16; %oversampling factor
EbN0lin = 10.^(EbN0dB/10); %converting dB values to linear scale
BER = zeros(length(EbN0dB),1); %SER values for each Eb/N0

%-----------------Transmitter---------------------
a = rand(N,1)>0.5; %random symbols for modulation
[s,t,s_complex] = gmsk_mod(a,Fc,L,BT);%GMSK modulation

for i=1:length(EbN0dB)
    
Eb = sum(abs(s_complex).^2)/(length(s_complex)); %compute Energy
N0 = Eb/EbN0lin(i); %required noise spectral density from Eb/N0
n = sqrt(N0/2)*(randn(size(s_complex))+1i*randn(size(s_complex))); %randn: gaussian distributed rand numbers
r = s_complex + n ; %noise added baseband GMSK signal

%---------------Receiver--------------------
a_cap = gmsk_demod(r,L);%Baseband GMSK demodulation
BER(i) = sum(a~=a_cap)/N;%Bit Error Rate Computation
end

figure;%Plot performance curves
semilogy(EbN0dB,BER,'g*-','LineWidth',1.5); hold on;%simulated BER
title('Probability of Bit Error for GMSK modulation');
xlabel('E_b/N_0 (dB)');ylabel('Probability of Bit Error - P_b');

