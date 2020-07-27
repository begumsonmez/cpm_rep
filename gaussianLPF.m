function [h,t] = gaussianLPF(BT,Tb,L,k)
%Function to generate impulse response of a Gaussian low pass filter
%% Parameters
%BT - BT product - Bandwidth x bit period
%Tb - bit period
%L - oversampling factor (number of samples per bit)

%k - span length of the pulse (bit interval): normally the gaussian filter's impulse response is of infinite length, 
%but this obviously is not applicable to a digital modulation scenario, so we truncate h(t) to k*bit_duration

%h - impulse response of the Gaussian pulse
%t - generated time base

B = BT/Tb;%bandwidth of the filter
t=-k*Tb:Tb/L:k*Tb; %truncated time limits for the filter
h = sqrt(2*pi*B^2/(log(2)))*exp(-t.^2*2*pi^2*B^2/(log(2)));
h=h/sum(h);


