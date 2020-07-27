%% Simulation Command Script
close all; clear all; clc;
EbNo = 0:10;
EbNolin = 10.^(EbNo/10);

open_system('Sim_CPM_GMSK');
in = Simulink.SimulationInput('Sim_CPM_GMSK');

BER_GMSK = zeros(length(EbNo),1);
BER_MSK = zeros(length(EbNo),1);
BER_CPMg = zeros(length(EbNo),1);
BER_CPMrc = zeros(length(EbNo),1);

ModelParameterNames = get_param('Sim_CPM_GMSK/AWGN Channel','ObjectParameters');

for i=1:length(EbNo)
value = EbNo(i);
in = in.setBlockParameter('Sim_CPM_GMSK/AWGN Channel','EbNodB','value');
in = in.setBlockParameter('Sim_CPM_GMSK/AWGN Channel1','EbNodB','value');
in = in.setBlockParameter('Sim_CPM_GMSK/AWGN Channel2','EbNodB','value');
in = in.setBlockParameter('Sim_CPM_GMSK/AWGN Channel3','EbNodB','value');
in.applyToModel
sim('Sim_CPM_GMSK')
BER_GMSK(i) = ans.berGMSK(end,1);
BER_MSK(i) = ans.berMSK(end,1);
BER_CPMg(i) = ans.berCPMg(end,1);
BER_CPMrc(i) = ans.berCPMrc(end,1);
end

% Theoretical Calculation
berGMSK = zeros(length(EbNo),1);
berMSK = zeros(length(EbNo),1);

alpha_GMSK = 0.68; 
for j=1:length(EbNo)
berGMSK(j) = 0.5*erfc(sqrt(EbNolin(j)*alpha_GMSK));
end

bermsk = berawgn(EbNo,'msk','off','coherent');

modIndex = 0.5; k=1; M=2;
berCPMrc = zeros(length(EbNo),1);
for t=1:length(EbNo)
berCPMrc(t) = (1-(1/M))*erfc(pi*modIndex*sqrt(EbNolin(t)*(log2(M)/2.957)));
end

figure;
subplot(221)
semilogy(EbNo,BER_MSK,'b*-','LineWidth', 1.5); hold on; semilogy(EbNo,bermsk,'r-','LineWidth', 1.5);
title('BER of MSK'); xlabel('Eb/No (dB)'); ylabel('Probability of error P_e');
legend('simulation', 'theoretical');

subplot(222)
semilogy(EbNo,BER_GMSK,'b*-','LineWidth', 1.5); hold on; semilogy(EbNo,berGMSK,'r-','LineWidth', 1.5);
title('BER of GMSK'); xlabel('Eb/No (dB)'); ylabel('Probability of error P_e');
legend('simulation', 'theoretical');

subplot(223)
semilogy(EbNo,BER_CPMg,'b*-','LineWidth', 1.5); hold on; semilogy(EbNo,berCPMrc,'r-','LineWidth', 1.5); 
title('BER of CPM with Gaussian Pulse'); xlabel('Eb/No (dB)'); ylabel('Probability of error P_e');
legend('simulation', 'theoretical');

subplot(224)
semilogy(EbNo,BER_CPMrc,'b*-','LineWidth', 1.5); hold on; semilogy(EbNo,berCPMrc,'r-','LineWidth', 1.5); 
title('BER of CPM with Raised Cosine Pulse'); xlabel('Eb/No (dB)'); ylabel('Probability of error P_e');
legend('simulation', 'theoretical');









