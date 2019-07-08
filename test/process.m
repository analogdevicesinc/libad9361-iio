chipA = VarName2;
chipB = VarName3;
freq = VarName1(1:10:end);

chipA_mean = mean(reshape(chipA,10,numel(chipA)/10),1);
chipA_var = var(reshape(chipA,10,numel(chipA)/10),1);
chipB_mean = mean(reshape(chipB,10,numel(chipB)/10),1);
chipB_var = var(reshape(chipB,10,numel(chipB)/10),1);

% error
% % plot(freq,chipA_mean,freq,chipB_mean)
errorbar(freq/1e6,chipA_mean,chipA_var); hold on;
errorbar(freq/1e6,chipB_mean,chipB_var); hold off;
xlabel('LO Frequency (MHz)');
ylabel('Phase Difference (Degrees)');
grid on;
legend('Same Chip','Across Chips','Location','best');