load('monkeydata_training.mat')
training_data = trial;
%% Tune parameters for SVM-Gaussian
window = (50:25:500)';
stdev = (50:25:500)';

kfold_SVM = nan(size(window,1),size(stdev,1));
time_SVM = nan(size(window,1),size(stdev,1));

for win= 1:size(window,1)
    fprintf('%i-',win)
    for sts=1:size(stdev,1)
        tic
        kfold_SVM(win,sts) = SVM_gauss_kfold(training_data,window(win),stdev(sts));
        time_SVM(win,sts) = toc;
    end
end

%% Tune parameters for SVM and CNB with dell-cells
ths = (0.1:0.1:5)'./1000;

kfold_CNB = nan(size(ths,1),50);
time_CNB = nan(size(ths,1),50);

kfold_dSVM = nan(size(ths,1),50);
time_dSVM = nan(size(ths,1),50);


for th=1:size(ths,1)
    fprintf('%i-',th)
    for iter =1:size(time_dSVM,2)
        tic
        kfold_dSVM(th,iter) = SVM_del_kfold(training_data,ths(th));
        time_dSVM(th,iter) = toc;

        tic
        kfold_CNB(th,iter) = CNB_kfold(training_data,ths(th));
        time_CNB(th,iter) = toc;
    end
end
% Save
save('time_CNB.mat', 'time_CNB')
save('time_dSVM.mat', 'time_dSVM')
save('kfold_CNB.mat', 'kfold_CNB')
save('kfold_dSVM.mat', 'kfold_dSVM')

% Plot

figure
hold on
yyaxis left


plotFill1(ths*1000, 100*kfold_CNB,[.1 .6 .8], 'sem',50,'-')
plotFill1(ths*1000, 100*kfold_dSVM,[.1 .6 .8], 'sem',50,'--')
xlabel('Firing Rate Threshold [Hz]', 'FontSize',16)
ylabel('Classification Error [%]', 'FontSize',16)

L(1) = plot(nan, nan, 'k-');
L(2) = plot(nan, nan, 'k--');
yyaxis right

plotFill1(ths*1000, time_CNB,[1 .75 .75], 'sem',50,'-')
plotFill1(ths*1000, time_dSVM,[1 .75 .75], 'sem',50,'--')
ylabel('Time [s]', 'FontSize',16)

yyaxis left
h = legend(L, {'CNB', 'dSVM'});
h.String(1:2) = [];
set(gca,'FontSize',20)
set(gca,'TickDir','out')
set(gca,'Box','off')


        