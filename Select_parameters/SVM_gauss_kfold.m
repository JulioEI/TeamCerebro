function [kfold] = SVM_gauss_kfold(training_data,window, stdv)  
    % define parameters
    n = size(training_data,1); %number of trials for each direction
    k = size(training_data,2); %number of directions
    c = size(training_data(1,1).spikes,1); %number of cells
    t_interval = 320;
    entry = 0;
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %             FIT MULTICLASS NAIVE BAYESIAN AND SVM MODEL             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % window and stdv is defined based on parameter tuning results
    modelParameters.window = window;
    modelParameters.stdv = stdv;
    X = zeros(k*n, c*t_interval);
    Y = zeros(1, k*n);
    for direc_num = 1:k
        for trial_num = 1:n
            index = (direc_num-1)*n + trial_num;
            entry = entry+1;
            filtered_neuron = gaussian_filter(training_data(trial_num,direc_num).spikes(:,1:t_interval), modelParameters.window, modelParameters.stdv);
            filtered_neuron_all = reshape(filtered_neuron, [1, size(filtered_neuron,1)*t_interval]);
            X(index,:,:) = filtered_neuron_all;
            Y(index) = direc_num;        
        end
    end
    %% TRAIN AND SAVE
    % Fit model
    modelParameters.svmTemplate = ...,
            templateSVM('Standardize',1, ...,
            'KernelFunction','linear');
    modelfit = fitcecoc(X,Y,'Learners',modelParameters.svmTemplate);%this works for predict

    % Save parameters
    kfold = kfoldLoss(crossval(modelfit));
end