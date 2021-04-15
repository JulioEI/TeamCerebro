function [kfold] = SVM_del_kfold(training_data,th) 
    % define parameters
    n = size(training_data,1); %number of trials for each direction
    k = size(training_data,2); %number of directions
    c = size(training_data(1,1).spikes,1); %number of cells
    t_interval = 320;
    entry = 0;
    del_cell = zeros(1,c); %cells to exclude from the model due to low rate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %             FIT MULTICLASS NAIVE BAYESIAN AND SVM MODEL             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% DATA CLEANING AND PREPARATION
        X = nan(n*k,c); 
        Y = nan(n*k,1);
        for direc_num = 1:k
            st_dir = entry+1;
            for trial_num = 1:n
                entry = entry+1;
                X(entry,:) = mean(training_data(trial_num,direc_num).spikes(:,1:t_interval),2)'; %compute mean activity of each cell during that interval
                Y(entry,1) = direc_num;    
            end
                del_cell(1,prctile(X(st_dir:entry,:),80)<=th) = 1; %delete cells with an 80 percentile activity lower than 1Hz
        end  
        X(:,del_cell==1) = []; %delete cells
    
    %% TRAIN AND SAVE
    % Fit model
    modelParameters.svmTemplate = ...,
            templateSVM('Standardize',1, ...,
            'KernelFunction','linear');
    modelfit = fitcecoc(X,Y,'Learners',modelParameters.svmTemplate);%this works for predict
    % Save parameters
    kfold = kfoldLoss(crossval(modelfit));
end