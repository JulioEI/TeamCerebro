% Team Cerebro: Angelica Parra, Zihao Wang, Julio Esparza, Maria Royo
function [modelParameters] = positionEstimatorTraining(training_data,classifier, predictor)
    % define parameters
    n = size(training_data,1); %number of trials for each direction
    k = size(training_data,2); %number of directions
    c = size(training_data(1,1).spikes,1); %number of cells
    t_interval = 320; %time interval to consider for the model (first 320 ms)

%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                              CLASSIFIER                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% WORKING MODE
    if strcmp(classifier, 'SVM_gauss')
        % window and stdv is defined based on parameter tuning results
        windowSize = 280;
        stdv = 500;
        alpha = (windowSize-1) / (2*stdv);
        w = gausswin(windowSize,alpha) / (stdv*sqrt(2*pi));
        X = arrayfun(@(s) reshape(conv2(1,w,s.spikes(:,1:320), 'same'), [1, c*320]), training_data, 'UniformOutput', false);
        X = reshape(X, [k*n,1]);
        X = cell2mat(X);
        Y = [repelem(1,n), repelem(2,n), repelem(3,n), repelem(4,n), repelem(5,n), repelem(6,n), repelem(7,n), repelem(8,n)];
        %% TRAIN AND SAVE
        % Fit model
        modelParameters.svmTemplate = ...,
                templateSVM('Standardize',1, ...,
                'KernelFunction','linear');
        modelfit = fitcecoc(X,Y,'Learners',modelParameters.svmTemplate);%this works for predict
        % Save parameters
        modelParameters.modelfit= modelfit;
        modelParameters.w = w;
    elseif strcmp(classifier, 'CNB')
        X = nan(n*k,c); 
        Y = nan(n*k,1);
        del_cell = zeros(1,c); %cells to exclude from the model due to low rate
        entry = 0;
        for direc_num = 1:k
            st_dir = entry+1;
            for trial_num = 1:n
                entry = entry+1;
                X(entry,:) = mean(training_data(trial_num,direc_num).spikes(:,1:t_interval),2)'; %compute mean activity of each cell during that interval
                Y(entry,1) = direc_num;
            end
            del_cell(1,prctile(X(st_dir:entry,:),80)<0.001) = 1; %delete cells with an 80 percentile activity lower than 1Hz
        end
        X(:,del_cell==1) = []; %delete cells
        %% TRAIN AND SAVE
        % Fit model
        modelfit = fitcnb(X,Y);%this works for predict
        % Save parameters
        modelParameters.modelfit= modelfit;
        modelParameters.del_cell = del_cell;
        modelParameters.t_interval = t_interval;
    elseif strcmp(classifier,'SVM_del')
        entry = 0;
        del_cell = zeros(1,c);
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
                    del_cell(1,prctile(X(st_dir:entry,:),80)<=0.001) = 1; %delete cells with an 80 percentile activity lower than 1Hz
            end  
            X(:,del_cell==1) = []; %delete cells

        %% TRAIN AND SAVE
        % Fit model
        modelParameters.svmTemplate = ...,
                templateSVM('Standardize',1, ...,
                'KernelFunction','linear');
        modelfit = fitcecoc(X,Y,'Learners',modelParameters.svmTemplate);%this works for predict
        % Save parameters
        modelParameters.modelfit= modelfit;
        modelParameters.del_cell = del_cell;
        modelParameters.t_interval = t_interval;
    end
      
 %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                             PREDICTOR                         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(predictor, 'MeanTraj')
        %% ADJUST DATA
        % match size
        state_motion = arrayfun(@(s) [s.handPos(1:2, 320:20:end), s.handPos(1:2,end)], training_data,'UniformOutput',false);
        mean_state_motion = cell(8,1);
        for dir = 1:k 
            max_size = max(cellfun(@(data) ceil(size(data,2)), state_motion(:,dir)))+1;
            state_motion(:,dir) = cellfun(@(data) [data, repmat(data(1:2,end),1,max_size - size(data,2))], state_motion(:,dir), 'UniformOutput', false);
            mean_state_motion{dir,1} = mean(cat(3,state_motion{:,dir}),3);
        end
        modelParameters.mean_state_motion = mean_state_motion;
    elseif strcmp(predictor, 'Kalman')
        %% ADJUST DATA
        % create cell structures, and implement delay into spikes
        lag_idx_delta = 50; % time_lag_options = 0:10:200; Optimal lag is 150~170ms, tested with M1 data
        target_on_idx = 300; %start of motion
        state_motion = arrayfun(@(s) s.handPos(1:2,target_on_idx:end-1), training_data, 'UniformOutput', false);
        spike_motion = arrayfun(@(s) s.spikes(:,target_on_idx-lag_idx_delta:end-1-lag_idx_delta),training_data,'UniformOutput',false);
        % down sample to bins of bin_size and compute frequency
        bin_size = 10;
        comp_len = @(s,bin_size) bin_size*fix(size(s,2)/bin_size);
        state_motion = cellfun(@(data) (squeeze(mean(reshape(data(:, 1:comp_len(data, bin_size)), ... 
            [2, bin_size, comp_len(data, bin_size)/bin_size]),2))'), state_motion,'UniformOutput',false);
        spike_motion = cellfun(@(data) sqrt(squeeze(sum(reshape(data(:, 1:comp_len(data, bin_size)), ...
            [c, bin_size, comp_len(data, bin_size)/bin_size]),2))'), spike_motion,'UniformOutput',false);
        % add velocity to handPos 
        state_motion = cellfun(@(data) [data, cat(1, diff(data,1,1),[0,0])], state_motion, 'UniformOutput', false);
        % Reduce spike count to 50 dimensions, explaining approx. 85% of variability
        PCA_coeff = cell(k,1);
        for dir= 1:k
            PCA_coeff{dir,1} = pca(vertcat(spike_motion{:,dir}), 'NumComponents', 50); 
            spike_motion(:,dir) = cellfun(@(data) data*PCA_coeff{dir,1}, spike_motion(:,dir), 'UniformOutput', false);
        end

        A = cell(k,1);
        Q = cell(k,1);
        pi_0 = cell(k,1);
        V = cell(k,1);
        C = cell(k,1);
        R = cell(k,1);
        %% TRAIN AND SAVE
        % Train Kalman filter
        for dir = 1:k
            total_num_of_timestamps = sum(cellfun(@(x) size(x,1), state_motion(:,dir)));
            %compute autocorrelation for each state on each trial (except for the last entry due for dimensionality consistency with state_onestep_corr
            state_auto_corr = cellfun(@(data) data(1:end-1,:)'*data(1:end-1,:), state_motion(:,dir), 'UniformOutput', false);
            %compute autocorrelation between past/present for each state on each trial 
            state_onestep_corr = cellfun(@(data) data(2:end,:)'*data(1:end-1,:), state_motion(:,dir), 'UniformOutput', false);
            %compute full autocorrelation for each state on each trial(including last entry)
            state_auto_corr_full = cellfun(@(data) data'*data, state_motion(:,dir), 'UniformOutput', false);
            %compute full cross-correlation for each state and its spike activity on each trial
            state_rate_crosscorr = cellfun(@(data_state, data_spike) data_spike'*data_state, state_motion(:,dir), spike_motion(:,dir), 'UniformOutput', false);
            %get initial state for all trials
            starting_state = cellfun(@(data) data(1,:), state_motion(:,dir), 'UniformOutput', false);
            %compute A (transformation matrix of state model, state past to present) using the sum of state_auto_corr of all trials
            A(dir,1) = {sum(cat(3,state_onestep_corr{:}),3)/sum(cat(3,state_auto_corr{:}),3)};
            %compute C (transformation matrix of spike to state model) using the sum of cross_correlations across trials 
            C(dir,1) = {sum(cat(3, state_rate_crosscorr{:}),3)/sum(cat(3, state_auto_corr_full{:}),3)}; 
            %compute mean starting position
            pi_0(dir,1) = {mean(cat(3,starting_state{:}),3)'};
            %compute V (covariance of initial state)
            V(dir,1) = {cov(squeeze(cat(3,starting_state{:}))',1)};
            %compute Q (noise covariance of state) as the mean of the covariance of the difference between the actual state and the predicted one across trials
            Q_temp = cellfun(@(data) (data(2:end,:)-data(1:end-1,:)*A{dir,1}')'*(data(2:end,:)-data(1:end-1,:)*A{dir,1}'), state_motion(:,dir), 'UniformOutput', false);
            Q(dir,1) = {sum(cat(3,Q_temp{:}), 3)./(total_num_of_timestamps-size(state_motion,1))};
            %compute R (noise covariance of neural activity)
            R_temp = cellfun(@(data_state, data_spike) (data_spike-data_state*C{dir,1}')'*(data_spike-data_state*C{dir,1}'), state_motion(:,dir), spike_motion(:,dir),'UniformOutput',false);
            R(dir,1) = {sum(cat(3,R_temp{:}), 3)./total_num_of_timestamps};
        end
        modelParameters.Kalman.A = A;
        modelParameters.Kalman.Q = Q;
        modelParameters.Kalman.pi_0 = pi_0;
        modelParameters.Kalman.V = V;
        modelParameters.Kalman.C = C;
        modelParameters.Kalman.R = R;
        modelParameters.Kalman.PCA_coeff = PCA_coeff;
        modelParameters.Kalman.bin_size = bin_size;
        modelParameters.Kalman.lag_idx_delta = lag_idx_delta;
    end
end