% Team Cerebro: Angelica Parra, Zihao Wang, Julio Esparza, Maria Royo
function [x, y, modelParameters] = positionEstimator(test_data, modelParameters, varargin)

    sj = inputParser;
    addParameter(sj,'classifier','SVM_del',@isstr)
    addParameter(sj, 'predictor', 'MeanTraj', @isstr)

    parse(sj,varargin{:})

    classifier = sj.Results.classifier;
    predictor = sj.Results.predictor;
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                             CLASSIFY                                %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    t = size(test_data.spikes,2) - 300;
    if (t==20)
        if strcmp(classifier, 'SVM_gauss')
            filtered_spike = conv2(1,modelParameters.w,test_data.spikes(:,1:320),'same');
            filtered_spike = reshape(filtered_spike,[1, size(test_data.spikes,1)*320]);
            dir = predict(modelParameters.modelfit, filtered_spike);
            modelParameters.dir = dir;
        else
            state_prev = mean(test_data.spikes(:,1:modelParameters.t_interval),2)';
            state_prev(:,modelParameters.del_cell==1) = [];
            dir = predict(modelParameters.modelfit,state_prev); 
            modelParameters.dir = dir;
        end
    end
    dir = modelParameters.dir;
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                              PREDICT                             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(predictor, 'MeanTraj')
        %% ADJUST DATA
        mean_trajectory = modelParameters.mean_state_motion{dir,1};
        t = min([t/20,size(mean_trajectory,2)]);
        x = mean_trajectory(1,t);
        y = mean_trajectory(2,t);
    elseif strcmp(predictor, 'Kalman')
        c = size(test_data.spikes,1); %number of cells
        %% ADJUST DATA
        % create cell structures, and implement delay into spikes
        lag_idx_delta = modelParameters.Kalman.lag_idx_delta; % time_lag_options = 0:10:200; Optimal lag is 150~170ms, tested with M1 data
        target_on_idx = 300; %start of motion
        spike_motion = test_data.spikes(:,target_on_idx-lag_idx_delta:end-1-lag_idx_delta);
        % down sample to bins of bin_size and compute frequency
        bin_size = modelParameters.Kalman.bin_size;
        len = bin_size*fix(size(spike_motion,2)/bin_size);
        spike_motion = sqrt(squeeze(sum(reshape(spike_motion(:, 1:len),[c, bin_size, len/bin_size]),2))');
        % project data
        spike_motion = spike_motion*modelParameters.Kalman.PCA_coeff{dir,1};
        % load kalman parameters
        A = modelParameters.Kalman.A{dir,1};
        Q = modelParameters.Kalman.Q{dir,1};
        C = modelParameters.Kalman.C{dir,1};
        V = modelParameters.Kalman.V{dir,1};
        R = modelParameters.Kalman.R{dir,1};
        pi_0 = modelParameters.Kalman.pi_0{dir,1};
        %get previous state
        state_prev = [test_data.startHandPos, test_data.decodedHandPos];
        if size(state_prev,2)==1
            state_prev = [state_prev;pi_0([3,4])];
            V_prev = zeros(size(V));
        else
            state_prev = [state_prev(:,end); diff(state_prev(:,end-1:end),1,2)./10];
            V_prev = modelParameters.Kalman.V_prev;
        end
        nsteps = floor(20/bin_size);
        for step = 1:nsteps
            Z = spike_motion(end-nsteps+step,:)';
            state_estimate = 1.01*A*state_prev;
            state_onestep_cov = A*V_prev*A' + Q;

            K = (state_onestep_cov*C')/(C*state_onestep_cov*C'+R);
            state_estimate = state_estimate + K*(Z - C*state_estimate);

            %update V and state
            V_prev = (eye(size(V_prev)) - K*C)*state_onestep_cov;
            state_prev = state_estimate;
        end   
        modelParameters.Kalman.V_prev = V_prev;   
        x = state_estimate(1);
        y = state_estimate(2);
        %if dir==1 && (norm(state_estimate)<50)
        %    x = x*cosd(-5) - y*sind(-5);
        %   y = x*sind(-5) + y*cosd(-5);
        %end
    end
end
