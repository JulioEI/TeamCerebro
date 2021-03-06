%load monkeydata0.mat
load monkeydata_training.mat

iters = (1:50:1000)';
%define algorithm combinations
protcs = [{'CNB'}, {'CNB'}, {'SVM_del'}, {'SVM_del'}, {'SVM_gauss'}, {'SVM_gauss'};...
    {'MeanTraj'}, {'Kalman'},{'MeanTraj'}, {'Kalman'},{'MeanTraj'}, {'Kalman'}];
kfold = nan(length(iters), length(protcs));
tim = nan(length(iters), length(protcs));
for iter = 1:length(iters)  %iterate over each train/test data split
    fprintf('%i-',iter)
    rng(iters(iter)); %set random seed for reproducibility
    ix = randperm(length(trial));
    % Split into training and testing data
    trainingData = trial(ix(1:60),:);
    testData = trial(ix(61:end),:);
    
    for ver = 1:length(protcs) %iterate over each algorithm combination as defined in protcs 
        classifier = protcs{1,ver}; %select classifier of the algorithm combination
        predictor = protcs{2,ver}; %select predictor of the algorithm combination
        meanSqError = 0;
        n_predictions = 0;  
        tic % keep track of running time 
        % Train Model
        modelParameters = positionEstimatorTraining(trainingData, 'classifier', classifier, 'predictor',  predictor); %train the model for the algorithm combination 
        for tr=1:size(testData,1) %iterate over each 8 trial batch in the test data
            for direc=randperm(8) %iterate over each direction in each test data batch
                decodedHandPos = [];

                times=320:20:size(testData(tr,direc).spikes,2); %start 20ms after movement onset in 20ms steps

                for t=times %for each time bin decode hand position
                    past_current_trial.trialId = testData(tr,direc).trialId;
                    past_current_trial.spikes = testData(tr,direc).spikes(:,1:t); 
                    past_current_trial.decodedHandPos = decodedHandPos;

                    past_current_trial.startHandPos = testData(tr,direc).handPos(1:2,1); 

                    [decodedPosX, decodedPosY, newParameters] = positionEstimator(past_current_trial, modelParameters,'classifier',classifier, 'predictor',predictor);
                    modelParameters = newParameters;

                    decodedPos = [decodedPosX; decodedPosY];
                    decodedHandPos = [decodedHandPos decodedPos];

                    meanSqError = meanSqError + norm(testData(tr,direc).handPos(1:2,t) - decodedPos)^2;

                end
                n_predictions = n_predictions+length(times);
            end
        end
        RMSE = sqrt(meanSqError/n_predictions); 
        t = toc;
        kfold(iter, ver) = RMSE;
        tim(iter,ver) =t;
    end
end

%% PLOT RESULTS

group_label = [{'CNB+Mean'}, {'CNB+Kalman'}, {'dSVM+Mean'}, {'dSVM+Kalman'}, {'gSVM+Mean'},{'gSVM+Kalman'}];
groups = repmat(group_label, 20, 1);

groups = groups(:);
X = [tim(:),kfold(:)];

figure 
gscatter(X(:,1),X(:,2),groups); 
h = gca; 
cxlim = h.XLim; 
cylim = h.YLim; 
hold on 

Mu = zeros(6,2); % Extract the means 
Sigma = zeros(2,2,6); 
col = 0.5*[1,0,0; 1,1,0;0,1,0;0,1,1;0,0,1;1,0,1];
txt_pos = [4.5,13;7,28;1,12;4,27;19,10;20,26];
for j = 1:6     
    Mu(j,1) = mean(tim(:,j));
    Mu(j,2) = mean(kfold(:,j));
    
    Sigma(:,:,j) = cov(tim(:,j), kfold(:,j))^2;
    xlim = Mu(j,1) + 4*[-1 1]*sqrt(Sigma(1,1,j));     
    ylim = Mu(j,2) + 4*[-1 1]*sqrt(Sigma(2,2,j));     
    f = @(x1,x2)reshape(mvnpdf([x1(:),x2(:)],Mu(j,:),Sigma(:,:,j)),size(x1));     
    fcontour(f,[xlim ylim],'--','LineColor', col(j,:)) % Draw contours for the multivariate normal distributions  
    text(txt_pos(j,1),txt_pos(j,2),group_label{j},'FontSize',14,'Color', col(j,:))
end

h.XLim = cxlim; 
h.YLim = cylim; 
xlabel('Running Time (s)') 
ylabel('RMSE') 
legend 'off'
set(gca,'FontSize',20)
set(gca,'TickDir','out')
set(gca,'Box','off')

