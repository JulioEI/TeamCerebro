%load monkeydata0.mat
load monkeydata_training.mat

iters = (1:50:1000)';
protcs = [{'CNB'}, {'CNB'}, {'SVM_del'}, {'SVM_del'}, {'SVM_gauss'}, {'SVM_gauss'};...
    {'MeanTraj'}, {'Kalman'},{'MeanTraj'}, {'Kalman'},{'MeanTraj'}, {'Kalman'}];
kfold = nan(length(iters), length(protcs));
tim = nan(length(iters), length(protcs));
for iter = 1:length(iters)
    fprintf('%i-',iter)
    rng(iters(iter));
    ix = randperm(length(trial));
    % Select training and testing data (you can choose to split your data in a different way if you wish)
    trainingData = trial(ix(1:60),:);
    testData = trial(ix(61:end),:);
    
    for ver = 1:length(protcs)
        classifier = protcs{1,ver};
        predictor = protcs{2,ver};
        meanSqError = 0;
        n_predictions = 0;  
        tic
        % Train Model
        modelParameters = positionEstimatorTraining(trainingData, classifier, predictor);
        for tr=1:size(testData,1)
            pause(0.001)
            for direc=randperm(8) 
                decodedHandPos = [];

                times=320:20:size(testData(tr,direc).spikes,2);

                for t=times
                    past_current_trial.trialId = testData(tr,direc).trialId;
                    past_current_trial.spikes = testData(tr,direc).spikes(:,1:t); 
                    past_current_trial.decodedHandPos = decodedHandPos;

                    past_current_trial.startHandPos = testData(tr,direc).handPos(1:2,1); 

                    [decodedPosX, decodedPosY, newParameters] = positionEstimator(past_current_trial, modelParameters,classifier, predictor);
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

