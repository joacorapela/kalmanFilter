
% PLDS toolbox example, with external input
%
% Lars Buesing, Jakob H Macke, 2014
%

clear all
close all

stateInputMemorySecs = 0;
analysisStartTimeSecs = 180;
trainDurSecs = 180;
codeRoot = '~/dev/research/gatsby-swc/gatsby/PLDS/matlabCode/pop_spike_dyn';
% timeSeriesFilename = '../../results/MPV18_2/task_2019-02-06_21-36-35_preprocessing_2019_04_05_15_04_39_ks2_timeSeries_10msBin.mat';
timeSeriesFilename = '../../results/MPV18_2/task_2019-02-06_21-36-35_preprocessing_2019_04_05_15_04_39_ks2_timeSeries.mat';
% resultsFilename = '../../results/MPV18_2/PLDSresults_10msBin.mat';
resultsFilename = '../../results/MPV18_2/PLDSresults.mat';

timeSeries = load(timeSeriesFilename);
timeSeries = timeSeries.timeSeries;
sRate = timeSeries.sRate;
trainSamples = (analysisStartTimeSecs*sRate)+(1:(trainDurSecs*sRate));

oldFolder = cd(codeRoot);
set_path
cd(oldFolder)

dbstop if error

spikeCounts = timeSeries.v1Shaft1SpikeCounts;

stateInputs = buildGoNogoVisualAndLaserInputs(timeSeries.goStim, timeSeries.nogoStim, timeSeries.laserStim, stateInputMemorySecs*sRate);

y = spikeCounts(:,trainSamples);
u = stateInputs(:,trainSamples);

uDim    = size(u,1);
xDim    = 9;
yDim    = size(y, 1);
T       = size(y, 2);
Trials  = 1;
maxIter = 100;
doff    = 0.0;

seqOrig.y = y;
seqOrig.u = u;
seqOrig.T = T;

fprintf('Max spike count:    %i \n', max(vec([seqOrig.y])))
fprintf('Mean spike count:   %d \n', mean(vec([seqOrig.y])))
fprintf('Freq non-zero bin:  %d \n', mean(vec([seqOrig.y])>0.5))

%%%% generate data

% trueparams = PLDSgenerateExample('xDim',xDim,'yDim',yDim,'doff',doff,'uDim',uDim);
% seqOrig    = PLDSsample(trueparams,T,Trials,'yMax',10);
% tp         = trueparams;

%%% fit model

seq    = seqOrig;
params = [];
% important: set flag to use external input
if uDim>0;params.model.notes.useB = true;end

params = PLDSInitialize(seq, xDim, 'NucNormMin', params);
% fprintf('Initial subspace angle:  %d \n', subspace(tp.model.C,params.model.C))

params.model.inferenceHandle = @PLDSLaplaceInference;
params.opts.algorithmic.EMIterations.maxIter     = maxIter;
params.opts.algorithmic.EMIterations.maxCPUTime  = inf;
tic; [params seq varBound EStepTimes MStepTimes] = PopSpikeEM(params,seq); toc
% fprintf('Final subspace angle:  %d \n', subspace(tp.model.C,params.model.C))

%%% save true and estimated models

y = seqOrig.y;
u = seqOrig.u;
save(resultsFilename, 'y', 'u', 'params', 'seq', 'varBound');

%%% check learned parameters

% figure
% plot(vec(tp.model.C*tp.model.B),vec(params.model.C*params.model.B),'xr')

%%% test model by predicting future spike trains 

% Tpred = 200;
% seqPred = PLDSsample(trueparams,Tpred,1);  % sample a test data set

% condRange = [1:100];                       % the time interval to condition on
% predRange = [101:200];                     % the time interval to predict

% predict with learned parameters
% tic; [ypred xpred xpredCov seqInf] = PLDSPredictRange(params,seqPred(1).y,condRange,predRange,'u',seqPred(1).u); toc

% figure;
% subplot(2,1,1)
% imagesc(seqPred.y(:,predRange))
% title('true data')
% subplot(2,1,2)
% imagesc(ypred)
% title('prediction')

cd(oldFolder)

