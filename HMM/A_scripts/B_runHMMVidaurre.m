%% prepare data

% add HMM-MAR toolbox
addpath(genpath('/Users/kosciessa/BrainHack/HMM/T_tools/HMM-MAR-master/'))

pn.dataDir = '/Users/kosciessa/BrainHack/HMM/B_data/';

fileNames = dir([pn.dataDir, '10subjects_mat/*.mat']);
fileNames = {fileNames(:).name}';

X=[]; T=[]; % data = X
for indData = 1:numel(fileNames)
    loadData = load([pn.dataDir, fileNames{indData}]);
    X{indData,1} = zscore(loadData.ts,[],1); % time x channels (extracted regions); z-score (see Vidaurre et al., 2017 PNAS)
    T{indData,1} = [size(loadData.ts,1)]; 
end

%% set up HMM model

no_states = 12; % the number of states depends a lot on the question at hand
Hz = 1; % the frequency of the data
stochastic_inference = 0; % set to 1 if a normal run is too computationally expensive (memory or time)
N = length(T); % number of subjects

ndim = size(X,2); 

% Setting the options

options = struct();
options.K = no_states;
options.standardise = 1;
options.verbose = 1;
options.Fs = Hz;

options.useParallel=0;

if iscell(T), sumT = 0; for j = 1:N, sumT = sumT + sum(T{j}); end
else, sumT = sum(T); 
end

% Gaussian observation model for fMRI
options.order = 0;
options.zeromean = 0;
options.covtype = 'uniquefull';     

% stochastic options
if stochastic_inference
    options.BIGNbatch = max(round(N/30),5);
    options.BIGtol = 1e-7;
    options.BIGcyc = 500;
    options.BIGundertol_tostop = 5;
    options.BIGforgetrate = 0.7;
    options.BIGbase_weights = 0.9;
end

%% run HMM model

% HMM computation
[hmm, Gamma] = hmmmar(data,T,options);

%% plot results

% Some useful information about the dynamics
maxFO = getMaxFractionalOccupancy(Gamma,T,options); % useful to diagnose if the HMM 
            % is capturing dynamics or grand between-subject 
            % differences (see Wiki)
FO = getFractionalOccupancy (Gamma,T,options); % state fractional occupancies per session
LifeTimes = getStateLifeTimes (Gamma,T,options); % state life times
Intervals = getStateIntervalTimes (Gamma,T,options); % interval times between state visits
SwitchingRate =  getSwitchingRate(Gamma,T,options); % rate of switching between stats

%% save outputs

save([pn.dataDir, 'A_HMMoutput.mat'], 'options', 'hmm', 'Gamma', 'maxFO', 'FO', 'LifeTimes', 'Intervals', 'SwitchingRate')

%% plot results

h = figure; 
subplot(2,2,[1,2]); imagesc(FO)
hold on;
for indSub = 0:4:40
    line([0 12.5], [indSub+.5 indSub+.5],  'Color','k', 'LineWidth', 2)
end
xlabel('Subject'); ylabel('Scan'); title('State Occupancy'); colorbar;
subplot(2,2,3); plot(SwitchingRate)
xlabel('Scan'); ylabel('Switching Rate'); title('Switching Rate');
subplot(2,2,4); plot(maxFO)
xlabel('Scan'); ylabel('maximum Fractional Occupancy'); title('maximum Fractional Occupancy');

pn.plotFolder = '/Users/kosciessa/BrainHack/C_figures/';
figureName = 'A_HMMoutput';

saveas(h, [pn.plotFolder, figureName], 'fig');
saveas(h, [pn.plotFolder, figureName], 'epsc');
saveas(h, [pn.plotFolder, figureName], 'png');