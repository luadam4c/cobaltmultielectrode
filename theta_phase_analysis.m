function [thetaPhaseVec, troughIndices, ...
                spikeIndices, thetaPhases, thetaPhaseCounts] = ...
                theta_phase_analysis(dataTheta, dataHighFreq, ...
                       samplingRate, troughThreshold, spikeThreshold, ...
                       outFolder, fileBase)
%% Analyze the theta phase of high frequency spikes
% File History:
%   2018-08-25 Created

%% Hard-coded parameters
nBins = 32;                 % number of bins in the histogram

%% Preparation
%% Create an output folder for cross correlation figures
figFolder = fullfile(outFolder, 'thetaphase');
if exist(figFolder, 'dir') ~= 7
    mkdir(figFolder);
end

% Count the number of samples
nSamples = size(dataTheta, 1);

% Count the number of channels
nChannels = size(dataTheta, 2);

% Create a time vector in seconds for the entire channel
timeVector = (1:nSamples)' * 1/samplingRate;

% Create a file name for the histogram
histFigName = fullfile(figFolder, [fileBase, '_theta_phase_distribution']);

%% Create a theta phase vector for each channel
thetaPhaseVec = zeros(nSamples, nChannels);
troughIndices = cell(nChannels, 1);
parfor iChannel = 1:nChannels
    % Get the theta-filtered data for this channel
    dataThetaThis = dataTheta(:, iChannel);

    % Negate the data
    dataThetaUpsideDown = -dataThetaThis;
    troughThresholdUpsideDown = -troughThreshold;

    % Find all the troughs of the theta-filtered data
    %   that has amplitude below troughThreshold
    [~, troughIndicesThis] = ...
        findpeaks(dataThetaUpsideDown, ...
                    'MinPeakHeight', troughThresholdUpsideDown);

    % Count the number of troughs
    nTroughs = length(troughIndicesThis);
                
    % Initialize a theta phase vector for all time points
    thetaPhaseThis = zeros(nSamples, 1);
    
    % Assign a theta phase to each time point  
    thetaPhaseThis(1:troughIndicesThis(1)-1) = NaN;
    for iTrough = 1:nTroughs
        % Get the starting index of this period
        idxStart = troughIndicesThis(iTrough);
        
        % Get the ending index of this period
        if iTrough < nTroughs
            idxEnd = troughIndicesThis(iTrough + 1);
        else
            idxEnd = nSamples + 1;
        end

        % Count the number of samples within this period
        nSamplesDiff = idxEnd - idxStart;

        % Assign the theta phase, with the trough being
        %   0 and the next trough being 2*pi
        thetaPhaseThis(idxStart:idxEnd - 1) = ...
            (0:(nSamplesDiff-1)) * 2 * pi / nSamplesDiff;
    end

    % Save in thetaPhaseVec
    thetaPhaseVec(:, iChannel) = thetaPhaseThis;
    troughIndices{iChannel} = troughIndicesThis;
end

%% Detect and assign high frequency spikes for each channel
spikeIndices = cell(nChannels, 1);
thetaPhases = cell(nChannels, 1);
thetaPhaseCounts = cell(nChannels, 1);
parfor iChannel = 1:nChannels
    % Get the high frequency-filtered data for this channel
    dataHighFreqThis = dataHighFreq(:, iChannel);
    
    % Get the theta phase vector for this channel
    thetaPhaseVecThis = thetaPhaseVec(:, iChannel);
    
    % Find the spike in the high frequency data
    [~, spikeIndicesThis] = ...
        findpeaks(dataHighFreqThis, 'MinPeakHeight', spikeThreshold);
    
    % Assign the theta phase to each spike
    thetaPhasesThis = thetaPhaseVecThis(spikeIndicesThis);
    
    % Create histogram counts
    thetaPhaseCountsThis = histcounts(thetaPhasesThis, nBins);

    % Store in cell arrays
    spikeIndices{iChannel} = spikeIndicesThis;
    thetaPhases{iChannel} = thetaPhasesThis;
    thetaPhaseCounts{iChannel} = thetaPhaseCountsThis;
end

%% Plot a histogram
% Set the number of rows and columns
nCols = 2;
nRows = ceil(nChannels / nCols);

% Create histogram plots
h = figure(10101);
clf(h);
for iChannel = 1:nChannels
    axes(iChannel) = subplot(nRows, nCols, iChannel);
    histogram(thetaPhases{iChannel}, nBins, ...
            'Displayname', ['Channel ', num2str(iChannel)], ...
            'Normalization', 'probability');
    legend('location', 'northwest');
    
    % Fix x limits and show only some x tick labels
    xlim([0, 2*pi]);
    xticks([0, pi, 2*pi]);
    xticklabels({'0', '\pi', '2\pi'});
end

% Link the y axes
linkaxes(axes, 'y');

% Create a title
suplabel(['Theta phase distribution for ', fileBase], 't');

% Save the figure
saveas(h, histFigName, 'jpg');

%% Plot the detection
parfor iChannel = 1:nChannels
    % Extract vectors for this channel
    dataThetaThis = dataTheta(:, iChannel);
    dataHighFreqThis = dataHighFreq(:, iChannel);
    thetaPhaseVecThis = thetaPhaseVec(:, iChannel);
    troughIndicesThis = troughIndices{iChannel};
    spikeIndicesThis = spikeIndices{iChannel};
    
    % Get a string for the current channel
    channelStr = [fileBase, ' Channel', num2str(iChannel)];

    % Plot the troughs on the data for visualization
    figName = fullfile(figFolder, [channelStr, '_theta_phase_detection']);
    h = figure('Visible', 'off', 'PaperPosition', [0, 0, 80, 5]);
    clf(h)
    subplot(3, 1, 1); hold on
    plot(timeVector, dataThetaThis, 'b', ...
            'DisplayName', 'Theta-filtered data');
    plot(timeVector(troughIndicesThis), ...
            dataThetaThis(troughIndicesThis), 'r.', ...
            'DisplayName', 'Detected troughs');
    ylabel('Voltage (mV)');

    subplot(3, 1, 2); hold on
    plot(timeVector, dataHighFreqThis, 'b');
    plot(timeVector(spikeIndicesThis), ...
        dataHighFreqThis(spikeIndicesThis), 'g.');
    ylabel('Voltage (mV)');

    subplot(3, 1, 3); hold on
    plot(timeVector, thetaPhaseVecThis, 'b');
    plot(timeVector(spikeIndicesThis), ...
        thetaPhaseVecThis(spikeIndicesThis), 'g.');
    xlabel('Time (s)');
    ylabel('Phase (rad)');
    
    suplabel(['Theta phase assignment for ', channelStr], 't');
%    saveas(h, figName, 'jpg');
%    print(h, figName, '-dpdf', '-fillpage');
    print(h, figName, '-djpeg');
    close(h)
end
