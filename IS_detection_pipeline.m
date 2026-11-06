%% INTERICTAL SPIKE DETECTION PIPELINE
%   PURPOSE:
%       This script serves as a guide/pipline for tuning the interictal
%       spike detector, as used in Yi et al. J. Neuro 2025. The .mat file
%       iteritively saved in each day folder adds information throughout
%       the process used in this pipeline. Each section can be run
%       individually to move step-wise through the process. The sections do
%       the following:
%           (1) Identify 4 random 3 minute segments from each monitoring
%               session. Calculate ripple power across all channels for
%               those time blocks. Select channel with highest power.
%           (2) Manually label ISs using the SignalLabeler GUI within the
%               time segments.
%           (3) Loop through and plot the manually labeled ISs to visually
%               check that everything was done and saved correctly.
%           (4) Loop and plot data from one of the four time bins at
%               different high pass filters. Input a starting point for a
%               detection threshold.
%           (5) Loop high pass filters and thresholds. Calculate TPs, FPs,
%               FNs, Precision, Recall, F05, and F1 for each combination.
%           (6) Loop and check if the highest or lowest threshold used
%               yielded the highest F05 value. Continue with more
%               thresholds until F05 maxes out, if necessary.
%           (7) Detect ISs from the entire monitoring sessions using the
%               optimized high pass frequency and threshold.
%                   (7.5) Detect from combined .bin files for days when
%                         monitoring was split into two folders.
%           (8) Visualize randomly selected time bins with detected events
%               as a check.
%                   (8.5) After checking, two days had to be re-run with a
%                   manually set skew direction (i.e., it was attempting to
%                   detect positive-going interictal spikes even though the
%                   spikes were negative-going, or vice versa). I then went
%                   back to re-run steps 5-8 for those days.
% 
% MMD
% Ewell Lab
% 08-11/2026

%% SET UP - MUST RUN THIS FIRST!

pyenv(Version = "C:\Users\megmd\miniconda3\envs\MATLAB_NEO\python.exe"); % EDIT THIS AS NEEDED - THIS IS THE PYTHON ENVIRONMENT YOU HAVE CREATED TO RUN THIS IN

%% (1) IDENTIFY RANDOM SEGMENTS AND GET HIGHEST RIPPLE POWER

clear

%parameters for reading in data
Fs = 3000;
ds = 30000 ./ Fs; % fs = 30 kHz - ds factor of 10 -> 3000 Hz
channelNums = -1; %-1 for all channels; otherwise  py.list([channel numbers]), starting at 0! not 1
% channelNums = py.list([0, 1]);
startInd = 0; %start frame, in terms of fs/ds sampling rate; start at time 0
% endTm = 5*60; %time in seconds
% endInd = fs .* endTm;

%filter
fMin = 100; %Hz
fMax = 400;

s = RandStream('mt19937ar');

dataDir = 'D:\AmygdalaMice';
neoDir = 'C:\Users\megmd\MATLAB_NEO'; %THIS IS THE CLONED FOLDER FROM GITHUB
mouseLog = readtable([dataDir '\mouse_masterlog.xlsx']);
kaInds = strmatch('KA', mouseLog.Group);
cd(neoDir) %need to be for python reasons

for m = kaInds'
    fprintf('%s\n', mouseLog.Number{m})
    % cd([dataDir '\' mouseLog.Number{m}])
    delDirs = dir([dataDir '\' mouseLog.Number{m} '\Delay*']); %some mice did not complete all 5 days (died /:)

    for d = 1:length(delDirs)
        fprintf('\t%s\n', delDirs(d).name)
        monDirs = dir([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name '\*onitor*']);
        monDirs(strmatch('Post', vertcat({monDirs(:).name}))) = []; %get rid of any post-behavior monitoring to keep things consistent
        if length(monDirs) > 1
            if m == 2 && d == 1 || m == 2 && d == 2 || m == 6 && d == 1 || m == 6 && d == 2 || m == 6 && d == 3
                if isfile([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                        '\' monDirs(1).name '\ISdetectWorkFile.mat'])
                    continue %to the next day
                end %if file already exists

                % have to do this in a slightly different order
                dur = zeros(1,2);
                for mon = 1:2
                    vidObj = VideoReader([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                        '\' monDirs(mon).name '\video.avi']); %for checking monitoring duration
                    dur(mon) = vidObj.Duration;
                    clear vidObj
                end %part

                durBins = floor(dur ./ (3*60));
                durBinsCum = [durBins(1) durBins(1)+durBins(2)]; %cumulative
                potBins = sum(durBins); %number of 3 min bins
                useBins = sort(randsample(s, potBins, 4));

                avgPow = zeros(4, 64); %for storing
                for b = 1:4
                    mon = find(useBins(b)<durBinsCum, 1, 'First');
                    binFile = dir([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                        '\' monDirs(mon).name '\ampData*.bin']);
                    binFile = [dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                        '\' monDirs(mon).name '\' binFile.name]; %convert to full

                    startTm = (useBins(b)-(mon-1)*durBinsCum(1) -1)*3*60; %in second
                    startInd = startTm * Fs; %in samples
                    endTm = startTm + 3*60;
                    endInd = endTm * Fs; %in samples

                    r = pyrunfile("python_code\load_bonsai_binary.py", "traces", ...
                        binary_file = binFile, ... %'path\to\file.bin', ...
                        ds = ds, ... % downsampling factor
                        channels = -1, ... % which channels to get, if you want all set channels = -1. NOTE, starts AT ZERO, ends at 64
                        start_frame = startInd, ... % start frame, in terms of fs/ds sampling rate
                        stop_frame = endInd); % end frame, in terms of fs/ds sampling rate


                    r_use_allCh = double(r); % the actual data to use
                    ts = ((startTm:size(r_use_allCh,1))/(30000/ds))';

                    for ch = 1:64
                        bpLfp = bandpass(r_use_allCh(:,ch), [fMin, fMax], Fs);
                        mhLfp = abs(hilbert(bpLfp)); %magnitude of the hilbert transform
                        avgPow(b,ch) = mean(mhLfp); %average - using same time bins for each channel
                    end %64 channels

                end %bin

                [~, maxCh] = max(mean(avgPow));

                for mon = 1:2
                    save([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                        '\' monDirs(mon).name '\ISdetectWorkFile.mat'], 'useBins', 'maxCh'); %save the time bins used and the channel w max ripple pow
                end %mon

            else
                keyboard
            end %monitoring is split into two on this day
        else %only one monitoring folder
            mon = 1;

            if isfile([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                    '\' monDirs(mon).name '\ISdetectWorkFile.mat'])
                continue %to the next day
            end %if file already exists

            binFile = dir([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                '\' monDirs(mon).name '\ampData*.bin']);
            binFile(vertcat(binFile(:).bytes)==0) = [];
            if length(binFile) > 1
                keyboard
            end %check
            binFile = [dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                '\' monDirs(mon).name '\' binFile.name]; %convert to full

            vidObj = VideoReader([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                '\' monDirs(mon).name '\video.avi']); %for checking monitoring duration
            dur = vidObj.Duration;
            clear vidObj

            potBins = floor(dur ./ (3*60)); %number of 3 min bins
            try
                useBins = sort(randsample(s, potBins, 4));
            catch
                if m == 13 && d == 2
                    dur = 2*60*60;
                    potBins = floor(dur ./ (3*60)); %number of 3 min bins
                    useBins = sort(randsample(s, potBins, 4));
                else
                    keyboard
                end
            end
            avgPow = zeros(4, 64); %for storing
            for b = 1:4
                startTm = (useBins(b)-1)*3*60; %in second
                startInd = startTm * Fs; %in samples
                endTm = startTm + 3*60;
                endInd = endTm * Fs; %in samples

                r = pyrunfile("python_code\load_bonsai_binary.py", "traces", ...
                    binary_file = binFile, ... %'path\to\file.bin', ...
                    ds = ds, ... % downsampling factor
                    channels = channelNums, ... % which channels to get, if you want all set channels = -1. NOTE, starts AT ZERO, ends at 64
                    start_frame = startInd, ... % start frame, in terms of fs/ds sampling rate
                    stop_frame = endInd); % end frame, in terms of fs/ds sampling rate

                r_use_allCh = double(r); % the actual data to use
                ts = ((startTm:size(r_use_allCh,1))/(30000/ds))';

                for ch = 1:64
                    bpLfp = bandpass(r_use_allCh(:,ch), [fMin, fMax], Fs);
                    mhLfp = abs(hilbert(bpLfp)); %magnitude of the hilbert transform
                    avgPow(b,ch) = mean(mhLfp); %average - using same time bins for each channel
                end %64 channels

            end %bins

            [~, maxCh] = max(mean(avgPow));

            % visual inspection
            % figure;
            % addAmt = 0:2000:63000*2;
            % reData = r_use_allCh + addAmt;
            % plot(reData, 'k')
            % hold on
            % plot(reData(:,maxCh), 'r')

            save([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                '\' monDirs(mon).name '\ISdetectWorkFile.mat'], 'useBins', 'maxCh'); %save the time bins used and the channel w max ripple pow
        end %more than one monitoring folder
    end %day

end %ka mice


%% (2) GET SIGNALS FOR SIGNAL LABELER

clear

dataDir = 'D:\AmygdalaMice';
neoDir = 'C:\Users\megmd\MATLAB_NEO';
mouseLog = readtable([dataDir '\mouse_masterlog.xlsx']);
kaInds = strmatch('KA', mouseLog.Group);

%parameters for reading in data
Fs = 3000;
ds = 30000 ./ Fs; % fs = 30 kHz - ds factor of 10 -> 3000 Hz
cd(neoDir) %need to be for python reasons

message = 'Please load "rawLfp" into signal labeler, label interictal spikes, and export labeled signals as "manIS"';
title = '';

for m = kaInds'
    fprintf('%s\n', mouseLog.Number{m})
    % cd([dataDir '\' mouseLog.Number{m}])
    delDirs = dir([dataDir '\' mouseLog.Number{m} '\Delay*']); %some mice did not complete all 5 days (died /:)

    for d = 1:length(delDirs)
        fprintf('\t%s\n', delDirs(d).name)
        monDirs = dir([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name '\*onitor*']);
        monDirs(strmatch('Post', vertcat({monDirs(:).name}))) = []; %get rid of any post-behavior monitoring to keep things consistent

        if m == 2 && d == 1 || m == 2 && d == 2 || m == 6 && d == 1

            load([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                '\' monDirs(1).name '\ISdetectWorkFile.mat']); %load in the whole thing so we can check if we already did this part!

            if exist('labTms')
                clear labTms
                continue
            end %we already did this day

            % have to do this in a slightly different order
            dur = zeros(1,2);
            for mon = 1:2
                vidObj = VideoReader([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                    '\' monDirs(mon).name '\video.avi']); %for checking monitoring duration
                dur(mon) = vidObj.Duration;
                clear vidObj
            end %part

            durBins = floor(dur ./ (3*60));
            durBinsCum = [durBins(1) durBins(1)+durBins(2)]; %cumulative

            binLfp = [];
            labInds = cell(4,1);
            labTms = cell(4,1);
            for b = 1:4
                mon = find(useBins(b)<durBinsCum, 1, 'First');
                binFile = dir([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                    '\' monDirs(mon).name '\ampData*.bin']);
                binFile = [dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                    '\' monDirs(mon).name '\' binFile.name]; %convert to full

                startTm = (useBins(b)-(mon-1)*durBinsCum(1) -1)*3*60; %in second
                startInd = startTm * Fs; %in samples
                endTm = startTm + 3*60;
                endInd = endTm * Fs; %in samples

                r = pyrunfile("python_code\load_bonsai_binary.py", "traces", ...
                    binary_file = binFile, ... %'path\to\file.bin', ...
                    ds = ds, ... % downsampling factor
                    channels = -1, ... % which channels to get, if you want all set channels = -1. NOTE, starts AT ZERO, ends at 64
                    start_frame = startInd, ... % start frame, in terms of fs/ds sampling rate
                    stop_frame = endInd); % end frame, in terms of fs/ds sampling rate

                r_use = double(r); % the actual data to use
                comAvg = mean(r_use,2);
                comAvgLfp = r_use(:,maxCh) - comAvg;
                comAvgAll = r_use - comAvg;
                rawLfp = r_use(:,maxCh);
                binLfp = cat(2, binLfp, rawLfp);

                % visual inspection
                % figure;
                % addAmt = 0:2000:63000*2;
                % reData = r_use + addAmt;
                % plot(reData, 'k')
                % hold on
                % plot(reData(:,maxCh), 'r')
                % hold off

                % label ISs, export ROIs as "manIS"
                signalLabeler
                % f = msgbox(message, title, 'help');
                keyboard
                
                if ~isempty( manIS.Labels.('manIS'){1}.ROILimits)
                    labInds{b} = manIS.Labels.('manIS'){1}.ROILimits;
                    labTms{b} =  manIS.Labels.('manIS'){1}.ROILimits ./ Fs;
                end %anything there
            end %bin

            for mon = 1:2
                save([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                    '\' monDirs(mon).name '\ISdetectWorkFile.mat'], 'useBins', 'maxCh', 'binLfp', 'labInds', 'labTms'); %save the time bins used and the channel w max ripple pow
                clear labInds labTms
            end %mon - save in both

        else
            mon = 1;
                load([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                    '\' monDirs(mon).name '\ISdetectWorkFile.mat']); %load in the whole thing so we can check if we already did this part!
            if exist('labTms')
                clear labTms
                continue
            end %we already did this day

            binFile = dir([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                '\' monDirs(mon).name '\ampData*.bin']);
            binFile(vertcat(binFile(:).bytes)==0) = [];
            if length(binFile) > 1
                keyboard
            end %check
            binFile = [dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                '\' monDirs(mon).name '\' binFile.name]; %convert to full

            binLfp = [];
            labInds = cell(4,1);
            labTms = cell(4,1);
            for b = 1:4
                startTm = (useBins(b)-1)*3*60; %in second
                startInd = startTm * Fs; %in samples
                endTm = startTm + 3*60;
                endInd = endTm * Fs; %in samples

                r = pyrunfile("python_code\load_bonsai_binary.py", "traces", ...
                    binary_file = binFile, ... %'path\to\file.bin', ...
                    ds = ds, ... % downsampling factor
                    channels = -1, ... % which channels to get, if you want all set channels = -1. NOTE, starts AT ZERO, ends at 64
                    start_frame = startInd, ... % start frame, in terms of fs/ds sampling rate
                    stop_frame = endInd); % end frame, in terms of fs/ds sampling rate

                r_use = double(r); % the actual data to use
                comAvg = mean(r_use,2);
                comAvgLfp = r_use(:,maxCh) - comAvg;
                rawLfp = r_use(:,maxCh);
                binLfp = cat(2, binLfp, rawLfp);

                % visual inspection
                % addAmt = 0:2000:63000*2;
                % reData = r_use + addAmt;
                % plot(reData, 'k')
                % hold on
                % plot(reData(:,maxCh), 'r')
                % hold off

                % label ISs, export ROIs as "manIS"
                signalLabeler
                % f = msgbox(message, title, 'help');
                keyboard

                if ~isempty( manIS.Labels.('manIS'){1}.ROILimits)
                    labInds{b} = manIS.Labels.('manIS'){1}.ROILimits;
                    labTms{b} =  manIS.Labels.('manIS'){1}.ROILimits ./ Fs;
                end %anything there
            end %bins

            save([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                '\' monDirs(mon).name '\ISdetectWorkFile.mat'], 'useBins', 'maxCh', 'binLfp', 'labInds', 'labTms'); %save the time bins used and the channel w max ripple pow
            clear labTms labInds

        end %monitoring file is split into two or not

    end %day
end %mouse

%% (3) CHECK MANUALLY LABELED EVENTS
% go through and quick check to make sure everything saved correctly

clear
close all
figure('Position', [76 128 1102 815])

dataDir = 'D:\AmygdalaMice';
neoDir = 'C:\Users\megmd\MATLAB_NEO';
mouseLog = readtable([dataDir '\mouse_masterlog.xlsx']);
kaInds = strmatch('KA', mouseLog.Group);

cd(neoDir) %need to be for python reasons

for m = kaInds'
    fprintf('%s\n', mouseLog.Number{m})
    % cd([dataDir '\' mouseLog.Number{m}])
    delDirs = dir([dataDir '\' mouseLog.Number{m} '\Delay*']); %some mice did not complete all 5 days (died /:)

    for d = 1:length(delDirs)
        fprintf('\t%s\n', delDirs(d).name)
        monDirs = dir([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name '\*onitor*']);
        monDirs(strmatch('Post', vertcat({monDirs(:).name}))) = []; %get rid of any post-behavior monitoring to keep things consistent

        load([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
            '\' monDirs(1).name '\ISdetectWorkFile.mat']); %load in the whole thing so we can check if we already did this part!
        clf
        for b = 1:4
            subplot(4,1,b)
            if ~isempty(labInds{b})
                xline(mean(labInds{b},2), 'r--')
            end
            hold on
            plot(binLfp(:,b), 'k')
            xlim([0 5.4*10^5])

        end %bin
        same_axes;
keyboard
    end %day
end %mouse


%% (4) EXPERIMENT A LITTLE WITH FILTERING THRESHOLDS AND PROMPT FOR THRESHOLD START

clear
close all
figure('Position', [76 128 1102 815])

dataDir = 'D:\AmygdalaMice';
neoDir = 'C:\Users\megmd\MATLAB_NEO';
mouseLog = readtable([dataDir '\mouse_masterlog.xlsx']);
kaInds = strmatch('KA', mouseLog.Group);

fLow = 2;
fHigh = [60 80 100];
Fs = 3000;

prompt = 'Input threshold start';
dlgtitle = 'Input threshold';
fieldsize = [1 45];

cd(neoDir) %need to be for python reasons

% for m = kaInds'
for m = 5
    fprintf('%s\n', mouseLog.Number{m})
    % cd([dataDir '\' mouseLog.Number{m}])
    delDirs = dir([dataDir '\' mouseLog.Number{m} '\Delay*']); %some mice did not complete all 5 days (died /:)

    for d = 1:length(delDirs)
        fprintf('\t%s\n', delDirs(d).name)
        monDirs = dir([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name '\*onitor*']);
        monDirs(strmatch('Post', vertcat({monDirs(:).name}))) = []; %get rid of any post-behavior monitoring to keep things consistent

        load([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
            '\' monDirs(1).name '\ISdetectWorkFile.mat']); %load in the whole thing so we can check if we already did this part!
        if exist('inThresh')
            clear inThresh
            continue
        end %we already did this day
        clf

        [~, b] = max(cellfun(@length, labInds)); %use the bin with the most events as a quick catch

        subplot(4,1,1)
        hold on
        plot(binLfp(:,b), 'k')
        yline([-5000:2500:5000], 'Color', rgb('LightGray'), 'Layer', 'bottom')
        xline(mean(labInds{b},2), 'r--', 'Layer', 'bottom')

        xlim([0 5.4*10^5])
        title('Original')

        for f = 1:3
            bpLfp = bandpass(binLfp(:,b), [fLow, fHigh(f)], Fs);

            subplot(4,1,f+1)

            hold on
            plot(bpLfp, 'k')
            xline(mean(labInds{b},2), 'r--', 'Layer', 'bottom')
            yline([-5000:2500:5000], 'Color', rgb('LightGray'), 'Layer', 'bottom')
            xlim([0 5.4*10^5])
            title(['High pass = ' num2str(fHigh(f)) ' Hz'])
        end %freq band
        same_axes;

        % inThresh = input(prompt);
        inThresh = inputdlg(prompt, dlgtitle, fieldsize);
        inThresh = str2num(inThresh{1});

        % save([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
        %     '\' monDirs(1).name '\ISdetectWorkFile.mat'], 'useBins', 'maxCh', 'binLfp', 'labInds', 'labTms', 'inThresh'); %save the time bins used and the channel w max ripple pow
        clear inThresh
    end %day
end %mouse

%% (5) START LOOPING

clear
cln

dataDir = 'D:\AmygdalaMice';
neoDir = 'C:\Users\megmd\MATLAB_NEO';
mouseLog = readtable([dataDir '\mouse_masterlog.xlsx']);
kaInds = strmatch('KA', mouseLog.Group);

fLow = 2;
fHigh = [60 80 100];
Fs = 3000;
threshAdds = [-750:250:750];
numCombos = length(fHigh)*length(threshAdds);

maxWidth = 0.1; %100 ms

cd(neoDir) %need to be for python reasons

for m = kaInds'
    % fprintf('%s\n', mouseLog.Number{m})
    % cd([dataDir '\' mouseLog.Number{m}])
    delDirs = dir([dataDir '\' mouseLog.Number{m} '\Delay*']); %some mice did not complete all 5 days (died /:)

    for d = 1:length(delDirs)
        % fprintf('\t%s\n', delDirs(d).name)
        monDirs = dir([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name '\*onitor*']);
        monDirs(strmatch('Post', vertcat({monDirs(:).name}))) = []; %get rid of any post-behavior monitoring to keep things consistent

        detArr = []; %initialize the detector array to become the detector table

        load([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
            '\' monDirs(1).name '\ISdetectWorkFile.mat']); %load in the whole thing so we can check if we already did this part!

        if exist('detTbl')
            clear detTbl
            continue
        end %we already did this day

        useThreshs = threshAdds + inThresh; %thresholds we'll try
        w = waitbar(0, 'Looping combinations...', 'Name', [mouseLog.Number{m} ' - ' delDirs(d).name]);
        cCntr = 0;
        for t = 1:length(useThreshs)
            if useThreshs(t) <= 0
                continue
            end %too low
            for f = 1:3 %this is not the fastest way but it's fine - consistency
                cCntr = cCntr + 1;
                waitbar(cCntr/numCombos, w, sprintf('Combo %d of %d', cCntr, numCombos));

                useDets = []; %initialize
                useLabs = [];

                for b = 1:4
                    useLabs = cat(1, useLabs, labInds{b} + (b-1)*Fs*3*60); %add time based on b so we can put everything together in the PercSummary
                    if exist('skewDir')
                        [spkTms, spkInds] = spike_detector_v2(binLfp(:,b), 0:1/Fs:3*60*Fs, Fs, fLow, fHigh(f), useThreshs(t), maxWidth, skewDir);
                    else
                        [spkTms, spkInds] = spike_detector_v2(binLfp(:,b), 0:1/Fs:3*60*Fs, Fs, fLow, fHigh(f), useThreshs(t), maxWidth);
                    end %do we need to manually reset the skew direction
                    useDets = cat(1, useDets, spkInds + (b-1)*Fs*3*60);

                end %bin

                [TP, FP, FN, Precision, Recall, F1, F05] = PrecSummary_v2(useDets, useLabs);

                tmpArr = [fLow fHigh(f) useThreshs(t) TP FP FN Precision Recall F1 F05]; %compile everything
                detArr = [detArr; tmpArr]; %add it in

            end %freqeuncies
        end %thresholds

        detTbl = array2table(detArr, 'VariableNames', {'fLow', 'fHigh', 'thresh', 'TP', 'FP', 'FN', 'Precision', 'Recall', 'F1', 'F05'});
        if exist('skewDir')
            save([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                '\' monDirs(1).name '\ISdetectWorkFile.mat'], 'useBins', 'maxCh', 'binLfp', 'labInds', 'labTms', 'inThresh', 'detTbl', 'skewDir');
            clear skewDir
        else
            save([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                '\' monDirs(1).name '\ISdetectWorkFile.mat'], 'useBins', 'maxCh', 'binLfp', 'labInds', 'labTms', 'inThresh', 'detTbl');
        end %if there is a re-set skew direction
        clear detTbl
        close (w)
    end %day
end %mouse

%% (6) GO BACK AND SEE IF WE ARE ON THE ENDS WITH OUR THRESHOLDS - DO MORE?

clear
cln

dataDir = 'D:\AmygdalaMice';
neoDir = 'C:\Users\megmd\MATLAB_NEO';
mouseLog = readtable([dataDir '\mouse_masterlog.xlsx']);
kaInds = strmatch('KA', mouseLog.Group);

fLow = 2;
fHigh = [60 80 100];
Fs = 3000;

maxWidth = 0.1; %100 ms

cd(neoDir) %need to be for python reasons

for m = kaInds'
    fprintf('%s\n', mouseLog.Number{m})
    % cd([dataDir '\' mouseLog.Number{m}])
    delDirs = dir([dataDir '\' mouseLog.Number{m} '\Delay*']); %some mice did not complete all 5 days (died /:)

    % for d = 1:length(delDirs)
    for d = 2
        fprintf('\t%s\n', delDirs(d).name)
        monDirs = dir([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name '\*onitor*']);
        monDirs(strmatch('Post', vertcat({monDirs(:).name}))) = []; %get rid of any post-behavior monitoring to keep things consistent

        load([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
            '\' monDirs(1).name '\ISdetectWorkFile.mat']); %load in the whole thing so we can check if we already did this part!

        [maxF05, maxInd] = max(detTbl.F05);

        while detTbl.thresh(maxInd) == max(detTbl.thresh)
            detArr = [];
            for f = 1:3 %this is not the fastest way but it's fine - consistency
                useDets = []; %initialize
                useLabs = [];

                for b = 1:4
                    useLabs = cat(1, useLabs, labInds{b} + (b-1)*Fs*3*60); %add time based on b so we can put everything together in the PercSummary
                    if exist('skewDir')
                        [spkTms, spkInds] = spike_detector_v2(binLfp(:,b), 0:1/Fs:3*60*Fs, Fs, fLow, fHigh(f), max(detTbl.thresh)+250, maxWidth, skewDir);
                    else
                        [spkTms, spkInds] = spike_detector_v2(binLfp(:,b), 0:1/Fs:3*60*Fs, Fs, fLow, fHigh(f), max(detTbl.thresh)+250, maxWidth);
                    end %do we need to manually reset the skew direction
                    useDets = cat(1, useDets, spkInds + (b-1)*Fs*3*60);

                end %bin

                [TP, FP, FN, Precision, Recall, F1, F05] = PrecSummary_v2(useDets, useLabs);

                tmpArr = [fLow fHigh(f) max(detTbl.thresh)+250 TP FP FN Precision Recall F1 F05]; %compile everything
                detArr = [detArr; tmpArr]; %add it in

            end %freq

            detTbl = cat(1, detTbl, array2table(detArr, 'VariableNames', {'fLow', 'fHigh', 'thresh', 'TP', 'FP', 'FN', 'Precision', 'Recall', 'F1', 'F05'}));
            [maxF05, maxInd] = max(detTbl.F05);
        end %while max

        while detTbl.thresh(maxInd) == min(detTbl.thresh) && min(detTbl.thresh) > 250
            detArr = [];
            for f = 1:3 %this is not the fastest way but it's fine - consistency
                useDets = []; %initialize
                useLabs = [];

                for b = 1:4
                    useLabs = cat(1, useLabs, labInds{b} + (b-1)*Fs*3*60); %add time based on b so we can put everything together in the PercSummary
                    if exist('skewDir')
                        [spkTms, spkInds] = spike_detector_v2(binLfp(:,b), 0:1/Fs:3*60*Fs, Fs, fLow, fHigh(f), min(detTbl.thresh)-250, maxWidth, skewDir);
                    else
                        [spkTms, spkInds] = spike_detector_v2(binLfp(:,b), 0:1/Fs:3*60*Fs, Fs, fLow, fHigh(f), min(detTbl.thresh)-250, maxWidth);
                    end %set skew dir or not
                    useDets = cat(1, useDets, spkInds + (b-1)*Fs*3*60);

                end %bin

                [TP, FP, FN, Precision, Recall, F1, F05] = PrecSummary_v2(useDets, useLabs);

                tmpArr = [fLow fHigh(f) min(detTbl.thresh)-250 TP FP FN Precision Recall F1 F05]; %compile everything
                detArr = [detArr; tmpArr]; %add it in

            end %freq

            detTbl = cat(1, detTbl, array2table(detArr, 'VariableNames', {'fLow', 'fHigh', 'thresh', 'TP', 'FP', 'FN', 'Precision', 'Recall', 'F1', 'F05'}));
            [maxF05, maxInd] = max(detTbl.F05);
        end %while

        if exist('skewDir')
            save([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                '\' monDirs(1).name '\ISdetectWorkFile.mat'], 'useBins', 'maxCh', 'binLfp', 'labInds', 'labTms', 'inThresh', 'detTbl', 'skewDir');
            clear skewDir
        else
            save([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                '\' monDirs(1).name '\ISdetectWorkFile.mat'], 'useBins', 'maxCh', 'binLfp', 'labInds', 'labTms', 'inThresh', 'detTbl');
        end %if there is a re-set skew direction
        clear detTbl

    end %day
end %mouse

%% (7) DETECT FROM WHOLE SESSIONS

clear
cln

pyenv(Version = "C:\Users\megmd\miniconda3\envs\MATLAB_NEO\python.exe"); % EDIT THIS AS NEEDED
cd(neoDir) %need to be for python reasons

dataDir = 'D:\AmygdalaMice';
neoDir = 'C:\Users\megmd\MATLAB_NEO';
mouseLog = readtable([dataDir '\mouse_masterlog.xlsx']);
kaInds = strmatch('KA', mouseLog.Group);

%parameters for reading in data
Fs = 3000;
ds = 30000 ./ Fs; % fs = 30 kHz - ds factor of 10 -> 3000 Hz
channelNums = -1; %-1 for all channels; otherwise  py.list([channel numbers]), starting at 0! not 1
% channelNums = py.list([0, 1]);
startInd = 0; %start frame, in terms of fs/ds sampling rate; start at time 0

maxWidth = 0.1; %100 ms

cd(neoDir) %need to be for python reasons

for m = kaInds'
    fprintf('%s\n', mouseLog.Number{m})

    delDirs = dir([dataDir '\' mouseLog.Number{m} '\Delay*']); %some mice did not complete all 5 days (died /:)

    for d = 1:length(delDirs)
        fprintf('\t%s\n', delDirs(d).name)
        monDirs = dir([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name '\*onitor*']);
        monDirs(strmatch('Post', vertcat({monDirs(:).name}))) = []; %get rid of any post-behavior monitoring to keep things consistent

        load([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
            '\' monDirs(1).name '\ISdetectWorkFile.mat']);

        for mon = 1:length(monDirs)
            if isfile([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                    '\' monDirs(mon).name '\interictalSpikes.mat'])
                continue %to the next day
            end %if file already exists

            binFile = dir([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                '\' monDirs(mon).name '\ampData*.bin']);
            binFile(vertcat(binFile(:).bytes)==0) = []; %get rid of any erroneous files - started bonsai with no camera input or something
            binFile = [dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                '\' monDirs(mon).name '\' binFile.name]; %convert to full

            vidObj = VideoReader([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                '\' monDirs(mon).name '\video.avi']); %for checking monitoring duration
            dur = vidObj.Duration;
             if m == 13 && d == 2
                    dur = 2*60*60;
             end
            % dur = min([dur 2*60*60]); %max out at 2 hours so we aren't getting excessive amounts of data and to keep things sort of consistent
            clear vidObj

            channelNums = py.list([maxCh-1 0]);

            startTm = 0;
            startInd = startTm * Fs + 1; %in samples
            endTm = dur;
            endInd = endTm*Fs;
            worked = 0; %initialize

            while worked == 0
                try
                    r = pyrunfile("python_code\load_bonsai_binary.py", "traces", ...
                        binary_file = binFile, ... %'path\to\file.bin', ...
                        ds = ds, ... % downsampling factor
                        channels = channelNums, ... % which channels to get, if you want all set channels = -1. NOTE, starts AT ZERO, ends at 64
                        start_frame = startInd, ... % start frame, in terms of fs/ds sampling rate
                        stop_frame = endInd); % end frame, in terms of fs/ds sampling rate
                    worked = 1;
                catch
                    endInd = endInd - Fs; %take off a second at a time so we don't accidentally go past the end of the lfp file; different sampling rates/data for lfp and vid
                end %try
            end %while loop

            r_use = double(r);
            lfpData = r_use(:,1);
            ts = ((startTm:size(r_use,1))/(30000/ds))';
            monDur = ts(end); %to be used for calculating rates in the future
            clear r r_use

            [maxF05, maxInd] = max(detTbl.F05);
            if exist('skewDir')
                [spkTms, spkInds] = spike_detector_v2(lfpData, ts, Fs, detTbl.fLow(maxInd), detTbl.fHigh(maxInd), detTbl.thresh(maxInd), maxWidth, skewDir);
                clear skewDir
            else
                [spkTms, spkInds] = spike_detector_v2(lfpData, ts, Fs, detTbl.fLow(maxInd), detTbl.fHigh(maxInd), detTbl.thresh(maxInd), maxWidth);
            end %do we need to manually reset the skew direction
            clear lfpData %goodbye to this massive variable

            isis = diff(spkTms); %get inter-spike intervals
            potBirds = find(isis<2);
            bursts = {};
                bursts{1} =[potBirds(1) potBirds(1)+1]; %initialize

            for s = 2:length(potBirds)
                if potBirds(s) == bursts{end}(end)
                    bursts{end} = [bursts{end} potBirds(s)+1];
                else
                    bursts{length(bursts)+1} = [potBirds(s) potBirds(s)+1];
                end %add
            end %iss

            birds = {};
            seiz = {};
            for b = 1:length(bursts)
                tmpLen = spkTms(bursts{b}(end)) - spkTms(bursts{b}(1)); %length of this burst
                if tmpLen < 10
                    birds = cat(2, birds, bursts{b});
                else
                    seiz = cat(2, seiz, bursts{b});
                end %bird or scs
            end %bird

            save([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                '\' monDirs(mon).name '\interictalSpikes.mat'], 'spkTms', 'spkInds', 'birds', 'seiz', 'isis', 'monDur'); %save - including monDur tells us how long our lfp data is for this day

        end %monitoring folder - some are split in two

    end %day
end %mouse

%% (7.5) DETECT FROM WHOLE SESSIONS - GO BACK AND CHECK DAYS WHERE MONITORING WAS SPLIT INTO TWO

clear
cln

dataDir = 'D:\AmygdalaMice';
neoDir = 'C:\Users\megmd\MATLAB_NEO';
mouseLog = readtable([dataDir '\mouse_masterlog.xlsx']);
kaInds = strmatch('KA', mouseLog.Group);

%parameters for reading in data
Fs = 3000;
ds = 30000 ./ Fs; % fs = 30 kHz - ds factor of 10 -> 3000 Hz
channelNums = -1; %-1 for all channels; otherwise  py.list([channel numbers]), starting at 0! not 1
% channelNums = py.list([0, 1]);
startInd = 0; %start frame, in terms of fs/ds sampling rate; start at time 0

maxWidth = 0.1; %100 ms

cd(neoDir) %need to be for python reasons

for m = kaInds'
    fprintf('%s\n', mouseLog.Number{m})

    delDirs = dir([dataDir '\' mouseLog.Number{m} '\Delay*']); %some mice did not complete all 5 days (died /:)

    for d = 1:length(delDirs)
        fprintf('\t%s\n', delDirs(d).name)
        monDirs = dir([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name '\*onitor*']);
        monDirs(strmatch('Post', vertcat({monDirs(:).name}))) = []; %get rid of any post-behavior monitoring to keep things consistent

        if length(monDirs) > 1
            mon = strmatch('Combined', vertcat({monDirs(:).name}));
            if isempty(mon)
                keyboard
            end

            %get the approximate length of the combined recordings
            dur = 0;
            loadData = 0;
            for checkMon = setxor(1:length(monDirs), mon)
                vidObj = VideoReader([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                    '\' monDirs(checkMon).name '\video.avi']); %for checking monitoring duration
                dur = dur + vidObj.Duration;
                if loadData == 0
                    load([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                        '\' monDirs(checkMon).name '\ISdetectWorkFile.mat']);
                    loadData = 1;
                end %load necessary data
            end %check for
              clear vidObj

            if isfile([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                    '\' monDirs(mon).name '\interictalSpikes.mat'])
                continue %to the next day
            end %if file already exists

            binFile = dir([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                '\' monDirs(mon).name '\combined_recording.bin']);
            if isempty(binFile)
                keyboard
            end

            channelNums = py.list([maxCh-1 0]);

            startTm = 0;
            startInd = startTm * Fs + 1; %in samples
            endTm = dur;
            endInd = endTm*Fs;
            worked = 0; %initialize

                 binFile = [dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                '\' monDirs(mon).name '\' binFile.name]; %convert to full

            while worked == 0
                try
                    r = pyrunfile("python_code\load_bonsai_binary.py", "traces", ...
                        binary_file = binFile, ... %'path\to\file.bin', ...
                        ds = ds, ... % downsampling factor
                        channels = channelNums, ... % which channels to get, if you want all set channels = -1. NOTE, starts AT ZERO, ends at 64
                        start_frame = startInd, ... % start frame, in terms of fs/ds sampling rate
                        stop_frame = endInd); % end frame, in terms of fs/ds sampling rate
                    worked = 1;
                catch
                    endInd = endInd - Fs; %take off a second at a time so we don't accidentally go past the end of the lfp file; different sampling rates/data for lfp and vid
                end %try
            end %while loop

               r_use = double(r);
            lfpData = r_use(:,1);
            ts = ((startTm:size(r_use,1))/(30000/ds))';
            monDur = ts(end); %to be used for calculating rates in the future
            clear r r_use

              [maxF05, maxInd] = max(detTbl.F05);
            if exist('skewDir')
                [spkTms, spkInds] = spike_detector_v2(lfpData, ts, Fs, detTbl.fLow(maxInd), detTbl.fHigh(maxInd), detTbl.thresh(maxInd), maxWidth, skewDir);
                clear skewDir
            else
                [spkTms, spkInds] = spike_detector_v2(lfpData, ts, Fs, detTbl.fLow(maxInd), detTbl.fHigh(maxInd), detTbl.thresh(maxInd), maxWidth);
            end %do we need to manually reset the skew direction
            clear lfpData %goodbye to this massive variable

            isis = diff(spkTms); %get inter-spike intervals
            potBirds = find(isis<2);
            bursts = {};
                bursts{1} =[potBirds(1) potBirds(1)+1]; %initialize

            for s = 2:length(potBirds)
                if potBirds(s) == bursts{end}(end)
                    bursts{end} = [bursts{end} potBirds(s)+1];
                else
                    bursts{length(bursts)+1} = [potBirds(s) potBirds(s)+1];
                end %add
            end %iss

            birds = {};
            seiz = {};
            for b = 1:length(bursts)
                tmpLen = spkTms(bursts{b}(end)) - spkTms(bursts{b}(1)); %length of this burst
                if tmpLen < 10
                    birds = cat(2, birds, bursts{b});
                else
                    seiz = cat(2, seiz, bursts{b});
                end %bird or scs
            end %bird

            save([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                '\' monDirs(mon).name '\interictalSpikes.mat'], 'spkTms', 'spkInds', 'birds', 'seiz', 'isis', 'monDur'); %save - including monDur tells us how long our lfp data is for this day

        end
    end %day
end %mouse


%% (8) VISUALLY CHECK SOME RANDOM TIME BINS

clear
cln

dataDir = 'D:\AmygdalaMice';
neoDir = 'C:\Users\megmd\MATLAB_NEO';
mouseLog = readtable([dataDir '\mouse_masterlog.xlsx']);
kaInds = strmatch('KA', mouseLog.Group);

%parameters for reading in data
Fs = 3000;
ds = 30000 ./ Fs; % fs = 30 kHz - ds factor of 10 -> 3000 Hz
% channelNums = py.list([0, 1]);
startInd = 0; %start frame, in terms of fs/ds sampling rate; start at time 0

figure;

cd(neoDir) %need to be for python reasons

for m = kaInds'
    fprintf('%s\n', mouseLog.Number{m})

    delDirs = dir([dataDir '\' mouseLog.Number{m} '\Delay*']); %some mice did not complete all 5 days (died /:)

    for d = 1:length(delDirs)
        fprintf('\t%s\n', delDirs(d).name)
        monDirs = dir([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name '\*onitor*']);
        monDirs(strmatch('Post', vertcat({monDirs(:).name}))) = []; %get rid of any post-behavior monitoring to keep things consistent
        mon = 1;

        load([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
            '\' monDirs(mon).name '\ISdetectWorkFile.mat'], 'maxCh');
        channelNums = py.list([maxCh-1 0]);

        binFile = dir([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
            '\' monDirs(mon).name '\ampData*.bin']);
        binFile(vertcat(binFile(:).bytes)==0) = []; %get rid of any erroneous files - started bonsai with no camera input or something
        binFile = [dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
            '\' monDirs(mon).name '\' binFile.name]; %convert to full

        load([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
            '\' monDirs(mon).name '\interictalSpikes.mat'], 'spkInds', 'spkTms', 'monDur');

        potBins = floor(monDur ./ (2*60)); %number of 2 min bins
        useBins = sort(randsample(potBins, 4));

        for b = 1:4
            startTm = (useBins(b)-1)*2*60; %in second
            startInd = startTm * Fs; %in samples
            endTm = startTm + 2*60;
            endInd = endTm * Fs; %in samples

            r = pyrunfile("python_code\load_bonsai_binary.py", "traces", ...
                binary_file = binFile, ... %'path\to\file.bin', ...
                ds = ds, ... % downsampling factor
                channels = channelNums, ... % which channels to get, if you want all set channels = -1. NOTE, starts AT ZERO, ends at 64
                start_frame = startInd, ... % start frame, in terms of fs/ds sampling rate
                stop_frame = endInd); % end frame, in terms of fs/ds sampling rate
            r_use = double(r);
            lfpData = r_use(:,1);

            subplot(4,1,b)
            plot(startTm:1/Fs:endTm-(1/Fs), lfpData, 'k')

            hold on
            xline(spkTms, 'r--', 'Layer', 'bottom')
            xlim([startTm endTm])
        end %time bin
        keyboard

        clf
    end %day
end %mouse

%% (8.5) NEED TO MANUALLY ADD FOR TWO DAYS
% After checking events, two days have the signal skewed negatively but
% events are being detected positively (M467 D1 and D4).

% load('D:\AmygdalaMice\M467\Delay1\Pre-Monitor\ISdetectWorkFile.mat')
% skewDir = -1;
% save('D:\AmygdalaMice\M467\Delay1\Pre-Monitor\ISdetectWorkFile.mat',...
%     'useBins', 'maxCh', 'binLfp', 'labInds', 'labTms', 'inThresh', 'skewDir');
% delete('D:\AmygdalaMice\M467\Delay1\Pre-Monitor\interictalSpikes.mat')
%
% load('D:\AmygdalaMice\M467\Delay3\Premonitor\ISdetectWorkFile.mat')
% skewDir = -1;
% save('D:\AmygdalaMice\M467\Delay3\Premonitor\ISdetectWorkFile.mat',...
%     'useBins', 'maxCh', 'binLfp', 'labInds', 'labTms', 'inThresh', 'skewDir');
% delete('D:\AmygdalaMice\M467\Delay3\Premonitor\interictalSpikes.mat')

% Now, need to go back and repeat section 5-8



%% VISUALIZE PIPELINE

clear
cln

%SET THESE TO CHANGE THE MOUSE AND DAY USED FOR THE EXAMPLES
m = 2;
d = 3;

saveDir = 'D:\RESULTS\pipeline';

dataDir = 'D:\AmygdalaMice';
neoDir = 'C:\Users\megmd\MATLAB_NEO';
mouseLog = readtable([dataDir '\mouse_masterlog.xlsx']);

mouseInfo = [mouseLog.Number{m} ' - Delay ' num2str(d)]; %annotated on the bottom of the figure - to keep track of where the data is coming from
saveDir = [saveDir '\' mouseLog.Number{m} '_Delay' num2str(d)];
if ~isfolder(saveDir)
    mkdir(saveDir)
end %check if we have a subfolder for this mouse

monDir = dir([dataDir '\' mouseLog.Number{m} '\Delay' num2str(d) '\*onitor*']);
if length(monDir) > 1
    warning('Monitoring split into multiple subfolders')
end %check
load([dataDir '\' mouseLog.Number{m} '\Delay' num2str(d) '\' monDir.name '\ISdetectWorkFile.mat']);
load([dataDir '\' mouseLog.Number{m} '\Delay' num2str(d)...
    '\' monDir.name '\interictalSpikes.mat'], 'spkTms', 'spkInds', 'monDur', 'birds', 'seiz'); %load

%parameters for reading in data
Fs = 3000;
ds = 30000 ./ Fs; % fs = 30 kHz - ds factor of 10 -> 3000 Hz
channelNums = -1; %-1 for all channels; otherwise  py.list([channel numbers]), starting at 0! not 1
startInd = 0; %start frame, in terms of fs/ds sampling rate; start at time 0
cd(neoDir) %need to be for python reasons

binFile = dir([dataDir '\' mouseLog.Number{m} '\Delay' num2str(d)...
    '\' monDir.name '\ampData*.bin']);
binFile = [dataDir '\' mouseLog.Number{m} '\Delay' num2str(d)...
    '\' monDir.name '\' binFile.name]; %convert to full

b = 3; %see how it looks
startTm = (useBins(b)-1)*3*60; %in second
startInd = startTm * Fs; %in samples
endTm = startTm + 3*60;
endInd = endTm * Fs; %in samples

r = pyrunfile("python_code\load_bonsai_binary.py", "traces", ...
    binary_file = binFile, ... %'path\to\file.bin', ...
    ds = ds, ... % downsampling factor
    channels = -1, ... % which channels to get, if you want all set channels = -1. NOTE, starts AT ZERO, ends at 64
    start_frame = startInd, ... % start frame, in terms of fs/ds sampling rate
    stop_frame = endInd); % end frame, in terms of fs/ds sampling rate

r_use_allCh = double(r); % the actual data to use
ts = startTm:1/Fs:endTm-(1/Fs);

figtitle = '1_detect_maxRipPowerCh';
figure('Name', figtitle, 'Position', [530 166 793 642])
addAmt = 0:2000:63000*2;
reData = r_use_allCh + addAmt;
plot(ts, reData, 'k')
hold on
plot(ts, reData(:,maxCh), 'r')

%adding scale bars
line([startTm+60 startTm+70], [-1e4 -1e4], 'Color', 'k')
text(startTm+60, -1.25e4, '10 s')
line([startTm+60 startTm+60], [-1e4 -1e4+5000], 'Color', 'k')
text(startTm+61, -0.75e4, '5 mV')
axis off

% line([startTm+20 startTm+25], [-1e4 -1e4], 'Color', 'k')
% text(startTm+20, -1.25e4, '5 s')
% line([startTm+20 startTm+20], [-1e4 -1e4+1000], 'Color', 'k')
% text(startTm+21, -0.75e4, '1 mV')

annotation('textbox', [0 0 0.1 0.1], 'string', mouseInfo, 'EdgeColor', 'None')

saveas(gcf, [saveDir '\' figtitle], 'epsc')
saveas(gcf, [saveDir '\' figtitle], 'fig')

figtitle = '2_labelIS';
figure('Name', figtitle, 'Position', [530 166 793 642])

for b = 1:4
    subplot(4,1,b)

    for lb = 1:size(labInds{b},1)
        x = [labInds{b}(lb,:) flip(labInds{b}(lb,:))];
        y = [-6000 -6000 6000 6000];
        patch(x, y, rgb('Orange'), 'EdgeColor', rgb('Orange'))
    end %labeled

    hold on
    startTm = (useBins(b)-1)*3*60; %in second
    endTm = startTm + 3*60;
    plot(binLfp(:,b), 'k')

    xlim([0 6e5])

    %scale bars
    line([5.55e5 5.55e6+(Fs*5)], [0 0], 'Color', 'k')
    text(5.55e5 , -1000, '5 s')
    line([5.55e5 5.55e5], [0 5000], 'Color', 'k')
    text(5.575e5, 1000, '5 mV')
    axis off

end %b

same_axes;

annotation('textbox', [0 0 0.1 0.1], 'string', mouseInfo, 'EdgeColor', 'None')

saveas(gcf, [saveDir '\' figtitle], 'epsc')
saveas(gcf, [saveDir '\' figtitle], 'fig')

figtitle = '3_test_filtersThresholds';
figure('Name', figtitle, 'Position', [530 166 793 642])

b = 4;

subplot(4,1,1)

for lb = 1:size(labInds{b},1)
    x = [labInds{b}(lb,:) flip(labInds{b}(lb,:))];
    y = [-6000 -6000 6000 6000];
    patch(x, y, rgb('Orange'), 'EdgeColor', rgb('Orange'))
end %labeled

hold on
startTm = (useBins(b)-1)*3*60; %in second
endTm = startTm + 3*60;
plot(binLfp(:,b), 'k')
xlim([0 5.4*10^5])
xticks(0:Fs*60:5.4e5)
xticklabels(0:3)
xlabel('Time (s)')
ylabel('Voltage (\muV)')
title('Unfiltered signal')

%scale bars
% line([5.55e5 5.55e6+(Fs*5)], [0 0], 'Color', 'k')
% text(5.55e5 , -1000, '5 s')
% line([5.55e5 5.55e5], [0 5000], 'Color', 'k')
% text(5.575e5, 1000, '5 mV')
% axis off

fLow = 2;
fHigh = [60 80 100];
yLines = cat(1, unique(detTbl.thresh(mod(detTbl.thresh(:),500)==0)),...
    -1*unique(detTbl.thresh(mod(detTbl.thresh(:),500)==0)));
for f = 1:3
    bpLfp = bandpass(binLfp(:,b), [fLow, fHigh(f)], Fs);
    skewDir = sign(skewness(bpLfp));

    subplot(4,1,f+1)

    hold on
    plot(bpLfp, 'k')
    % xline(mean(labInds{b},2), 'r--', 'Layer', 'bottom')
    yline(yLines, 'Color', rgb('LightGray'), 'Layer', 'bottom')
    % ylim([-2200 2200])
    ylim([-max(detTbl.thresh)-250 max(detTbl.thresh)+250])
    xlim([0 5.4*10^5])
    xticks(0:Fs*60:5.4e5)
    xticklabels(0:3)
    xlabel('Time (s)')
    ylabel('Voltage (\muV)')
    title(['High pass = ' num2str(fHigh(f)) ' Hz'])

end %frequency

[maxF05, maxInd] = max(detTbl.F05);

subplot(4,1,match(detTbl.fHigh(maxInd), fHigh)+1)
yline(skewDir*detTbl.thresh(maxInd), 'Color', rgb('Crimson'), 'Layer', 'bottom')
yLoc = (match(detTbl.fHigh(maxInd), flip(fHigh)))*0.25;
if skewDir
    yLoc = yLoc + 0.05;
else
    yLoc = yLoc - 0.1;
end
annotation('textbox', [0.9, yLoc-0.05, 0.5, 0], 'string', ['F_1_/_2 = ' num2str(round(maxF05,2))],...
    'EdgeColor', 'None')

annotation('textbox', [0 0 0.1 0.1], 'string', [mouseInfo ' - bin ' num2str(b)], 'EdgeColor', 'None')

saveas(gcf, [saveDir '\' figtitle], 'epsc')
saveas(gcf, [saveDir '\' figtitle], 'fig')

figtitle = '4_detect_categorize_events';
figure('Name', figtitle, 'Position', [530 166 793 642])

for b = 1:4
    subplot(4,1,b)

    startTm = (useBins(b)-1)*3*60; %in second
    endTm = startTm + 3*60;
    plot(binLfp(:,b), 'k')

    hold on
    tmpSpks = spkInds(spkTms >= startTm & spkTms <= endTm);
    tmpSpkTms =  spkTms(spkTms >= startTm & spkTms <= endTm);

    xline(tmpSpks - startTm*Fs + 1,  'r', 'Layer', 'bottom')

    xlim([0 6e5])

    %scale bars
    line([5.55e5 5.55e6+(Fs*5)], [0 0], 'Color', 'k')
    text(5.55e5 , -1000, '5 s')
    line([5.55e5 5.55e5], [0 5000], 'Color', 'k')
    text(5.575e5, 1000, '5 mV')
    axis off

    %find first and last spk in this time bin that is part of a bird
    minSpk = min(intersect(horzcat(birds{:}), match(tmpSpks(1), spkInds):match(tmpSpks(end), spkInds)));
    maxSpk = max(intersect(horzcat(birds{:}), match(tmpSpks(1), spkInds):match(tmpSpks(end), spkInds)));
    if ~isempty(minSpk)
        tmpBirds = find(cellfun(@(x)(ismember(minSpk, x)), birds)):find(cellfun(@(x)(ismember(maxSpk, x)), birds)); %birds in this window

        ax = gca;
        yMax = ax.YLim(end);
        for bd = 1:length(tmpBirds)
            xLims = [spkInds(birds{tmpBirds(bd)}(1)) - startTm*Fs + 1 spkInds(birds{tmpBirds(bd)}(end)) - startTm*Fs + 1];
            yLims = [yMax yMax];
            line([xLims], [yLims], 'Color', rgb('DarkGreen'), 'LineWidth', 10)
        end %birds
    end %anything to plot

    %okay now scs
    minSpk = min(intersect(horzcat(seiz{:}), match(tmpSpks(1), spkInds):match(tmpSpks(end), spkInds)));
    maxSpk = max(intersect(horzcat(seiz{:}), match(tmpSpks(1), spkInds):match(tmpSpks(end), spkInds)));
    if ~isempty(minSpk)
        tmpSeiz = find(cellfun(@(x)(ismember(minSpk, x)), seiz)):find(cellfun(@(x)(ismember(maxSpk, x)), seiz)); %birds in this window

        ax = gca;
        yMax = ax.YLim(end);
        for sz = 1:length(tmpSeiz)
            xLims = [spkInds(seiz{tmpSeiz(sz)}(1)) - startTm*Fs + 1 spkInds(seiz{tmpSeiz(sz)}(end)) - startTm*Fs + 1];
            yLims = [yMax yMax];
            line([xLims], [yLims], 'Color', rgb('DeepPink'), 'LineWidth', 10)
        end %birds
    end %there is anything to plot
end %b

annotation('textbox', [0 0 0.1 0.1], 'string', mouseInfo, 'EdgeColor', 'None')

saveas(gcf, [saveDir '\' figtitle], 'epsc')
saveas(gcf, [saveDir '\' figtitle], 'fig')