function kaMon_IS_triggered_avg
%
% PURPOSE:
%   Plot the interictal spike-triggered average of putative interneuron and
%   pyramidal cell firing by animal and summed. Separate out plots for all
%   interictal spikes, the initiating spike in a BIRD, and an initiating
%   spike in a seizure.

%% INITIALIZE

mouseLog = readtable('D:\AmygdalaMice\mouse_masterlog.xlsx');
kaInds = strmatch('KA', mouseLog.Group);

mice2remove = {'M471'}; 
keepMice = setxor(mouseLog.Number(kaInds), mice2remove);
[~, keepMiceInd] = intersect(mouseLog.Number, keepMice);

staPyr = []; %initialize
staInt = [];

minSpkBird = 3; %minimum number of spikes to be considered a bird; Laura has suggested that we maybe change this from 2 to higher
birdFunc = @(x)(length(x) >= minSpkBird);
pullBirdSpkFunc = @(x)(x(1));
btaPyr = [];
btaInt = [];

sztaPyr = [];
sztaInt = [];

Fs = 1000; %Hz
timeAround = 0.5; %time in seconds - and + spike time
pullFunc = @(x)(x - timeAround*Fs:x + timeAround*Fs);
ttlNames = {'Pyramidal cells', 'Interneurons'};

saveDir = 'F:\RESULTS';
saveOrNot = 0;

%% GATHER DATA AND PLOT BY MOUSE

for m = kaInds'
    if ~ismember(m, keepMiceInd)
        continue %o next mouse
    end %not a kept mouse
    delDirs = dir(['D:\AmygdalaMice\' mouseLog.Number{m} '\Delay*']); %some mice did not complete all 5 days (died /:)

    mouseStaPyr = []; %initialize in mouse
    mouseStaInt = []; %initialize in mouse

    mouseBtaPyr = []; %initialize in mouse
    mouseBtaInt = []; %initialize in mouse

    mouseSztaPyr = []; %initialize in mouse
    mouseSztaInt = []; %initialize in mouse

    for d = 1:length(delDirs)
        monDirs = dir(['D:\AmygdalaMice\' mouseLog.Number{m} '\' delDirs(d).name '\*onitor*']);
        monDirs(strmatch('Post', vertcat({monDirs(:).name}))) = []; %get rid of any post-behavior monitoring to keep things consistent

        if length(monDirs) > 1
            mon = strmatch('Combined', vertcat({monDirs(:).name}));
            if isempty(mon)
                keyboard
            end
        else
            mon = 1;
        end %set which monitoring folder to go into

        load(['D:\AmygdalaMice\' mouseLog.Number{m} '\' delDirs(d).name '\' monDirs(mon).name...
            '\interictalSpikes.mat'], 'spkTms', 'birds', 'seiz', 'monDur'); %load
        unit = read_in_abr_phy(['D:\AmygdalaMice\' mouseLog.Number{m} '\' delDirs(d).name...
            '\' monDirs(mon).name '\SpikeSort\Pipeline\Phy']);
        if isempty(unit)
            continue %to next day
        end %if no good units

        [spkRstr, timeMat] = make_spike_raster(vertcat({unit(:).spkTms}), [0 monDur], Fs);

        tmpFr = cellfun(@length, vertcat({unit(:).spkTms})) ./ monDur; %for putative pyr vs int identification

        %do this for all interictal spikes
        matchInds = match(spkTms, timeMat);
        pullSBins = arrayfun(pullFunc, matchInds, 'UniformOutput', false);
        pullSBins = vertcat(pullSBins{:});
        [row, ~] = find(pullSBins<1);
        pullSBins(unique(row),:) = []; %delete what would be negative time bins
        [row, ~] = find(pullSBins>size(spkRstr,2));
        pullSBins(unique(row),:) = []; %delete what would be time bins past the recording end

        %now do it for birds
        keepBirds = cellfun(birdFunc, birds); %see what matches the minimum spk set above
        birds = birds(keepBirds);
        matchBInds = match(spkTms(cellfun(pullBirdSpkFunc, birds)), timeMat);
        pullBBins = arrayfun(pullFunc, matchBInds, 'UniformOutput', false);
        pullBBins = vertcat(pullBBins{:});
        [row, ~] = find(pullBBins<1);
        pullBBins(unique(row),:) = []; %delete what would be negative time bins
        [row, ~] = find(pullBBins>size(spkRstr,2));
        pullBBins(unique(row),:) = []; %delete what would be time bins past the recording end

        %now do it for seizures
        if isempty(seiz)
            pullSzBins = [];
        else
            matchSzInds = match(spkTms(cellfun(pullBirdSpkFunc, seiz)), timeMat);
            pullSzBins = arrayfun(pullFunc, matchSzInds, 'UniformOutput', false);
            pullSzBins = vertcat(pullSzBins{:});
            [row, ~] = find(pullSzBins<1);
            pullSzBins(unique(row),:) = []; %delete what would be negative time bins
            [row, ~] = find(pullSzBins>size(spkRstr,2));
            pullSzBins(unique(row),:) = []; %delete what would be time bins past the recording end
        end %if there are seizures

        for u = 1:length(unit)
            uRstr = spkRstr(u,:);
            uSBins = full(uRstr(pullSBins)) * Fs; %convert to Hz
            uSMeanBin = mean(uSBins,1);

            uBBins = full(uRstr(pullBBins)) * Fs; %convert to Hz
            uBMeanBin = mean(uBBins,1);

            uSzBins = full(uRstr(pullSzBins)) * Fs; %convert to Hz
            uSzMeanBin = mean(uSzBins,1);
            if tmpFr(u) > 10
                mouseStaInt = [mouseStaInt; uSMeanBin];
                mouseBtaInt = [mouseBtaInt; uBMeanBin];
                mouseSztaInt = [mouseSztaInt; uSzMeanBin]; %initialize in mouse
            else
                mouseStaPyr = [mouseStaPyr; uSMeanBin];
                mouseBtaPyr = [mouseBtaPyr; uBMeanBin];
                mouseSztaPyr = [mouseSztaPyr; uSzMeanBin]; %initialize in mouse
            end %where to put
        end %unit

    end %DAY

    figtitle = [mouseLog.Number{m} '_InterictalSpikeTriggeredAverages'];
    figure('Name', figtitle, 'Position', [202 344 1142 429])
    for sp = 1:2
        subplot(1,2,sp)
        if sp == 1
            useData = mouseStaPyr;
        else
            useData = mouseStaInt;
        end %grab data

        plot(-timeAround:1/Fs:timeAround, mean(useData))
        xlabel('Time (s)')
        ylabel('Firing rate (Hz)')
        axis square
        title([ttlNames{sp} ' (n = ' num2str(size(useData,1)) ')'])
    end %subplot

    if saveOrNot == 1
        saveas(gcf, [saveDir '\' figtitle], 'epsc')
        saveas(gcf, [saveDir '\' figtitle], 'fig')
        saveas(gcf, [saveDir '\' figtitle], 'png')
    end % save
    close all

    staPyr = [staPyr; mouseStaPyr];
    staInt = [staInt; mouseStaInt];

    btaPyr = [btaPyr; mouseBtaPyr];
    btaInt = [btaInt; mouseBtaInt];

    sztaPyr = [sztaPyr; mouseSztaPyr];
    sztaInt = [sztaInt; mouseSztaInt];

end %mouse


keyboard

%% OVERALL FIG

figtitle = 'InterictalSpikeTriggeredAverages';
    figure('Name', figtitle, 'Position', [202 344 1142 429])
for sp = 1:2
    subplot(1,2,sp)
    if sp == 1
        useData = staPyr;
    else
        useData = staInt;
    end %grab data

    plot(-timeAround:1/Fs:timeAround, mean(useData))
    xlabel('Time (s)')
    ylabel('Firing rate (Hz)')
    title([ttlNames{sp} ' (n = ' num2str(size(useData,1)) ')'])
end %subplot

if saveOrNot == 1
    saveas(gcf, [saveDir '\' figtitle], 'epsc')
    saveas(gcf, [saveDir '\' figtitle], 'fig')
    saveas(gcf, [saveDir '\' figtitle], 'png')
end % save

figtitle = 'BIRDTriggeredAverages';
    figure('Name', figtitle, 'Position', [202 344 1142 429])
for sp = 1:2
    subplot(1,2,sp)
    if sp == 1
        useData = btaPyr;
    else
        useData = btaInt;
    end %grab data

    plot(-timeAround:1/Fs:timeAround, mean(useData))
    xlabel('Time (s)')
    ylabel('Firing rate (Hz)')
    title([ttlNames{sp} ' (n = ' num2str(size(useData,1)) ')'])
end %subplot

if saveOrNot == 1
    saveas(gcf, [saveDir '\' figtitle], 'epsc')
    saveas(gcf, [saveDir '\' figtitle], 'fig')
    saveas(gcf, [saveDir '\' figtitle], 'png')
end % save

figtitle = 'SeizureTriggeredAverages';
    figure('Name', figtitle, 'Position', [202 344 1142 429])
for sp = 1:2
    subplot(1,2,sp)
    if sp == 1
        useData = sztaPyr;
    else
        useData = sztaInt;
    end %grab data

    plot(-timeAround:1/Fs:timeAround, mean(useData))
    xlabel('Time (s)')
    ylabel('Firing rate (Hz)')
    title([ttlNames{sp} ' (n = ' num2str(size(useData,1)) ')'])
end %subplot

if saveOrNot == 1
    saveas(gcf, [saveDir '\' figtitle], 'epsc')
    saveas(gcf, [saveDir '\' figtitle], 'fig')
    saveas(gcf, [saveDir '\' figtitle], 'png')
end % save

end %fxn