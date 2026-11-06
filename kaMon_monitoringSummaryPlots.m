function kaMon_monitoringSummaryPlots
% function kaMon_monitoringSummaryPlots
% 
% PURPOSE:
%   Creates summary plots for each mouse showing the interictal spike rate,
%   BIRD rate, and number of detected seizures. These figures are saved as
%   a .fig and .eps in the raw data folder for that mouse. Also creates a
%   seizure table for each mouse, which is saved in the raw data folder as
%   a .mat file.
% 
% MMD
% 07/2026


%% MAKE PLOTS BY MOUSE

clear
cln

dataDir = 'D:\AmygdalaMice';
mouseLog = readtable([dataDir '\mouse_masterlog.xlsx']);
kaInds = strmatch('KA', mouseLog.Group);

cd(neoDir) %need to be for python reasons

for m = kaInds'
    % fprintf('%s\n', mouseLog.Number{m})

    figtitle =  [mouseLog.Number{m} '_monitoringSummary'];
    figure('Name', figtitle, 'Position', [1 1 1142 429])
    mouseArr = []; %initialize
    mouseSeizArr = [];

    delDirs = dir([dataDir '\' mouseLog.Number{m} '\Delay*']); %some mice did not complete all 5 days (died /:)

    for d = 1:length(delDirs)
        % fprintf('\t%s\n', delDirs(d).name)
        monDirs = dir([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name '\*onitor*']);
        monDirs(strmatch('Post', vertcat({monDirs(:).name}))) = []; %get rid of any post-behavior monitoring to keep things consistent

        allSpks = 0;
        allBirds = 0;
        allSeiz = 0;
        allDur = 0;
        for mon = 1:length(monDirs)

            load([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                '\' monDirs(mon).name '\interictalSpikes.mat'], 'spkTms', 'monDur', 'birds', 'seiz'); %load

            allSpks = allSpks + length(spkTms);
            allBirds = allBirds + length(birds);
            allSeiz = allSeiz + length(seiz);
            allDur = allDur + monDur;

            for s = 1:length(seiz)
                tmpArr = [d mon spkTms(seiz{s}(1)) spkTms(seiz{s}(end))-spkTms(seiz{s}(1))];
                mouseSeizArr = [mouseSeizArr; tmpArr];
            end %seiz
        end %mon dir

        subplot(1,3,1)
        hold on
        plot(d, allSpks/allDur, 'ko')

        subplot(1,3,2)
        hold on
        plot(d, allBirds/allDur, 'ko')

        subplot(1,3,3)
        hold on
        plot(d, allSeiz, 'ko')

        tmpArr = [d monDur allSpks allSpks/allDur allBirds allBirds/allDur allSeiz];
        mouseArr = [mouseArr; tmpArr];

    end %day

    subplot(1,3,1)
    xlim([0.5 length(delDirs)+0.5])
    xlabel('Day')
    ax = gca;
    ylim([0 ax.YLim(2)])
    ylabel('Interictal spike rate (Hz)')
    axis square

    subplot(1,3,2)
    xlim([0.5 length(delDirs)+0.5])
    xlabel('Day')
    ax = gca;
    ylim([0 ax.YLim(2)])
    ylabel('BIRD rate (Hz)')
    axis square

    subplot(1,3,3)
    xlim([0.5 length(delDirs)+0.5])
    xlabel('Day')
    ax = gca;
    ylim([-0.5 ax.YLim(2)])
    ylabel('Number of detected seizures')
    axis square
    
    mouseTbl = array2table(mouseArr, 'VariableNames', {'day', 'monDur', 'numSpks', 'spkRate', 'numBirds', 'birdRate', 'numSeiz'});
    if isempty(mouseSeizArr)
        mouseSeizTbl = table();
    else
        mouseSeizTbl = array2table(mouseSeizArr, 'VariableNames', {'day', 'monFolder', 'onsetTm', 'seizDur'});
    end %this mouse had seizures or not
    save([dataDir '\' mouseLog.Number{m} '\' mouseLog.Number{m} '_monitorSummary.mat'], 'mouseTbl', 'mouseSeizTbl')
    saveas(gcf, [dataDir '\' mouseLog.Number{m} '\' figtitle], 'epsc')
    saveas(gcf, [dataDir '\' mouseLog.Number{m} '\' figtitle], 'fig')

end %mouse

%% OVERALL SUMMARY PLOTS

clear
cln

saveDir = 'D:\RESULTS\summary_examples';

dataDir = 'D:\AmygdalaMice';
mouseLog = readtable([dataDir '\mouse_masterlog.xlsx']);
kaInds = strmatch('KA', mouseLog.Group);

figtitle = 'amygMice_summary';
figure('Name', figtitle, 'Position', [228 413 1099 435])

for m = kaInds'
    % fprintf('%s\n', mouseLog.Number{m})

    delDirs = dir([dataDir '\' mouseLog.Number{m} '\Delay*']); %some mice did not complete all 5 days (died /:)

    for d = 1:length(delDirs)
        % fprintf('\t%s\n', delDirs(d).name)
        monDirs = dir([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name '\*onitor*']);
        monDirs(strmatch('Post', vertcat({monDirs(:).name}))) = []; %get rid of any post-behavior monitoring to keep things consistent

        allSpks = 0;
        allBirds = 0;
        allSeiz = 0;
        allDur = 0;
        for mon = 1:length(monDirs)

            load([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                '\' monDirs(mon).name '\interictalSpikes.mat'], 'spkTms', 'monDur', 'birds', 'seiz'); %load

            allSpks = allSpks + length(spkTms);
            allBirds = allBirds + length(birds);
            allSeiz = allSeiz + length(seiz);
            allDur = allDur + monDur;
        end %mon dir

        subplot(1,3,1)
        hold on
        scatter(match(m,kaInds), allSpks/allDur, 'MarkerEdgeColor', rgb('White'), 'MarkerFaceColor', rgb('Red'))

        subplot(1,3,2)
        hold on
        scatter(match(m,kaInds), allBirds/allDur, 'MarkerEdgeColor', rgb('White'), 'MarkerFaceColor', rgb('DarkGreen'))

        subplot(1,3,3)
        hold on
        scatter(match(m,kaInds), allSeiz, 'MarkerEdgeColor', rgb('White'), 'MarkerFaceColor', rgb('DeepPink'))

    end %day

end %mouse

subplot(1,3,1)
xlim([0.5 length(kaInds)+0.5])
% xlabel('Day')
xticks(1:length(kaInds))
xticklabels(mouseLog.Number(kaInds))
ax = gca;
ylim([0 ax.YLim(2)]) %make sure we start at 0
ylabel('Interictal spike rate (Hz)')
axis square

subplot(1,3,2)
xlim([0.5 length(kaInds)+0.5])
xticks(1:length(kaInds))
xticklabels(mouseLog.Number(kaInds))
ax = gca;
ylim([0 ax.YLim(2)])
ylabel('BIRD rate (Hz)')
axis square

subplot(1,3,3)
xlim([0.5 length(kaInds)+0.5])
xticks(1:length(kaInds))
xticklabels(mouseLog.Number(kaInds))
ax = gca;
ylim([-0.5 ax.YLim(2)])
ylabel('Number of detected seizures')
axis square

saveas(gcf, [saveDir '\' figtitle], 'epsc')
saveas(gcf, [saveDir '\' figtitle], 'fig')

%% OVERALL SUMMARY PLOTS - EDITED FOR LAURA'S GRANT

clear
cln

saveDir = 'D:\RESULTS\summary_examples';

dataDir = 'D:\AmygdalaMice';
neoDir = 'C:\Users\megmd\MATLAB_NEO';
mouseLog = readtable([dataDir '\mouse_masterlog.xlsx']);
kaInds = strmatch('KA', mouseLog.Group);

mice2remove = {'M454', 'M471', 'M493'}; %as requested by Laura
keepMice = setxor(mouseLog.Number(kaInds), mice2remove);
[plotMice, plotMiceInd] = intersect(mouseLog.Number, keepMice);

figtitle = 'amygMice_summary_LauraGrantVersion';
figure('Name', figtitle, 'Position', [228 413 1099 435])

nboot = 1000;

pullRates = [];

for m = plotMiceInd'
    % fprintf('%s\n', mouseLog.Number{m})

    delDirs = dir([dataDir '\' mouseLog.Number{m} '\Delay*']); %some mice did not complete all 5 days (died /:)

    tmpSData = [];
    tmpBData = [];
    tmpSzData = [];

    for d = 1:length(delDirs)
        % fprintf('\t%s\n', delDirs(d).name)
        monDirs = dir([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name '\*onitor*']);
        monDirs(strmatch('Post', vertcat({monDirs(:).name}))) = []; %get rid of any post-behavior monitoring to keep things consistent

        allSpks = 0;
        allBirds = 0;
        allSeiz = 0;
        allDur = 0;
        for mon = 1:length(monDirs)

            load([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name...
                '\' monDirs(mon).name '\interictalSpikes.mat'], 'spkTms', 'monDur', 'birds', 'seiz'); %load

            allSpks = allSpks + length(spkTms);
            allBirds = allBirds + length(birds);
            allSeiz = allSeiz + length(seiz);
            allDur = allDur + monDur;
        end %mon dir

        subplot(1,3,1)
        hold on
        scatter(match(m,plotMiceInd), allSpks/allDur, 'MarkerEdgeColor', rgb('White'), 'MarkerFaceColor', rgb('Red'), 'MarkerFaceAlpha', 0.5)
        tmpSData = [tmpSData; allSpks/allDur];
        pullRates = [pullRates; allSpks/allDur];

        subplot(1,3,2)
        hold on
        scatter(match(m,plotMiceInd), allBirds/allDur, 'MarkerEdgeColor', rgb('White'), 'MarkerFaceColor', rgb('DarkGreen'), 'MarkerFaceAlpha', 0.5)
        tmpBData = [tmpBData; allBirds/allDur];

        subplot(1,3,3)
        hold on
        scatter(match(m,plotMiceInd), allSeiz, 'MarkerEdgeColor', rgb('White'), 'MarkerFaceColor', rgb('DeepPink'), 'MarkerFaceAlpha', 0.5)
        tmpSzData = [tmpSzData; allSeiz];

    end %day

    subplot(1,3,1)
    tmpMean = mean(tmpSData);
    tmpCI = bootci(nboot, {@mean, tmpSData}, 'type', 'per');

    hold on
    line([match(m, plotMiceInd) match(m, plotMiceInd)], [tmpCI], 'Color', rgb('Red'), 'LineWidth', 0.5)
    line([match(m, plotMiceInd)-0.2 match(m, plotMiceInd)+0.2], [tmpMean tmpMean], 'Color', rgb('Red'), 'LineWidth', 3)

    subplot(1,3,2)
    tmpMean = mean(tmpBData);
    tmpCI = bootci(nboot, {@mean, tmpBData}, 'type', 'per');

    hold on
    line([match(m, plotMiceInd) match(m, plotMiceInd)], [tmpCI], 'Color', rgb('DarkGreen'), 'LineWidth', 0.5)
    line([match(m, plotMiceInd)-0.2 match(m, plotMiceInd)+0.2], [tmpMean tmpMean], 'Color', rgb('DarkGreen'), 'LineWidth', 3)

    subplot(1,3,3)
    tmpMean = mean(tmpSzData);
    tmpCI = bootci(nboot, {@mean, tmpSzData}, 'type', 'per');

    hold on
    line([match(m, plotMiceInd) match(m, plotMiceInd)], [tmpCI], 'Color', rgb('DeepPink'), 'LineWidth', 0.5)
    line([match(m, plotMiceInd)-0.2 match(m, plotMiceInd)+0.2], [tmpMean tmpMean], 'Color', rgb('DeepPink'), 'LineWidth', 3)
end %mouse

subplot(1,3,1)
xlim([0.5 length(plotMiceInd)+0.5])
% xlabel('Day')
xticks(1:length(plotMiceInd))
xticklabels(plotMice)
ax = gca;
ylim([0 ax.YLim(2)]) %make sure we start at 0
ylabel('Interictal spike rate (Hz)')
axis square

subplot(1,3,2)
xlim([0.5 length(plotMiceInd)+0.5])
% xlabel('Day')
xticks(1:length(plotMiceInd))
xticklabels(plotMice)
ax = gca;
ylim([0 ax.YLim(2)]) %make sure we start at 0
ylabel('BIRD rate (Hz)')
axis square

subplot(1,3,3)
xlim([0.5 length(plotMiceInd)+0.5])
% xlabel('Day')
xticks(1:length(plotMiceInd))
xticklabels(plotMice)
ax = gca;
ylim([0 ax.YLim(2)]) %make sure we start at 0
ylabel('Number of detected seizures')
axis square

% set(gcf, 'Renderer', 'Painters')
% saveas(gcf, [saveDir '\' figtitle], 'epsc')
% saveas(gcf, [saveDir '\' figtitle], 'fig')


end %function