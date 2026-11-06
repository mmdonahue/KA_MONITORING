function kaMon_visualize_ISunits

%% EXAMPLE DAY
% use initial day where I have detected units to visualize simulatenous
% spiking and LFP

saveDir = 'D:\RESULTS\summary_examples';
dayFolder = 'D:\AmygdalaMice\M454\Delay3\Pre_monitor';

Fs = 3000; %for lfp
ds = 30000 ./ Fs; % fs = 30 kHz - ds factor of 10 -> 3000 Hz

load([dayFolder '\ISdetectWorkFile.mat']); %load
load([dayFolder '\interictalSpikes.mat'], 'spkTms', 'spkInds', 'monDur', 'birds', 'seiz'); %load

strInd1 = strfind(dayFolder, 'M4');
strInd2 = strfind(dayFolder, 'Delay');
ttl = dayFolder(strInd1:strInd2+5);
ttl(5) = '_';

for b = 1:4
    figtitle = [ttl '_bin' num2str(b) '_spikeRstr'];
    figure('Name', figtitle, 'Position', [1 1 1287 857])
    subplot(6,1,1)
    startTm = (useBins(b)-1)*3*60; %in second
    endTm = startTm + 3*60;
    plot(binLfp(:,b), 'k')
    xlim([0 6e5])

    hold on
    tmpSpks = spkInds(spkTms >= startTm & spkTms <= endTm);
    tmpSpkTms =  spkTms(spkTms >= startTm & spkTms <= endTm);

    if ~isempty(tmpSpks)
        xline(tmpSpks - startTm*Fs + 1,  'r', 'Layer', 'bottom')
    end

    xlim([0 6e5])

    %scale bars
    line([5.55e5 5.55e6+(Fs*5)], [0 0], 'Color', 'k')
    text(5.55e5 , -1000, '5 s')
    line([5.55e5 5.55e5], [0 5000], 'Color', 'k')
    text(5.575e5, 1000, '5 mV')
    axis off

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

    subplot(6,1,2:6)

    unit = read_in_abr_phy([dayFolder '\SpikeSort\Pipeline\Phy']);

    colorMap = jet(length(unit));
    for u = 1:length(unit)
        hold on
        tmpSpkTms =  unit(u).spkTms(unit(u).spkTms >= startTm & unit(u).spkTms <= endTm);

        plot(tmpSpkTms - startTm, repmat(u, length(tmpSpkTms), 1), '|', 'Color', colorMap(u,:))
    end %unit

    xlim([0 6e5./Fs])
    ylim([0 length(unit)+1])
    axis off

    % saveas(gcf, [saveDir '\' figtitle], 'epsc')
    % saveas(gcf, [saveDir '\' figtitle], 'fig')
end %begin



end %funciton