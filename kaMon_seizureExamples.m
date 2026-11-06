function kaMon_seizureExamples

%% BEHAVIORAL SEIZURE EXAMPLES

clear
cln

saveDir = 'D:\RESULTS\summary_examples';

dataDir = 'D:\AmygdalaMice';
neoDir = 'C:\Users\megmd\MATLAB_NEO';
mouseLog = readtable([dataDir '\mouse_masterlog.xlsx']);
kaInds = strmatch('KA', mouseLog.Group);

Fs = 3000;
ds = 30000 ./ Fs; % fs = 30 kHz - ds factor of 10 -> 3000 Hz
cd(neoDir) %need to be for python reasons

m = 2;
d = 1;
mon = 1;
mouseInfo = [mouseLog.Number{m} ' - Delay ' num2str(d)]; %annotated on the bottom of the figure - to keep track of where the data is coming from
monDirs = dir([dataDir '\' mouseLog.Number{m} '\Delay' num2str(d) '\*onitor*']);

load([dataDir '\' mouseLog.Number{m} '\Delay' num2str(d) '\' monDirs(mon).name '\ISdetectWorkFile.mat']);
load([dataDir '\' mouseLog.Number{m} '\Delay' num2str(d)...
    '\' monDirs(mon).name '\interictalSpikes.mat'], 'spkTms', 'spkInds', 'monDur', 'birds', 'seiz'); %load

binFile = dir([dataDir '\' mouseLog.Number{m} '\Delay' num2str(d)...
    '\' monDirs(mon).name '\ampData*.bin']);
binFile = [dataDir '\' mouseLog.Number{m} '\Delay' num2str(d)...
    '\' monDirs(mon).name '\' binFile.name]; %convert to full

startTm = 0;
startInd = 0; %start frame, in terms of fs/ds sampling rate; start at time 0
endTm = monDur;
endInd = monDur*Fs;

r = pyrunfile("python_code\load_bonsai_binary.py", "traces", ...
    binary_file = binFile, ... %'path\to\file.bin', ...
    ds = ds, ... % downsampling factor
    channels = py.list([0 maxCh-1]), ... % which channels to get, if you want all set channels = -1. NOTE, starts AT ZERO, ends at 64
    start_frame = startInd, ... % start frame, in terms of fs/ds sampling rate
    stop_frame = endInd); % end frame, in terms of fs/ds sampling rate

r_use = double(r);
lfpData = r_use(:,2);

%two behavioral seizures in this monitoring window - sz = 7 and sz = 9, both ~1 min long
useSz = [7 9]; %so we can loop!

for sz = useSz
    figtitle = [mouseLog.Number{m} '_Delay' num2str(d) '_behavSeiz' num2str(match(sz,useSz))];
    figure('Name', figtitle)

    startTm = spkTms(seiz{sz}(1));
    endTm = spkTms(seiz{sz}(end));

    subplot(3,1,1)
    plot(-60:1/Fs:5*60, lfpData((startTm-60)*Fs:(startTm+5*60)*Fs), 'k')
    xlim([-60 5.5*60])

    %identify the seizure
    ax = gca;
    yMax = ax.YLim(end);
    xLims = [0 spkTms(seiz{sz}(end)) - spkTms(seiz{sz}(1))];
    yLims = [yMax yMax];
    line([xLims], [yLims], 'Color', rgb('DeepPink'), 'LineWidth', 5)

    %indicate that the next subplot will be
    yMin = ax.YLim(1);
    xLims = [-5 90];
    yLims = [yMin yMin];
    line([xLims], [yLims], 'Color', rgb('Gray'), 'LineWidth', 5)

    %scale bar
    line([310 320], [0 0], 'Color', 'k')
    text(310, -700, '5 s')
    line([310 310], [0 1000], 'Color', 'k')
    text(312, 605, '1 mV')
    axis off

    subplot(3,1,2)
    plot(-5:1/Fs:90, lfpData((startTm-5)*Fs:(startTm+90)*Fs), 'k')
    xlim([-5 100])

    %identify the seizure
    ax = gca;
    yMax = ax.YLim(end);
    xLims = [0 spkTms(seiz{sz}(end)) - spkTms(seiz{sz}(1))];
    yLims = [yMax yMax];
    line([xLims], [yLims], 'Color', rgb('DeepPink'), 'LineWidth', 5)

    %indicate that the next subplot will be
    yMin = ax.YLim(1);
    xLims = [-2 60];
    yLims = [yMin yMin];
    line([xLims], [yLims], 'Color', rgb('Gray'), 'LineWidth', 5)

    %scale bar
    line([92 97], [0 0], 'Color', 'k')
    text(92, -700, '5 s')
    line([92 92], [0 1000], 'Color', 'k')
    text(92, 605, '1 mV')
    axis off

    subplot(3,1,3)
    plot(-2:1/Fs:60, lfpData((startTm-2)*Fs:(startTm+60)*Fs), 'k')
    xlim([-2 65])

    %identify the seizure
    ax = gca;
    yMax = ax.YLim(end);
    xLims = [0 spkTms(seiz{sz}(end)) - spkTms(seiz{sz}(1))];
    yLims = [yMax yMax];
    line([xLims], [yLims], 'Color', rgb('DeepPink'), 'LineWidth', 5)

    %scale bar
    line([61 62], [0 0], 'Color', 'k')
    text(61, -700, '1 s')
    line([61 61], [0 1000], 'Color', 'k')
    text(61, 605, '1 mV')
    axis off

    annotation('textbox', [0 0 0.1 0.1], 'string', mouseInfo, 'EdgeColor', 'None')

    saveas(gcf, [saveDir '\' figtitle], 'epsc')
    saveas(gcf, [saveDir '\' figtitle], 'fig')
end %sz

%% SUBCLINICAL SEIZURE EXAMPLES

clear
cln

saveDir = 'D:\RESULTS\summary_examples';

dataDir = 'D:\AmygdalaMice';
neoDir = 'C:\Users\megmd\MATLAB_NEO';
mouseLog = readtable([dataDir '\mouse_masterlog.xlsx']);
kaInds = strmatch('KA', mouseLog.Group);

Fs = 3000;
ds = 30000 ./ Fs; % fs = 30 kHz - ds factor of 10 -> 3000 Hz
cd(neoDir) %need to be for python reasons

m = 5;
d = 2;
mon = 1;
mouseInfo = [mouseLog.Number{m} ' - Delay ' num2str(d)]; %annotated on the bottom of the figure - to keep track of where the data is coming from
monDirs = dir([dataDir '\' mouseLog.Number{m} '\Delay' num2str(d) '\*onitor*']);

load([dataDir '\' mouseLog.Number{m} '\Delay' num2str(d) '\' monDirs(mon).name '\ISdetectWorkFile.mat']);
load([dataDir '\' mouseLog.Number{m} '\Delay' num2str(d)...
    '\' monDirs(mon).name '\interictalSpikes.mat'], 'spkTms', 'spkInds', 'monDur', 'birds', 'seiz'); %load

binFile = dir([dataDir '\' mouseLog.Number{m} '\Delay' num2str(d)...
    '\' monDirs(mon).name '\ampData*.bin']);
binFile = [dataDir '\' mouseLog.Number{m} '\Delay' num2str(d)...
    '\' monDirs(mon).name '\' binFile.name]; %convert to full

startTm = 0;
startInd = 0; %start frame, in terms of fs/ds sampling rate; start at time 0
endTm = monDur;
endInd = monDur*Fs;

r = pyrunfile("python_code\load_bonsai_binary.py", "traces", ...
    binary_file = binFile, ... %'path\to\file.bin', ...
    ds = ds, ... % downsampling factor
    channels = py.list([0 maxCh-1]), ... % which channels to get, if you want all set channels = -1. NOTE, starts AT ZERO, ends at 64
    start_frame = startInd, ... % start frame, in terms of fs/ds sampling rate
    stop_frame = endInd); % end frame, in terms of fs/ds sampling rate

r_use = double(r);
lfpData = r_use(:,2);

%two behavioral seizures in this monitoring window - sz = 7 and sz = 9, both ~1 min long
useSz = [1 2]; %so we can loop!

for sz = useSz
    figtitle = [mouseLog.Number{m} '_Delay' num2str(d) '_subClinSeiz' num2str(sz)];
    figure('Name', figtitle)

    startTm = spkTms(seiz{sz}(1));
    endTm = spkTms(seiz{sz}(end));

    subplot(3,1,1)
    plot(-60:1/Fs:5*60, lfpData((startTm-60)*Fs:(startTm+5*60)*Fs), 'k')
    xlim([-60 5.5*60])

    %identify the seizure
    ax = gca;
    yMax = ax.YLim(end);
    xLims = [0 spkTms(seiz{sz}(end)) - spkTms(seiz{sz}(1))];
    yLims = [yMax yMax];
    line([xLims], [yLims], 'Color', rgb('DeepPink'), 'LineWidth', 5)

    %indicate that the next subplot will be
    yMin = ax.YLim(1);
    xLims = [-10 85];
    yLims = [yMin yMin];
    line([xLims], [yLims], 'Color', rgb('Gray'), 'LineWidth', 5)

    %scale bar
    line([310 320], [0 0], 'Color', 'k')
    text(310, -700, '5 s')
    line([310 310], [0 1000], 'Color', 'k')
    text(312, 605, '1 mV')
    axis off

    subplot(3,1,2)
    plot(-10:1/Fs:85, lfpData((startTm-10)*Fs:(startTm+85)*Fs), 'k')
    xlim([-10 95])

    %identify the seizure
    ax = gca;
    yMax = ax.YLim(end);
    xLims = [0 spkTms(seiz{sz}(end)) - spkTms(seiz{sz}(1))];
    yLims = [yMax yMax];
    line([xLims], [yLims], 'Color', rgb('DeepPink'), 'LineWidth', 5)

    %indicate that the next subplot will be
    yMin = ax.YLim(1);
    xLims = [-2 60];
    yLims = [yMin yMin];
    line([xLims], [yLims], 'Color', rgb('Gray'), 'LineWidth', 5)

    %scale bar
    line([87 92], [0 0], 'Color', 'k')
    text(87, -700, '5 s')
    line([87 87], [0 1000], 'Color', 'k')
    text(87, 605, '1 mV')
    axis off

    subplot(3,1,3)
    plot(-2:1/Fs:60, lfpData((startTm-2)*Fs:(startTm+60)*Fs), 'k')
    xlim([-2 65])

    %identify the seizure
    ax = gca;
    yMax = ax.YLim(end);
    xLims = [0 spkTms(seiz{sz}(end)) - spkTms(seiz{sz}(1))];
    yLims = [yMax yMax];
    line([xLims], [yLims], 'Color', rgb('DeepPink'), 'LineWidth', 5)

    %scale bar
    line([61 62], [0 0], 'Color', 'k')
    text(61, -700, '1 s')
    line([61 61], [0 1000], 'Color', 'k')
    text(61, 605, '1 mV')
    axis off

    annotation('textbox', [0 0 0.1 0.1], 'string', mouseInfo, 'EdgeColor', 'None')

    saveas(gcf, [saveDir '\' figtitle], 'epsc')
    saveas(gcf, [saveDir '\' figtitle], 'fig')
end %sz









end %function