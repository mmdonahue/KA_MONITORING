function kaMon_checkCells
% function kaMon_checkCells
% 
% PURPOSE:
%   This function goes into every folder and checkes that I have cut the
%   cells for a given day. Also displays the number of clusters and the
%   clusters identified as "good" cells.
% 
% MMD
% Ewell Lab
% 10/2026

%% CHECK THAT WE HAVE THE CELLS DONE

clear
cln

dataDir = 'D:\AmygdalaMice';
mouseLog = readtable([dataDir '\mouse_masterlog.xlsx']);
kaInds = strmatch('KA', mouseLog.Group);

for m = kaInds'
    fprintf('%s\n', mouseLog.Number{m})

    delDirs = dir([dataDir '\' mouseLog.Number{m} '\Delay*']); %some mice did not complete all 5 days (died /:)

    for d = 1:length(delDirs)
        fprintf('\t%s\n', delDirs(d).name)
        monDirs = dir([dataDir '\' mouseLog.Number{m} '\' delDirs(d).name '\*onitor*']);
        monDirs(strmatch('Post', vertcat({monDirs(:).name}))) = []; %get rid of any post-behavior monitoring to keep things consistent
        if length(monDirs) == 1
            mon = 1;
        else
            mon = strmatch('Combined', vertcat({monDirs(:).name}));
        end
        try
            cluster_info = readtable([monDirs(mon).folder '\' monDirs(mon).name...
                '\SpikeSort\Pipeline\Phy\cluster_info.tsv'], 'FileType', 'text', 'Delimiter', '\t');
            goodUnits = find(strcmp('good', cluster_info.group));

            fprintf('\t\t%d clusters - %d good units\n', size(cluster_info,1), size(goodUnits,1))
        catch
            fprintf('\t\tCheck cluster detection\n')
        end %try to see if there is something there

    end %day

end %mouse

end %function