function units = read_in_abr_phy(path, varargin)
% function read_in_phy_data(path)
%
% PURPOSE:
%   Standardize a function to read in and format the data from phy after
%   manual curation. Only reads in the bare bones for the good clusters
%   (unit id, spike times, spike inds).
%
% INPUT:
%   path =  folder where data is
%   Fs = sampling frequency (optional, default = 30000)
%
% OUTPUT:
%   units = struct with info about each unit
%
% MMD
% Ewell Lab
% 08/2026

if nargin == 1
    Fs = 30000;
end %check for fs

spike_times = double(readNPY([path '\spike_times.npy']));
spike_clusters = readNPY([path '\spike_clusters.npy']);
cluster_info = readtable([path '\cluster_info.tsv'], 'FileType', 'text', 'Delimiter', '\t');

goodUnits = find(strcmp('good', cluster_info.group));

if isempty(goodUnits)
    units = [];
else
    for u = 1:length(goodUnits)
        cl = cluster_info.cluster_id(goodUnits(u)); %cluster id
        units(u).id = cl;
        findSpks = spike_clusters==cl;
        units(u).spkInds = spike_times(findSpks);
        units(u).spkTms = spike_times(findSpks) ./ Fs; %convert to s

    end %unit
end %if good units

end %funciton