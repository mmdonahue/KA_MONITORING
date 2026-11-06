



monDirs = dir(['D:\AmygdalaMice\' mouseLog.Number{m} '\' delDirs(d).name '\*onitor*']); %finds all of the folders with "onitor" in them
%   this is so we don't have to worry about the capitalization of the M,
%   and will include all of the potential names
monDirs(strmatch('Post', vertcat({monDirs(:).name}))) = []; %get rid of any post-behavior monitoring folders

if length(monDirs) > 1 %if there are still more than one folder in our list, it is mostly likely due to monitoring being split in two
    mon = strmatch('Combined', vertcat({monDirs(:).name})); %so we check to find the "combined_monitoring" folder that I made
%       this is the folder that contains the merged .bin file with all of
%       the ephys data
    if isempty(mon) %if we now have NO folder left
        keyboard %this is a check that will pause the code if we have somehow eliminated very folder
    end %if none
else %if there is only one folder remaining for us to do into
    mon = 1;
end %set which monitoring folder to go into

mouseLog = readtable('D:\AmygdalaMice\mouse_masterlog.xlsx'); %change to the destination of the file on your computer
kaInds = strmatch('KA', mouseLog.Group); %this pulls out the indicies of the KA mice