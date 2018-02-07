function secog_packEEG(subjNum)
%% Uses the PathInfo excel file. Loads the .mat data files extracted using Spike2.8, that's saved as individual channels 
% packs all the labeled channles inot one data structure and save into packed folder
%% load the first recording
subjname = {'P2'};
mainDir = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/' subjname{subjNum} , '/'] ;

%% load up the path info file from the patients file 
[~, ~, PathInfo] = xlsread([mainDir , 'PathInfo.xlsx'],'Sheet1');
PathInfo(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),PathInfo)) = {''};

idx = cellfun(@ischar, PathInfo);
PathInfo(idx) = cellfun(@(x) string(x), PathInfo(idx), 'UniformOutput', false);

clearvars idx;
%% load the .mat up, pack individual channle to a single data structure and save 
% only consider channles that have a channel lable, or else is the marker. for every patient this can be differnet. so far marker (TTL) always comes out as CH141
% load labeld channels
for i = 5:size(PathInfo , 1)
    Data = [];
    cd(char(strcat(mainDir , PathInfo{i,1})));
    filename = PathInfo{i,2};
    % get all the filds (channles) stored in the filename
    vars = whos('-file',char(filename));
    for chans = 1:length(vars)
        D = load(filename , vars(chans).name);
        eval(['D = D.' , vars(chans).name , ';']);
        id = 3;
        while isempty(str2num(vars(chans).name(end-id:end)))
            id = id - 1;
        end
        D.ChannelNumber = str2num(vars(chans).name(end-id:end));
        if isfield(D , 'values')
            D.values = D.values';
            D.label{1} =D.title;
            fields = {'comment' , 'title'}; % fields title and comment are the same
            D = rmfield(D , fields);
            if D.ChannelNumber==141
                D.label{1} = 'TTL';
            end
            % make sure that the channle has a separate label and is not null
            if ~strcmp(vars(chans).name([53,55:end]) , D.label{1}) || D.ChannelNumber==141
                % channle numbers matter! get them right
                Data = addstruct(Data , D);
            end
        end
        disp(['Channle ' , num2str(D.ChannelNumber) , ' Read'])
    end
    filename = char(filename);
    saveName = [mainDir , 'Packed/' , strcat(filename(1:end-4) ,  '_packed.mat')];
    disp(['Saving ' , saveName , ' ...'])
    save(saveName , 'Data' , '-v7.3')
    
end





