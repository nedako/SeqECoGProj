function seqeeg_packEEG(what , subjNum, varargin)

c = 1;
while(c<=length(varargin))
    switch(varargin{c})
        case {'Channels'}
            % channels of interest Default : everythig
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        otherwise
            error(sprintf('Unknown option: %s',varargin{c}));
    end
end

subjname = {'P2' , 'P4' , 'P5'};
saveDir = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/' subjname{subjNum} , '/Packed/'] ;
mainDir = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/' subjname{subjNum} , '/'] ;
cd(mainDir)
%% Import PathInfo
[~, ~, PathInfo] = xlsread(['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/' ,subjname{subjNum} ,'/PathInfo.xlsx'],'Sheet1');
PathInfo(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),PathInfo)) = {''};

idx = cellfun(@ischar, PathInfo);
PathInfo(idx) = cellfun(@(x) string(x), PathInfo(idx), 'UniformOutput', false);

%% Clear temporary variables
clearvars idx;
%% chop up the EEG data into trials, and filter out the power line noise
switch what
    case 'PackEEG'
        for p = 1:size(PathInfo , 1)
            Data = [];
            cd([mainDir ,char(PathInfo{p , 1})])
            allInfo = whos('-file' ,  char(PathInfo{p,2}));
            for ch = 1:length(allInfo)
                name = allInfo(ch).name;
                eeg = load(PathInfo{p,2} , name);
                eval(['eeg = eeg.' , name , ';']);
                
                if (isfield(eeg , 'values') & isfield(eeg , 'interval')) & ~sum(strcmp(eeg.title(1) , {'C'})) & ~sum(strcmp(eeg.title(1:2) , {'Tr'}))
                    disp([char(PathInfo{p,2}) , ' Channel ' , eeg.title , ' read.'])
                    chnum = 0;
                    while ~isempty(str2num(name(end-chnum:end)))
                        chnum = chnum +1;
                    end
                    chnum = chnum -1;
                    D.ChannelNumber = str2num(name(end-chnum:end));
                    D.interval   = eeg.interval;
                    D.length     = eeg.length;
                    D.offset     = eeg.offset;
                    D.scale      = eeg.scale;
                    D.start      = eeg.start;
                    D.units      = eeg.units;
                    D.label      = {eeg.title};
                    D.values     = eeg.values';
                    Data         = addstruct(Data , D);
                    
                elseif str2num(name(end-2:end)) == 141
                    chnum = 0;
                    while ~isempty(str2num(name(end-chnum:end)))
                        chnum = chnum +1;
                    end
                    chnum = chnum -1;
                    D.ChannelNumber = str2num(name(end-chnum:end));
                    D.interval   = eeg.interval;
                    D.length     = eeg.length;
                    D.offset     = eeg.offset;
                    D.scale      = eeg.scale;
                    D.start      = eeg.start;
                    D.units      = eeg.units;
                    D.label      = {eeg.title};
                    D.values     = eeg.values';
                    % TTL on DC 13 always shows up on C144
                    D.label = {'TTL'};
                    disp([char(PathInfo{p,2}) , ' Channel ' , eeg.title , ' read.'])
                    Data         = addstruct(Data , D);
                end
                
            end
            savename = char(PathInfo{p,2});
            Marker = getrow(Data, ismember(Data.label , 'TTL'));
            figure;
            plot(Marker.values)
            title('Trigger Channel')
            disp('*************************************************************************')
            disp('Complete the BlockInfo.xlsx file by entering the start and end of blocks')
            disp('***If a block is not marked put -1 for beginnign and end of the block.***')
            disp('*************************************************************************')
            keyboard
            save(strcat([saveDir , savename(1:end-4)] , '_packed.mat') , 'Data' , '-v7.3')
        end
        ChanLabels = Data.label;
        save([mainDir , 'ChanLabels.mat'] , 'ChanLabels')
        
end


