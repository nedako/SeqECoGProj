secog_makeEEG(subjNum , fileInfo)

%% load the first recording
subjname = {'P2'};
mainDir = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/' subjname{subjNum} , '/'] ;

cd([mainDir , 'iEEG_P2_Jan16/rec1'])
filename = 'P2_Jan16-rec1.mat';
prefix = 'Thissen__Peter_01d34644_b566_4679_84d6_d3f415b203f6_Ch';
n = 1;
fName = [prefix num2str(n)];
Data = load(filename , fName);
cmdstr = ['Data = Data.' , fName];
eval(cmdstr);
Data.ChannelNumber = n;
Data.values = Data.values';
fields = {'comment'}; % fields title and comment are the same
Data = rmfield(Data , fields);


for n = 2:200
    fName = [prefix num2str(n)];
    temp = load(filename , fName);
    temp.ChannelNumber = n;
    cmdstr = ['temp = temp.' , fName];
    eval(cmdstr);
    temp.ChannelNumber = n;
    temp.values = temp.values';
    temp = rmfield(temp , fields);
    
    Data = addstruct(Data , temp);
    clear temp
end
%% load the first channel of the second recording to figure out the shared amount of data
filename = 'P1_sep14-rec2.mat';
prefix = 'Gibson__Darah_0db15236_15fa_4bf5_bb87_704564c8b16c_Ch';
n = 1;
fName = [prefix num2str(n)];
temp = load(filename , fName);
cmdstr = ['temp = temp.' , fName];
eval(cmdstr);
A = temp.values;
% figure out how many samples at the end of recording 1, comes in the beggining of recording 1
i = 1:length(A);
cnt = 1;
while ~ismember((A(1:i(cnt)))' , Data.values(1,end-(i(cnt)-1):end) , 'rows')
    cnt = cnt + 1
end

disp(['The shared amount of data is ' , num2str(cnt) , ' samples or ' , num2str(cnt/(1024*60)) , ' minutes.'])

%%

filename = 'P1_sep14-rec2.mat';
prefix = 'Gibson__Darah_0db15236_15fa_4bf5_bb87_704564c8b16c_Ch';
n = 1;
fName = [prefix num2str(n)];
Data2 = load(filename , fName);
cmdstr = ['Data2 = Data2.' , fName];
eval(cmdstr);
Data2.ChannelNumber = n;
id  = ones(length(Data2.values) , 1);
id(1:cnt) = 0;
Data2.values = Data2.values(logical(id))';
fields = {'title' , 'comment'};
Data2 = rmfield(Data2 , fields);



for n = 2:200
    fName = [prefix num2str(n)];
    temp = load(filename , fName);
    temp.ChannelNumber = n;
    cmdstr = ['temp = temp.' , fName];
    eval(cmdstr);
    temp.ChannelNumber = n;
    temp.values = temp.values(logical(id))';
    temp = rmfield(temp , fields);
    
    Data2 = addstruct(Data2 , temp);
    clear temp
end

Data.values = [Data.values Data2.values];
