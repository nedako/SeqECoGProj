for i = 3:44

    
    saveName = ['Raw_PSD_B',num2str(i) ,'.mat'];
    load(saveName);
    A = whos('-file', saveName);
    for tr = 1:length(A)
        eval(['Pall.PSD',num2str(tr),' = 10*log10(PSD',num2str(tr) , ');']);
    end
    disp(['PSD calculation for block ' , num2str(i) ,' completed'])
    save(saveName,'-struct','Pall', '-v7.3');
    clear Pall

end
%% Calculate the PSDs
load('AllData_Behav.mat')
chans = [3:15 , 99:109 , 121:129 ,36];

Dout  = secog_parseEEG_PSD('ParseEEG-freqDecomp' , Dall, 1 , 'Channels' , chans);


Dout  = secog_parseEEG_PSD('ParseEEG-calc_norm_PSD' , Dall, 1 , 'Channels' , chans);
% time warp
Dout  = secog_parseEEG_PSD('TimeWarpPSD_Raw_Binned' , Events, 1);
Dout  = secog_parseEEG_PSD('TimeWarpPSD_Raw_Binned_seqType' , Events, 1);

cd('/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/P2')
load('AllData_PSD_Warped.mat')
secog_BlockGroupAverage(Pall,1,'raw_BlockGroup', 'Channels' , chans)

secog_BlockGroupAverage(Pall,1,'binned_BlockGroup', 'Channels' , chans)
secog_BlockGroupAverage([],1,'raw_AvgPower_BlockGroup', 'Channels' , chans)
secog_BlockGroupAverage([],1,'raw_SeqType', 'Channels' , chans)


%% plotting case summary
%      'binned_SingleTrial'
%      'binned_SingleTrial_AvgChann'
%      'binned_BlockGroup_AvgChann'
%      'binned_BlockGroup'
%      'raw_SingleTrial'
%      'raw_BlockGroup'
%      'raw_BlockGroup_AvgChann'
blockGroupNames = {'SingleFingNat' , 'SingleFingSlow' , 'SingleFingFast' , 'Intermixed1' , 'Intermixed2' , 'ChunkDay1' , 'Intermixed3' , 'Intermixed4' , 'Intermixed5',...
    'ChunkDay2' , 'Intermixed6' , 'Intermixed7' , 'Intermixed8', 'ChunkDay3', 'Intermixed9'}';


secog_visualizePSD(Pall,1,'raw_BlockGroup' , 'BlockGroup' , 'SingleFingSlow','Chan2Plot' , [20 :22, 24, 12, 14],'Channels' , chans) % PMd
secog_visualizePSD(Pall,1,'raw_BlockGroup' , 'BlockGroup' , 'SingleFingSlow','Chan2Plot' , [27:30],'Channels' , chans) % SMA
secog_visualizePSD(Pall,1,'binned_BlockGroup' , 'BlockGroup' , 'SingleFingFast','Chan2Plot' , [20 :22, 24, 12, 14],'Channels' , chans); % PMd
secog_visualizePSD(Pall,1,'binned_BlockGroup' , 'BlockGroup' , 'SingleFingFast','Chan2Plot' , [27:30],'Channels' , chans); % SMA



secog_visualizePSD(Pall,1,'raw_SingleTrial' , 'TBNum', [10 1], 'Chan2Plot' , [23 8],'Channels' , chans)

load('AllData_PSD_StimNorm.mat')
secog_visualizePSD(Pall,1,'binned_SingleTrial' , 'TBNum', [10 4], 'Chan2Plot' , [104:107 , 109])
secog_visualizePSD(Pall,1,'binned_SingleTrial_AvgChann' , 'TBNum', [10 40], 'Chan2Plot' , [104:107 , 109])
secog_visualizePSD(Pall,1,'raw_SingleTrial' , 'TBNum', [1 1], 'Chan2Plot' , [104:107 , 109])


%%
load('AllData_PSD_Warped.mat')
% PMD layer 2
secog_visualizePSD(Pall,1,'binned_BlockGroup_AvgChann' , 'BlockGroup' , 'SingleFingSlow','Chan2Plot' , [104:107 , 109])
secog_visualizePSD(Pall,1,'binned_BlockGroup' , 'BlockGroup' , 'SingleFingSlow','Chan2Plot' , [104:107 , 109])
% PMD layer 1
secog_visualizePSD(Pall,1,'binned_BlockGroup' , 'BlockGroup' , 'SingleFingSlow','Chan2Plot' , [12 14])
% PMv
secog_visualizePSD(Pall,1,'binned_BlockGroup' , 'BlockGroup' , 'SingleFingSlow','Chan2Plot' , [6])
% SMA
secog_visualizePSD(Pall,1,'binned_BlockGroup' , 'BlockGroup' , 'SingleFingSlow','Chan2Plot' , [123:126])
secog_visualizePSD(Pall,1,'binned_BlockGroup' , 'BlockGroup' , blockGroupNames{3},'Chan2Plot' , [123:126])


