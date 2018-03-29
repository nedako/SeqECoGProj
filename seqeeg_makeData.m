function seqeeg_makeData(subjNum , chans)
%% This function generates and saves all the data structures needed for further analysis
  % Provide subject number
  % Provide channels of interest
%%
subjname = {'P2' , 'P4' , 'P5'};
mainDir = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/' subjname{subjNum} ,'/'];

Dall = seqeeg_analyze ('sing_subj' , 'subjCode' , subjname{subjNum});
save([mainDir , 'AllData_Behav_Force.mat'] , 'Dall')
Dall = rmfield(Dall , 'F');
save([mainDir , 'AllData_Behav.mat'] , 'Dall')

%% Pack EEG
% seqeeg_packEEG('PackEEG' , subjNum);
%% Calculate the PSDs
% cd(['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/' , subjname{subjNum}]);
load([mainDir,'AllData_Behav.mat'])

% Dout  = seqeeg_parseEEG_PSD('ParseEEG-freqDecomp' , Dall, subjNum , 'Channels' , chans);
% Dout  = seqeeg_parseEEG_PSD('ParseEEG-calc_norm_PSD' , Dall, subjNum , 'Channels' , chans);
load([mainDir,'AllData_Events.mat'])
% time warp
% Dout  = seqeeg_parseEEG_PSD('TimeWarpPSD_Raw_Binned' , Events, subjNum);
% Dout  = seqeeg_parseEEG_PSD('TimeWarpPSD_Raw_Binned_seqType' , Events, subjNum);
% Dout  = seqeeg_parseEEG_PSD('AlignEvents_SeqType' , Events, subjNum);
% Dout  = seqeeg_parseEEG_PSD('AlignEvents_SeqType_Warped' , Events, subjNum);


%% average pattern calculation
cd(['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/' , subjname{subjNum}]);
load('AllData_PSD_Warped.mat')
seqeeg_BlockGroup(Pall,subjNum,'binned_BlockGroup', 'Channels' , chans)
seqeeg_BlockGroup(Pall,subjNum,'raw_BlockGroup', 'Channels' , chans)

load('AllData_PSD_Warped_SeqType.mat')
seqeeg_BlockGroup(Pall,subjNum,'raw_SeqType', 'Channels' , chans)  
seqeeg_BlockGroup(Pall,subjNum,'binned_SeqType', 'Channels' , chans)  
seqeeg_BlockGroup([],subjNum,'raw_AvgPower_BlockGroup', 'Channels' , chans)
seqeeg_BlockGroup([],subjNum,'raw_AvgPower_SeqType', 'Channels' , chans)
seqeeg_BlockGroup([],subjNum,'raw_Power_SeqType', 'Channels' , chans)
seqeeg_BlockGroup([],subjNum,'binned_Power_SeqType', 'Channels' , chans)
seqeeg_BlockGroup([],subjNum,'Aligned_SeqType', 'Channels' , chans)
seqeeg_BlockGroup([],subjNum,'AlignedWarped_SeqType', 'Channels' , chans)
