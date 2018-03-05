#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 
Created on Wed Feb 21 09:56:22 2018

@author: nkordjazi
"""
import tensorflow  
import statsmodels
import os
import numpy as np 
import os
import scipy.io as sio
from numpy.linalg import inv  # Matrix invesion
from numpy import dot         # matrix multiplication
import itertools
import matplotlib.pylab as plt
from numpy import linalg as LA
import copy 
import pandas as pd
from xlrd import open_workbook


def secog_parseEEG_PSD(what , subjNum, DownsampleRate = 10 , NumWarpSampFast = 200 ,
                       NumWarpSampSlow = 500 , TimeDelay = 0.5, FreqRange = [2 ,180], 
                       numFreqBins = 90 , Channels)

## reads the all channels-packed EEG data, uses the BlockInfo file to parse the EEG into single trials and ammend the bihavioral data structure > 
## It also calculates the PSD on the whole block and then parses up the PSD inot trials. This is mainly to avoid any window effect on single trials
## setup the defaults and deal with the varargin
subjname = ['P2' , 'P4'];
mainDir = '/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/' + subjname[subjNum-1] + '/Packed/' ;
saveDir = '/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/' + subjname[subjNum-1] + '/' ;

temp = sio.loadmat(saveDir + 'ChanLabels.mat');
# get the field of the dictionary that contains the data
temp = temp['ChanLabels'];
ChanLabels=list();#np.zeros((1 ,len(temp)));
for i in range(0,len(temp)):
    ChanLabels.append(list(list(temp[i])[0]))
    ChanLabels[i]
# load behavioral data
temp = sio.loadmat(saveDir + 'AllData_Behav.mat')# , squeeze_me=True,struct_as_record=False);
Dall = temp['Dall'];
mdtype = Dall.dtype;
# turn Dall into a paython dictonary
Dall = {n: Dall[n][0, 0] for n in mdtype.names}
# turn Dall into a Pandas dataframe
Dall = pd.Series(Dall).to_frame()
# get the list of keys in the Dall dictionary
keyList = list()
for i in Dall.keys():
    keyList.append(i)



# if Channels has not been set, take everything
try:
    Channels
except NameError:
    Channels = range(0 , ChanLabels.size);
    
##  control for too short IPIs that the keys get accidentally pressed
if (subjNum==1):
    for i in range(len(Dall.TN)):
        if (sum(Dall.IPI[i,:]<120)):
            Dall.isError[i] =1;      

## HousKeeping 
min_freq =  FreqRange[0];
max_freq = FreqRange[1];
frex = np.linspace(min_freq, max_freq,numFreqBins);
BandInfo = {'bandsLab' : ['Delta <4Hz' , 'Theta 5-8Hz' , 'Alpha 9-16Hz' , 'Beta 17-36Hz' , 
                          'L-Gamma 37-70Hz' , 'H-Gamma 70-130' , 'NoBandLand 130-180'],
    'bands':[[0, 4], [5, 8], [9, 16], [17, 36], [37, 70], [70, 110], [110 ,179]] , 
    'bandid':[]};
# find the indiceis is frex that correspond to each band
for b in range(len(BandInfo['bands'])):
    Fid = next(x[0] for x in enumerate(frex) if x[1]>BandInfo['bands'][b][0]);
    Lid = next(x[0] for x in enumerate(frex) if x[1]>BandInfo['bands'][b][1])-1;
    BandInfo['bandid'].append ([Fid , Lid]);

# Import the BLock Info
BlockInfo = open_workbook(mainDir+'BlockInfo.xlsx')
df = pd.read_excel(open(mainDir+'BlockInfo.xlsx','rb'))
[~, ~, BlockInfo] = xlsread([mainDir , 'BlockInfo.xlsx'],'Sheet1');
BlockInfo = BlockInfo(2:end,:);
BlockInfo(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),BlockInfo)) = {''};

idx = cellfun(@ischar, BlockInfo);
BlockInfo(idx) = cellfun(@(x) string(x), BlockInfo(idx), 'UniformOutput', false);
Fs = 1024;
Fs_ds = Fs/DownsampleRate;
clearvars idx;
Dout = [];


## chop up the EEG data into trials, and filter out the power line noise
switch what
    case 'ParseEEG-freqDecomp'
        
        ## preprocess EEG and filtering
        
        ChanLabels = ChanLabels(Channels);
        # Q : quality factor is the center frequency divided by the bandwidth.
        # Q = 35;
        # BW = Fo/(Fs/2);
        # [b,a] = iircomb(10,BW,'notch');
        Fo = 60;
        Fs = 1024;
        wo = Fo/(Fs/2);  bw = wo/35; # notch to eliminate 60 Hz
        [b,a] = iirnotch(wo,bw);
        
        wo = 2*Fo/(Fs/2);  bw = wo/60;# notch to eliminate the first harmonic
        [c,d] = iirnotch(wo,bw);
        # definitions, selections...
        
        ## multiple blocks are sotred in the same file. so avoid loading them up multiple times.
        fName1 =  BlockInfo{1,4};
        tn = 1;
        load(fName1);
        
        
        ##
        Events = [];
        for i = 1:size(BlockInfo , 1)
            clear Pall
            
            D = getrow(Dall , Dall.BN == i);
            fName = BlockInfo{i,4};
            if ~strcmp(fName , fName1)
                load(fName);
                # filter the power line noise out of the whole data and then chopp it up
                
                fName1 = fName;
            end
            #     extract the data for the block at hand
            BlockRang = [BlockInfo{i,2} : BlockInfo{i,3}];
            if ~ismember(BlockRang , -1)
                Beeg = getrow(Data , Channels);
                Beeg.values = Beeg.values(:,BlockRang);
                # get the indecies for starts ans ends of the trials
                marker = Beeg.values(find(strcmp(Beeg.label , 'TTL')) , :);
                marker = [0 diff(marker <-2*10^6)];
                for ch = 1:size(Beeg.values , 1)
                    B = Beeg.values(ch , :);
                    A = filter(b,a , B);
                    Beeg.values(ch , :) = A;
                end
                start_tr = find(marker == 1); # starts of trials
                # right now the end TTl pulse is being sent by the releas eof
                # the last finger, so is not aligned to the last press time. so
                # better define it this way for trials with presses
                end_tr = find(marker == -1);  # ends of trials
                for tn = 3:length(start_tr)
                    tn
                    if ~isnan(D.AllPressTimes(tn , D.seqlength(tn)))
                        end_tr(tn) = start_tr(tn) + Fs*(D.AllPressTimes(tn , D.seqlength(tn))/1000)';
                    end
                end
                start_tr = floor(start_tr / DownsampleRate);
                end_tr = floor(end_tr / DownsampleRate);
                for ch = 1:size(Beeg.values , 1)
                    Chname{ch} = ['RawEEGpower',num2str(ch)];
                    [REG, BandInfo] = secog_waveletPSD(Beeg.values(ch , :) , Fs , 'DownsampleRate' , DownsampleRate);
                    # normalize each trial to baseline : TimeDelay ms before the stim  onset
                    for tr = 1:length(start_tr)
                        [ch tr]
                        X = nanmean(REG(:,start_tr(tr)-floor(Fs_ds*TimeDelay):start_tr(tr)) , 2);
                        eval(['D.decompBL{tr,1}(ch , :,:)  = X;']);
                        
                        X = REG(:,start_tr(tr) : end_tr(tr));
                        eval(['D.decompTR{tr,1}(ch , :,:)  = X;']);
                        
                        X = REG(:,start_tr(tr)-floor(Fs_ds*TimeDelay) : start_tr(tr)-1);
                        eval(['D.decompBefTR{tr,1}(ch , :,:)  = X;']);
                        
                        X = REG(:,end_tr(tr)+1 : end_tr(tr)+floor(Fs_ds*2*TimeDelay));
                        eval(['D.decompAftTR{tr,1}(ch , :,:)  = X;']);
                    end
                    disp(['PSD calculation for block ' , num2str(i) , ', channel ' , num2str(ch) , ' completed'])
                    clear REG
                end
                if length(start_tr)<length(D.TN)
                    missing = abs(length(start_tr)-length(D.TN));
                    missingTR = [length(D.TN)-(missing-1) : length(D.TN)];
                    for mtr = missingTR
                        eval(['D.decompBL{mtr,1}  = NaN']);
                        eval(['D.decompTR{mtr,1}  = NaN;']);
                        eval(['D.decompBefTR{mtr,1}= NaN;']);
                        eval(['D.decompAftTR{mtr,1} = NaN;']);
                    end
                end
                
                
                ## save the unbinned data in separate blocks for managability in size
                #             save the raw unbinned psd
                for tn = 1 :length(D.TN)
                    # prepare the individual trials to be saved as the
                    # fileds of a structre to make loading easier
                    eval(['Pall.DEC',num2str(tn),'.decompBL    = D.decompBL{', num2str(tn), '};']);
                    eval(['Pall.DEC',num2str(tn),'.decompTR    = D.decompTR{', num2str(tn), '};']);
                    eval(['Pall.DEC',num2str(tn),'.decompBefTR = D.decompBefTR{', num2str(tn), '};']);
                    eval(['Pall.DEC',num2str(tn),'.decompAftTR = D.decompAftTR{', num2str(tn), '};']);
                end
            else
                for tn = 1 :length(D.TN)
                    # prepare the individual trials to be saved as the
                    # fileds of a structre to make loading easier
                    eval(['Pall.DEC',num2str(tn),'.decompBL    = NaN;']);
                    eval(['Pall.DEC',num2str(tn),'.decompTR    = NaN;']);
                    eval(['Pall.DEC',num2str(tn),'.decompBefTR = NaN;']);
                    eval(['Pall.DEC',num2str(tn),'.decompAftTR = NaN;']);
                end
                
            end
            saveName = [saveDir,'Raw_Decomp_B',num2str(i) ,'.mat'];
            save(saveName,'-struct','Pall', '-v7.3');
            
        end
        
        ## Filter Visualization
        # h = fvtool(b,a);
        # h.Fs = Fs;
        # h.FrequencyRange='[-Fs/2, Fs/2)';
        # zplane(b,a)
        #
        # h = fvtool(c,d);
        # h.Fs = Fs;
        # h.FrequencyRange='[-Fs/2, Fs/2)';
        # zplane(c,d)
        # filt_A = Dall.EEG{tn}(10,:);
        # t = 0:(1/Fs):(length(A)/Fs) - (1/Fs);
        # figure('color' , 'white')
        # subplot(2,1,1)
        # periodogram(A(1,:),[],length(A(1,:)),Fs,'power')
        # subplot(2,1,2)
        # periodogram(filt_A(1,:),[],length(filt_A(1,:)),Fs,'power')
    case 'ParseEEG-calc_norm_PSD'
        ## preprocess EEG and filtering
        load([saveDir , 'ChanLabels.mat']);
        ChanLabels = ChanLabels(Channels);
        # Q : quality factor is the center frequency divided by the bandwidth.
        # Q = 35;
        # BW = Fo/(Fs/2);
        # [b,a] = iircomb(10,BW,'notch');
        Fo = 60;
        Fs = 1024;
        wo = Fo/(Fs/2);  bw = wo/35; # notch to eliminate 60 Hz
        [b,a] = iirnotch(wo,bw);
        Dall  = secog_addEventMarker(Dall, subjNum, Fs_ds , 'addEvent' , 'NumWarpSampFast' , NumWarpSampFast, 'NumWarpSampSlow'  ,NumWarpSampSlow)';
        Events = Dall;
        save([saveDir , 'AllData_Events.mat'] , 'Events');
        wo = 2*Fo/(Fs/2);  bw = wo/60;# notch to eliminate the first harmonic
        [c,d] = iirnotch(wo,bw);
        # definitions, selections...
        
        ## multiple blocks are sotred in the same file. so avoid loading them up multiple times.
        fName1 =  BlockInfo{1,4};
        tn = 1;
        load(fName1);
        
        
        ##
        for i = 1:size(BlockInfo , 1)
            clear Pall
            
            D = getrow(Dall , Dall.BN == i);
            fName = BlockInfo{i,4};
            if ~strcmp(fName , fName1)
                load(fName);
                # filter the power line noise out of the whole data and then chopp it up
                
                fName1 = fName;
            end
            #     extract the data for the block at hand
            BlockRang = [BlockInfo{i,2} : BlockInfo{i,3}];
            if ~ismember(BlockRang , -1)
                Beeg = getrow(Data , Channels);
                Beeg.values = Beeg.values(:,BlockRang);
                # get the indecies for starts ans ends of the trials
                marker = Beeg.values(find(strcmp(Beeg.label , 'TTL')) , :);
                marker = [0 diff(marker <-2*10^6)];
                for ch = 1:size(Beeg.values , 1)
                    B = Beeg.values(ch , :);
                    A = filter(b,a , B);
                    Beeg.values(ch , :) = A;
                end
                start_tr = find(marker == 1); # starts of trials
                # right now the end TTl pulse is being sent by the releas eof
                # the last finger, so is not aligned to the last press time. so
                # better define it this way for trials with presses
                end_tr = find(marker == -1);  # ends of trials
                for tn = 3:length(start_tr)
                    if ~isnan(D.AllPressTimes(tn , D.seqlength(tn)))
                        end_tr(tn) = start_tr(tn) + Fs*(D.AllPressTimes(tn , D.seqlength(tn))/1000)';
                    end
                end
                start_tr = floor(start_tr / DownsampleRate);
                end_tr = floor(end_tr / DownsampleRate);
                for ch = 1:size(Beeg.values , 1)
                    Chname{ch} = ['RawEEGpower',num2str(ch)];
                    [REG, BandInfo] = secog_waveletPSD(Beeg.values(ch , :) , Fs , 'DownsampleRate' , DownsampleRate);
                    REG = 10*log10(abs(REG).^2);
                    # normalize each trial to baseline : TimeDelay ms before the stim  onset
                    for tr = 1:length(start_tr)
                        baseline = nanmean(REG(:,start_tr(tr)-floor(Fs_ds*TimeDelay):start_tr(tr)) , 2);
                        X = REG(:,start_tr(tr)-floor(Fs_ds*TimeDelay) : end_tr(tr)+floor(Fs_ds*2*TimeDelay));
                        X = (X - repmat(baseline , 1,size(X,2)))./repmat(baseline , 1,size(X,2));
                        eval(['D.PSD{tr,1}(ch , :,:)  = X;']);
                    end
                    disp(['PSD calculation for block ' , num2str(i) , ', channel ' , num2str(ch) , ' completed'])
                    clear REG
                end
                # in case some trials at the end didnt get recorded on the marker
                if length(start_tr)<length(D.TN)
                    missing = abs(length(start_tr)-length(D.TN));
                    missingTR = [length(D.TN)-(missing-1) : length(D.TN)];
                    for mtr = missingTR
                        eval(['D.PSD{mtr,1}  = NaN;']);
                    end
                end
            else
                missing = length(D.TN);
                missingTR = [length(D.TN)-(missing-1) : length(D.TN)];
                for mtr = missingTR
                    eval(['D.PSD{mtr,1}  = NaN;']);
                end
            end
            
            ## find the event markers and normalize the power to TimeDelay ms before the stimulus came on - or press
            # complete the structure with behavior again
            D1 = getrow(Dall , Dall.BN == i);
            D1.Pow_Norm_stim = D.PSD;
            D = D1; clear D1

            ## save the unbinned data in separate blocks for managability in size
            #             save the raw unbinned psd
            for tn = 1 :length(D.TN)
                # prepare the individual trials to be saved as the
                # fileds of a structre to make loading easier
                eval(['Pall.PSD',num2str(tn),' = D.Pow_Norm_stim{', num2str(tn), '};']);
            end
            
            saveName = [saveDir,'Raw_PSD_B',num2str(i) ,'.mat'];
            save(saveName,'-struct','Pall', '-v7.3');
            
            # then band-average and stack up
            for tn = 1 :length(D.TN)
                if ~isnan(D.Pow_Norm_stim{tn})
                    temp = D.Pow_Norm_stim{tn};
                    for b =1:length(BandInfo.bandid)
                        P.Pow_Norm_stim{tn,1}(:,b, :) =  nanmean(temp(:,BandInfo.bandid{b}(1) : BandInfo.bandid{b}(2),:) , 2);
                    end
                else
                    P.Pow_Norm_stim{tn,1} = NaN;
                end
            end
            Dout = addstruct(Dout , P);
            clear D P
        end
        load([saveDir , 'AllData_Behav.mat'])
        Dall.Pow_Norm_stim = Dout.Pow_Norm_stim;
        Pall = Dall;
        clear Dall
        saveName = [saveDir,'AllData_PSD_StimNorm.mat'];
        save(saveName , 'Pall' , '-v7.3');
        
    case 'TimeWarpPSD_Raw_Binned'
        # this case loads up it's own input so the Dall structure can be left blac
        load([saveDir,'AllData_Events.mat']);
        Fs = 1024;
        Fs_ds = floor(Fs/DownsampleRate);
        # find average event markers
        Events  = secog_addEventMarker(Dall,subjNum, Fs_ds, 'addEvent' , 'NumWarpSampFast' , NumWarpSampFast, 'NumWarpSampSlow'  ,NumWarpSampSlow);
        E  = secog_addEventMarker(Events, subjNum, Fs_ds , 'CalcAveragePattern' , 'NumWarpSampFast' , NumWarpSampFast, 'NumWarpSampSlow'  ,NumWarpSampSlow)';
        save([saveDir , 'AllData_AvgMarker.mat'] , 'E');
        BN = unique(Events.BN);
        
        Pall_binned = [];
        for bn = 1:length(BN)
            bg = 1;
            mem = 0;
            while mem == 0
                if ismember(BN(bn) , E.blockGroups{bg})
                    mem = 1;
                    BG = bg;
                end
                bg = bg+1;
            end
            
            clear D P
            D = getrow(Events , Events.BN == BN(bn));
            Bname = [saveDir,'Raw_PSD_B',num2str(bn) ,'.mat'];
            P = load(Bname);
            for tn = 1 :length(D.TN)
                # prepare the individual trials to be saved as the
                # fileds of a structre to make loading easier
                eval(['D.Pow_Norm_stim{',num2str(tn),'} = P.PSD', num2str(tn), ';']);
            end
            # this loads up the data strutre where the raw PSD for that block is already stored
            # 'ParseEEG-calcPSD'
            # Use the average patterns of block groups to warp them
            nPSD = D.Pow_Norm_stim;
            
            # set the sequence length for "*******" trilas to one
            D.seqlength(D.seqNumb==5) = 1;
            D.NumWarpSamp = NumWarpSampSlow* ones(size(D.TN));
            D.NumWarpSamp(logical(D.Fast)) = NumWarpSampFast;
            for tn = 1:length(D.TN)
                if length(find(D.NormEventMarker{tn})) == D.seqlength(tn) +1 & ~D.isError(tn) & ...
                        length(D.EventMarker{tn})>D.NumWarpSamp(tn) & ~isnan(nPSD{tn})
                    # check which block group this block falls into
                    
                    A = nPSD{tn};
                    # the normalized event markers in trial tn
                    idx = [0 find(D.NormEventMarker{tn}) D.NumWarpSamp(tn)];
                    
                    # find the row number corresponding to the seqNumb
                    sn = find(E.SN{BG}==D.seqNumb(tn));
                    # average normalized time stamps for the sn , BG
                    if tn>2
                        idn = [0 E.NEM{BG}(sn , ~isnan(E.NEM{BG}(sn ,:))) D.NumWarpSamp(tn)];
                        diffNEM = diff(idn);
                    else
                        # for "* * * * * * * *" trials, the average pattern is the same as the trial since there is no variability
                        idn = idx;
                        diffNEM = diff(idn);
                    end
                    
                    for e = 2:length(idn)
                        # make sure that each
                        idd   = linspace(1 , [idx(e)+1 - (idx(e-1) + 1)] , diffNEM(e-1));
                        for ch = 1:size(A,1)
                            for fi=1:size(A,2)
                                #                                 idd = floor(linspace(1, length([1:idx(e) - idx(e-1)]) , length(idd)));
                                D.PSD_stim{tn ,1}(ch,fi,idn(e-1)+1:idn(e)) = interp1([1:idx(e) - idx(e-1)] , squeeze(A(ch,fi,idx(e-1)+1:idx(e))) , idd);
                            end
                        end
                    end
                else
                    D.PSD_stim{tn ,1} = nan(size(A,1) , length(frex) , D.NumWarpSamp(tn));
                end
            end
            for tn = 1:length(D.PSD_stim)
                trialName = ['Pall.PSD',num2str(tn)];
                eval([trialName, ' = D.PSD_stim{' , num2str(tn) , '};'])
            end
            
            saveName = [saveDir,'warped_PSD_B',num2str(bn) ,'.mat'];
            save(saveName,'-struct','Pall', '-v7.3');
            # stace the binned PSDs
            
            
            D = getrow(Events , Events.BN == BN(bn));
            
            for tn = 1 :length(D.TN)
                PSD = eval(['Pall.PSD' , num2str(tn)]);
                
                for b = 1:length(BandInfo.bands)
                    D.Pow_Norm_stim{tn,1}(:,b, :) =  nanmean(PSD(:,BandInfo.bandid{b}(1) : BandInfo.bandid{b}(2),:) , 2);
                end
                
            end
            Pall_binned = addstruct(Pall_binned, D);
            clear D
            disp(['PSD warping for block ' , num2str(bn) , ' complete'])
        end
        Pall = Pall_binned;
        saveName = [saveDir,'AllData_PSD_Warped.mat'];
        save(saveName,'Pall', '-v7.3');
        Dout = Pall;

    case 'TimeWarpPSD_Raw_Binned_seqType'
        # this case loads up it's own input so the Dall structure can be left blac
        ## the goal here is to get average time pattern for general sequence types regardless of finger
        # so all the slow single fingers, fast singel finger, random, structured, triplets, quadruples
        # so we will change the SeqNumb and average the average times of the SeqNumbs
        # within the same type
        
        
        load([saveDir,'AllData_Events.mat']);
        
        
        Fs = 1024;
        Fs_ds = floor(Fs/DownsampleRate);
        # find average event markers
        Events  = secog_addEventMarker(Dall,subjNum, Fs_ds, 'addEvent' , 'NumWarpSampFast' , NumWarpSampFast, 'NumWarpSampSlow'  ,NumWarpSampSlow);
        E  = secog_addEventMarker(Events, subjNum, Fs_ds , 'CalcAveragePattern_seqType' , 'NumWarpSampFast' , NumWarpSampFast, 'NumWarpSampSlow'  ,NumWarpSampSlow)';
        save([saveDir , 'AllData_AvgMarker_SeqType.mat'] , 'E');
        
        # Define sequence numbers and their transformations:
        SeqTrans = [5 11 22 33 44 55 0 1 2 3 4 103 104 203 204;...
            100 10 10 10 10 10 20 30 30 30 30 40 50 40 50];
        for sn = 1:length(SeqTrans)
            id  = Events.seqNumb == SeqTrans(1 , sn);
            Events.seqNumb(id) = SeqTrans(2 , sn);
        end
        
        BN = unique(Events.BN);
        Pall_binned = [];
        for bn = 1:length(BN)
            bg = 1;
            mem = 0;
            while mem == 0
                if ismember(BN(bn) , E.blockGroups{bg})
                    mem = 1;
                    BG = bg;
                end
                bg = bg+1;
            end
            clear D P
            D = getrow(Events , Events.BN == BN(bn));
            Bname = [saveDir,'Raw_PSD_B',num2str(bn) ,'.mat'];
            P = load(Bname);
            for tn = 1 :length(D.TN)
                # prepare the individual trials to be saved as the
                # fileds of a structre to make loading easier
                eval(['D.Pow_Norm_stim{',num2str(tn),'} = P.PSD', num2str(tn), ';']);
            end
            # this loads up the data strutre where the raw PSD for that block is already stored
            # 'ParseEEG-calcPSD'
            # Use the average patterns of block groups to warp them
             nPSD = D.Pow_Norm_stim;
            
            # set the sequence length for "*******" trilas to one
            D.seqlength(D.seqNumb==100) = 1;
            D.NumWarpSamp = NumWarpSampSlow* ones(size(D.TN));
            D.NumWarpSamp(logical(D.Fast)) = NumWarpSampFast;
            for tn = 1:length(D.TN)
                if length(find(D.NormEventMarker{tn})) == D.seqlength(tn) +1 & ~D.isError(tn) & ...
                        length(D.EventMarker{tn})>D.NumWarpSamp(tn) & ~isnan(nPSD{tn})
                    # check which block group this block falls into
                    
                    A = nPSD{tn};
                    # the normalized event markers in trial tn
                    idx = [0 find(D.NormEventMarker{tn}) D.NumWarpSamp(tn)];
                    
                    # find the row number corresponding to the seqNumb
                    sn = find(E.SN{BG}==D.seqNumb(tn));
                    # average normalized time stamps for the sn , BG
                    if tn>2
                        idn = [0 E.NEM{BG}(sn , ~isnan(E.NEM{BG}(sn ,:))) D.NumWarpSamp(tn)];
                        diffNEM = diff(idn);
                    else
                        # for "* * * * * * * *" trials, the average pattern is the same as the trial since there is no variability
                        idn = idx;
                        diffNEM = diff(idn);
                    end
                    
                    for e = 2:length(idn)
                        # make sure that each
                        idd   = linspace(1 , [idx(e)+1 - (idx(e-1) + 1)] , diffNEM(e-1));
                        for ch = 1:size(A,1)
                            for fi=1:size(A,2)
                                #                                 idd = floor(linspace(1, length([1:idx(e) - idx(e-1)]) , length(idd)));
                                D.PSD_stim{tn ,1}(ch,fi,idn(e-1)+1:idn(e)) = interp1([1:idx(e) - idx(e-1)] , squeeze(A(ch,fi,idx(e-1)+1:idx(e))) , idd);
                            end
                        end
                    end
                else
                    D.PSD_stim{tn ,1} = nan(size(A,1) , length(frex) , D.NumWarpSamp(tn));
                end
            end
            for tn = 1:length(D.PSD_stim)
                trialName = ['Pall.PSD',num2str(tn)];
                eval([trialName, ' = D.PSD_stim{' , num2str(tn) , '};'])
            end
            
            saveName = [saveDir,'warped_PSD_B_SeqType',num2str(bn) ,'.mat'];
            save(saveName,'-struct','Pall', '-v7.3');
            # stace the binned PSDs
            
            
            D = getrow(Events , Events.BN == BN(bn));
            
            for tn = 1 :length(D.TN)
                PSD = eval(['Pall.PSD' , num2str(tn)]);
                
                for b = 1:length(BandInfo.bands)
                    D.Pow_Norm_stim{tn,1}(:,b, :) =  nanmean(PSD(:,BandInfo.bandid{b}(1) : BandInfo.bandid{b}(2),:) , 2);
                end
                
            end
            Pall_binned = addstruct(Pall_binned, D);
            clear D
            disp(['PSD warping for block ' , num2str(bn) , ' complete'])
        end
        
        Pall = Pall_binned;
        saveName = [saveDir,'AllData_PSD_Warped_SeqType.mat'];
        save(saveName,'Pall', '-v7.3');
        Dout = Pall;
    case 'AlignEvents_SeqType'
        load([saveDir,'AllData_Events.mat']);
        Fs = 1024;
        Fs_ds = floor(Fs/DownsampleRate);
        # find average event markers
        Events  = secog_addEventMarker(Dall,subjNum, Fs_ds, 'addEvent' , 'NumWarpSampFast' , NumWarpSampFast, 'NumWarpSampSlow'  ,NumWarpSampSlow);
        E  = secog_addEventMarker(Events, subjNum, Fs_ds , 'CalcAveragePattern' , 'NumWarpSampFast' , NumWarpSampFast, 'NumWarpSampSlow'  ,NumWarpSampSlow);
        BN = unique(Events.BN);
        
        Pall_binned = [];
        for bn = 1:length(BN)
            bg = 1;
            mem = 0;
            while mem == 0
                if ismember(BN(bn) , E.blockGroups{bg})
                    mem = 1;
                    BG = bg;
                end
                bg = bg+1;
            end
            
            clear D P
            D = getrow(Events , Events.BN == BN(bn));
            Bname = [saveDir,'Raw_Decomp_B',num2str(bn) ,'.mat'];

        
            # set the sequence length for "*******" trilas to one
            D.seqlength(D.seqNumb==100) = 1;
            D.NumWarpSamp = NumWarpSampSlow* ones(size(D.TN));
            D.NumWarpSamp(logical(D.Fast)) = NumWarpSampFast;
            numEvents = max(D.seqlength) + 1;
            for tn = 1:2 # Null Trails
                for e = 1:length(numEvents )
                    PSD_Aligned{tn,e} = [];
                end
            end
            for tn = 3:length(D.TN)

                [bn tn]
                trialName = ['DEC',num2str(tn)];
                A = load(Bname , trialName);
                temptr = eval(['10*log10(abs(A.',trialName,'.decompTR).^2);']);
                AvgPowTR = squeeze(mean(temptr  , 3));
                tempbl = eval(['10*log10(abs(A.',trialName,'.decompBefTR).^2);']);
                AvgPowBL = squeeze(mean(tempbl  , 3));
                tempaf = eval(['10*log10(abs(A.',trialName,'.decompAftTR).^2);']);
                AvgPowAF = squeeze(mean(tempaf  , 3));
                if  ~D.isError(tn) & sum(sum(~isnan(D.EventMarker{tn}))) & sum(sum(~isnan(AvgPowTR)))
                    # check which block group this block falls into
                    tempAll = cat(3 , tempbl , temptr , tempaf);
                    A  = tempAll  - repmat(AvgPowBL , 1,1,size(tempAll , 3));
                    # the normalized event markers in trial tn
                    idx = find(D.EventMarker{tn});
                    idx = [idx , idx(end) + 30]; # take the last event with 300 ms after the last press
                    for e = 1:length(idx)-1
                        PSD_Aligned{tn,e} = A(:,:,idx(e)-10:idx(e+1) - 10);
                    end
                else
                    PSD_Aligned(tn,:) = PSD_Aligned(1,:);
                end
            end
            for tn = 1:length(PSD_Aligned)
                trialName = ['Pall.PSD_Aligned',num2str(tn)];
                eval([trialName, ' = PSD_Aligned(' , num2str(tn) , ',:);'])
            end
            saveName = [saveDir,'EventAligned_PSD_B',num2str(bn) ,'.mat'];
            save(saveName,'-struct','Pall', '-v7.3');
            # stace the binned PSDs
            Dout = Pall;
            clear PSD_Aligned Pall
            
            disp(['PSD Alignment for block ' , num2str(bn) , ' complete'])
        end        
    case 'AlignEvents_SeqType_Warped'
        load([saveDir,'AllData_Events.mat']);
        Fs = 1024;
        Fs_ds = floor(Fs/DownsampleRate);
        # find average event markers
        Events  = secog_addEventMarker(Dall,subjNum, Fs_ds, 'addEvent' , 'NumWarpSampFast' , NumWarpSampFast, 'NumWarpSampSlow'  ,NumWarpSampSlow);
        E  = secog_addEventMarker(Events, subjNum, Fs_ds , 'CalcAveragePattern_seqType' , 'NumWarpSampFast' , NumWarpSampFast, 'NumWarpSampSlow'  ,NumWarpSampSlow)';
        save([saveDir , 'AllData_AvgMarker_SeqType.mat'], 'E');
        BN = unique(Events.BN);
        
        
        SeqTrans = [5 11 22 33 44 55 0 1 2 3 4 103 104 203 204;...
            100 10 10 10 10 10 20 30 30 30 30 40 50 40 50];
        for sn = 1:length(SeqTrans)
            id  = Events.seqNumb == SeqTrans(1 , sn);
            Events.seqNumb(id) = SeqTrans(2 , sn);
        end
        
        BN = unique(Events.BN);
       
        
        Pall_binned = [];
        for bn = 1:length(BN)
            bg = 1;
            mem = 0;
            while mem == 0
                if ismember(BN(bn) , E.blockGroups{bg})
                    mem = 1;
                    BG = bg;
                end
                bg = bg+1;
            end
            E1 = getrow(E , BG);
            clear D P
            D = getrow(Events , Events.BN == BN(bn));
            Bname = [saveDir,'warped_PSD_B',num2str(bn) ,'.mat'];

        
            # set the sequence length for "*******" trilas to one
            D.seqlength(D.seqNumb==5) = 1;
            D.NumWarpSamp = NumWarpSampSlow* ones(size(D.TN));
            D.NumWarpSamp(logical(D.Fast)) = NumWarpSampFast;
            numEvents = max(D.seqlength) + 1;
            for tn = 1:2 # Null Trails
                for e = 1:numEvents
                    PSD_Aligned{tn,e} = [];
                end
            end
            for tn = 3:length(D.TN)
                sn = find(E1.SN{1} == D.seqNumb(tn));
                trialName = ['PSD',num2str(tn)];
                A = load(Bname , trialName);
                eval(['A = A.' , trialName , ';']);
                if   ~D.isError(tn) & sum(sum(~isnan(D.EventMarker{tn}))) & sum(sum(~isnan(A)))
                    # check which block group this block falls into
                    # the normalized event markers in trial tn
                    idx = E1.NEM{1}(sn , 1:D.seqlength(tn)+1);
                    idx = [idx , idx(end) + 10]; # take the last event with 300 ms after the last press
                    for e = 1:length(idx)-1
                        PSD_Aligned{tn,e} = A(:,:,idx(e)-5:idx(e+1) - 5);
                    end
                else
                    PSD_Aligned(tn,:) = PSD_Aligned(1,:);
                end
            end
            for tn = 1:length(PSD_Aligned)
                trialName = ['Pall.PSD_Aligned',num2str(tn)];
                eval([trialName, ' = PSD_Aligned(' , num2str(tn) , ',:);'])
            end
            saveName = [saveDir,'EventAligned_WarpedPSD_B',num2str(bn) ,'.mat'];
            save(saveName,'-struct','Pall', '-v7.3');
            # stace the binned PSDs
            Dout = Pall;
            clear PSD_Aligned Pall
            
            disp(['PSD Alignment for block ' , num2str(bn) , ' complete'])
        end


end
