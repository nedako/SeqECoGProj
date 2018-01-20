function out  = secog_visualize(Dall , subjnum, what, distance, calc , day , rep, GroupCode)
%%  distances:
% 'euclidean'	      Euclidean distance (default).
% 'squaredeuclidean'  Squared Euclidean distance. (This option is provided for efficiency only. It does not satisfy the triangle inequality.)
% 'seuclidean'	      Standardized Euclidean distance. Each coordinate difference between rows in X and Y is scaled by dividing by the corresponding element of the standard deviation computed from X, S=nanstd(X). To specify another value for S, use D = PDIST2(X,Y,'seuclidean',S).
% 'cityblock'	      City block metric.
% 'minkowski'         Minkowski distance. The default exponent is 2. To compute the distance with a different exponent, use D = pdist2(X,Y,'minkowski',P), where the exponent P is a scalar positive value.
% 'Chebychev'	      Chebychev distance (maximum coordinate difference).
% 'mahalanobis'	      Mahalanobis distance, using the sample covariance of X as computed by nancov. To compute the distance with a different covariance, use D = pdist2(X,Y,'mahalanobis',C) where the matrix C is symmetric and positive definite.
% 'cosine'	          One minus the cosine of the included angle between points (treated as vectors).
% 'correlation'	      One minus the sample correlation between points (treated as sequences of values).
% 'spearman'	      One minus the sample Spearman's rank correlation between observations, treated as sequences of values.
% 'hamming'           Hamming distance, the percentage of coordinates that differ.
% 'jaccard'           One minus the Jaccard coefficient, the percentage of nonzero coordinates that differ.%     case 'chunk_est_instance'
%%  Cases
%     case 'chunk_est_instance'
%     case 'eye_vel_instance'
%     case 'eye_mainSequence'
%     case 'eye_vs_press_instance'
%     case 'eyePos_eyeVel'
%     case 'IPI_dist'
%     case 'IPI_ttest_rand'
%     case 'chunk_dist'
%     case 'Mychunk_dist'
%     case 'Avg_pattern_sh3'
%     case 'eye_vel_chunkplace'
%     case 'eye_vel_seqplace'
%     case 'eye_pos_seqplace'


%     case 'eyepress_pos_traces'
%     case 'eyepress_pos_distances'
%     case 'eyepress_pos_avg_distances'


%     case 'eyepress_vel_traces'
%     case 'eyepress_vel_distances'
%     case 'eyepress_vel_avg_distances'


%     case 'run_f-tests_mt'
%     case 'dtw'
%     case 'eyeEndPos_ChunkLength'
%     case 'eyeFrstPos_ChunkLength'


%     case 'crossvaldist_pos'
%     case 'crossvaldist_vel'


%     case 'crossvaldist_pos_presstim'
%     case 'crossvaldist_vel_presstime'
%     case 'crossval_IPI_dist'
%     case 'crossvaldist_chunk'
%     case 'perveiw_benefit'
%     case 'saccades'

%%

baseDir = '/Users/nkordjazi/Documents/SeqECoG/analyze';
% baseDir = '/Users/nkordjazi/Documents/SeqEye/se1/SeqEye1/se1_data/analyze';
% subj_name = {'XW' , 'ML' , 'DS' , 'BM' , 'HK' , 'BW' 'XX'};
subj_name = {'P2', 'XX'};
if subjnum == length(subj_name)
      subjnum = 1:length(subj_name)-1;
end


load([baseDir , '/CMB.mat'])
%load([baseDir , '/se1_all.mat'])
days = {1 ,2 ,3 ,4 ,[2:4] ,[3:4] [1:4]};
% subjnum = 2;

tid = {[1:250] [251 :750] [751:1000] [1:1000]};
switch what
   case 'IPI_ttest_rand'
        
        %% within and between IPI ttest with randomization
        %% IPI per day
%         
        ANA = getrow(Dall , Dall.isgood & ~Dall.isError & ismember(Dall.seqNumb , [0:4]) & Dall.Group == 1);
        ANA.seqNumb(ANA.seqNumb > 0) = 1;
        BN = unique(ANA.BN);        
        for bn = 1:length(BN)
            temp = getrow(ANA , ANA.BN == BN(bn));
            boxplot(reshape(temp.IPI , numel(temp.IPI),1) , reshape(temp.IPIarrangement  , numel(temp.IPIarrangement),1))
            hold on
        end

        
        ANA1 = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.seqNumb , [1:4])  & ~Dall.isError &  ismember(Dall.Rep , rep) );
        ANA0 = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.seqNumb , [0])  & ~Dall.isError &  ismember(Dall.Rep , rep)); 
        
       
        IPInums = 1:6;
        BN = unique(ANA1.BN);
        clear IPIwithin IPIbetween
        IPIw = [];
        IPIb = [];
        ANA1.MT = ANA1.AllPressTimes(:,7) - ANA1.AllPressTimes(:,1);
        ANA0.MT = ANA0.AllPressTimes(:,7) - ANA0.AllPressTimes(:,1);
        for bn = 1:length(BN)
            IPIwithin{bn} = [];
            IPIbetween{bn} = [];
            ANA = getrow(ANA1 , ANA1.BN == BN(bn));
            for tn = 1:length(ANA.TN)
                Arr = ANA.IPIarrangement(tn , :);
                IPIwithin{bn}  = [IPIwithin{bn} ;[ANA.IPI(tn,Arr == 2)' BN(bn)*ones(size(ANA.IPI(tn,Arr == 2)')) ANA.SN(tn)*ones(size(ANA.IPI(tn,Arr == 2)')) IPInums(Arr == 2)']];
                IPIbetween{bn} = [IPIbetween{bn} ;[ANA.IPI(tn,Arr == 1)' BN(bn)*ones(size(ANA.IPI(tn,Arr == 1)')) ANA.SN(tn)*ones(size(ANA.IPI(tn,Arr == 1)')) IPInums(Arr == 1)']];
            end
            IPIw = [IPIw ;  IPIwithin{bn}];
            IPIb = [IPIb ; IPIbetween{bn}];
%             D(bn) = ANA.Day(1);
        end
        
        
%         out.IPI = anovaMixed([IPIw(:,1) ; IPIb(:,1)] , [IPIw(:,3) ; IPIb(:,3)],'within',[0*IPIw(:,1) ; 1+0*IPIb(:,1)],{'Within/Between'},'intercept',1)  ;
        idw = ~ismember(IPIw(:,end),[1 5]);
        idb = ~ismember(IPIb(:,end),[1 5]);
        out.IPI = anovan([IPIw(idw,1) ; IPIb(idb,1)] , [[0*IPIw(idw,1) ; 1+0*IPIb(idb,1)]],'varnames',{'Within/Between'} , 'display' , 'off' , 'model' , 'full');
        
        A = reshape(ANA0.IPI(:,2:6) , numel(ANA0.IPI(:,2:6)),1);
        
        
        out.IPI_WithinVsRand = anovan([A ; IPIw(idw,1)] , [[0*A ; 1+0*IPIw(idw,1)]],'varnames',{'Within/Between'} , 'display' , 'off' , 'model' , 'full');
        out.IPI_betweenVsRand = anovan([A ; IPIb(idb,1)] , [[0*A ; 1+0*IPIb(idb,1)]],'varnames',{'Within/Between'} , 'display' , 'off' , 'model' , 'full');
        
        h1 = figure('color' , 'white');
        [xb,ePLOTb,ERRORb] = lineplot(IPIb(idb,2) , IPIb(idb,1));
        hold on 
        [xw,ePLOTw,ERRORw] = lineplot(IPIw(idw,2) , IPIw(idw,1));
        A = repmat(ANA0.BN , 1,size(ANA0.IPI(:,2:6) , 2));
        [xr,ePLOTr,ERRORr] = lineplot(reshape(A , numel(A) , 1) , reshape(ANA0.IPI(:,2:6) , numel(ANA0.IPI(:,2:6)) , 1));
        close(h1);
        
        figure('color' , 'white')
        hb = plotshade(xb',ePLOTb,ERRORb,'transp' , .2 , 'patchcolor' , 'b' , 'linecolor' , 'b' , 'linewidth' , 3 , 'linestyle' , ':')
        hold on
        hw = plotshade(xw',ePLOTw , ERRORw,'transp' , .2 , 'patchcolor' , 'm' , 'linecolor' , 'm' , 'linewidth' , 3 , 'linestyle' , ':')  
        hr = plotshade(xr',ePLOTr , ERRORr,'transp' , .2 , 'patchcolor' , 'r' , 'linecolor' , 'r' , 'linewidth' , 3 , 'linestyle' , ':')  
        ax = gca;
        ax.FontSize = 20;
        legend([hb,hw,hr] , {'Between chunk IPIs' , 'Within chunk IPIs' , 'Random IPIs'} , 'FontSize' , 20)
        title(['Within chunk vs Between chunk Chunk IPIs , p(W/B) = ' , num2str(out.IPI) , '   p(W/R) = ' , num2str(out.IPI_WithinVsRand) , '   p(B/R) = ' , num2str(out.IPI_betweenVsRand)],'FontSize' , 16)
        xlabel('Training Blocks', 'FontSize' , 20)
        ylabel('msec'  ,'FontSize' , 20)
        grid on
        
        
%         out.MT = anovaMixed([ANA1.MT ; ANA0.MT] , [ANA1.SN ;ANA0.SN],'within',[0*ANA1.MT ; 1+0*ANA0.MT],{'Random/Chunked'},'intercept',1)  ;
        out.MT = anovan([ANA1.MT ; ANA0.MT] , [0*ANA1.MT ; 1+0*ANA0.MT] , 'display' , 'off' , 'varnames' , 'Rand/Chunked');
        h1 = figure('color' , 'white');
        [xs,PLOTs,ERRORs] = lineplot(ANA1.BN,ANA1.MT);
        hold on 
        [xr,PLOTr,ERRORr] = lineplot(ANA0.BN,ANA0.MT);
        close(h1);
        
        figure('color' , 'white')
        hs = plotshade(xs',PLOTs,ERRORs,'transp' , .2 , 'patchcolor' , 'black' , 'linecolor' , 'k' , 'linewidth' , 3 , 'linestyle' , ':')
        hold on
        hr = plotshade(xr',PLOTr,ERRORr,'transp' , .2 , 'patchcolor' , 'r' , 'linecolor' , 'r' , 'linewidth' , 3 , 'linestyle' , ':')  
        hold on
        ax = gca;
        ax.FontSize = 20;
        title(['Movement Time, p = ' , num2str(out.MT)] ,'FontSize' , 20)
        ylabel('msec'  ,'FontSize' , 20)
        xlabel('Training Blocks', 'FontSize' , 20)
%         ax.XTick = [1:4];
%         ax.XTickLabel = {'Day 0 ' ,'Day 1' , 'Day 2' , 'Day 3'};
        legend([hs ,hr] , {'Chunked' , 'Random'})
        grid on        
    case 'Errors'
        GroupCode = 1;
        errorrate0 = [];
        errorrate1 = [];
        ANA1 = getrow(Dall , ismember(Dall.SN , subjnum) & ismember(Dall.seqNumb , [1:2]) & Dall.isgood & ismember(Dall.Rep , rep) & ismember(Dall.Day , days{day}) & ismember(Dall.Group , GroupCode));
        ANA0 = getrow(Dall , ismember(Dall.SN , subjnum) & ismember(Dall.seqNumb , 0) & Dall.isgood & ismember(Dall.Rep , rep) & ismember(Dall.Day , days{day}) & ismember(Dall.Group , GroupCode));
        uSub = unique(ANA1.SN);
        for sub = 1 : length(uSub)
            ANA1sn = getrow(ANA1 , ismember(ANA1.SN,uSub(sub)));
            ANA0sn = getrow(ANA0 , ismember(ANA0.SN,uSub(sub)));
            
            errorrate0 = [errorrate0 ; [100 * sum(ANA0sn.isError)/length(ANA0sn.TN) ,  uSub(sub) ,0]];
            errorrate1 = [errorrate1 ; [100 * sum(ANA1sn.isError)/length(ANA1sn.TN) ,  uSub(sub) ,1]];
        end
        errorrate_g1 = [[errorrate1 ; errorrate0] ones(length([errorrate1 ; errorrate0]) , 1)];
        
        GroupCode = 2;
        errorrate0 = [];
        errorrate1 = [];
        ANA1 = getrow(Dall , ismember(Dall.SN , subjnum) & ismember(Dall.seqNumb , [1:2]) & Dall.isgood & ismember(Dall.Rep , rep) & ismember(Dall.Day , days{day}) & ismember(Dall.Group , GroupCode));
        ANA0 = getrow(Dall , ismember(Dall.SN , subjnum) & ismember(Dall.seqNumb , 0) & Dall.isgood & ismember(Dall.Rep , rep) & ismember(Dall.Day , days{day}) & ismember(Dall.Group , GroupCode));
        uSub = unique(ANA1.SN);
        for sub = 1 : length(uSub)
            ANA1sn = getrow(ANA1 , ismember(ANA1.SN,uSub(sub)));
            ANA0sn = getrow(ANA0 , ismember(ANA0.SN,uSub(sub)));
            
            errorrate0 = [errorrate0 ; [100 * sum(ANA0sn.isError)/length(ANA0sn.TN) ,  uSub(sub) ,0]];
            errorrate1 = [errorrate1 ; [100 * sum(ANA1sn.isError)/length(ANA1sn.TN) ,  uSub(sub) ,1]];
        end
        errorrate_g2 = [[errorrate1 ; errorrate0] 2*ones(length([errorrate1 ; errorrate0]) , 1)];
        
        error_all = [errorrate_g1 ; errorrate_g2];
        
        out.Err1 = anovaMixed(errorrate_g1(:,1) , errorrate_g1(:,2),'within',errorrate_g1(:,3),{'Random/Chunked'},'intercept',1)  ;
        out.Err2 = anovaMixed(errorrate_g2(:,1) , errorrate_g2(:,2),'within',errorrate_g2(:,3),{'Random/Chunked'},'intercept',1)  ;
      
        h1 = figure('color' , 'white');
        [xcoord1,PLOT1,ERROR1] = lineplot(errorrate_g1(:,3),errorrate_g1(:,1));
        hold on 
        [xcoord2,PLOT2,ERROR2] = lineplot(errorrate_g2(:,3),errorrate_g2(:,1));
        close (h1)
       
        
        
        figure('color' , 'white')
        errorbar(xcoord1,PLOT1,ERROR1, 'LineWidth' , 3);
        hold on
        errorbar(xcoord2,PLOT2,ERROR2, 'LineWidth' , 3);
        grid on
        ylabel('Percent')
        title(['Error rate , p1 = ' ,  num2str(out.Err1.eff(2).p) , '    p2 =', num2str(out.Err2.eff(2).p)])
        hold on
        ax = gca;
        ax.XTick = [0:1];
        ax.XTickLabel = {'Random' , 'Chunked'};
        ax.FontSize = 20;
        legend({'1-Explicite Chunking' , '2-No Explicite chunking'});        
    case 'Error_chunk'
        for tn = 1:length(Dall.TN)
            Dall.MT(tn , 1) = Dall.AllPressTimes(tn , Dall.seqlength(tn)) - Dall.AllPressTimes(tn , 1); 
        end
        err = [];
        err0 = [];
        proberr = [];
        proberr0 = [];
        
        
        err = [];
        err0 = [];
        for GroupCode = 1:2
            errhs = [];
            errhs0 = [];
            ANA  = getrow(Dall ,  ismember(Dall.seqNumb , [1:6]) & Dall.isError & ismember(Dall.Rep , rep) & ismember(Dall.Day , days{day}) & ismember(Dall.Group , GroupCode));
            ANA0 = getrow(Dall ,  ismember(Dall.seqNumb , 0) & Dall.isError & ismember(Dall.Rep , rep) & ismember(Dall.Day , days{day}) & ismember(Dall.Group , GroupCode));
            
            for tn  = 1:length(ANA.TN)
                ANA.ChunkBndry(tn , :) = diff(ANA.ChnkArrang(tn,:));
                a = find(ANA.ChunkBndry(tn , :));
                ANA.ChunkBndry(tn , a-1) = 3;
                ANA.ChunkBndry(tn , end) = 3;
                ANA.ChunkBndry(tn , ANA.ChunkBndry(tn , :) == 0) = 2;
                %                 ANA.ChunkBndry(tn , 1:2) = [-1 -1];  % dont account for the first and last sseqeuce presses
                %                 ANA.ChunkBndry(tn , end-1:end) = [-1 -1];% dont account for the first and last sseqeuce presses
                a = find(ANA.AllPress(tn,:) ~= ANA.AllResponse(tn,:));
                plc = [1 ANA.ChunkBndry(tn , :)];
                err   = [err   ; [a(1) plc(a(1)) , ANA.SN(tn) GroupCode]];
            end
            
            
            for tn  = 1:length(ANA0.TN)
                a = find(ANA0.AllPress(tn,:) ~= ANA0.AllResponse(tn,:));
                err0   = [err0   ; [a(1) 0 , ANA0.SN(tn) GroupCode]];
            end
            uSub = unique(ANA.SN);
            for sub = 1:length(uSub)
                for p = 1:7
                    proberr  = [proberr  ; 100*sum(err(:,1) == p)/size(err,1)   sub GroupCode p];
                    proberr0 = [proberr0 ; 100*sum(err0(:,1) == p)/size(err0,1) sub GroupCode p];
                end
            end
            %                 subplot(2,4,GroupCode)
        end
        figure('color' , 'white')
       
        subplot(2,1,1)
        ax = gca;
        hold on
        histogram(err(err(:,end) == 1 ,2) , 'Normalization' , 'probability');
        title('Error distribution for Group 1')
        ylabel('Probability');
        xlabel('Chunk placement')
        title(['Errors in different chunk placements, GroupCode = ' , num2str(GroupCode)])
        ax.YLim =[0: 1];
        ax.XTick = [1:3];
        ax.XTickLabel = {'First' , 'Middle' , 'Last'};
        
        subplot(2,1,2)
        ax = gca;
        hold on
        histogram(err(err(:,end) == 2 ,2) , 'Normalization' , 'probability');
        title('Error distribution for Group 2')
        ylabel('Probability');
        xlabel('Chunk placement')
        title(['Errors in different chunk placements, GroupCode = ' , num2str(GroupCode)])
        ax.YLim =[0: 1];
        ax.XTick = [1:3];
        ax.XTickLabel = {'First' , 'Middle' , 'Last'};
    case 'IPI_image'
        
        ANA = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.seqNumb , [0:2])  & ~Dall.isError &  ismember(Dall.Rep , rep));
        
        IPI1 = [];
        IPI2 = [];
        for ipi = 1:6
            a = pivottable(ANA.seqNumb, ANA.Group , ANA.IPI(:,ipi),'nanmean');
            IPI1 = [IPI1 a(:,1)];
            IPI2 = [IPI2 a(:,2)];
        end
        figure('color' , 'white')
        subplot(1,2,1)
        imagesc(IPI1 , [200,900])
        colorbar
        hold on
        ax = gca;
        title('IPI distribution for Group 1')
        ax.YTick = [1:3];
        ax.YTickLabel = {'Random' , 'Structure 1' , 'Structure 2'};
        xlabel('IPI')
        
        subplot(1,2,2)
        imagesc(IPI2, [200,900])
        colorbar
        hold on
        ax = gca;
        title('IPI distribution for Group 2')
        ax.YTick = [1:3];
        ax.YTickLabel = {'Random' , 'Structure 1' , 'Structure 2'};
         xlabel('IPI')
        
        out = [];    
       
    case 'chunk_dist' % CityBlock
        %% chunk distances
        ANA = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.Rep , rep) & ismember(Dall.Day , days{day}));
        for tn = 1:length(ANA.TN)
            n = (ANA.AllPressIdx(tn , sum(~isnan(ANA.AllPressIdx(tn , :))))  - ANA.AllPressIdx(tn , 1)) / 1000;
            nIdx(tn , :) = (ANA.AllPressIdx(tn , :) - ANA.AllPressIdx(tn , 1))/n;
            IPI(tn , :) = detrend(diff(nIdx(tn ,:) , 1 , 2) , 'linear' , 2);
        end
        
        for tn = 1:length(ANA.TN)
            thresh = .3 * std(IPI(tn , :));
            [dum , estChnkBndry] = findpeaks(IPI(tn , :));% ,'MinPeakProminence', thresh); % slope of IPIs
            
            if ~isempty(estChnkBndry)
                goodpeak = ones(1, length(estChnkBndry));
                for cb = 1:length(estChnkBndry)
                    if IPI(tn , estChnkBndry(cb)) < nanmean(IPI(tn  ,:)) + thresh
                        goodpeak(cb) = 0;
                    end
                end
                if sum(goodpeak)
                    estChnkBndry = estChnkBndry(logical(goodpeak));
                else
                    estChnkBndry = [];
                end
            end
            
            
            ANA.estChnkBndry(tn , :) = zeros(1, 14);  % first presses of chunks will be 1
            ANA.estChnkBndry(tn , estChnkBndry+1) = 1;
            ANA.estChnkBndry(tn , 1) = 1;
        end
        
        if GroupCode == 1
            ChnkArrang = CMB.DesiredCnhkargmnt(CMB.Group1 , :);
            for j = 1:size(ChnkArrang , 1)
                temp = [];
                temp1 = [];
                for k = 1:length(find(ChnkArrang(j , :)))
                    temp = [temp k*ones(1,ChnkArrang(j,k))];
                    temp1 = [temp1 1:ChnkArrang(j,k)];
                end
                ChnkArrng(j , :) = temp;
                ChnkPlcmnt(j,:)  = temp1;
            end
            IPIarrangement = diff(ChnkArrng , 1 , 2);   % between IPIs will have the index of 1
            IPIarrangement(~IPIarrangement) = 2;         % within IPIs will have the index of 2
        elseif GroupCode == 2
            ChnkArrang = CMB.DesiredCnhkargmnt(CMB.Group2 , :);
            for j = 1:size(ChnkArrang , 1)
                temp = [];
                temp1 = [];
                for k = 1:length(find(ChnkArrang(j , :)))
                    temp = [temp k*ones(1,ChnkArrang(j,k))];
                    temp1 = [temp1 1:ChnkArrang(j,k)];
                end
                ChnkArrng(j , :) = temp;
                ChnkPlcmnt(j,:)  = temp1;
            end
            IPIarrangement = diff(ChnkArrng , 1 , 2); % between IPIs will have the index of 1
            IPIarrangement(~IPIarrangement) = 2;       % within IPIs will have the index of 2
        end
        
        
        
        %% the boundry distance
        temp = diff(ChnkPlcmnt,1,2);
        temp(temp<0) = 0;
        chbndry = [ones(6,1) ~temp]; % all the first presses are one
        
        clear meanbnd allmeanbnd CB est_dist act_est_dist
        
        for s = 0:6
            A = getrow(ANA , ANA.seqNumb == s);
            CBD{s+1} = A.estChnkBndry;
        end
        
        for s = 1:7
            for s1 = 1:7
                est_dist(s,s1) = nanmean(nanmean(pdist2(CBD{s} , CBD{s1} , distance)));
            end
            for s2 = 1:6
                act_est_dist(s,s2) = nanmean(nanmean(pdist2(CBD{s} , chbndry(s2 , :) , distance)));
            end
        end
        
        figure('color' , [1 1 1])
        subplot(1,3,1)
        
        imagesc(est_dist);
        title('Estimated Chunking Structure Dissimilarity Matrix' , 'FontSize' , 20)
        hold on
        ax = gca;
        ax.XTick = [1:7];
        ax.YTick = [1:7];
        ax.XTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
        ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
        ax.XTickLabelRotation = 45;
        ax.FontSize = 20;
        hold on
        line([1.5 7.5] , [1.5 1.5] ,  'LineWidth' , 3 , 'color' , [0 0 0])
        line([1.5 1.5] ,[1.5 7.5] ,  'LineWidth' , 3 , 'color' , [0 0 0])
        axis square
        colorbar
        
        subplot(1,3,2)
        
        imagesc(act_est_dist);
        title('Estimated vs. Expected Chunking Structure Dissimilarity Matrix'  , 'FontSize' , 20)
        hold on
        ax = gca;
        line([.5 7.5] , [1.5 1.5] ,  'LineWidth' , 3 , 'color' , [0 0 0])
        ax.XTick = [1:6];
        ax.YTick = [1:7];
        ax.XTickLabel = {'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
        ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
        xlabel('Expected')
        ylabel('Estimated')
        ax.XTickLabelRotation = 45;
        ax.FontSize = 20;
        axis square
        colorbar
        
        subplot(1,3,3)
        imagesc(squareform(pdist(chbndry)));
        title('Chunking Structure Dissimilarity Matrix'  , 'FontSize' , 20)
        hold on
        ax = gca;
        ax.XTick = [1:6];
        ax.YTick = [1:6];
        ax.XTickLabel = {'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
        ax.YTickLabel = {'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
        ax.XTickLabelRotation = 45;
        ax.FontSize = 20;
        axis square
        colorbar
        out = [];
    case 'Mychunk_dist'
        %% chunk distances
        ANA = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood);
        ANA = getrow(ANA ,  ismember(ANA.Rep , rep) & ismember(ANA.Day , days{day}));
        for tn = 1:length(ANA.TN)
            n = (ANA.AllPressIdx(tn , sum(~isnan(ANA.AllPressIdx(tn , :))))  - ANA.AllPressIdx(tn , 1)) / 1000;
            nIdx(tn , :) = (ANA.AllPressIdx(tn , :) - ANA.AllPressIdx(tn , 1))/n;
            IPI(tn , :) = detrend(diff(nIdx(tn ,:) , 1 , 2) , 'linear' , 2);
        end
        
        for tn = 1:length(ANA.TN)
            thresh = .3 * std(IPI(tn , :));
            [dum , estChnkBndry] = findpeaks(IPI(tn , :));% ,'MinPeakProminence', thresh); % slope of IPIs
            
            if ~isempty(estChnkBndry)
                goodpeak = ones(1, length(estChnkBndry));
                for cb = 1:length(estChnkBndry)
                    if IPI(tn , estChnkBndry(cb)) < nanmean(IPI(tn  ,:)) + thresh
                        goodpeak(cb) = 0;
                    end
                end
                if sum(goodpeak)
                    estChnkBndry = estChnkBndry(logical(goodpeak));
                else
                    estChnkBndry = [];
                end
            end
            
            
            ANA.estChnkBndry(tn , :) = zeros(1, 14);  % first presses of chunks will be 1
            ANA.estChnkBndry(tn , estChnkBndry+1) = 1;
            ANA.estChnkBndry(tn , 1) = 1;
        end
        
        if GroupCode == 1
            ChnkArrang = CMB.DesiredCnhkargmnt(CMB.Group1 , :);
            for j = 1:size(ChnkArrang , 1)
                temp = [];
                temp1 = [];
                for k = 1:length(find(ChnkArrang(j , :)))
                    temp = [temp k*ones(1,ChnkArrang(j,k))];
                    temp1 = [temp1 1:ChnkArrang(j,k)];
                end
                ChnkArrng(j , :) = temp;
                ChnkPlcmnt(j,:)  = temp1;
            end
            IPIarrangement = diff(ChnkArrng , 1 , 2);   % between IPIs will have the index of 1
            IPIarrangement(~IPIarrangement) = 2;         % within IPIs will have the index of 2
        elseif GroupCode == 2
            ChnkArrang = CMB.DesiredCnhkargmnt(CMB.Group2 , :);
            for j = 1:size(ChnkArrang , 1)
                temp = [];
                temp1 = [];
                for k = 1:length(find(ChnkArrang(j , :)))
                    temp = [temp k*ones(1,ChnkArrang(j,k))];
                    temp1 = [temp1 1:ChnkArrang(j,k)];
                end
                ChnkArrng(j , :) = temp;
                ChnkPlcmnt(j,:)  = temp1;
            end
            IPIarrangement = diff(ChnkArrng , 1 , 2); % between IPIs will have the index of 1
            IPIarrangement(~IPIarrangement) = 2;       % within IPIs will have the index of 2
        end
        
        
        for d = 2:4
            for s = 0:6
                A = getrow(ANA , ANA.seqNumb == s & ~ANA.isError & ismember(ANA.Day , d) & ismember(ANA.Rep , rep) & ANA.isgood);
                CPD{d-1 , s+1} = A.estChnkPlcmnt;
            end
        end
        for d= 1:3
            for s = 1:7
                for s1 = 1:7
                    cp_est_dist(d,s,s1) = nanmean(nanmean(pdist2(CPD{d,s} , CPD{d,s1} , distance)));
                end
                for s2 = 1:6
                    act_cp_est_dist(d,s,s2) = nanmean(nanmean(pdist2(CPD{d,s} , ChnkPlcmnt(s2 , :) , distance)));
                end
            end
        end
        
        figure('color' , [1 1 1])
        subplot(1,2,1)
        
        imagesc(squeeze(mean(cp_est_dist , 1)));
        title('Estimated Chunking Placement Dissimilarity Matrix' , 'FontSize' , 20)
        hold on
        ax = gca;
        ax.XTick = [1:7];
        ax.YTick = [1:7];
        ax.XTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
        ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
        ax.XTickLabelRotation = 45;
        ax.FontSize = 20;
        hold on
        line([1.5 7.5] , [1.5 1.5] ,  'LineWidth' , 3 , 'color' , [0 0 0])
        line([1.5 1.5] ,[1.5 7.5] ,  'LineWidth' , 3 , 'color' , [0 0 0])
        axis square
        colorbar
        
        subplot(1,2,2)
        imagesc(squareform(pdist(ChnkPlcmnt)));
        title('Chunking Placement Dissimilarity Matrix'  , 'FontSize' , 20)
        hold on
        ax = gca;
        ax.XTick = [1:6];
        ax.YTick = [1:6];
        ax.XTickLabel = {'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'};
        ax.YTickLabel = {'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'};
        ax.XTickLabelRotation = 45;
        ax.FontSize = 20;
        axis square
        colorbar
        out = [];
        
        
        figure('color' , [1 1 1])
        subplot(1,2,1)
        
        imagesc(squeeze(mean(act_cp_est_dist , 1)));
        title('Estimated vs. Expected Chunking Placement Dissimilarity Matrix'  , 'FontSize' , 20)
        hold on
        ax = gca;
        line([.5 7.5] , [1.5 1.5] ,  'LineWidth' , 3 , 'color' , [0 0 0])
        ax.XTick = [1:6];
        ax.YTick = [1:7];
        ax.XTickLabel = {'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
        ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
        xlabel('Expected')
        ylabel('Estimated')
        ax.XTickLabelRotation = 45;
        ax.FontSize = 20;
        axis square
        colorbar
        
        subplot(1,2,2)
        imagesc(squareform(pdist(ChnkPlcmnt)));
        title('Chunking Placement Dissimilarity Matrix'  , 'FontSize' , 20)
        hold on
        ax = gca;
        ax.XTick = [1:6];
        ax.YTick = [1:6];
        ax.XTickLabel = {'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
        ax.YTickLabel = {'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
        ax.XTickLabelRotation = 45;
        ax.FontSize = 20;
        axis square
        colorbar
        
        
        out = [];
    case 'Avg_pattern_sh3'
        D = load('/Users/nedakordjazi/Documents/SeqEye/SequenceHierarchical/Analysis/sh3_avrgPattern.mat');
        MT = D.MT(126:end);
        sq = D.Sequence(126:end , 1:2);
        for i = 1:length(sq)
            AvgMT(sq(i,1),sq(i,2)) = D.MT(i);
        end
        figure('color' , [1 1 1])
        imagesc(AvgMT)
        axis square
        hold on
        ax = gca;
        ax.XTick = [1:5];
        ax.YTick = [1:5];
        colorbar
        title('Double Average Patterns');
        xlabel('First press')
        ylabel('Second press')
        min(MT)/max(MT);
        
        out = [];
    case 'eye_vel_chunkplace'
        % eye velocity as a function of chunk placement
        ANA = getrow(Dall , ismember(Dall.SN , subjnum) & ~Dall.isError & Dall.isgood & ismember(Dall.Rep ,rep) & ismember(Dall.Day , days{day}));
        figure('color' , [1 1 1])
        subplot(1,2,1)
        col = {'b' , 'r' ,'c' ,'g'};
        for cp = 1:4
            A{cp} = ANA.xEyeVelCnkPlcmnt(~isnan(ANA.xEyeVelCnkPlcmnt(:,cp)),cp);
            %histogram(A{cp} ,'EdgeColor' , col{cp} , 'FaceColor' , col{cp} , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
            plot(cp*ones(size(A{cp})) ,A{cp} , '*' , 'color' , col{cp})
            hold on
        end
        legend({'First chunk place' , 'Second chunk place' , 'Third chunk place' , 'Forth chunk place'})
        %         for cp = 1:4
        %             line([nanmean(A{cp}) nanmean(A{cp})] , [1 250] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , col{cp});
        %         end
        xlabel('degree/s')
        grid on
        
        
        All = NaN*ones(max([length(A{1}) , length(A{2}) , length(A{3}) , length(A{4})]) , 4);
        for cp  = 1:4
            All(1:length(A{cp}) , cp) = A{cp};
        end
        
        subplot(1,2,2)
        hold on
        boxplot(All)
        title('Eye velocity per chunk placement')
        xlabel('Chunk place')
        
        
        
        figure('color' , [1 1 1])
        subplot(1,2,1)
        col = {'b' , 'r' ,'c' ,'g'};
        for cp = 1:4
            A{cp} = ANA.xEyeVelEstCnkPlcmnt(~isnan(ANA.xEyeVelEstCnkPlcmnt(:,cp)),cp);
            histogram(A{cp} ,'EdgeColor' , col{cp} , 'FaceColor' , col{cp} , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
            hold on
        end
        legend({'First chunk place' , 'Second chunk place' , 'Third chunk place' , 'Forth chunk place'})
        for cp = 1:4
            line([nanmean(A{cp}) nanmean(A{cp})] , [1 250] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , col{cp});
        end
        xlabel('degree/s')
        grid on
        
        
        All = NaN*ones(max([length(A{1}) , length(A{2}) , length(A{3}) , length(A{4})]) , 4);
        for cp  = 1:4
            All(1:length(A{cp}) , cp) = A{cp};
        end
        
        subplot(1,2,2)
        hold on
        boxplot(All)
        title('Eye Velocity Per Estimated Chunk Placement')
        xlabel('Chunk place')
        out =[];
    case 'eye_vel_seqplace'
        %  velocity per sequence quarter --------CLA
        ANA = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood);
        clear A
        
        id  = ismember(ANA.seqNumb , [1:6]) & ~ANA.isError;
        A{1} = nanmean(ANA.xEyePressAngVelocity(id,1:4) , 2);
        A{2} = nanmean(ANA.xEyePressAngVelocity(id,5:7) , 2);
        A{3} = nanmean(ANA.xEyePressAngVelocity(id,8:11) , 2);
        A{4} = nanmean(ANA.xEyePressAngVelocity(id,12:14) , 2);
        
        out.seqQuarter_vel = A;
        
        figure('color' , [1 1 1])
        subplot(1,2,1)
        
        col = {'b' , 'r' ,'c' ,'g'};
        for cp = 1:4
            histogram(A{cp} ,'EdgeColor' , col{cp} , 'FaceColor' , col{cp} , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
            hold on
        end
        
        legend({'First quarter' , 'Second quarter' , 'Third quarter' , 'Forth quarter'})
        for cp = 1:4
            line([nanmean(A{cp}) nanmean(A{cp})] , [0 250] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , col{cp});
        end
        grid on
        xlabel('deg/sec')
        title('Average eye velocity over all trials in 4 quarters in CLA sequences')
        
        
        All = NaN*ones(max([length(A{1}) , length(A{2}) , length(A{3}) , length(A{4})]) , 4);
        for sq  = 1:4
            All(1:length(A{sq}) , sq) = A{sq};
        end
        
        subplot(1,2,2)
        hold on
        boxplot(All)
        title('Eye velocity per sequence quarter -- CLA sequences')
        xlabel('Seq Quarter')
    case 'eye_pos_seqplace'
        %  velocity per sequence quarter --------Random
        ANA = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood);
        clear A
        
        ANAcla = getrow(ANA , ismember(ANA.seqNumb , [1:6]) & ~ANA.isError & ismember(ANA.Rep , rep));
        ANAr   = getrow(ANA , ismember(ANA.seqNumb , [0]) & ~ANA.isError & ismember(ANA.Day , [1:4]) & ismember(ANA.Rep , rep));
        subplot(2,1,1)
        boxplot(ANAcla.EyePressTimePos);
        title('Median eye position in a 50 ms vicinity of the press times - CLA')
        ylabel('cm')
        ylim([1 14])
        grid on
        
        subplot(2,1,2)
        boxplot(ANAr.EyePressTimePos);
        title('Median eye position in a 50 ms vicinity of the press times - Random')
        ylabel('cm')
        grid on
        ylim([1 14])
    case 'eyepress_pos_traces'
        clear len EyePatterne pressPattern eye press chunk chunkPattern estchunkPattern
        if isequal(subjnum , [1:length(subj_name)-1])
            subjnum = 13;
        end
        %         calc = 0;
        if calc
            N = se1_timeNorm('single_sub' ,subjnum,day);
            out.eye = N.norm.eye;
            out.press= N.norm.press;
            out.eyeveloc= N.norm.eyeveloc;
            out.prsveloc= N.norm.prsveloc;
            for s = 0:6
                out.EyePattern(s+1, :)  = nanmean(out.eye{s+1});
                out.pressPattern(s+1,:) = nanmean(out.press{s+1});
                out.eyevelocPattern(s+1,:) = nanmean(out.eyeveloc{s+1});
                out.prsvelocPattern(s+1,:) = nanmean(out.prsveloc{s+1});
            end
        else
            switch rep
                case 1
                    load([baseDir , '/se1_tnorm_r1.mat'])
                case 2
                    load([baseDir , '/se1_tnorm_r2.mat'])
                otherwise
                    load([baseDir , '/se1_tnorm.mat'])
            end
            
            
            
            out.eye = N.norm(day,subjnum).eye;
            out.press= N.norm(day,subjnum).press;
            out.eyeveloc= N.norm(day,subjnum).eyeveloc;
            out.prsveloc= N.norm(day,subjnum).prsveloc;
            
            out.EyePattern  = N.norm(day,subjnum).eyePattern;
            out.pressPattern = N.norm(day,subjnum).pressPattern;
            out.eyevelocPattern = N.norm(day,subjnum).eyevelocPattern;
            out.prsvelocPattern = N.norm(day,subjnum).prsvelocPattern;
            
        end
        
        figure('color' , [1 1 1])
        
        for s = 2:7
            plot(out.pressPattern(s,:) , out.EyePattern(s,:) , 'LineWidth' , 3);
            hold on
        end
        plot(out.pressPattern(1,:) , out.EyePattern(1,:)  , 'LineWidth' , 5 , 'color' , 'black');
        grid on
        xlabel('Finger Press' , 'FontSize' , 20)
        ylabel('Eye Position' , 'FontSize' , 20)
        title(['The Average Time Normalized Eye Position vs. Finger Press Time series - Subject ' , num2str(subjnum)] , 'FontSize' , 20)
        legend({ 'Structure 1','Structure 2','Structure 3','Structure 4','Structure 5','Structure 6' 'Random'} , 'FontSize' , 20)
        hold on
        ax = gca;
        ax.YTick = [1:14];
        ax.YLim = [1 14];
        ax.XTick = [1:14];
        ax.XLim = [1 14];
        ax.FontSize = 20;
        ax.Box = 'off';
        axis square
        
        
        figure('color' , [1 1 1])
        subplot(2,1,2)
        plot(out.pressPattern(2:end , :)' , 'LineWidth' , 3)
        grid on
        hold on
        plot(out.pressPattern(1 , :)' , 'LineWidth' , 5 , 'color' , 'black')
        xlabel('Normalized time' , 'FontSize' , 20)
        ylabel('Press Position' , 'FontSize' , 20)
        title(['The Average Time Normalized Finger Press Time series - Subject ' , num2str(subjnum)] , 'FontSize' , 20)
%         legend({ 'Structure 1','Structure 2','Structure 3','Structure 4','Structure 5','Structure 6' 'Random'} , 'FontSize' , 20)
        
        ax = gca;
        ax.YTick = [1:14];
        ax.YLim = [1 14];
        ax.FontSize = 20;
        ax.Box = 'off';
        
        
        subplot(2,1,1)
        plot(out.EyePattern(2:end , :)' , 'LineWidth' , 3)
        hold on
        plot(out.EyePattern(1 , :)' , 'LineWidth' , 5 , 'color' , 'black')
        grid on
        xlabel('Normalized time' , 'FontSize' , 20)
        ylabel('Eye Position' , 'FontSize' , 20)
        title(['The Average Time Normalized Eye Position Time series - Subject ' , num2str(subjnum)] , 'FontSize' , 20)
%         legend({ 'Structure 1','Structure 2','Structure 3','Structure 4','Structure 5','Structure 6' 'Random'} , 'FontSize' , 20)
        
        ax = gca;
        ax.YTick = [1:14];
        ax.YLim = [1 max(max(out.EyePattern))];
        ax.FontSize = 20;
        ax.Box = 'off';
        
        figure('color' , [1 1 1])
        subplot(2,1,1)
        plot(out.EyePattern(2:end , 901:end)' , 'LineWidth' , 3)
        hold on
        plot(out.EyePattern(1 , 901:end)' , 'LineWidth' , 5 , 'color' , 'black')
        grid on
        xlabel('Normalized time' , 'FontSize' , 20)
        ylabel('Eye Position' , 'FontSize' , 20)
        title(['The Average Last 10% Eye Position Time series - Subject ' , num2str(subjnum)] , 'FontSize' , 20)
        ax = gca;
        ax.FontSize = 20;
        ax.Box = 'off';
        ax.XLim = [1 100];
        ax.YLim = [min(min(out.EyePattern(: , 901:end))) max(max(out.EyePattern(: , 901:end)))];
        
        subplot(2,1,2)
        plot(out.EyePattern(2:end , 1:100)' , 'LineWidth' , 3)
        hold on
        plot(out.EyePattern(1 , 1:100)' , 'LineWidth' , 5 , 'color' , 'black')
        grid on
        xlabel('Normalized time' , 'FontSize' , 20)
        ylabel('Eye Position' , 'FontSize' , 20)
        title(['The Average First 10% Eye Position Time series - Subject ' , num2str(subjnum)] , 'FontSize' , 20)
        legend({ 'Structure 1','Structure 2','Structure 3','Structure 4','Structure 5','Structure 6' 'Random'} , 'FontSize' , 20)        
        ax = gca;
        ax.XLim = [1 100];
        ax.YLim = [min(min(out.EyePattern(: , 1:100))) max(max(out.EyePattern(: , 1:100)))];
        ax.FontSize = 20;
        ax.Box = 'off';
        
        out = [];
    case 'eyepress_pos_distances'
        clear len EyePatterne pressPattern eye press chunk chunkPattern estchunkPattern
        %         calc = 0;
        if isequal(subjnum , [1 3 4 5 6 7 8 9 10] , 'rows')
            subjnum = 11;
        end
        
        titleAdd = {'First Third' , 'Middle Third' , 'Last Third' , 'The Entire Time Series'};
        if calc
            N = se1_timeNorm('single_sub' ,subjnum,day);
            out.eye = N.norm.eye;
            out.press= N.norm.press;
            out.eyeveloc= N.norm.eyeveloc;
            out.prsveloc= N.norm.prsveloc;
            for s = 0:6
                out.EyePattern(s+1, :)  = nanmean(out.eye{s+1});
                out.pressPattern(s+1,:) = nanmean(out.press{s+1});
                out.eyevelocPattern(s+1,:) = nanmean(out.eyeveloc{s+1});
                out.prsvelocPattern(s+1,:) = nanmean(out.prsveloc{s+1});
            end
        else
            switch rep
                case 1
                    load([baseDir , '/se1_tnorm_r1.mat'])
                case 2
                    load([baseDir , '/se1_tnorm_r2.mat'])
                otherwise
                    load([baseDir , '/se1_tnorm.mat'])
            end
            % Use the random sequences in day 1 as the baseline. On day one we only have random sequences and  no CLAs
            
            
            
            out.eye = N.norm(day,subjnum).eye;
            out.press= N.norm(day,subjnum).press;
            out.eyeveloc= N.norm(day,subjnum).eyeveloc;
            out.prsveloc= N.norm(day,subjnum).prsveloc;
            out.eyePattern  = N.norm(day,subjnum).eyePattern;
            out.pressPattern = N.norm(day,subjnum).pressPattern;
            out.eyevelocPattern = N.norm(day,subjnum).eyevelocPattern;
            out.prsvelocPattern = N.norm(day,subjnum).prsvelocPattern;
            
        end
        
        
        figure('color' , [1 1 1])
        for t = 1:length(tid)
            for s = 1:length(out.eye)
                for s1 = 1:length(out.eye)
                    out.prsDist{t}(s,s1) = nanmean(nanmean(pdist2(out.press{s}(:,tid{t}) , out.press{s1}(:,tid{t}), distance)));
                end
            end
            subplot(1,length(tid),t)
            imagesc(out.prsDist{t})
            colorbar
            hold on
            ax = gca;
            ax.XTick = [1:7];
            ax.YTick = [1:7]
            ax.XTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6' };
            ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6' };
            ax.XTickLabelRotation = 45;
            axis square
            title(['Time-Normalized Press Position Time-Series Distances -  ' , titleAdd{t} , ' - Subject ' , num2str(subjnum)])
        end
        
        figure('color' , [1 1 1])
        for t = 1:length(tid)
            for s = 1:length(out.eye)
                for s1 = 1:length(out.eye)
                    out.eyeDist{t}(s,s1) = nanmean(nanmean(pdist2(out.eye{s}(:,tid{t}) , out.eye{s1}(:,tid{t}), distance)));
                end
            end
            subplot(1,length(tid),t)
            imagesc(out.eyeDist{t})
            colorbar
            hold on
            ax = gca;
            ax.XTick = [1:7];
            ax.YTick = [1:7]
            ax.XTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6' };
            ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6' };
            ax.XTickLabelRotation = 45;
            axis square
            title(['Time-Normalized Eye Position Time-Series Distances -  ' , titleAdd{t} , ' - Subject ' , num2str(subjnum)])
        end
    case 'eyepress_pos_avg_distances'
        
        clear len eyePatterne pressPattern eye press chunk chunkPattern estchunkPattern
        if isequal(subjnum , [1 3 4 5 6 7 8 9 10] , 'rows')
            subjnum = 11;
        end
        if calc
            N = se1_timeNorm('single_sub' ,subjnum,day);
            out.eye = N.norm.eye;
            out.press= N.norm.press;
            out.eyeveloc= N.norm.eyeveloc;
            out.prsveloc= N.norm.prsveloc;
            for s = 0:6
                out.eyePattern(s+1, :)  = nanmean(out.eye{s+1});
                out.pressPattern(s+1,:) = nanmean(out.press{s+1});
                out.eyevelocPattern(s+1,:) = nanmean(out.eyeveloc{s+1});
                out.prsvelocPattern(s+1,:) = nanmean(out.prsveloc{s+1});
            end
        else
            switch rep
                case 1
                    load([baseDir , '/se1_tnorm_r1.mat'])
                case 2
                    load([baseDir , '/se1_tnorm_r2.mat'])
                otherwise
                    load([baseDir , '/se1_tnorm.mat'])
            end
            out.eye = N.norm(day,subjnum).eye;
            out.press= N.norm(day,subjnum).press;
            out.eyeveloc= N.norm(day,subjnum).eyeveloc;
            out.prsveloc= N.norm(day,subjnum).prsveloc;
            out.eyePattern  = N.norm(day,subjnum).eyePattern;
            out.pressPattern = N.norm(day,subjnum).pressPattern;
            out.eyevelocPattern = N.norm(day,subjnum).eyevelocPattern;
            out.prsvelocPattern = N.norm(day,subjnum).prsvelocPattern;
        end
        
        titleAdd = {'First Third' , 'Middle Third' , 'Last Third' , 'The Entire Time Series'};
        
        figure('color' , [1 1 1])
        for t = 1:length(tid)-1
            out.prseyeDist{t} = pdist2(out.pressPattern(:,tid{t}) , out.eyePattern(:,tid{t}) , distance);
            subplot(1,length(tid),t)
            imagesc(out.prseyeDist{t})
            colorbar
            hold on
            ax = gca;
            ax.XTick = [1:7];
            ax.YTick = [1:7];
            ax.XTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6' }
            ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6' }
            xlabel('Eye')
            ylabel('Press')
            ax.XTickLabelRotation = 45;
            axis square
            ax.FontSize = 20;
            title(['Finger Press/Eye Position Dissimilarity Matrix - ' ,titleAdd{t}] , 'FontSize' , 20)
        end
        
        figure('color' , [1 1 1])
        for t = length(tid)
            out.prseyeDist{t} = pdist2(out.pressPattern(:,tid{t}) , out.eyePattern(:,tid{t}) , distance);
            subplot(1,length(tid),t)
            imagesc(out.prseyeDist{t})
            colorbar
            hold on
            ax = gca;
            ax.XTick = [1:7];
            ax.YTick = [1:7];
            ax.XTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6' }
            ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6' }
            xlabel('Eye')
            ylabel('Press')
            ax.XTickLabelRotation = 45;
            axis square
            ax.FontSize = 20;
            title('Finger Press/Eye Position Dissimilarity Matrix' , 'FontSize' , 20)
        end
    case 'eyepress_vel_traces'
        clear len eyePatterne pressPattern eye press chunk chunkPattern estchunkPattern
        if isequal(subjnum , [1 3 4 5 6 7 8 9 10])
            subjnum = 11;
        end
        %         calc = 0;
        if calc
            N = se1_timeNorm('single_sub' ,subjnum,day);
            out.eye = N.norm.eye;
            out.press= N.norm.press;
            out.eyeveloc= N.norm.eyeveloc;
            out.prsveloc= N.norm.prsveloc;
            out.eyePattern  = N.norm.eyePattern;
            out.pressPattern = N.norm.pressPattern;
            out.eyevelocPattern = N.norm.eyevelocPattern;
            out.prsvelocPattern = N.norm.prsvelocPattern;
            
        else
            switch rep
                case 1
                    load([baseDir , '/se1_tnorm_r1.mat'])
                case 2
                    load([baseDir , '/se1_tnorm_r2.mat'])
                otherwise
                    load([baseDir , '/se1_tnorm.mat'])
            end
            
            out.eye = N.norm(day,subjnum).eye;
            out.press= N.norm(day,subjnum).press;
            out.eyeveloc= N.norm(day,subjnum).eyeveloc;
            out.prsveloc= N.norm(day,subjnum).prsveloc;
            out.eyePattern  = N.norm(day,subjnum).eyePattern;
            out.pressPattern = N.norm(day,subjnum).pressPattern;
            out.eyevelocPattern = N.norm(day,subjnum).eyevelocPattern;
            out.prsvelocPattern = N.norm(day,subjnum).prsvelocPattern;
        end
        
        figure('color' , [1 1 1])
        
        hold on
        for s = 2:7
            plot(smooth(out.prsvelocPattern(s,:) , 10) ,smooth(out.eyevelocPattern(s,:) ,10), 'LineWidth' , 1);
            hold on
        end
        plot(smooth(out.prsvelocPattern(1,:) ,10),smooth(out.eyevelocPattern(1,:) ,10), 'LineWidth' , 3);
        legend({'Structure 1','Structure 2','Structure 3','Structure 4','Structure 5','Structure 6' 'Random'})
        axis square
        xlabel('Press/sec')
        ylabel('deg/sec')
        title(['Press Velocity vs. Eye Position Time Series, Subject ' ,num2str(subjnum)])
        ax = gca;
        ax.FontSize = 20;
        grid on
        
        figure('color' , [1 1 1])
        subplot(2,1,1)
        for m  = 2:7
            plot(smooth(out.prsvelocPattern(m,:) ,10) , 'LineWidth' , 1);
            hold on
        end
        plot(smooth(out.prsvelocPattern(1,:)' ,10)   , 'LineWidth' , 3);
        xlabel('Normalized Time')
        ylabel('Press/sec')
        title(['Finger Press Velocity, Subject  ' , num2str(subjnum)])
        grid on
        ax = gca;
        ax.FontSize = 20;
        legend({'Structure 1','Structure 2','Structure 3','Structure 4','Structure 5','Structure 6' 'Random' })
        
        
        subplot(2,1,2)
        for m = 2:7
            plot(smooth(out.eyevelocPattern(m , :) , 10), 'LineWidth' , 1)
            hold on
        end
        plot(smooth(out.eyevelocPattern(1, :)' , 10), 'LineWidth' , 3)
        xlabel('Normalized Time')
        ylabel('deg/sec')
        title(['Eye Angular Velocity, Subject ' , num2str(subjnum)])
        grid on
        ax = gca;
        ax.FontSize = 20;
        legend({ 'Structure 1','Structure 2','Structure 3','Structure 4','Structure 5','Structure 6' 'Random'})
    case 'eyepress_vel_distances'
        if isequal(subjnum , [1 3 4 5 6 7 8 9 10])
            subjnum = 11;
        end
        
        clear len eyePatterne pressPattern eye press chunk chunkPattern estchunkPattern
        %         calc = 0;
        if calc
            N = se1_timeNorm('single_sub' ,subjnum,day);
            out.eye = N.norm.eye;
            out.press= N.norm.press;
            out.eyeveloc= N.norm.eyeveloc;
            out.prsveloc= N.norm.prsveloc;
            for s = 0:6
                out.eyePattern(s+1, :)  = nanmean(out.eye{s+1});
                out.pressPattern(s+1,:) = nanmean(out.press{s+1});
                out.eyevelocPattern(s+1,:) = nanmean(out.eyeveloc{s+1});
                out.prsvelocPattern(s+1,:) = nanmean(out.prsveloc{s+1});
            end
        else
            switch rep
                case 1
                    load([baseDir , '/se1_tnorm_r1.mat'])
                case 2
                    load([baseDir , '/se1_tnorm_r2.mat'])
                otherwise
                    load([baseDir , '/se1_tnorm.mat'])
            end
            out.eye = N.norm(day,subjnum).eye;
            out.press= N.norm(day,subjnum).press;
            out.eyeveloc= N.norm(day,subjnum).eyeveloc;
            out.prsveloc= N.norm(day,subjnum).prsveloc;
            out.eyePattern  = N.norm(day,subjnum).eyePattern;
            out.pressPattern = N.norm(day,subjnum).pressPattern;
            out.eyevelocPattern = N.norm(day,subjnum).eyevelocPattern;
            out.prsvelocPattern = N.norm(day,subjnum).prsvelocPattern;
        end
        
        eyeid = [1 2 3 4];
        prsid = [5 6 7 8];
        titleAdd = {'First Quarter' , 'Second Quarter' , 'Third Quarter' , 'Fourth Quarter' , 'The Entire Time Series'};
        figure('color' , [1 1 1])
        for t = 1:length(tid)-1
            for s = 1:length(out.eye)
                for s1 = 1:length(out.eye)
                    out.prsDist{t}(s,s1) = nanmedian(nanmedian(pdist2(out.prsveloc{s}(:,tid{t}),out.prsveloc{s1}(:,tid{t}), distance)));
                end
            end
            subplot(1,length(tid),t)
            imagesc(out.prsDist{t})
            colorbar
            hold on
            ax = gca;
            ax.XTick = [1:7];
            ax.YTick = [1:7];
            ax.XTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6' }
            ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6' }
            xlabel('Eye')
            ylabel('Press')
            ax.XTickLabelRotation = 45;
            axis square
            title(['Dissimilarity Between Average Time-Normalized Eye Angular Velocities - ' , titleAdd{t}])
        end
        
        figure('color' , [1 1 1])
        for t = 1:length(tid)-1
            for s = 1:length(out.eye)
                for s1 = 1:length(out.eye)
                    out.eyeDist{t}(s,s1) = nanmedian(nanmedian(pdist2(out.eyeveloc{s}(:,tid{t}),out.eyeveloc{s1}(:,tid{t}), distance)));
                end
            end
            subplot(1,length(tid),t)
            imagesc(out.eyeDist{t})
            colorbar
            hold on
            ax = gca;
            ax.XTick = [1:7];
            ax.YTick = [1:7];
            ax.XTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6' }
            ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6' }
            xlabel('Eye')
            ylabel('Press')
            ax.XTickLabelRotation = 45;
            axis square
            title(['Dissimilarity Between Average Time-Normalized Eye Angular Velocities - ' , titleAdd{t}])
        end
    case 'eyepress_vel_avg_distances'
        %%
        clear len eyePatterne pressPattern eye press chunk chunkPattern estchunkPattern
        if isequal(subjnum , [1 3 4 5 6 7 8 9 10])
            subjnum = 11;
        end
        if calc
            N = se1_timeNorm('single_sub' ,subjnum,day);
            out.eye = N.norm.eye;
            out.press= N.norm.press;
            out.eyeveloc= N.norm.eyeveloc;
            out.prsveloc= N.norm.prsveloc;
            for s = 0:6
                out.eyePattern(s+1, :)  = nanmean(out.eye{s+1});
                out.pressPattern(s+1,:) = nanmean(out.press{s+1});
                out.eyevelocPattern(s+1,:) = nanmean(out.eyeveloc{s+1});
                out.prsvelocPattern(s+1,:) = nanmean(out.prsveloc{s+1});
            end
        else
            switch rep
                case 1
                    load([baseDir , '/se1_tnorm_r1.mat'])
                case 2
                    load([baseDir , '/se1_tnorm_r2.mat'])
                otherwise
                    load([baseDir , '/se1_tnorm.mat'])
            end
            out.eye_bl{1}      = 1 + N.norm(1,subjnum).eye{1};
            out.press_bl{1}    = N.norm(1,subjnum).press{1};
            out.eyeveloc_bl{1} = N.norm(1,subjnum).eyeveloc{1};
            out.prsveloc_bl{1} = N.norm(1,subjnum).prsveloc{1};
            
            
            out.eyePattern_bl      = nanmean(out.eye_bl{1});
            out.pressPattern_bl    = nanmean(out.press_bl{1});
            out.eyevelocPattern_bl = nanmean(out.eyeveloc_bl{1});
            out.prsvelocPattern_bl = nanmean(out.prsveloc_bl{1});
            
            
            out.eye = 1 + N.norm(day,subjnum).eye;
            out.press= N.norm(day,subjnum).press;
            out.eyeveloc= N.norm(day,subjnum).eyeveloc;
            out.prsveloc= N.norm(day,subjnum).prsveloc;
            for s = 0:6
                out.eyePattern(s+1, :)  = nanmean(out.eye{s+1});
                out.pressPattern(s+1,:) = nanmean(out.press{s+1});
                out.eyevelocPattern(s+1,:) = nanmean(out.eyeveloc{s+1});
                out.prsvelocPattern(s+1,:) = nanmean(out.prsveloc{s+1});
            end
        end
        
        
        eyeid = [1 2 3 4];
        prsid = [5 6 7 8];
        
        
        figure('color' , [1 1 1])
        for t = 1:length(tid)-1
            out.prseyeDist{t} = pdist2(out.prsvelocPattern(2:7,tid{t}) , out.eyevelocPattern(2:7,tid{t}) , distance);
            subplot(1,length(tid)-1,t)
            imagesc(out.prseyeDist{t})
            colorbar
            hold on
            ax = gca;
            ax.XTick = [1:7];
            ax.YTick = [1:7];
            ax.XTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6' }
            ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6' }
            xlabel('Eye')
            ylabel('Press')
            ax.XTickLabelRotation = 45;
            axis square
            title(['Dissimilarity between time-normalized press and eye velocity - quarter ' , num2str(t)])
        end
        
        figure('color' , [1 1 1])
        for t = length(tid)
            out.prseyeDist{t} = pdist2(out.prsvelocPattern(2:7,tid{t}) , out.eyevelocPattern(2:7,tid{t}) , distance);
            imagesc(out.prseyeDist{t})
            colorbar
            hold on
            ax = gca;
            ax.XTick = [1:7];
            ax.YTick = [1:7];
            ax.XTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6' }
            ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6' }
            xlabel('Eye')
            ylabel('Press')
            ax.XTickLabelRotation = 45;
            axis square
            title('Dissimilarity between time-normalized press and eye velocity')
        end
    case 'run_f-tests_mt'
        % /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\movement time of Structres vs. Random/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
        %//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\
        %\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//
        ANA1 = getrow(Dall , ismember(Dall.seqNumb,[0:6]) &  ~Dall.isError & Dall.isgood);
        ANA1.tempSeqNum = ANA1.seqNumb;
        ANA1.tempSeqNum(ismember(ANA1.tempSeqNum , [1:6])) = 1;
        ANA1.mt = ANA1.AllPressTimes(:,14) - ANA1.AllPressTimes(:,1);
        
        id1 = ANA1.tempSeqNum == 1;
        id2 = ANA1.tempSeqNum == 0;
        out.CLAvsRand_mt = anovaMixed(ANA1.mt , ANA1.SN,'within',[ANA1.tempSeqNum],{'Sequence type'},'subset',~isnan(ANA1.mt),'intercept',1)  ;
        out.CLAvsRand_rt = anovaMixed(ANA1.AllPressTimes(:,1) , ANA1.SN,'within',[ANA1.tempSeqNum],{'Sequence type'},'subset',~isnan(ANA1.AllPressTimes(:,1)),'intercept',1)  ;
        
        figure('color' , [1 1 1])
        subplot(1,2,1)
        histogram(ANA1.mt(id1) ,'EdgeColor' , 'b' , 'FaceColor' , 'b' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        hold on
        
        histogram(ANA1.mt(id2) ,'EdgeColor' , 'r' , 'FaceColor' , 'r' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        legend({'Structured sequences' , 'Random sequences'})
        
        line([nanmean(ANA1.mt(id1)) nanmean(ANA1.mt(id1))] , [1 300] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'b');
        line([nanmean(ANA1.mt(id2)) nanmean(ANA1.mt(id2))] , [1 300] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'r');
        title(['Distribution of movement time, p = ' , num2str( out.CLAvsRand_mt.eff(2).p)])
        grid on
        
        subplot(1,2,2)
        histogram(ANA1.AllPressTimes(id1,1) ,'EdgeColor' , 'b' , 'FaceColor' , 'b' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        hold on
        id2 = ANA1.tempSeqNum == 0;
        histogram(ANA1.AllPressTimes(id2,1) ,'EdgeColor' , 'r' , 'FaceColor' , 'r' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        legend({'Structured sequences' , 'Random sequences'})
        
        line([nanmean(ANA1.AllPressTimes(id1,1)) nanmean(ANA1.AllPressTimes(id1,1))] , [1 600] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'b');
        line([nanmean(ANA1.AllPressTimes(id2,1)) nanmean(ANA1.AllPressTimes(id2,1))] , [1 600] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'r');
        title(['Distribution of reaction time, p = ' , num2str(out.CLAvsRand_mt.eff(2).p)])
        grid on
        
        
        out.StrucvsStruct = anovaMixed(ANA1.mt , ANA1.SN,'within',[ANA1.seqNumb],{'Sertucture number'},'subset',ismember(ANA1.seqNumb ,[1:6]) & ~isnan(ANA1.mt),'intercept',1);
        color = {'b' , 'r' ,'g' ,'c' ,'y' ,'black'};
        figure
        for s= 1:6
            histogram(ANA1.mt(ANA1.seqNumb == s) ,'EdgeColor' , color{s} , 'FaceColor' , color{s} , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
            hold on
        end
        legend({'Structured 1' , 'Structured 2' ,'Structured 3' , 'Structured 4' , 'Structured 5' , 'Structured 6'})
        for s= 1:6
            line([nanmean(ANA1.mt(ANA1.seqNumb == s)) nanmean(ANA1.mt(ANA1.seqNumb == s))] , [0 150] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' ,  color{s});
            hold on
        end
        grid on
        title(['Distribution of movement time, p = ' , num2str( out.StrucvsStruct.eff(2).p)])
        
        
        
        
        % /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\movement time of repetition 1 vs repetition 2/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
        %//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\
        %\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//
        
        id1 = ANA1.Rep == 1 & ANA1.tempSeqNum == 1;
        id2 = ANA1.Rep == 2 & ANA1.tempSeqNum == 1;
        out.S_R1vsR2_mt = anovaMixed(ANA1.mt , ANA1.SN,'within',[ANA1.Rep],{'Repetition'},'subset',~isnan(ANA1.mt) & ANA1.tempSeqNum == 1,'intercept',1);
        out.R_R1vsR2_mt = anovaMixed(ANA1.mt , ANA1.SN,'within',[ANA1.Rep],{'Repetition'},'subset',~isnan(ANA1.mt) & ANA1.tempSeqNum == 0,'intercept',1);
        
        
        figure('color' , [1 1 1])
        subplot(1,2,1)
        histogram(ANA1.mt(id1) ,'EdgeColor' , 'b' , 'FaceColor' , 'b' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        hold on
        
        histogram(ANA1.mt(id2) ,'EdgeColor' , 'r' , 'FaceColor' , 'r' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        legend({'Repetition 1' , 'Repetition 2'})
        
        line([nanmean(ANA1.mt(id1)) nanmean(ANA1.mt(id1))] , [1 250] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'b');
        line([nanmean(ANA1.mt(id2)) nanmean(ANA1.mt(id2))] , [1 250] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'r');
        title(['Distribution of movement time over across all structured sequences, p = ' , num2str(out.S_R1vsR2_mt.eff(2).p)])
        grid on
        
        
        id1 = ANA1.Rep == 1 & ANA1.tempSeqNum == 0;
        id2 = ANA1.Rep == 2 & ANA1.tempSeqNum == 0;
        subplot(1,2,2)
        histogram(ANA1.mt(id1) ,'EdgeColor' , 'b' , 'FaceColor' , 'b' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        hold on
        id2 = ANA1.tempSeqNum == 0;
        histogram(ANA1.mt(id2) ,'EdgeColor' , 'r' , 'FaceColor' , 'r' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        legend({'Repetition 1' , 'Repetition 2'})
        
        line([nanmean(ANA1.mt(id1)) nanmean(ANA1.mt(id1))] , [1 250] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'b');
        line([nanmean(ANA1.mt(id2)) nanmean(ANA1.mt(id2))] , [1 250] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'r');
        title(['Distribution of reaction time over across all random sequences, p = ' , num2str(out.R_R1vsR2_mt.eff(2).p)])
        grid on
        
        % /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\only structure blocks vs intermixed/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
        %//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\
        %\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//
        
        % intermixed --> BT = 2
        % CLAT --> BT = 3
        id1 = ANA1.BT == 2;
        id2 = ANA1.BT == 3;
        out.SvsI_mt = anovaMixed(ANA1.mt , ANA1.SN,'within',[ANA1.BT],{'Block Type'},'subset',~isnan(ANA1.mt) & ANA1.tempSeqNum == 1,'intercept',1);
        out.SvsI_rt = anovaMixed(ANA1.AllPressTimes(:,1) , ANA1.SN,'within',[ANA1.BT],{'Block Type'},'subset',~isnan(ANA1.AllPressTimes(:,1)) & ANA1.tempSeqNum == 1,'intercept',1);
        
        figure('color' , [1 1 1])
        subplot(1,2,1)
        histogram(ANA1.mt(id1) ,'EdgeColor' , 'b' , 'FaceColor' , 'b' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        hold on
        
        histogram(ANA1.mt(id2) ,'EdgeColor' , 'r' , 'FaceColor' , 'r' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        legend({'Intermixed block' , 'Structured block'})
        
        line([nanmean(ANA1.mt(id1)) nanmean(ANA1.mt(id1))] , [1 250] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'b');
        line([nanmean(ANA1.mt(id2)) nanmean(ANA1.mt(id2))] , [1 250] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'r');
        title(['Distribution of movement time over across all structured trails, p = ' , num2str(out.SvsI_mt.eff(2).p)])
        grid on
        
        subplot(1,2,2)
        histogram(ANA1.AllPressTimes(id1,1) ,'EdgeColor' , 'b' , 'FaceColor' , 'b' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        hold on
        id2 = ANA1.tempSeqNum == 0;
        histogram(ANA1.AllPressTimes(id2,1) ,'EdgeColor' , 'r' , 'FaceColor' , 'r' , 'FaceAlpha' , .5 , 'EdgeAlpha' , .5);
        legend({'Intermixed block' , 'Structured block'})
        
        line([nanmean(ANA1.AllPressTimes(id1,1)) nanmean(ANA1.AllPressTimes(id1,1))] , [1 500] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'b');
        line([nanmean(ANA1.AllPressTimes(id2,1)) nanmean(ANA1.AllPressTimes(id2,1))] , [1 500] , 'LineWidth' , 3 ,'LineStyle' , ':' , 'color' , 'r');
        title(['Distribution of reaction time over across all structured trails, p = ' , num2str(out.SvsI_rt.eff(2).p)])
        grid on
    case 'dtw'
        
        %% Dynamic time warping
        ANA = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood);
        ANA1 = getrow(ANA , ismember(ANA.seqNumb , [0:6]) & ~ANA.isError);
        if calc
            N = se1_timeNorm('single_sub' ,subjnum,day);
            out.eye = N.norm.eye;
            out.press= N.norm.press;
            out.eyeveloc= N.norm.eyeveloc;
            out.prsveloc= N.norm.prsveloc;
            for s = 0:6
                out.eyePattern(s+1, :)  = nanmean(out.eye{s+1});
                out.pressPattern(s+1,:) = nanmean(out.press{s+1});
                out.eyevelocPattern(s+1,:) = nanmean(out.eyeveloc{s+1});
                out.prsvelocPattern(s+1,:) = nanmean(out.prsveloc{s+1});
            end
        else
            switch rep
                case 1
                    load([baseDir , '/se1_tnorm_r1.mat'])
                case 2
                    load([baseDir , '/se1_tnorm_r2.mat'])
                otherwise
                    load([baseDir , '/se1_tnorm.mat'])
            end
            out.eye_bl{1}      = 1 + N.norm(1,subjnum).eye{1};
            out.press_bl{1}    = N.norm(1,subjnum).press{1};
            out.eyeveloc_bl{1} = N.norm(1,subjnum).eyeveloc{1};
            out.prsveloc_bl{1} = N.norm(1,subjnum).prsveloc{1};
            
            
            out.eyePattern_bl      = nanmean(out.eye_bl{1});
            out.pressPattern_bl    = nanmean(out.press_bl{1});
            out.eyevelocPattern_bl = nanmean(out.eyeveloc_bl{1});
            out.prsvelocPattern_bl = nanmean(out.prsveloc_bl{1});
            
            
            out.eye = 1 + N.norm(day,subjnum).eye;
            out.press= N.norm(day,subjnum).press;
            out.eyeveloc= N.norm(day,subjnum).eyeveloc;
            out.prsveloc= N.norm(day,subjnum).prsveloc;
            for s = 0:6
                out.eyePattern(s+1, :)  = nanmean(out.eye{s+1});
                out.pressPattern(s+1,:) = nanmean(out.press{s+1});
                out.eyevelocPattern(s+1,:) = nanmean(out.eyeveloc{s+1});
                out.prsvelocPattern(s+1,:) = nanmean(out.prsveloc{s+1});
            end
        end
        
        for tn = 1:length(ANA1.TN)
            if ANA1.AllPressIdx(tn,14) < length(ANA1.xEyePosDigit{tn})
                ANA1.xEyePosDigit{tn}          = ANA1.xEyePosDigit{tn}(ANA1.AllPressIdx(tn,1) : ANA1.AllPressIdx(tn,14));
                ANA1.PressTimeSeries{tn}    = ANA1.PressTimeSeries{tn}(ANA1.AllPressIdx(tn,1) : ANA1.AllPressIdx(tn,14));
                %                 ANA1.ChunkTimeSeries{tn}    = ANA1.ChunkTimeSeries{tn}(ANA1.AllPressIdx(tn,1) : ANA1.AllPressIdx(tn,14));
                %                 ANA1.estChunkTimeSeries{tn} = ANA1.estChunkTimeSeries{tn}(ANA1.AllPressIdx(tn,1) : ANA1.AllPressIdx(tn,14));
            end
        end
        
        dxeye  = zeros(length(ANA1.TN) , length(ANA1.TN));
        dpress = zeros(length(ANA1.TN) , length(ANA1.TN));
        figure('color' , [1 1 1])
        for tn1 = 1:length(ANA1.TN)
            for tn2 = tn1+1:length(ANA1.TN)
                [dxeye(tn1 , tn2) , ~] = dtw(ANA1.xEyePosDigit{tn1},ANA1.xEyePosDigit{tn2});
                dxeye(tn2 , tn1)  = dxeye(tn1 , tn2);
                
                [dpress(tn1 , tn2) , ~] = dtw(ANA1.PressTimeSeries{tn1},ANA1.PressTimeSeries{tn2});
                dpress(tn2 , tn1)  = dpress(tn1 , tn2);
                [tn1 tn2]
                %                 imagesc(dxeye)
                %                 drawnow()
            end
            
        end
        out.seqNumb = ANA1.seqNumb;
        out.dxeye = dxeye;
        out.dpress = dpress;
    case 'eyeEndPos_ChunkLength'
        %%
        ANA = getrow(Dall , ismember(Dall.seqNumb , [1:6]) & Dall.Group == 1 & ismember(Dall.Rep , rep) & ~Dall.isError & Dall.isgood & ismember(Dall.Day , days{day}));
        ANA.EndChunkLength = NaN*ones(length(ANA.TN) , 1);
        ANA.EyeLastPos = NaN*ones(length(ANA.TN) , 1);
        ANA.numChunks = NaN*ones(length(ANA.TN) , 1);
        ChnkArrang = CMB.DesiredCnhkargmnt(CMB.Group1 , :);
        for j = 1:size(ChnkArrang , 1)
            temp = [];
            temp1 = [];
            for k = 1:length(find(ChnkArrang(j , :)))
                temp = [temp k*ones(1,ChnkArrang(j,k))];
                temp1 = [temp1 1:ChnkArrang(j,k)];
            end
            ChnkArrng(j , :) = temp;
            ChnkPlcmnt(j,:)  = temp1;
        end
        IPIarrangement = diff(ChnkArrng , 1 , 2);   % between IPIs will have the index of 1
        IPIarrangement(~IPIarrangement) = 2;         % within IPIs will have the index of 2
        id = find(ismember(ANA.seqNumb , [1:6]));
        
        for tn = 1:length(id)
            ANA.EndChunkLength(id(tn)) = ChnkPlcmnt(ANA.seqNumb(id(tn)) , end);
            ANA.EyeLastPos(id(tn) ,1) = ANA.EyePressTimePos(id(tn) , end);
            ANA.numChunks(id(tn) ,1) = length(find(ANA.ChnkPlcmnt(tn,:) == 1));
        end

        out.endEffect1 = anovaMixed(ANA.EyeLastPos , ANA.SN,'within',[ANA.EndChunkLength],{'Last Chunk Length'},'subset',~isnan(ANA.EyeLastPos),'intercept',1)  ;
        out.endEffect2 = anovaMixed(ANA.EyeLastPos , ANA.SN,'within',[ANA.numChunks],{'Number of Chunks'},'subset',~isnan(ANA.EyeLastPos),'intercept',1)  ;
    case 'eyeFrstPos_ChunkLength'
        ANA = getrow(Dall , Dall.Group == 1 & ismember(Dall.Rep , rep) & ~Dall.isError & Dall.isgood & ismember(Dall.Day , days{day}));
        
        ANA.FrstChunkLength = NaN*ones(length(ANA.TN) , 1);
        ANA.EyeFrstPos = NaN*ones(length(ANA.TN) , 1);
        ANA.numChunks = NaN*ones(length(ANA.TN) , 1);
        ChnkArrang = CMB.DesiredCnhkargmnt(CMB.Group1 , :);
        for j = 1:size(ChnkArrang , 1)
            temp = [];
            temp1 = [];
            for k = 1:length(find(ChnkArrang(j , :)))
                temp = [temp k*ones(1,ChnkArrang(j,k))];
                temp1 = [temp1 1:ChnkArrang(j,k)];
            end
            ChnkArrng(j , :) = temp;
            ChnkPlcmnt(j,:)  = temp1;
        end
        IPIarrangement = diff(ChnkArrng , 1 , 2);   % between IPIs will have the index of 1
        IPIarrangement(~IPIarrangement) = 2;         % within IPIs will have the index of 2
        id = find(ismember(ANA.seqNumb , [1:6]));
        for tn = 1:length(id)
            [a,b] = unique(ChnkArrng(ANA.seqNumb(id(tn)),:));
            ANA.FrstChunkLength(id(tn)) = b(2)-1;
            ANA.EyeFrstPos(id(tn) ,1) = ANA.EyePressTimePos(id(tn) , 1);
            ANA.numChunks(id(tn) ,1) = length(find(ANA.ChnkPlcmnt(tn,:) == 1));
        end
        A1 = ANA;
        
        
        out.strtEffect1 = anovaMixed(ANA.EyeFrstPos , ANA.SN,'within',[ANA.FrstChunkLength],{'Last Chunk Length'},'subset',~isnan(ANA.EyeFrstPos),'intercept',1)  ;
        out.strtEffect2 = anovaMixed(ANA.EyeFrstPos , ANA.SN,'within',[ANA.numChunks],{'Last Chunk Length'},'subset',~isnan(ANA.EyeFrstPos),'intercept',1)  ;
    case 'eye10p_crossvaldist_pos'
%         %         if isequal(subjnum , [1 3 4 5 6 7 8 9 10])
        %             subjnum = 11;
        %         end
        calcDist = 1;
        if calcDist
            switch rep
                case 1
                    load([baseDir , '/se1_tnorm_r1.mat'])
                case 2
                    load([baseDir , '/se1_tnorm_r2.mat'])
                otherwise
                    load([baseDir , '/se1_tnorm.mat'])
            end
            
            %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EYE
            seqs = [1:7];
            for d = 2:5
                dwithclass{d-1} = [];
                dbetclass{d-1} = [];
                for SubN = 1:length(subjnum)
                    out(SubN, d-1).first10p_eye= N.norm(d,subjnum(SubN)).eye;
                    X = [];clear dis
                    lab = [];
                    for s = 1:length(seqs)
                        X = [X ; out(SubN, d-1).first10p_eye{seqs(s)}(:,1:100)];
                        %                         X = [X ; reshape(smooth(out(SubN, d-1).eyeveloc{s}' , 10) ,size(out(SubN, d-1).eyeveloc{s}'))'];
                        lab = [lab ; seqs(s)*ones(size(out(SubN, d-1).first10p_eye{seqs(s)}(:,1:100) , 1) ,1)];
                    end
                    lab(lab>1) = 2;
                    seqs = unique(lab);
                    for i = 1:size(X,1)
                        id = ones(size(X , 1) , 1);
                        id(i) = 0;
                        Y = inpaint_nans(X(~id , :)); % The left out trial
                        X1 = X(id==1 , :);     % all the other trials
                        lab1 = lab(id==1);
                        s1 = lab(id ==0);
                        
                        SS = seqs(seqs~= s1);
                        for s = SS
                            m = nanmean(X1(lab1==s , :));
                            dbetclass{d-1} = [dbetclass{d-1} ;[pdist([Y;m], distance) s1 s SubN d-1]];%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                        end
                        m = nanmean(X1(lab1==s1 , :));
                        dwithclass{d-1} = [dwithclass{d-1} ;[pdist([Y;m], distance) s1 s1 SubN d-1]];%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                    end
                    alldist = [dwithclass{d-1}(dwithclass{d-1}(:,4) == SubN , :) ;dbetclass{d-1}(dbetclass{d-1}(:,4) == SubN , :)];
                    for l =  1:length(seqs)
                        for l1 =  1:length(seqs)
                            out(SubN, d-1).first10_D_eye(l,l1) = nanmean(alldist(alldist(:,2) == seqs(l) & alldist(:,3)==seqs(l1)) , 1);
                        end
                    end
                    [SubN d]
                end
                alldist = [dwithclass{d-1} ;dbetclass{d-1}];
                alldist = [alldist [ones(length(dwithclass{d-1}),1);zeros(length(dbetclass{d-1}),1)]];
                for l =  1:length(seqs)
                    for l1 =  1:length(seqs)
                        out(SubN+1, d-1).D_eyfirst10_D_eyee(l,l1) = nanmean(alldist(alldist(:,2) == l & alldist(:,3)==l1,1));
                    end
                end
                first10_D_eyeSig(d-1).veldist = anovaMixed(alldist(:,1) ,alldist(:,4) ,'within',alldist(:,6),{'Within/Between'},'intercept',1)  ;
            end
            allday_dwithclass = [];
            allday_dbetclass = [];
            for d = 1:4
                allday_dwithclass = [allday_dwithclass ; dwithclass{d}];
                allday_dbetclass = [allday_dbetclass ; dbetclass{d}];
            end
            h1 = figure('color' , 'white');
            for SubN = 1:length(subjnum)
                [xcoord,f10_ePLOTw(SubN , :),f10_eERRORw(SubN , :)] = lineplot(allday_dwithclass(allday_dwithclass(:,4) == SubN , 5) , allday_dwithclass(allday_dwithclass(:,4) == SubN , 1));
                hold on
                [xcoord,f10_ePLOTb(SubN , :),f10_eERRORb(SubN , :)] = lineplot(allday_dbetclass(allday_dbetclass(:,4) == SubN , 5) , allday_dbetclass(allday_dbetclass(:,4) == SubN , 1));
            end
            [xcoord,f10_ePLOTw(SubN+1 , :),f10_eERRORw(SubN+1 , :)] = lineplot(allday_dwithclass(: , 5) , allday_dwithclass(: , 1));
            hold on
            [xcoord,f10_ePLOTb(SubN+1 , :),f10_eERRORb(SubN+1 , :)] = lineplot(allday_dbetclass(: , 5) , allday_dbetclass(: , 1));
            close(h1);
            
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRESS
            seqs = [1:7];
            for d = 2:5
                dwithclass{d-1} = [];
                dbetclass{d-1} = [];
                for SubN = 1:length(subjnum)
                    out(SubN, d-1).last10p_eye= N.norm(d,subjnum(SubN)).eye;
                
                    X = [];clear dis
                    lab = [];
                    for s = 1:length(seqs)
                        X = [X ; out(SubN, d-1).last10p_eye{seqs(s)}(:,901:1000)];
                        lab = [lab ; seqs(s)*ones(size(out(SubN, d-1).last10p_eye{seqs(s)}(:,1:100) , 1) ,1)];
                    end
                    lab(lab>1) = 2;
                    seqs = unique(lab);
                    for i = 1:size(X,1)
                        id = ones(size(X , 1) , 1);
                        id(i) = 0;
                        Y = inpaint_nans(X(~id , :)); % The left out trial
                        X1 = X(id==1 , :);     % all the other trials
                        lab1 = lab(id==1);
                        s1 = lab(id ==0);
                        
                        SS = seqs(seqs~= s1);
                        for s = SS
                            m = nanmean(X1(lab1==s , :));
                            dbetclass{d-1} = [dbetclass{d-1} ;[pdist([Y;m], distance) s1 s SubN d-1]];%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                        end
                        m = nanmean(X1(lab1==s1 , :));
                        dwithclass{d-1} = [dwithclass{d-1} ;[pdist([Y;m], distance) s1 s1 SubN d-1]];%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                    end
                    alldist = [dwithclass{d-1}(dwithclass{d-1}(:,4) == SubN , :) ;dbetclass{d-1}(dbetclass{d-1}(:,4) == SubN , :)];
                    for l =  1:length(seqs)
                        for l1 =  1:length(seqs)
                            out(SubN, d-1).last10_D_eye(l,l1) = nanmean(alldist(alldist(:,2) == seqs(l) & alldist(:,3)==seqs(l1)) , 1);
                        end
                    end
                    [SubN d]
                end
                alldist = [dwithclass{d-1} ;dbetclass{d-1}];
                alldist = [alldist [ones(length(dwithclass{d-1}),1);zeros(length(dbetclass{d-1}),1)]];
                for l =  1:length(seqs)
                    for l1 =  1:length(seqs)
                        out(SubN+1, d-1).D_eylast10_D_eyee(l,l1) = nanmean(alldist(alldist(:,2) == l & alldist(:,3)==l1,1));
                    end
                end
                last10_D_eyeSig(d-1).veldist = anovaMixed(alldist(:,1) ,alldist(:,4) ,'within',alldist(:,6),{'Within/Between'},'intercept',1)  ;
            end
            allday_dwithclass = [];
            allday_dbetclass = [];
            for d = 1:4
                allday_dwithclass = [allday_dwithclass ; dwithclass{d}];
                allday_dbetclass = [allday_dbetclass ; dbetclass{d}];
            end
            h1 = figure('color' , 'white');
            for SubN = 1:length(subjnum)
                [xcoord,l10_ePLOTw(SubN , :),l10_eERRORw(SubN , :)] = lineplot(allday_dwithclass(allday_dwithclass(:,4) == SubN , 5) , allday_dwithclass(allday_dwithclass(:,4) == SubN , 1));
                hold on
                [xcoord,l10_ePLOTb(SubN , :),l10_eERRORb(SubN , :)] = lineplot(allday_dbetclass(allday_dbetclass(:,4) == SubN , 5) , allday_dbetclass(allday_dbetclass(:,4) == SubN , 1));
            end
            [xcoord,l10_ePLOTw(SubN+1 , :),l10_eERRORw(SubN+1 , :)] = lineplot(allday_dwithclass(: , 5) , allday_dwithclass(: , 1));
            hold on
            [xcoord,l10_ePLOTb(SubN+1 , :),l10_eERRORb(SubN+1 , :)] = lineplot(allday_dbetclass(: , 5) , allday_dbetclass(: , 1));
            close(h1);
            
        else
            seqs = [2:7];
            load([baseDir , '/eyepos_dbetclass.mat'])
            load([baseDir , '/eyepos_dwithclass.mat'])
            if isequall(distance , 'euclidean')
                load([baseDir , '/EUC_eyepos_dbetclass.mat'])
                load([baseDir , '/EUC_eyepos_dwithclass.mat'])
            end
            for d = 1:4
                alldist = [dwithclass{d} ;dbetclass{d}];
                alldist = [alldist [ones(length(dwithclass{d}),1);zeros(length(dbetclass{d}),1)]];
               
                eyesig(d).veldist = anovaMixed(alldist(:,1) ,alldist(:,4) ,'within',alldist(:,6),{'Within/Between'},'intercept',1)  ;
            end
            allday_dwithclass = [];
            allday_dbetclass = [];
            for d = 1:4
                allday_dwithclass = [allday_dwithclass ; dwithclass{d}];
                allday_dbetclass = [allday_dbetclass ; dbetclass{d}];
            end
            h1 = figure('color' , 'white');
            for SubN = 1:length(subjnum)
                [xcoord,ePLOTw(SubN , :),eERRORw(SubN , :)] = lineplot(allday_dwithclass(allday_dwithclass(:,4) == SubN , 5) , allday_dwithclass(allday_dwithclass(:,4) == SubN , 1));
                hold on
                [xcoord,ePLOTb(SubN , :),eERRORb(SubN , :)] = lineplot(allday_dbetclass(allday_dbetclass(:,4) == SubN , 5) , allday_dbetclass(allday_dbetclass(:,4) == SubN , 1));
            end
            [xcoord,ePLOTw(SubN+1 , :),eERRORw(SubN+1 , :)] = lineplot(allday_dwithclass(: , 5) , allday_dwithclass(: , 1));
            hold on
            [xcoord,ePLOTb(SubN+1 , :),eERRORb(SubN+1 , :)] = lineplot(allday_dbetclass(: , 5) , allday_dbetclass(: , 1));
            close(h1);
            
            
            for d = 1:4
                alldist = [dwithclass{d} ;dbetclass{d}];
                alldist = [alldist [ones(length(dwithclass{d}),1);zeros(length(dbetclass{d}),1)]];
                
                prssig(d).veldist = anovaMixed(alldist(:,1) ,alldist(:,4) ,'within',alldist(:,6),{'Within/Between'},'intercept',1)  ;
            end
            load([baseDir , '/prspos_dbetclass.mat'])
            load([baseDir , '/prspos_dwithclass.mat'])
            if isequall(distance , 'euclidean')
                load([baseDir , '/EUC_prspos_dbetclass.mat'])
                load([baseDir , '/EUC_prspos_dwithclass.mat'])
            end
            allday_dwithclass = [];
            allday_dbetclass = [];
            for d = 1:4
                allday_dwithclass = [allday_dwithclass ; dwithclass{d}];
                allday_dbetclass = [allday_dbetclass ; dbetclass{d}];
            end
            h1 = figure('color' , 'white');
            for SubN = 1:length(subjnum)
                [xcoord,pPLOTw(SubN , :),pERRORw(SubN , :)] = lineplot(allday_dwithclass(allday_dwithclass(:,4) == SubN , 5) , allday_dwithclass(allday_dwithclass(:,4) == SubN , 1));
                hold on
                [xcoord,pPLOTb(SubN , :),pERRORb(SubN , :)] = lineplot(allday_dbetclass(allday_dbetclass(:,4) == SubN , 5) , allday_dbetclass(allday_dbetclass(:,4) == SubN , 1));
            end
            [xcoord,pPLOTw(SubN+1 , :),pERRORw(SubN+1 , :)] = lineplot(allday_dwithclass(: , 5) , allday_dwithclass(: , 1));
            hold on
            [xcoord,pPLOTb(SubN+1 , :),pERRORb(SubN+1 , :)] = lineplot(allday_dbetclass(: , 5) , allday_dbetclass(: , 1));
            close(h1);
        end
       
        
            

        figure('color' , 'white')
        for SubN = subjnum
            subplot(4,3,SubN)
            errorbar(f10_ePLOTw(SubN , :),f10_eERRORw(SubN , :) , 'LineWidth' , 3 )
            grid on
            hold on
            errorbar(f10_ePLOTb(SubN , :),f10_eERRORb(SubN , :) , 'LineWidth' , 3 )
            title(['F10 Eye Velocity, Subject ' , num2str(SubN)])
            
            legend({'W' , 'B'})
            hold on
            ax = gca;
            ax.XTick = [1:4];
            ax.XTickLabel = {'Day 1' , 'Day 2' , 'Day 3', 'All Days'};
%             ax.FontSize = 16;
            ylabel('Average distance');
        end
        
        
        figure('color' , 'white')
        errorbar(f10_ePLOTw(SubN+1 , :),f10_eERRORw(SubN+1 , :) , 'LineWidth' , 3 )
        grid on
        hold on
        errorbar(f10_ePLOTb(SubN+1 , :),f10_eERRORb(SubN+1 , :) , 'LineWidth' , 3 )
        title(['f10 Subject Average Dissimilarity in Eye position,  p =' , num2str(first10_D_eyeSig(4).veldist.eff(2).p)])
        legend({'Within Class' , 'Between Class'})
        hold on
        ax = gca;
        ax.XTick = [1:4];
        ax.XTickLabel = {'Day 1' , 'Day 2' , 'Day 3', 'All Days'};
        ax.FontSize = 16;
        ylabel('Average distance');
        
        
        figure('color' , 'white')
        for SubN = subjnum
            subplot(4,3,SubN)
            errorbar(l10_ePLOTw(SubN , :),l10_eERRORw(SubN , :) , 'LineWidth' , 3 )
            grid on
            hold on
            errorbar(l10_ePLOTb(SubN , :),l10_eERRORb(SubN , :) , 'LineWidth' , 3 )
            title(['l10 eye Velocity, Subject ' , num2str(SubN)])
            
            legend({'W' , 'B'})
            hold on
            ax = gca;
            ax.XTick = [1:4];
            ax.XTickLabel = {'Day 1' , 'Day 2' , 'Day 3', 'All Days'};
%             ax.FontSize = 16;
            ylabel('Average distance');
        end
        
        
        figure('color' , 'white')
        errorbar(l10_ePLOTw(SubN+1 , :),l10_eERRORw(SubN+1 , :) , 'LineWidth' , 3 )
        grid on
        hold on
        errorbar(l10_ePLOTb(SubN+1 , :),l10_eERRORb(SubN+1 , :) , 'LineWidth' , 3 )
        title(['l10 eye Velocity,  p =' , num2str(last10_D_eyeSig(4).veldist.eff(2).p)])
        legend({'Within Class' , 'Between Class'})
        hold on
        ax = gca;
        ax.XTick = [1:4];
        ax.XTickLabel = {'Day 1' , 'Day 2' , 'Day 3', 'All Days'};
        ax.FontSize = 16;
        ylabel('Average distance');
        out = [];
    
    
    case 'crossvaldist_pos'
%         %         if isequal(subjnum , [1 3 4 5 6 7 8 9 10])
        %             subjnum = 11;
        %         end
        calcDist = 1;
        if calcDist
            switch rep
                case 1
                    load([baseDir , '/se1_tnorm_r1.mat'])
                case 2
                    load([baseDir , '/se1_tnorm_r2.mat'])
                otherwise
                    load([baseDir , '/se1_tnorm.mat'])
            end
            
            %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EYE
            seqs = [2:7];
            for d = 2:5
                dwithclass{d-1} = [];
                dbetclass{d-1} = [];
                for SubN = 1:length(subjnum)
                    out(SubN, d-1).eye= N.norm(d,subjnum(SubN)).eye;
                
                    X = [];clear dis
                    lab = [];
                    for s = 1:length(seqs)
                        X = [X ; out(SubN, d-1).eye{seqs(s)}];
                        lab = [lab ; s*ones(size(out(SubN, d-1).eye{seqs(s)} , 1) ,1)];
                    end
                    
                    for i = 1:size(X,1)
                        id = ones(size(X , 1) , 1);
                        id(i) = 0;
                        Y = inpaint_nans(X(~id , :)); % The left out trial
                        X1 = X(id==1 , :);     % all the other trials
                        lab1 = lab(id==1);
                        s1 = lab(id ==0);
                        
                        SS = seqs(seqs~= s1);
                        for s = SS
                            m = nanmean(X1(lab1==s , :));
                            dbetclass{d-1} = [dbetclass{d-1} ;[pdist([Y;m], distance) s1 s SubN d-1]];%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                        end
                        m = nanmean(X1(lab1==s1 , :));
                        dwithclass{d-1} = [dwithclass{d-1} ;[pdist([Y;m], distance) s1 s1 SubN d-1]];%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                    end
                    alldist = [dwithclass{d-1}(dwithclass{d-1}(:,4) == SubN , :) ;dbetclass{d-1}(dbetclass{d-1}(:,4) == SubN , :)];
                    for l =  1:length(seqs)
                        for l1 =  1:length(seqs)
                            out(SubN, d-1).D_eye(l,l1) = nanmean(alldist(alldist(:,2) == seqs(l) & alldist(:,3)==seqs(l1)) , 1);
                        end
                    end
                    [SubN d]
                end
                alldist = [dwithclass{d-1} ;dbetclass{d-1}];
                alldist = [alldist [ones(length(dwithclass{d-1}),1);zeros(length(dbetclass{d-1}),1)]];
                for l =  1:length(seqs)
                    for l1 =  1:length(seqs)
                        out(SubN+1, d-1).D_eye(l,l1) = nanmean(alldist(alldist(:,2) == l & alldist(:,3)==l1,1));
                    end
                end
                eyesig(d-1).veldist = anovaMixed(alldist(:,1) ,alldist(:,4) ,'within',alldist(:,6),{'Within/Between'},'intercept',1)  ;
            end
            allday_dwithclass = [];
            allday_dbetclass = [];
            for d = 1:4
                allday_dwithclass = [allday_dwithclass ; dwithclass{d}];
                allday_dbetclass = [allday_dbetclass ; dbetclass{d}];
            end
            h1 = figure('color' , 'white');
            for SubN = 1:length(subjnum)
                [xcoord,ePLOTw(SubN , :),eERRORw(SubN , :)] = lineplot(allday_dwithclass(allday_dwithclass(:,4) == SubN , 5) , allday_dwithclass(allday_dwithclass(:,4) == SubN , 1));
                hold on
                [xcoord,ePLOTb(SubN , :),eERRORb(SubN , :)] = lineplot(allday_dbetclass(allday_dbetclass(:,4) == SubN , 5) , allday_dbetclass(allday_dbetclass(:,4) == SubN , 1));
            end
            [xcoord,ePLOTw(SubN+1 , :),eERRORw(SubN+1 , :)] = lineplot(allday_dwithclass(: , 5) , allday_dwithclass(: , 1));
            hold on
            [xcoord,ePLOTb(SubN+1 , :),eERRORb(SubN+1 , :)] = lineplot(allday_dbetclass(: , 5) , allday_dbetclass(: , 1));
            close(h1);
            
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRESS
            
            for d = 2:5
                dwithclass{d-1} = [];
                dbetclass{d-1} = [];
                for SubN = 1:length(subjnum)
                    out(SubN, d-1).press= N.norm(d,subjnum(SubN)).press;

                    X = [];clear dis
                    lab = [];
                    for s = 1:length(out(SubN, d-1).press)
                        X = [X ; out(SubN, d-1).press{s}];
                        %                         X = [X ; reshape(smooth(out(SubN, d-1).eyeveloc{s}' , 10) ,size(out(SubN, d-1).eyeveloc{s}'))'];
                        lab = [lab ; s*ones(size(out(SubN, d-1).press{s} , 1) ,1)];
                    end
                    
                    for i = 1:size(X,1)
                        id = ones(size(X , 1) , 1);
                        id(i) = 0;
                        Y = inpaint_nans(X(~id , :)); % The left out trial
                        X1 = X(id==1 , :);     % all the other trials
                        lab1 = lab(id==1);
                        s1 = lab(id ==0);
                        
                        SS = seqs(seqs~= s1);
                        for s = SS
                            m = nanmean(X1(lab1==s , :));
                            dbetclass{d-1} = [dbetclass{d-1} ;[pdist([Y;m], distance) s1 s SubN d-1]];%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                        end
                        m = nanmean(X1(lab1==s1 , :));
                        dwithclass{d-1} = [dwithclass{d-1} ;[pdist([Y;m], distance) s1 s1 SubN d-1]];%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                    end
                    alldist = [dwithclass{d-1}(dwithclass{d-1}(:,4) == SubN , :) ;dbetclass{d-1}(dbetclass{d-1}(:,4) == SubN , :)];
                    for l =  1:length(seqs)
                        for l1 =  1:length(seqs)
                            out(SubN, d-1).D_eye(l,l1) = nanmean(alldist(alldist(:,2) == seqs(l) & alldist(:,3)==seqs(l1)) , 1);
                        end
                    end
                    [SubN d]
                end
                alldist = [dwithclass{d-1} ;dbetclass{d-1}];
                alldist = [alldist [ones(length(dwithclass{d-1}),1);zeros(length(dbetclass{d-1}),1)]];
                for l =  1:length(seqs)
                    for l1 =  1:length(seqs)
                        out(SubN+1, d-1).D_eye(l,l1) = nanmean(alldist(alldist(:,2) == l & alldist(:,3)==l1,1));
                    end
                end
                prssig(d-1).veldist = anovaMixed(alldist(:,1) ,alldist(:,4) ,'within',alldist(:,6),{'Within/Between'},'intercept',1)  ;
            end
            allday_dwithclass = [];
            allday_dbetclass = [];
            for d = 1:4
                allday_dwithclass = [allday_dwithclass ; dwithclass{d}];
                allday_dbetclass = [allday_dbetclass ; dbetclass{d}];
            end
            h1 = figure('color' , 'white');
            for SubN = 1:length(subjnum)
                [xcoord,pPLOTw(SubN , :),pERRORw(SubN , :)] = lineplot(allday_dwithclass(allday_dwithclass(:,4) == SubN , 5) , allday_dwithclass(allday_dwithclass(:,4) == SubN , 1));
                hold on
                [xcoord,pPLOTb(SubN , :),pERRORb(SubN , :)] = lineplot(allday_dbetclass(allday_dbetclass(:,4) == SubN , 5) , allday_dbetclass(allday_dbetclass(:,4) == SubN , 1));
            end
            [xcoord,pPLOTw(SubN+1 , :),pERRORw(SubN+1 , :)] = lineplot(allday_dwithclass(: , 5) , allday_dwithclass(: , 1));
            hold on
            [xcoord,pPLOTb(SubN+1 , :),pERRORb(SubN+1 , :)] = lineplot(allday_dbetclass(: , 5) , allday_dbetclass(: , 1));
            close(h1);
            
        else
            seqs = [2:7];
            load([baseDir , '/eyepos_dbetclass.mat'])
            load([baseDir , '/eyepos_dwithclass.mat'])
            if isequall(distance , 'euclidean')
                load([baseDir , '/EUC_eyepos_dbetclass.mat'])
                load([baseDir , '/EUC_eyepos_dwithclass.mat'])
            end
            for d = 1:4
                alldist = [dwithclass{d} ;dbetclass{d}];
                alldist = [alldist [ones(length(dwithclass{d}),1);zeros(length(dbetclass{d}),1)]];
               
                eyesig(d).veldist = anovaMixed(alldist(:,1) ,alldist(:,4) ,'within',alldist(:,6),{'Within/Between'},'intercept',1)  ;
            end
            allday_dwithclass = [];
            allday_dbetclass = [];
            for d = 1:4
                allday_dwithclass = [allday_dwithclass ; dwithclass{d}];
                allday_dbetclass = [allday_dbetclass ; dbetclass{d}];
            end
            h1 = figure('color' , 'white');
            for SubN = 1:length(subjnum)
                [xcoord,ePLOTw(SubN , :),eERRORw(SubN , :)] = lineplot(allday_dwithclass(allday_dwithclass(:,4) == SubN , 5) , allday_dwithclass(allday_dwithclass(:,4) == SubN , 1));
                hold on
                [xcoord,ePLOTb(SubN , :),eERRORb(SubN , :)] = lineplot(allday_dbetclass(allday_dbetclass(:,4) == SubN , 5) , allday_dbetclass(allday_dbetclass(:,4) == SubN , 1));
            end
            [xcoord,ePLOTw(SubN+1 , :),eERRORw(SubN+1 , :)] = lineplot(allday_dwithclass(: , 5) , allday_dwithclass(: , 1));
            hold on
            [xcoord,ePLOTb(SubN+1 , :),eERRORb(SubN+1 , :)] = lineplot(allday_dbetclass(: , 5) , allday_dbetclass(: , 1));
            close(h1);
            
            
            for d = 1:4
                alldist = [dwithclass{d} ;dbetclass{d}];
                alldist = [alldist [ones(length(dwithclass{d}),1);zeros(length(dbetclass{d}),1)]];
                
                prssig(d).veldist = anovaMixed(alldist(:,1) ,alldist(:,4) ,'within',alldist(:,6),{'Within/Between'},'intercept',1)  ;
            end
            load([baseDir , '/prspos_dbetclass.mat'])
            load([baseDir , '/prspos_dwithclass.mat'])
            if isequall(distance , 'euclidean')
                load([baseDir , '/EUC_prspos_dbetclass.mat'])
                load([baseDir , '/EUC_prspos_dwithclass.mat'])
            end
            allday_dwithclass = [];
            allday_dbetclass = [];
            for d = 1:4
                allday_dwithclass = [allday_dwithclass ; dwithclass{d}];
                allday_dbetclass = [allday_dbetclass ; dbetclass{d}];
            end
            h1 = figure('color' , 'white');
            for SubN = 1:length(subjnum)
                [xcoord,pPLOTw(SubN , :),pERRORw(SubN , :)] = lineplot(allday_dwithclass(allday_dwithclass(:,4) == SubN , 5) , allday_dwithclass(allday_dwithclass(:,4) == SubN , 1));
                hold on
                [xcoord,pPLOTb(SubN , :),pERRORb(SubN , :)] = lineplot(allday_dbetclass(allday_dbetclass(:,4) == SubN , 5) , allday_dbetclass(allday_dbetclass(:,4) == SubN , 1));
            end
            [xcoord,pPLOTw(SubN+1 , :),pERRORw(SubN+1 , :)] = lineplot(allday_dwithclass(: , 5) , allday_dwithclass(: , 1));
            hold on
            [xcoord,pPLOTb(SubN+1 , :),pERRORb(SubN+1 , :)] = lineplot(allday_dbetclass(: , 5) , allday_dbetclass(: , 1));
            close(h1);
            
            
            
        end
       
        
            

        figure('color' , 'white')
        for SubN = subjnum
            subplot(4,3,SubN)
            errorbar(ePLOTw(SubN , :),eERRORw(SubN , :) , 'LineWidth' , 3 )
            grid on
            hold on
            errorbar(ePLOTb(SubN , :),eERRORb(SubN , :) , 'LineWidth' , 3 )
            title(['Eye Velocity, Subject ' , num2str(SubN)])
            
            legend({'W' , 'B'})
            hold on
            ax = gca;
            ax.XTick = [1:4];
            ax.XTickLabel = {'Day 1' , 'Day 2' , 'Day 3', 'All Days'};
%             ax.FontSize = 16;
            ylabel('Average distance');
        end
        
        
        figure('color' , 'white')
        errorbar(ePLOTw(SubN+1 , :),eERRORw(SubN+1 , :) , 'LineWidth' , 3 )
        grid on
        hold on
        errorbar(ePLOTb(SubN+1 , :),eERRORb(SubN+1 , :) , 'LineWidth' , 3 )
        title(['Subject Average Dissimilarity in Eye Velocity,  p =' , num2str(eyesig(4).veldist.eff(2).p)])
        legend({'Within Class' , 'Between Class'})
        hold on
        ax = gca;
        ax.XTick = [1:4];
        ax.XTickLabel = {'Day 1' , 'Day 2' , 'Day 3', 'All Days'};
        ax.FontSize = 16;
        ylabel('Average distance');
        
        
        figure('color' , 'white')
        for SubN = subjnum
            subplot(4,3,SubN)
            errorbar(pPLOTw(SubN , :),pERRORw(SubN , :) , 'LineWidth' , 3 )
            grid on
            hold on
            errorbar(pPLOTb(SubN , :),pERRORb(SubN , :) , 'LineWidth' , 3 )
            title(['Press Velocity, Subject ' , num2str(SubN)])
            
            legend({'W' , 'B'})
            hold on
            ax = gca;
            ax.XTick = [1:4];
            ax.XTickLabel = {'Day 1' , 'Day 2' , 'Day 3', 'All Days'};
%             ax.FontSize = 16;
            ylabel('Average distance');
        end
        
        
        figure('color' , 'white')
        errorbar(pPLOTw(SubN+1 , :),pERRORw(SubN+1 , :) , 'LineWidth' , 3 )
        grid on
        hold on
        errorbar(pPLOTb(SubN+1 , :),pERRORb(SubN+1 , :) , 'LineWidth' , 3 )
        title(['Subject Average Dissimilarity in Press Velocity,  p =' , num2str(prssig(4).veldist.eff(2).p)])
        legend({'Within Class' , 'Between Class'})
        hold on
        ax = gca;
        ax.XTick = [1:4];
        ax.XTickLabel = {'Day 1' , 'Day 2' , 'Day 3', 'All Days'};
        ax.FontSize = 16;
        ylabel('Average distance');
        out = [];
    case 'crossval_IPI_dist'
        ANA = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.seqNumb , [0:6]) & ~Dall.isError & ismember(Dall.Rep , rep) & ismember(Dall.Day , days{day}));
        for tn = 1:length(ANA.TN)
            n = (ANA.AllPressIdx(tn , sum(~isnan(ANA.AllPressIdx(tn , :))))  - ANA.AllPressIdx(tn , 1)) / 1000;
            nIdx(tn , :) = (ANA.AllPressIdx(tn , :) - ANA.AllPressIdx(tn , 1))/n;
            ANA.IPI(tn , :) = diff(nIdx(tn ,:) , 1 , 2);
        end
        X = ANA.IPI;clear d
        lab = ANA.seqNumb +1;
        
        for i = 1:size(X,1)
            id = ones(size(X , 1) , 1);
            id(i) = 0;
            Y = inpaint_nans(X(~id , :));
            X1 = X(id==1 , :);
            lab1 = lab(id==1);
            for s = 1:7
                m = nanmean(X1(lab1==s , :));
                d(i , s) = pdist([Y;m], distance);%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
            end
            
        end
        seqs =  unique(lab);
        for l =  1:length(seqs)
            for l1 =  1:length(seqs)
                id = lab == l;s
                out.D_IPI(l,l1) = nanmean(d(id , l1));
            end
        end
        figure('color' , [1 1 1])
        imagesc(out.D_IPI);
        title('Crossvalidated IPI Dissimilarity Matrix' , 'FontSize' , 20)
        hold on
        ax = gca;
        ax.XTick = [1:7];
        ax.YTick = [1:7];
        ax.XTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
        ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
        ax.XTickLabelRotation = 45;
        ax.FontSize = 20;
        hold on
        line([1.5 7.5] , [1.5 1.5] ,  'LineWidth' , 3 , 'color' , [0 0 0])
        line([1.5 1.5] ,[1.5 7.5] ,  'LineWidth' , 3 , 'color' , [0 0 0])
        axis square
        colorbar
        
    case 'crossvaldist_vel'
        %         if isequal(subjnum , [1 3 4 5 6 7 8 9 10])
        %             subjnum = 11;
        %         end
        calcDist = 0;
        if calcDist
            switch rep
                case 1
                    load([baseDir , '/se1_tnorm_r1.mat'])
                case 2
                    load([baseDir , '/se1_tnorm_r2.mat'])
                otherwise
                    load([baseDir , '/se1_tnorm.mat'])
            end
            
            %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EYE
            seqs = [2:7];
            for d = 2:5
                dwithclass{d-1} = [];
                dbetclass{d-1} = [];
                for SubN = 1:length(subjnum)
                    out(SubN, d-1).eyeveloc= N.norm(d,subjnum(SubN)).eyeveloc;
                
                    X = [];clear dis
                    lab = [];
                    for s = 1:length(out(SubN, d-1).eyeveloc)
                        X = [X ; out(SubN, d-1).eyeveloc{s}];
                        %                         X = [X ; reshape(smooth(out(SubN, d-1).eyeveloc{s}' , 10) ,size(out(SubN, d-1).eyeveloc{s}'))'];
                        lab = [lab ; s*ones(size(out(SubN, d-1).eyeveloc{s} , 1) ,1)];
                    end
                    
                    for i = 1:size(X,1)
                        id = ones(size(X , 1) , 1);
                        id(i) = 0;
                        Y = inpaint_nans(X(~id , :)); % The left out trial
                        X1 = X(id==1 , :);     % all the other trials
                        lab1 = lab(id==1);
                        s1 = lab(id ==0);
                        
                        SS = seqs(seqs~= s1);
                        for s = SS
                            m = nanmean(X1(lab1==s , :));
                            dbetclass{d-1} = [dbetclass{d-1} ;[pdist([Y;m], distance) s1 s SubN d-1]];%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                        end
                        m = nanmean(X1(lab1==s1 , :));
                        dwithclass{d-1} = [dwithclass{d-1} ;[pdist([Y;m], distance) s1 s1 SubN d-1]];%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                    end
                    alldist = [dwithclass{d-1}(dwithclass{d-1}(:,4) == SubN , :) ;dbetclass{d-1}(dbetclass{d-1}(:,4) == SubN , :)];
                    for l =  1:length(seqs)
                        for l1 =  1:length(seqs)
                            out(SubN, d-1).D_eye(l,l1) = nanmean(alldist(alldist(:,2) == seqs(l) & alldist(:,3)==seqs(l1)) , 1);
                        end
                    end
                    [SubN d]
                end
                alldist = [dwithclass{d-1} ;dbetclass{d-1}];
                alldist = [alldist [ones(length(dwithclass{d-1}),1);zeros(length(dbetclass{d-1}),1)]];
                for l =  1:length(seqs)
                    for l1 =  1:length(seqs)
                        out(SubN+1, d-1).D_eye(l,l1) = nanmean(alldist(alldist(:,2) == l & alldist(:,3)==l1,1));
                    end
                end
                eyesig(d-1).veldist = anovaMixed(alldist(:,1) ,alldist(:,4) ,'within',alldist(:,6),{'Within/Between'},'intercept',1)  ;
            end
            allday_dwithclass = [];
            allday_dbetclass = [];
            for d = 1:4
                allday_dwithclass = [allday_dwithclass ; dwithclass{d}];
                allday_dbetclass = [allday_dbetclass ; dbetclass{d}];
            end
            h1 = figure('color' , 'white');
            for SubN = 1:length(subjnum)
                [xcoord,ePLOTw(SubN , :),eERRORw(SubN , :)] = lineplot(allday_dwithclass(allday_dwithclass(:,4) == SubN , 5) , allday_dwithclass(allday_dwithclass(:,4) == SubN , 1));
                hold on
                [xcoord,ePLOTb(SubN , :),eERRORb(SubN , :)] = lineplot(allday_dbetclass(allday_dbetclass(:,4) == SubN , 5) , allday_dbetclass(allday_dbetclass(:,4) == SubN , 1));
            end
            [xcoord,ePLOTw(SubN+1 , :),eERRORw(SubN+1 , :)] = lineplot(allday_dwithclass(: , 5) , allday_dwithclass(: , 1));
            hold on
            [xcoord,ePLOTb(SubN+1 , :),eERRORb(SubN+1 , :)] = lineplot(allday_dbetclass(: , 5) , allday_dbetclass(: , 1));
            close(h1);
            
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRESS
            
            for d = 2:5
                dwithclass{d-1} = [];
                dbetclass{d-1} = [];
                for SubN = 1:length(subjnum)
                    out(SubN, d-1).prsveloc= N.norm(d,subjnum(SubN)).prsveloc;

                    X = [];clear dis
                    lab = [];
                    for s = 1:length(out(SubN, d-1).prsveloc)
                        X = [X ; out(SubN, d-1).prsveloc{s}];
                        %                         X = [X ; reshape(smooth(out(SubN, d-1).eyeveloc{s}' , 10) ,size(out(SubN, d-1).eyeveloc{s}'))'];
                        lab = [lab ; s*ones(size(out(SubN, d-1).prsveloc{s} , 1) ,1)];
                    end
                    
                    for i = 1:size(X,1)
                        id = ones(size(X , 1) , 1);
                        id(i) = 0;
                        Y = inpaint_nans(X(~id , :)); % The left out trial
                        X1 = X(id==1 , :);     % all the other trials
                        lab1 = lab(id==1);
                        s1 = lab(id ==0);
                        
                        SS = seqs(seqs~= s1);
                        for s = SS
                            m = nanmean(X1(lab1==s , :));
                            dbetclass{d-1} = [dbetclass{d-1} ;[pdist([Y;m], distance) s1 s SubN d-1]];%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                        end
                        m = nanmean(X1(lab1==s1 , :));
                        dwithclass{d-1} = [dwithclass{d-1} ;[pdist([Y;m], distance) s1 s1 SubN d-1]];%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                    end
                    alldist = [dwithclass{d-1}(dwithclass{d-1}(:,4) == SubN , :) ;dbetclass{d-1}(dbetclass{d-1}(:,4) == SubN , :)];
                    for l =  1:length(seqs)
                        for l1 =  1:length(seqs)
                            out(SubN, d-1).D_eye(l,l1) = nanmean(alldist(alldist(:,2) == seqs(l) & alldist(:,3)==seqs(l1)) , 1);
                        end
                    end
                    [SubN d]
                end
                alldist = [dwithclass{d-1} ;dbetclass{d-1}];
                alldist = [alldist [ones(length(dwithclass{d-1}),1);zeros(length(dbetclass{d-1}),1)]];
                for l =  1:length(seqs)
                    for l1 =  1:length(seqs)
                        out(SubN+1, d-1).D_eye(l,l1) = nanmean(alldist(alldist(:,2) == l & alldist(:,3)==l1,1));
                    end
                end
                prssig(d-1).veldist = anovaMixed(alldist(:,1) ,alldist(:,4) ,'within',alldist(:,6),{'Within/Between'},'intercept',1)  ;
            end
            allday_dwithclass = [];
            allday_dbetclass = [];
            for d = 1:4
                allday_dwithclass = [allday_dwithclass ; dwithclass{d}];
                allday_dbetclass = [allday_dbetclass ; dbetclass{d}];
            end
            h1 = figure('color' , 'white');
            for SubN = 1:length(subjnum)
                [xcoord,pPLOTw(SubN , :),pERRORw(SubN , :)] = lineplot(allday_dwithclass(allday_dwithclass(:,4) == SubN , 5) , allday_dwithclass(allday_dwithclass(:,4) == SubN , 1));
                hold on
                [xcoord,pPLOTb(SubN , :),pERRORb(SubN , :)] = lineplot(allday_dbetclass(allday_dbetclass(:,4) == SubN , 5) , allday_dbetclass(allday_dbetclass(:,4) == SubN , 1));
            end
            [xcoord,pPLOTw(SubN+1 , :),pERRORw(SubN+1 , :)] = lineplot(allday_dwithclass(: , 5) , allday_dwithclass(: , 1));
            hold on
            [xcoord,pPLOTb(SubN+1 , :),pERRORb(SubN+1 , :)] = lineplot(allday_dbetclass(: , 5) , allday_dbetclass(: , 1));
            close(h1);
            
        else
            seqs = [2:7];
            load([baseDir , '/eye_dbetclass.mat'])
            load([baseDir , '/eye_dwithclass.mat'])
            for d = 1:4
                alldist = [dwithclass{d} ;dbetclass{d}];
                alldist = [alldist [ones(length(dwithclass{d}),1);zeros(length(dbetclass{d}),1)]];
               
                eyesig(d).veldist = anovaMixed(alldist(:,1) ,alldist(:,4) ,'within',alldist(:,6),{'Within/Between'},'intercept',1)  ;
            end
            allday_dwithclass = [];
            allday_dbetclass = [];
            for d = 1:4
                allday_dwithclass = [allday_dwithclass ; dwithclass{d}];
                allday_dbetclass = [allday_dbetclass ; dbetclass{d}];
            end
            h1 = figure('color' , 'white');
            for SubN = 1:length(subjnum)
                [xcoord,ePLOTw(SubN , :),eERRORw(SubN , :)] = lineplot(allday_dwithclass(allday_dwithclass(:,4) == SubN , 5) , allday_dwithclass(allday_dwithclass(:,4) == SubN , 1));
                hold on
                [xcoord,ePLOTb(SubN , :),eERRORb(SubN , :)] = lineplot(allday_dbetclass(allday_dbetclass(:,4) == SubN , 5) , allday_dbetclass(allday_dbetclass(:,4) == SubN , 1));
            end
            [xcoord,ePLOTw(SubN+1 , :),eERRORw(SubN+1 , :)] = lineplot(allday_dwithclass(: , 5) , allday_dwithclass(: , 1));
            hold on
            [xcoord,ePLOTb(SubN+1 , :),eERRORb(SubN+1 , :)] = lineplot(allday_dbetclass(: , 5) , allday_dbetclass(: , 1));
            close(h1);
            
            
            for d = 1:4
                alldist = [dwithclass{d} ;dbetclass{d}];
                alldist = [alldist [ones(length(dwithclass{d}),1);zeros(length(dbetclass{d}),1)]];
                
                prssig(d).veldist = anovaMixed(alldist(:,1) ,alldist(:,4) ,'within',alldist(:,6),{'Within/Between'},'intercept',1)  ;
            end
            load([baseDir , '/prs_dbetclass.mat'])
            load([baseDir , '/prs_dwithclass.mat'])
            allday_dwithclass = [];
            allday_dbetclass = [];
            for d = 1:4
                allday_dwithclass = [allday_dwithclass ; dwithclass{d}];
                allday_dbetclass = [allday_dbetclass ; dbetclass{d}];
            end
            h1 = figure('color' , 'white');
            for SubN = 1:length(subjnum)
                [xcoord,pPLOTw(SubN , :),pERRORw(SubN , :)] = lineplot(allday_dwithclass(allday_dwithclass(:,4) == SubN , 5) , allday_dwithclass(allday_dwithclass(:,4) == SubN , 1));
                hold on
                [xcoord,pPLOTb(SubN , :),pERRORb(SubN , :)] = lineplot(allday_dbetclass(allday_dbetclass(:,4) == SubN , 5) , allday_dbetclass(allday_dbetclass(:,4) == SubN , 1));
            end
            [xcoord,pPLOTw(SubN+1 , :),pERRORw(SubN+1 , :)] = lineplot(allday_dwithclass(: , 5) , allday_dwithclass(: , 1));
            hold on
            [xcoord,pPLOTb(SubN+1 , :),pERRORb(SubN+1 , :)] = lineplot(allday_dbetclass(: , 5) , allday_dbetclass(: , 1));
            close(h1);
            
            
            
        end
       
        
            

        figure('color' , 'white')
        for SubN = subjnum
            subplot(4,3,SubN)
            errorbar(ePLOTw(SubN , :),eERRORw(SubN , :) , 'LineWidth' , 3 )
            grid on
            hold on
            errorbar(ePLOTb(SubN , :),eERRORb(SubN , :) , 'LineWidth' , 3 )
            title(['Eye Velocity, Subject ' , num2str(SubN)])
            
            legend({'W' , 'B'})
            hold on
            ax = gca;
            ax.XTick = [1:4];
            ax.XTickLabel = {'Day 1' , 'Day 2' , 'Day 3', 'All Days'};
%             ax.FontSize = 16;
            ylabel('Average distance');
        end
        
        
        figure('color' , 'white')
        errorbar(ePLOTw(SubN+1 , :),eERRORw(SubN+1 , :) , 'LineWidth' , 3 )
        grid on
        hold on
        errorbar(ePLOTb(SubN+1 , :),eERRORb(SubN+1 , :) , 'LineWidth' , 3 )
        title(['Subject Average Dissimilarity in Eye Velocity,  p =' , num2str(eyesig(4).veldist.eff(2).p)])
        legend({'Within Class' , 'Between Class'})
        hold on
        ax = gca;
        ax.XTick = [1:4];
        ax.XTickLabel = {'Day 1' , 'Day 2' , 'Day 3', 'All Days'};
        ax.FontSize = 16;
        ylabel('Average distance');
        
        
        figure('color' , 'white')
        for SubN = subjnum
            subplot(4,3,SubN)
            errorbar(pPLOTw(SubN , :),pERRORw(SubN , :) , 'LineWidth' , 3 )
            grid on
            hold on
            errorbar(pPLOTb(SubN , :),pERRORb(SubN , :) , 'LineWidth' , 3 )
            title(['Press Velocity, Subject ' , num2str(SubN)])
            
            legend({'W' , 'B'})
            hold on
            ax = gca;
            ax.XTick = [1:4];
            ax.XTickLabel = {'Day 1' , 'Day 2' , 'Day 3', 'All Days'};
%             ax.FontSize = 16;
            ylabel('Average distance');
        end
        
        
        figure('color' , 'white')
        errorbar(pPLOTw(SubN+1 , :),pERRORw(SubN+1 , :) , 'LineWidth' , 3 )
        grid on
        hold on
        errorbar(pPLOTb(SubN+1 , :),pERRORb(SubN+1 , :) , 'LineWidth' , 3 )
        title(['Subject Average Dissimilarity in Press Velocity,  p =' , num2str(prssig(4).veldist.eff(2).p)])
        legend({'Within Class' , 'Between Class'})
        hold on
        ax = gca;
        ax.XTick = [1:4];
        ax.XTickLabel = {'Day 1' , 'Day 2' , 'Day 3', 'All Days'};
        ax.FontSize = 16;
        ylabel('Average distance');
        out = [];
    case 'crossval_IPI_dist'
        ANA = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.seqNumb , [0:6]) & ~Dall.isError & ismember(Dall.Rep , rep) & ismember(Dall.Day , days{day}));
        for tn = 1:length(ANA.TN)
            n = (ANA.AllPressIdx(tn , sum(~isnan(ANA.AllPressIdx(tn , :))))  - ANA.AllPressIdx(tn , 1)) / 1000;
            nIdx(tn , :) = (ANA.AllPressIdx(tn , :) - ANA.AllPressIdx(tn , 1))/n;
            ANA.IPI(tn , :) = diff(nIdx(tn ,:) , 1 , 2);
        end
        X = ANA.IPI;clear d
        lab = ANA.seqNumb +1;
        
        for i = 1:size(X,1)
            id = ones(size(X , 1) , 1);
            id(i) = 0;
            Y = inpaint_nans(X(~id , :));
            X1 = X(id==1 , :);
            lab1 = lab(id==1);
            for s = 1:7
                m = nanmean(X1(lab1==s , :));
                d(i , s) = pdist([Y;m], distance);%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
            end
            
        end
        seqs =  unique(lab);
        for l =  1:length(seqs)
            for l1 =  1:length(seqs)
                id = lab == l;s
                out.D_IPI(l,l1) = nanmean(d(id , l1));
            end
        end
        figure('color' , [1 1 1])
        imagesc(out.D_IPI);
        title('Crossvalidated IPI Dissimilarity Matrix' , 'FontSize' , 20)
        hold on
        ax = gca;
        ax.XTick = [1:7];
        ax.YTick = [1:7];
        ax.XTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
        ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
        ax.XTickLabelRotation = 45;
        ax.FontSize = 20;
        hold on
        line([1.5 7.5] , [1.5 1.5] ,  'LineWidth' , 3 , 'color' , [0 0 0])
        line([1.5 1.5] ,[1.5 7.5] ,  'LineWidth' , 3 , 'color' , [0 0 0])
        axis square
        colorbar
    case 'crossvaldist_chunk'
        %% chunk distances
        ANA = getrow(Dall ,ismember(Dall.seqNumb , [0:6]) & ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.Rep , rep));
        for tn = 1:length(ANA.TN)
            n = (ANA.AllPressIdx(tn , sum(~isnan(ANA.AllPressIdx(tn , :))))  - ANA.AllPressIdx(tn , 1)) / 1000;
            nIdx(tn , :) = (ANA.AllPressIdx(tn , :) - ANA.AllPressIdx(tn , 1))/n;
            IPI(tn , :) = detrend(diff(nIdx(tn ,:) , 1 , 2) , 'linear' , 2);
        end
        
        for tn = 1:length(ANA.TN)
            thresh = .3 * std(IPI(tn , :));
            [dum , estChnkBndry] = findpeaks(IPI(tn , :));% ,'MinPeakProminence', thresh); % slope of IPIs
            
            if ~isempty(estChnkBndry)
                goodpeak = ones(1, length(estChnkBndry));
                for cb = 1:length(estChnkBndry)
                    if IPI(tn , estChnkBndry(cb)) < nanmean(IPI(tn  ,:)) + thresh
                        goodpeak(cb) = 0;
                    end
                end
                if sum(goodpeak)
                    estChnkBndry = estChnkBndry(logical(goodpeak));
                else
                    estChnkBndry = [];
                end
            end
            
            
            ANA.estChnkBndry(tn , :) = zeros(1, 14);  % first presses of chunks will be 1
            ANA.estChnkBndry(tn , estChnkBndry+1) = 1;
            ANA.estChnkBndry(tn , 1) = 1;
        end
        
        if GroupCode == 1
            ChnkArrang = CMB.DesiredCnhkargmnt(CMB.Group1 , :);
            for j = 1:size(ChnkArrang , 1)
                temp = [];
                temp1 = [];
                for k = 1:length(find(ChnkArrang(j , :)))
                    temp = [temp k*ones(1,ChnkArrang(j,k))];
                    temp1 = [temp1 1:ChnkArrang(j,k)];
                end
                ChnkArrng(j , :) = temp;
                ChnkPlcmnt(j,:)  = temp1;
            end
            IPIarrangement = diff(ChnkArrng , 1 , 2);   % between IPIs will have the index of 1
            IPIarrangement(~IPIarrangement) = 2;         % within IPIs will have the index of 2
        elseif GroupCode == 2
            ChnkArrang = CMB.DesiredCnhkargmnt(CMB.Group2 , :);
            for j = 1:size(ChnkArrang , 1)
                temp = [];
                temp1 = [];
                for k = 1:length(find(ChnkArrang(j , :)))
                    temp = [temp k*ones(1,ChnkArrang(j,k))];
                    temp1 = [temp1 1:ChnkArrang(j,k)];
                end
                ChnkArrng(j , :) = temp;
                ChnkPlcmnt(j,:)  = temp1;
            end
            IPIarrangement = diff(ChnkArrng , 1 , 2); % between IPIs will have the index of 1
            IPIarrangement(~IPIarrangement) = 2;       % within IPIs will have the index of 2
        end
        
        
        
        %% the boundry distance
        temp = diff(ChnkPlcmnt,1,2);
        temp(temp<0) = 0;
        chbndry = [ones(6,1) ~temp]; % all the first presses are one
        
        
        
        for s = 0:6
            A = getrow(ANA , ANA.seqNumb == s);
            CBD{s+1} = A.estChnkBndry;
        end
        diag_IPI = [];
        diag_IPI_s = [];
        offdiag_IPI = [];
        offdiag_IPI_s = [];
        for d = 2:5
            clear meanbnd allmeanbnd CB est_dist act_est_dist dis lab lab1
            ANA1 = getrow(ANA , ismember(ANA.Day , days{d}));
            
            out(d-1).X = ANA1.estChnkBndry;
            lab = ANA1.seqNumb +1;
            
            for i = 1:size(out(d-1).X,1)
                id = ones(size(out(d-1).X , 1) , 1);
                id(i) = 0;
                Y = out(d-1).X(~id , :);
                X1 = out(d-1).X(id==1 , :);
                lab1 = lab(id==1);
                for s = 1:7
                    m = mode(X1(lab1==s , :));
                    out(d-1).dis(i , s) = pdist([Y;m], distance);%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                end
                
            end
            seqs =  unique(lab);
            for l =  1:length(seqs)
                for l1 =  1:length(seqs)
                    id = lab == l;
                    out(d-1).D_IPI(l,l1) = nanmean(out(d-1).dis(id , l1));
                end
            end
            
            for s = 1:7
                for s2 = 1:6
                    out(d-1).act_est_dist(s,s2) = nanmean(nanmean(pdist2(CBD{s} , chbndry(s2 , :) , distance)));
                end
            end
            
            out(d-1).offdiag_IPI = nanmean(nanmean(out(d-1).D_IPI(2:end ,2:end)-diag(NaN*ones(length(out(d-1).D_IPI(2:end ,2:end)),1))));
            out(d-1).offdiag_IPI_std  =nanstd(nanstd(out(d-1).D_IPI(2:end ,2:end)-diag(NaN*ones(length(out(d-1).D_IPI(2:end ,2:end)),1))));
            
            out(d-1).diag_IPI = mean(diag(out(d-1).D_IPI(2:end ,2:end)));
            out(d-1).diag_IPI_std  =std(diag(out(d-1).D_IPI(2:end ,2:end)));
            
            
            
            figure('color' , [1 1 1])
            
            imagesc(out(d-1).D_IPI);
            title('Crossvalidated Dissimilarity Matrix for Estimated Chunking Structures' , 'FontSize' , 20)
            hold on
            ax = gca;
            ax.XTick = [1:7];
            ax.YTick = [1:7];
            ax.XTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
            ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
            ax.XTickLabelRotation = 45;
            ax.FontSize = 20;
            hold on
            line([1.5 7.5] , [1.5 1.5] ,  'LineWidth' , 3 , 'color' , [0 0 0])
            line([1.5 1.5] ,[1.5 7.5] ,  'LineWidth' , 3 , 'color' , [0 0 0])
            axis square
            colorbar
            
            figure('color' , [1 1 1])
            
            imagesc(out(d-1).act_est_dist);
            title('Estimated vs. Expected Chunking Structure Dissimilarity Matrix'  , 'FontSize' , 20)
            hold on
            ax = gca;
            line([.5 7.5] , [1.5 1.5] ,  'LineWidth' , 3 , 'color' , [0 0 0])
            ax.XTick = [1:6];
            ax.YTick = [1:7];
            ax.XTickLabel = {'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
            ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
            xlabel('Expected')
            ylabel('Estimated')
            ax.XTickLabelRotation = 45;
            ax.FontSize = 20;
            axis square
            colorbar
            diag_IPI = [diag_IPI out(d-1).diag_IPI];
            diag_IPI_s = [diag_IPI_s out(d-1).diag_IPI_std];
            
            offdiag_IPI = [offdiag_IPI out(d-1).offdiag_IPI];
            offdiag_IPI_s = [offdiag_IPI_s out(d-1).offdiag_IPI_std];
        end
        
        
        
        
        figure('color' , 'white')
        errorbar(diag_IPI , diag_IPI_s , 'LineWidth' , 3 )
        grid on
        hold on
        errorbar(offdiag_IPI , offdiag_IPI_s , 'LineWidth' , 3 )
        title('Average Dissimilarity Between and Within Classes in City Block Chunk boundary distance')
        legend({'Within Class' , 'Between Class'})
        hold on
        ax = gca;
        ax.XTick = [1:3];
        ax.XTickLabel = {'Day 1' , 'Day 2' , 'Day 3', 'All Days'};
        ax.FontSize = 20;
        ylabel('Average distance');
 
    case 'saccs_all'
        h1 = figure;
        for d = 2:5
            ANA = getrow(Dall ,ismember(Dall.seqNumb , [0:6]) & ismember(Dall.SN , subjnum) & Dall.isgood & ~Dall.isError & ismember(Dall.Rep , rep) & ismember(Dall.Day , days{d}));
            ANA.seqNumb(ANA.seqNumb > 0) = 1;
            %%%%%%%%%%%%%%% Saccade rate
            [xcoord,PLOT_sp(d-1,:),ERROR_sp(d-1 , :)]=lineplot(ANA.seqNumb,ANA.SaccPerSec,'plotfcn','nanmean');
            hold on
            out.SP(d-1) = anovaMixed(ANA.SaccPerSec , ANA.SN ,'within',[ANA.seqNumb],{'Seq type'},'subset' , ~isnan(ANA.SaccPerSec) , 'intercept',1);
            
            %%%%%%%%%%%%%%% Number of Saccades
            [xcoord,PLOT_ns(d-1,:),ERROR_ns(d-1 , :)]=lineplot(ANA.seqNumb,ANA.NumSaccs,'plotfcn','nanmean');
            out.NS(d-1) = anovaMixed(ANA.NumSaccs , ANA.SN ,'within',[ANA.seqNumb],{'Seq type'},'subset' , ~isnan(ANA.NumSaccs) , 'intercept',1);
            
            %%%%%%%%%%%%%%% Fixation Durarion
            FD = [];
            seq = [];
            subs = [];
            for tn = 1:length(ANA.TN)
                FD  = [FD ; ANA.EyeFixDuration{tn}];
                seq = [seq ; ANA.seqNumb(tn) *ones(length(ANA.EyeFixDuration{tn}) ,1)];
                subs = [subs ; ANA.SN(tn)*ones(length(ANA.EyeFixDuration{tn}) ,1)];
            end
            out.FD(d-1) = anovaMixed(FD , subs ,'within',[seq],{'Seq type'},'subset' , ~isnan(FD) , 'intercept',1);
            [xcoord,PLOT_fd(d-1,:),ERROR_fd(d-1 , :)]=lineplot(seq,FD,'plotfcn','nanmean');

            
            %%%%%%%%%%%%%%% Saccade Durarion
            SD = [];
            seq = [];
            subs = [];
            for tn = 1:length(ANA.TN)
                SD  = [SD ; ANA.SaccDuration{tn}];
                seq = [seq ; ANA.seqNumb(tn) *ones(length(ANA.SaccDuration{tn}) ,1)];
                subs = [subs ; ANA.SN(tn)*ones(length(ANA.SaccDuration{tn}) ,1)];
            end
            out.SD(d-1) = anovaMixed(SD , subs ,'within',[seq],{'Seq type'},'subset' , ~isnan(SD) , 'intercept',1);
            [xcoord,PLOT_sd(d-1,:),ERROR_sd(d-1 , :)]=lineplot(seq,SD,'plotfcn','nanmean');
        
            %%%%%%%%%%%%%%% Saccade peak velocity
            for tn = 1:length(ANA.TN)
                SPV(tn,1) = nanmean(ANA.SaccPeakVel{tn});
                subs1(tn ,1) = ANA.SN(tn);
                seq1(tn ,1) = ANA.seqNumb(tn);
            end
            out.SPV(d-1) = anovaMixed(SPV , subs1 ,'within',[seq1],{'Seq type'},'subset' , ~isnan(SPV) , 'intercept',1);
            [xcoord,PLOT_spv(d-1,:),ERROR_spv(d-1 , :)]=lineplot(seq1,SPV,'plotfcn','nanmean');

            
            %%%%%%%%%%%%%%% Saccade Amplitude
            for tn = 1:length(ANA.TN)
                SA(tn,1) = nanmean(ANA.SaccAmplitude{tn});
            end
            out.SA(d-1) = anovaMixed(SA , subs1 ,'within',[seq1],{'Seq type'},'subset' , ~isnan(SA) , 'intercept',1);
            [xcoord,PLOT_sa(d-1,:),ERROR_sa(d-1 , :)]=lineplot(seq1,SA,'plotfcn','nanmean');

            
            %%%%%%%%%%%%%%% Preview effect
            for tn = 1:length(ANA.TN)
                perv_Ben(tn , :)        = nanmean([1:14] - ANA.EyePressTimePos(tn , :));
            end
            out.PB(d-1) = anovaMixed(perv_Ben , subs1 ,'within',[seq1],{'Seq type'},'subset' , ~isnan(perv_Ben) , 'intercept',1);
            [xcoord,PLOT_pb(d-1,:),ERROR_pb(d-1 , :)]=lineplot(seq1,perv_Ben,'plotfcn','nanmean');
            
            pb(d-1,1)   = nanmean(perv_Ben(ismember(ANA.seqNumb , 0)));
            pb(d-1,2)   = nanmean(perv_Ben(ismember(ANA.seqNumb , [1:6])));
            pb_s(d-1,1) = nanstd(perv_Ben(ismember(ANA.seqNumb , 0)));
            pb_s(d-1,2) = nanstd(perv_Ben(ismember(ANA.seqNumb , [1:6])));
        end
        close(h1)

        figure('color' , 'white')
        subplot(2,1,1)
        errorbar(PLOT_sp(:,2) , ERROR_sp(:,2), 'LineWidth' , 3);
        hold on
        errorbar(PLOT_sp(:,1) , ERROR_sp(:,1), 'LineWidth' , 3);
        grid on
        title(['Saccade rate - Allp = ',num2str(out.SP(4).eff(2).p)])
        ylabel('Saccades per Second')
        hold on
        ax = gca;
        ax.XTick = [1:3];
        ax.XTickLabel = {'Day 1' , 'Day2' , 'Day 3','AllDays'};
        legend({'Chunked' , 'Random'});
        disp(['Day 1 SaccPerSec p-val = ' , num2str(out.SP(1).eff(2).p)])
        disp(['Day 2 SaccPerSec p-val = ' , num2str(out.SP(2).eff(2).p)])
        disp(['Day 3 SaccPerSec p-val = ' , num2str(out.SP(3).eff(2).p)])
        
        
        
        subplot(2,1,2)
        errorbar(PLOT_ns(:,2) , ERROR_ns(:,2), 'LineWidth' , 3);
        hold on
        errorbar(PLOT_ns(:,1) , ERROR_ns(:,1), 'LineWidth' , 3);
        grid on
        title(['Number of Saccades per trial - Allp = ' , num2str(out.NS(4).eff(2).p)])
        ylabel('Saccades')
        hold on
        ax = gca;
        ax.XTick = [1:3];
        ax.XTickLabel = {'Day 1' , 'Day2' , 'Day 3','AllDays'};
        legend({'Chunked' , 'Random'});
        disp(['Day 1 NumSaccades p-val = ' , num2str(out.NS(1).eff(2).p)])
        disp(['Day 2 NumSaccades p-val = ' , num2str(out.NS(2).eff(2).p)])
        disp(['Day 3 NumSaccades p-val = ' , num2str(out.NS(3).eff(2).p)])
        
        
        
        
        
        figure('color' , 'white')
        errorbar(PLOT_fd(:,2) , ERROR_fd(:,2), 'LineWidth' , 3);
        hold on
        errorbar(PLOT_fd(:,1) , ERROR_fd(:,1), 'LineWidth' , 3);
        grid on
        ylabel('msec')
        title(['Average Fixation Duration - Allp = ' ,num2str(out.FD(4).eff(2).p)])
        hold on
        ax = gca;
        ax.XTick = [1:3];
        ax.XTickLabel = {'Day 1' , 'Day2' , 'Day 3','AllDays'};
        legend({'Chunked' , 'Random'});
        disp(['Day 1 FixDuration p-val = ' , num2str(out.FD(1).eff(2).p)])
        disp(['Day 2 FixDuration p-val = ' , num2str(out.FD(2).eff(2).p)])
        disp(['Day 3 FixDuration p-val = ' , num2str(out.FD(3).eff(2).p)])
        
        
        figure('color' , 'white')
        errorbar(PLOT_sd(:,2) , ERROR_sd(:,2), 'LineWidth' , 3);
        hold on
        errorbar(PLOT_sd(:,1) , ERROR_sd(:,1), 'LineWidth' , 3);
        grid on
        ylabel('msec')
        title(['Average Saccade Duration - Allp = ',num2str(out.SD(4).eff(2).p)])
        hold on
        ax = gca;
        ax.XTick = [1:3];
        ax.XTickLabel = {'Day 1' , 'Day2' , 'Day 3','AllDays'};
        legend({'Chunked' , 'Random'});
        disp(['Day 1 SacDuration p-val = ' , num2str(out.SD(1).eff(2).p)])
        disp(['Day 2 SacDuration p-val = ' , num2str(out.SD(2).eff(2).p)])
        disp(['Day 3 SacDuration p-val = ' , num2str(out.SD(3).eff(2).p)])
        
        
        figure('color' , 'white')
        subplot(2,1,1)
        errorbar(PLOT_sa(:,2) , ERROR_sa(:,2), 'LineWidth' , 3);
        hold on
        errorbar(PLOT_sa(:,1) , ERROR_sa(:,1), 'LineWidth' , 3);
        grid on
        ylabel('deg')
        title(['Average Saccade Amplitude - Allp = ' , num2str(out.SA(4).eff(2).p)])
        hold on
        ax = gca;
        ax.XTick = [1:3];
        ax.XTickLabel = {'Day 1' , 'Day2' , 'Day 3','AllDays'};
        legend({'Chunked' , 'Random'});
        disp(['Day 1 SacAmplitude p-val = ' , num2str(out.SA(1).eff(2).p)])
        disp(['Day 2 SacAmplitude p-val = ' , num2str(out.SA(2).eff(2).p)])
        disp(['Day 3 SacAmplitude p-val = ' , num2str(out.SA(3).eff(2).p)])
        
        
        subplot(2,1,2)
        errorbar(PLOT_spv(:,2) , ERROR_spv(:,2), 'LineWidth' , 3);
        hold on
        errorbar(PLOT_spv(:,1) , ERROR_spv(:,1), 'LineWidth' , 3);
        grid on
        ylabel('deg/sec')
        title(['Average Saccade Peak Velocity - Allp = ' , num2str(out.SPV(4).eff(2).p)])
        hold on
        ax = gca;
        ax.XTick = [1:3];
        ax.XTickLabel = {'Day 1' , 'Day2' , 'Day 3','AllDays'};
        legend({'Chunked' , 'Random'});
        disp(['Day 1 SacPeakVelocity p-val = ' , num2str(out.SPV(1).eff(2).p)])
        disp(['Day 2 SacPeakVelocity p-val = ' , num2str(out.SPV(2).eff(2).p)])
        disp(['Day 3 SacPeakVelocity p-val = ' , num2str(out.SPV(3).eff(2).p)])
        
        
        figure('color' , 'white')
        errorbar(-PLOT_pb(:,2) , ERROR_pb(:,2), 'LineWidth' , 3);
        hold on
        errorbar(-PLOT_pb(:,1) , ERROR_pb(:,1), 'LineWidth' , 3);
        grid on
        ylabel('Digits')
        title(['Average Preview Benefit in - Allp = ' , num2str(out.PB(4).eff(2).p)])
        hold on
        ax = gca;
        ax.XTick = [1:3];
        ax.XTickLabel = {'Day 1' , 'Day2' , 'Day 3','AllDays'};
        legend({'Chunked' , 'Random'});
        disp(['Day 1 PervBenefit p-val = ' , num2str(out.PB(1).eff(2).p)])
        disp(['Day 2 PervBenefit p-val = ' , num2str(out.PB(2).eff(2).p)])
        disp(['Day 3 PervBenefit p-val = ' , num2str(out.PB(3).eff(2).p)])
        

    case 'saccs_Singlesubj'
        SNu = [1 3 4 5 6 7 8 9 10];
        for subjnum = 1:length(SNu)
            for d = 2:4
                ANA = getrow(Dall ,ismember(Dall.seqNumb , [0:6]) & ismember(Dall.SN , SNu(subjnum)) & Dall.isgood & ~Dall.isError & ismember(Dall.Rep , rep) & ismember(Dall.Day , d));
                ANA.seqNumb(ANA.seqNumb > 0) = 1;
                sp(subjnum , d-1,1)   = nanmean(ANA.SaccPerSec(ismember(ANA.seqNumb , 0)));
                sp(subjnum , d-1,2)   = nanmean(ANA.SaccPerSec(ismember(ANA.seqNumb , [1:6])));
                sp_s(subjnum , d-1,1) = nanstd(ANA.SaccPerSec(ismember(ANA.seqNumb , 0)));
                sp_s(subjnum , d-1,2) = nanstd(ANA.SaccPerSec(ismember(ANA.seqNumb , [1:6])));
                
                [~ , out.SP(subjnum , d-1)] = ttest2(ANA.SaccPerSec(ANA.seqNumb == 1) , ANA.SaccPerSec(ANA.seqNumb == 0));
                [subjnum , d]
                %             out.SP(d-1) = anovaMixed(ANA.SaccPerSec , ANA.SN ,'within',[ANA.seqNumb],{'Seq type'},'subset' , ~isnan(ANA.SaccPerSec) , 'intercept',1);
                
                
                ns(subjnum , d-1,1)   = nanmean(ANA.NumSaccs(ismember(ANA.seqNumb , 0)));
                ns(subjnum , d-1,2)   = nanmean(ANA.NumSaccs(ismember(ANA.seqNumb , [1:6])));
                ns_s(subjnum , d-1,1) = nanstd(ANA.NumSaccs(ismember(ANA.seqNumb , 0)));
                ns_s(subjnum , d-1,2) = nanstd(ANA.NumSaccs(ismember(ANA.seqNumb , [1:6])));
                
                [~ , out.NS(subjnum , d-1)] = ttest2(ANA.NumSaccs(ANA.seqNumb == 1) , ANA.NumSaccs(ANA.seqNumb == 0));
                %             out.NS(d-1) = anovaMixed(ANA.NumSaccs , ANA.SN ,'within',[ANA.seqNumb],{'Seq type'},'subset' , ~isnan(ANA.NumSaccs) , 'intercept',1);
                
                FD = [];
                seq = [];
                subs = [];
                for tn = 1:length(ANA.TN)
                    FD  = [FD ; ANA.EyeFixDuration{tn}];
                    seq = [seq ; ANA.seqNumb(tn) *ones(length(ANA.EyeFixDuration{tn}) ,1)];
                    subs = [subs ; ANA.SN(tn)*ones(length(ANA.EyeFixDuration{tn}) ,1)];
                end
                [~ , out.FD(subjnum , d-1)] = ttest2(FD(seq == 1) , FD(seq == 0));
                %             out.FD(d-1) = anovaMixed(FD , subs ,'within',[seq],{'Seq type'},'subset' , ~isnan(FD) , 'intercept',1);
                
                
                fd(subjnum , d-1,1)   = nanmean(FD(ismember(seq , 0)));
                fd(subjnum , d-1,2)   = nanmean(FD(ismember(seq , [1:6])));
                fd_s(subjnum , d-1,1) = nanstd(FD(ismember(seq , 0)));
                fd_s(subjnum , d-1,2) = nanstd(FD(ismember(seq , [1:6])));
                SD = [];
                seq = [];
                subs = [];
                for tn = 1:length(ANA.TN)
                    SD  = [SD ; ANA.SaccDuration{tn}];
                    seq = [seq ; ANA.seqNumb(tn) *ones(length(ANA.SaccDuration{tn}) ,1)];
                    subs = [subs ; ANA.SN(tn)*ones(length(ANA.SaccDuration{tn}) ,1)];
                end
                %             out.SD(d-1) = anovaMixed(SD , subs ,'within',[seq],{'Seq type'},'subset' , ~isnan(SD) , 'intercept',1);
                [~ , out.SD(subjnum , d-1)] = ttest2(SD(seq == 1) , SD(seq == 0));
                
                sd(subjnum , d-1,1)   = nanmean(SD(ismember(seq , 0)));
                sd(subjnum , d-1,2)   = nanmean(SD(ismember(seq , [1:6])));
                sd_s(subjnum , d-1,1) = nanstd(SD(ismember(seq , 0)));
                sd_s(subjnum , d-1,2) = nanstd(SD(ismember(seq , [1:6])));
                
                for tn = 1:length(ANA.TN)
                    SPV(tn,1) = nanmean(ANA.SaccPeakVel{tn});
                    subs1(tn ,1) = ANA.SN(tn);
                    seq1(tn ,1) = ANA.seqNumb(tn);
                end
                [~ , out.SPV(subjnum , d-1)] = ttest2(SPV(seq1 == 1) , SPV(seq1 == 0));
                %             out.SPV(d-1) = anovaMixed(SPV , subs1 ,'within',[seq1],{'Seq type'},'subset' , ~isnan(SPV) , 'intercept',1);
                spv(subjnum , d-1,1)   = nanmean(SPV(ismember(ANA.seqNumb , 0)));
                spv(subjnum , d-1,2)   = nanmean(SPV(ismember(ANA.seqNumb , [1:6])));
                spv_s(subjnum , d-1,1) = nanstd(SPV(ismember(ANA.seqNumb , 0)));
                spv_s(subjnum , d-1,2) = nanstd(SPV(ismember(ANA.seqNumb , [1:6])));
                
                
                for tn = 1:length(ANA.TN)
                    SA(tn,1) = nanmean(ANA.SaccAmplitude{tn});
                end
                [~ , out.SA(subjnum , d-1)] = ttest2(SA(seq1 == 1) , SA(seq1 == 0));
                %             out.SA(d-1) = anovaMixed(SA , subs1 ,'within',[seq1],{'Seq type'},'subset' , ~isnan(SA) , 'intercept',1);
                sa(subjnum , d-1,1)   = nanmean(SA(ismember(ANA.seqNumb , 0)));
                sa(subjnum , d-1,2)   = nanmean(SA(ismember(ANA.seqNumb , [1:6])));
                sa_s(subjnum , d-1,1) = nanstd(SA(ismember(ANA.seqNumb , 0)));
                sa_s(subjnum , d-1,2) = nanstd(SA(ismember(ANA.seqNumb , [1:6])));
                
                for tn = 1:length(ANA.TN)
                    perv_Ben(tn , :)        = nanmean([1:14] - ANA.EyePressTimePos(tn , :));
                end
                %             out.PB(d-1) = anovaMixed(perv_Ben , subs1 ,'within',[seq1],{'Seq type'},'subset' , ~isnan(perv_Ben) , 'intercept',1);
                [~ , out.PB(subjnum , d-1)] = ttest2(perv_Ben(seq1 == 1) , perv_Ben(seq1 == 0));
                pb(subjnum , d-1,1)   = nanmean(perv_Ben(ismember(ANA.seqNumb , 0)));
                pb(subjnum , d-1,2)   = nanmean(perv_Ben(ismember(ANA.seqNumb , [1:6])));
                pb_s(subjnum , d-1,1) = nanstd(perv_Ben(ismember(ANA.seqNumb , 0)));
                pb_s(subjnum , d-1,2) = nanstd(perv_Ben(ismember(ANA.seqNumb , [1:6])));
            end
        end
        figure('color' , 'white')
        for subjnum = 1:length(SNu)
            subplot(3,3,subjnum)
            errorbar(sp(subjnum,:,2) , sp_s(subjnum,:,2), 'LineWidth' , 3);
            hold on
            errorbar(sp(subjnum,:,1) , sp_s(subjnum,:,1), 'LineWidth' , 3);
            grid on
            title(['Sub ' , num2str(SNu(subjnum)), ' Saccade rate - LDp = ',num2str(out.SP(subjnum , 3))])
            ylabel('Saccades per Second')
            hold on
            ax = gca;
            ax.XTick = [1:3];
            ax.XTickLabel = {'Day 1' , 'Day2' , 'Day 3'};
            legend({'Chunked' , 'Random'});
            disp(['Day 1 SaccPerSec p-val = ' , num2str(out.SP(subjnum , 1))])
            disp(['Day 2 SaccPerSec p-val = ' , num2str(out.SP(subjnum , 2))])
            disp(['Day 3 SaccPerSec p-val = ' , num2str(out.SP(subjnum , 3))])
        end
        
        figure('color' , 'white')
        for subjnum = 1:length(SNu)
            subplot(3,3,subjnum)
            errorbar(ns(subjnum,:,2) , ns_s(subjnum,:,2), 'LineWidth' , 3);
            hold on
            errorbar(ns(subjnum,:,1) , ns_s(subjnum,:,1), 'LineWidth' , 3);
            grid on
            title(['Sub ' , num2str(subjnum), ' Number of Saccades per trial - LDp = ' , num2str(out.NS(subjnum , 3))])
            ylabel('Saccades')
            hold on
            ax = gca;
            ax.XTick = [1:3];
            ax.XTickLabel = {'Day 1' , 'Day2' , 'Day 3'};
            legend({'Chunked' , 'Random'});
            disp(['Day 1 NumSaccades p-val = ' , num2str(out.NS(subjnum , 1))])
            disp(['Day 2 NumSaccades p-val = ' , num2str(out.NS(subjnum , 2))])
            disp(['Day 3 NumSaccades p-val = ' , num2str(out.NS(subjnum , 3))])
        end
        
        
        
        
        figure('color' , 'white')
        for subjnum = 1:length(SNu)
            subplot(3,3,subjnum)
            errorbar(fd(subjnum,:,2) , fd_s(subjnum,:,2), 'LineWidth' , 3);
            hold on
            errorbar(fd(subjnum,:,1) , fd_s(subjnum,:,1), 'LineWidth' , 3);
            grid on
            ylabel('msec')
            title(['Average Fixation Duration - LDp = ' ,num2str(out.FD(subjnum , 3))])
            hold on
            ax = gca;
            ax.XTick = [1:3];
            ax.XTickLabel = {'Day 1' , 'Day2' , 'Day 3'};
            legend({'Chunked' , 'Random'});
            disp(['Day 1 FixDuration p-val = ' , num2str(out.FD(subjnum , 1))])
            disp(['Day 2 FixDuration p-val = ' , num2str(out.FD(subjnum , 2))])
            disp(['Day 3 FixDuration p-val = ' , num2str(out.FD(subjnum , 3))])
        end

        
        figure('color' , 'white')
        for subjnum = 1:length(SNu)
            subplot(3,3,subjnum)
            errorbar(sd(subjnum,:,2) , sd_s(subjnum,:,2), 'LineWidth' , 3);
            hold on
            errorbar(sd(subjnum,:,1) , sd_s(subjnum,:,1), 'LineWidth' , 3);
            grid on
            ylabel('msec')
            title(['Average Saccade Duration - LDp = ',num2str(out.SD(subjnum , 3))])
            hold on
            ax = gca;
            ax.XTick = [1:3];
            ax.XTickLabel = {'Day 1' , 'Day2' , 'Day 3'};
            legend({'Chunked' , 'Random'});
            disp(['Day 1 SacDuration p-val = ' , num2str(out.SD(subjnum , 1))])
            disp(['Day 2 SacDuration p-val = ' , num2str(out.SD(subjnum , 2))])
            disp(['Day 3 SacDuration p-val = ' , num2str(out.SD(subjnum , 3))])
        end
        
        
        figure('color' , 'white')
        for subjnum = 1:length(SNu)
            subplot(3,3,subjnum)
            errorbar(sa(subjnum,:,2) , sa_s(subjnum,:,2), 'LineWidth' , 3);
            hold on
            errorbar(sa(subjnum,:,1) , sa_s(subjnum,:,1), 'LineWidth' , 3);
            grid on
            ylabel('deg')
            title(['Average Saccade Amplitude - LDp = ' , num2str(out.SA(subjnum , 3))])
            hold on
            ax = gca;
            ax.XTick = [1:3];
            ax.XTickLabel = {'Day 1' , 'Day2' , 'Day 3'};
            legend({'Chunked' , 'Random'});
            disp(['Day 1 SacAmplitude p-val = ' , num2str(out.SA(subjnum , 1))])
            disp(['Day 2 SacAmplitude p-val = ' , num2str(out.SA(subjnum , 2))])
            disp(['Day 3 SacAmplitude p-val = ' , num2str(out.SA(subjnum , 3))])
        end
        
        
        figure('color' , 'white')
        for subjnum = 1:length(SNu)
            subplot(3,3,subjnum)
            errorbar(spv(subjnum,:,2) , spv_s(subjnum,:,2), 'LineWidth' , 3);
            hold on
            errorbar(spv(subjnum,:,1) , spv_s(subjnum,:,1), 'LineWidth' , 3);
            grid on
            ylabel('deg/sec')
            title(['Average Saccade Peak Velocity - LDp = ' , num2str(out.SPV(subjnum , 3))])
            hold on
            ax = gca;
            ax.XTick = [1:3];
            ax.XTickLabel = {'Day 1' , 'Day2' , 'Day 3'};
            legend({'Chunked' , 'Random'});
            disp(['Day 1 SacPeakVelocity p-val = ' , num2str(out.SPV(subjnum , 1))])
            disp(['Day 2 SacPeakVelocity p-val = ' , num2str(out.SPV(subjnum , 2))])
            disp(['Day 3 SacPeakVelocity p-val = ' , num2str(out.SPV(subjnum , 3))])
        end
        
        
        figure('color' , 'white')
        for subjnum = 1:length(SNu)
            subplot(3,3,subjnum)
            errorbar(-pb(subjnum,:,2) , pb_s(subjnum,:,2), 'LineWidth' , 3);
            hold on
            errorbar(-pb(subjnum,:,1) , pb_s(subjnum,:,1), 'LineWidth' , 3);
            grid on
            ylabel('Digits')
            title(['Average Preview Benefit in - LDp = ' , num2str(out.PB(subjnum , 3))])
            hold on
            ax = gca;
            ax.XTick = [1:3];
            ax.XTickLabel = {'Day 1' , 'Day2' , 'Day 3'};
            legend({'Chunked' , 'Random'});
            disp(['Day 1 PervBenefit p-val = ' , num2str(out.PB(subjnum , 1))])
            disp(['Day 2 PervBenefit p-val = ' , num2str(out.PB(subjnum , 2))])
            disp(['Day 3 PervBenefit p-val = ' , num2str(out.PB(subjnum , 3))])
        end

        
    case 'sacc_Chunk_all'
        h1  = figure;
        for d = 2:5
            ANA = getrow(Dall ,ismember(Dall.seqNumb , [1:6]) & ismember(Dall.SN , subjnum) & Dall.isgood & ~Dall.isError & ismember(Dall.Rep , rep) & ismember(Dall.Day , days{d}));
            DigFix = [];
            PervBen = [];
            isSacc = [];
            for tn  = 1:length(ANA.TN)
                ANA.ChunkBndry(tn , :) = [1 diff(ANA.ChnkArrang(tn,:))];
                a = find(ANA.ChunkBndry(tn , :));
                ANA.ChunkBndry(tn , a(2:end)-1) = 2;
                ANA.ChunkBndry(tn , 1:3) = [-1 -1 -1];  % dont account for the first and last sseqeuce presses
                ANA.ChunkBndry(tn , end-2:end) = [-1 -1 -1];% dont account for the first and last sseqeuce presses
                ANA.DigFixDur(tn , :) = zeros(1 ,14);
                
                for p = 1:14
                    ANA.DigFixDur(tn , p) = 2*nansum(ANA.EyeFixDigit{tn} == p);
                end
                DigFix = [DigFix ;[ANA.DigFixDur(tn , ANA.ChunkBndry(tn , :) == 1)' ...
                    (d-1)*ones(size(ANA.DigFixDur(tn , ANA.ChunkBndry(tn , :) == 1)'))...
                    0*ones(size(ANA.DigFixDur(tn , ANA.ChunkBndry(tn , :) == 1)'))...
                    ANA.SN(tn)*ones(size(ANA.DigFixDur(tn , ANA.ChunkBndry(tn , :) == 1)'))...
                    ANA.seqNumb(tn)*ones(size(ANA.DigFixDur(tn , ANA.ChunkBndry(tn , :) == 1)'))]];
                
                DigFix = [DigFix ;[ANA.DigFixDur(tn , ANA.ChunkBndry(tn , :) == 0)' ...
                    (d-1)*ones(size(ANA.DigFixDur(tn , ANA.ChunkBndry(tn , :) == 0)'))...
                    ones(size(ANA.DigFixDur(tn , ANA.ChunkBndry(tn , :) == 0)'))...
                    ANA.SN(tn)*ones(size(ANA.DigFixDur(tn , ANA.ChunkBndry(tn , :) == 0)'))...
                    ANA.seqNumb(tn)*ones(size(ANA.DigFixDur(tn , ANA.ChunkBndry(tn , :) == 0)'))]];
                
                DigFix = [DigFix ;[ANA.DigFixDur(tn , ANA.ChunkBndry(tn , :) == 2)' ...
                    (d-1)*ones(size(ANA.DigFixDur(tn , ANA.ChunkBndry(tn , :) == 2)'))...
                    2*ones(size(ANA.DigFixDur(tn , ANA.ChunkBndry(tn , :) == 2)'))...
                    ANA.SN(tn)*ones(size(ANA.DigFixDur(tn , ANA.ChunkBndry(tn , :) == 2)'))...
                    ANA.seqNumb(tn)*ones(size(ANA.DigFixDur(tn , ANA.ChunkBndry(tn , :) == 2)'))]];
                
                
                perv_Ben        = [1:14] - ANA.EyePressTimePos(tn , :);
                
                PervBen = [PervBen ;[perv_Ben(ANA.ChunkBndry(tn , :) == 1)' ...
                    (d-1)*ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 1)'))...
                    0*ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 1)'))...
                    ANA.SN(tn)*ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 1)'))...
                    ANA.seqNumb(tn)*ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 1)'))]];
                
                PervBen = [PervBen ;[perv_Ben(ANA.ChunkBndry(tn , :) == 0)' ...
                    (d-1)*ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 0)'))...
                    ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 0)'))...
                    ANA.SN(tn)*ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 0)'))...
                    ANA.seqNumb(tn)*ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 0)'))]];
                
                PervBen = [PervBen ;[perv_Ben(ANA.ChunkBndry(tn , :) == 2)' ...
                    (d-1)*ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 2)'))...
                    2*ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 2)'))...
                    ANA.SN(tn)*ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 2)'))...
                    ANA.seqNumb(tn)*ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 2)'))]];
                
                
                isSacc = [isSacc ;[sum(ANA.isSaccWhilePress(tn , ANA.ChunkBndry(tn , :) == 0))/sum(ANA.ChunkBndry(tn , :) == 1) ...
                    d-1 0 ANA.SN(tn)  ANA.seqNumb(tn)]];
                
                isSacc = [isSacc ;[sum(ANA.isSaccWhilePress(tn , ANA.ChunkBndry(tn , :) == 1))/sum(ANA.ChunkBndry(tn , :) == 0) ...
                    d-1 1 ANA.SN(tn)  ANA.seqNumb(tn)]];
                
                isSacc = [isSacc ;[sum(ANA.isSaccWhilePress(tn , ANA.ChunkBndry(tn , :) == 2))/sum(ANA.ChunkBndry(tn , :) == 2) ...
                    d-1 2 ANA.SN(tn)  ANA.seqNumb(tn)]];
                                
            end
            
            
            
            out.FD(d-1) = anovaMixed(DigFix(:,1) , DigFix(:,4) ,'within',[DigFix(:,3)],{'Chunk place'},'subset' , DigFix(:,2) == d-1 & ~isnan(DigFix(:,1)), 'intercept',1);
            out.PB(d-1) = anovaMixed(PervBen(:,1) , PervBen(:,4) ,'within',[PervBen(:,3)],{'Chunk place'},'subset' , PervBen(:,2) == d-1 & ~isnan(PervBen(:,1)), 'intercept',1);
            out.iS(d-1) = anovaMixed(isSacc(:,1) , isSacc(:,4) ,'within',[isSacc(:,3)],{'Chunk place'},'subset' , isSacc(:,2) == d-1 &~isnan(isSacc(:,1)), 'intercept',1);
            
            [xcoord,PLOT_df(d-1,:),ERROR_df(d-1 , :)]=lineplot(DigFix(:,3),DigFix(:,1),'plotfcn','nanmean');
            hold on
            [xcoord,PLOT_pb(d-1,:),ERROR_pb(d-1 , :)]=lineplot(PervBen(:,3),PervBen(:,1),'plotfcn','nanmean');
            [xcoord,PLOT_is(d-1,:),ERROR_is(d-1 , :)]=lineplot(isSacc(:,3),100*isSacc(:,1),'plotfcn','nanmean');
 
        end
        close(h1)
        
        a = tapply(ANA , {'seqNumb' , 'SN' , 'Day'},{'DigFixDur','nanmean','name','meanDigFixDur' , 'subset' , ANA.SN ~=2} , {'ChunkBndry','nanmean','name','ChunkBndry' , 'subset' , ANA.SN ~=2});
        a.meanDigFixDur = a.meanDigFixDur ./ repmat(max(a.meanDigFixDur') , 14,1)';
        figure('color' , 'white')
        for seqnum = 1:length(unique(a.seqNumb))
            subplot(3,2,seqnum)
            imagesc([a.meanDigFixDur(a.seqNumb == seqnum & a.Day == 2 ,4:11) ; NaN*ones(5,8) ;...
                repmat(mean(a.meanDigFixDur(a.seqNumb == seqnum & a.Day == 4 , 4:11))/max(mean(a.meanDigFixDur(a.seqNumb == seqnum & a.Day == 4 , 4:11))) , 5, 1)]);
            title(['Normalized Fixation Duration - Structure ' , num2str(seqnum)])
            hold on 
            ax = gca;
            ax.YTick = [1:11,19];
            ax.YTickLabel = {'Sub 1' 'Sub 2' 'Sub 3' 'Sub 4' 'Sub 5' 'Sub 6' 'Sub 7' 'Sub 8' 'Sub 9' 'Sub 10' 'Sub 11' , 'Sub Avg'};
            ax.XTick = [1:8];
            ax.XTickLabel = mean(a.ChunkBndry(a.seqNumb == seqnum , 4:11));
        end
%         figure('color' , 'white')
%         for subj = 1:length(unique(a.SN))
%             subplot(3,4,subj)
%             imagesc(a.meanDigFixDur(a.SN == subj & a.Day == 4 , 3:12));
%         end
        
        figure('color' , 'white')
        errorbar(PLOT_df(1,:) , ERROR_df(1,:), 'LineWidth' , 3);
        hold on
        errorbar(PLOT_df(2,:) , ERROR_df(2,:), 'LineWidth' , 3);
        errorbar(PLOT_df(3,:) , ERROR_df(3,:), 'LineWidth' , 3);
        errorbar(PLOT_df(4,:) , ERROR_df(4,:), 'LineWidth' , 3);
        grid on
        ylabel('mssec')
        title('Digit Fixation per Chunk Placement')
        hold on
        ax = gca;
        ax.XTick = [1:3];
        legend({'Day 1' , 'Day2' , 'Day 3' , 'All Days'});
        ax.XTickLabel = {'First Chunk Presses' , 'Middle Chunk Presses' , 'Last Chunk Presses'};
        ax.XTickLabelRotation = 45;
        
        disp(['Day 1 FixDuration p-val = ' , num2str(out.FD(1).eff(2).p)])
        disp(['Day 2 FixDuration p-val = ' , num2str(out.FD(2).eff(2).p)])
        disp(['Day 3 FixDuration p-val = ' , num2str(out.FD(3).eff(2).p)])
        disp(['All Day FixDuration p-val = ' , num2str(out.FD(4).eff(2).p)])
        
        
        
        
        figure('color' , 'white')
        errorbar(PLOT_is(1,:) , ERROR_is(1,:), 'LineWidth' , 3);
        hold on
        errorbar(PLOT_is(2,:) , ERROR_is(2,:), 'LineWidth' , 3);
        errorbar(PLOT_is(3,:) , ERROR_is(3,:), 'LineWidth' , 3);
        errorbar(PLOT_is(4,:) , ERROR_is(4,:), 'LineWidth' , 3);
        ylabel('percent')
        title('Percentage of Presses Happend While Making a Saccade')
        hold on
        grid
        ax = gca;
        ax.XTick = [1:3];
        legend({'Day 1' , 'Day2' , 'Day 3' , 'All Days'});
        ax.XTickLabel = {'First Chunk Presses' , 'Middle Chunk Presses' , 'Last Chunk Presses'};
        ax.XTickLabelRotation = 45;
        disp(['Day 1 isSacc p-val = ' , num2str(out.iS(1).eff(2).p)])
        disp(['Day 2 isSacc p-val = ' , num2str(out.iS(2).eff(2).p)])
        disp(['Day 3 isSacc p-val = ' , num2str(out.iS(3).eff(2).p)])
        disp(['All Day isSacc p-val = ' , num2str(out.iS(4).eff(2).p)])
        
        
        
        figure('color' , 'white')
        errorbar(-PLOT_pb(1,:) , ERROR_pb(1,:), 'LineWidth' , 3);
        hold on
        errorbar(-PLOT_pb(2,:) , ERROR_pb(2,:), 'LineWidth' , 3);
        errorbar(-PLOT_pb(3,:) , ERROR_pb(3,:), 'LineWidth' , 3);
        errorbar(-PLOT_pb(4,:) , ERROR_pb(4,:), 'LineWidth' , 3);
        ylabel('Digits')
        hold on
        grid
        ax = gca;
        ax.XTick = [1:3];
        legend({'Day 1' , 'Day2' , 'Day 3' , 'All Days'});
        ax.XTickLabel = {'First Chunk Presses' , 'Middle Chunk Presses' , 'Last Chunk Presses'};
        ax.XTickLabelRotation = 45;
        grid on
        
        disp(['Day 1 PervBenefit p-val = ' , num2str(out.PB(1).eff(2).p)])
        disp(['Day 2 PervBenefit p-val = ' , num2str(out.PB(2).eff(2).p)])
        disp(['Day 3 PervBenefit p-val = ' , num2str(out.PB(3).eff(2).p)])
        disp(['All Day PervBenefit p-val = ' , num2str(out.PB(4).eff(2).p)])
        
        
    case 'sacc_Chunk_single'
        SNu = [1 3 4 5 6 7 8 9 10];
        for subjnum = 1:length(SNu)
            for d = 2:4
                ANA = getrow(Dall ,ismember(Dall.seqNumb , [1:6]) & ismember(Dall.SN , SNu(subjnum)) & Dall.isgood & ~Dall.isError & ismember(Dall.Rep , rep) & ismember(Dall.Day , d));
                iSwithin{subjnum , d-1} = [];
                iSfirst{subjnum , d-1} = [];
                iSlast{subjnum , d-1} = [];
                
                for tn  = 1:length(ANA.TN)
                    ANA.ChunkBndry(tn , :) = [1 diff(ANA.ChnkArrang(tn,:))];
                    a = find(ANA.ChunkBndry(tn , :));
                    ANA.ChunkBndry(tn , a(2:end)-1) = 2;
                    ANA.ChunkBndry(tn , 1) = -1;  % dont account for the first seqeuce presses
                    ANA.ChunkBndry(tn , end) = -1;% dont account for the last seqeuce presses
                    DigFixDur(tn , :) = zeros(1 ,14);
                    for p = 1:14
                        DigFixDur(tn , p) = 2*nansum(ANA.EyeFixDigit{tn} == p);
                    end
                    withinFix{subjnum}(tn,d-1) = nanmean(DigFixDur(tn , ANA.ChunkBndry(tn , :) == 0));
                    firstFix{subjnum}(tn,d-1)  = nanmean(DigFixDur(tn , ANA.ChunkBndry(tn , :) == 1));
                    lastFix{subjnum}(tn,d-1)   = nanmean(DigFixDur(tn , ANA.ChunkBndry(tn , :) == 2));
                    
                    perv_Ben        = [1:14] - ANA.EyePressTimePos(tn , :);
                    PBwithin{subjnum}(tn ,d-1) = nanmean(perv_Ben(ANA.ChunkBndry(tn , :) == 0));
                    PBfirst{subjnum}(tn , d-1) = nanmean(perv_Ben(ANA.ChunkBndry(tn , :) == 1));
                    PBlast{subjnum}(tn , d-1)  = nanmean(perv_Ben(ANA.ChunkBndry(tn , :) == 2));
                    subs{subjnum}(tn , d-1) = ANA.SN(tn);
                    ChunklabelW{subjnum}(tn , d-1) = 0;
                    ChunklabelF{subjnum}(tn , d-1) = 1;
                    ChunklabelL{subjnum}(tn , d-1) = 2;
                    
                    
                    iSwithin{subjnum , d-1} = [iSwithin{subjnum , d-1} ANA.isSaccWhilePress(tn ,ANA.ChunkBndry(tn , :) == 0)];
                    iSfirst{subjnum , d-1}  = [iSfirst{subjnum , d-1}  ANA.isSaccWhilePress(tn ,ANA.ChunkBndry(tn , :) == 1)];
                    iSlast{subjnum , d-1}   = [iSlast{subjnum , d-1}   ANA.isSaccWhilePress(tn ,ANA.ChunkBndry(tn , :) == 2)];
                    
                end
                PB = [PBwithin{subjnum}(: ,d-1) ; PBfirst{subjnum}(:,d-1) ; PBlast{subjnum}(:,d-1)];
                FD = [withinFix{subjnum}(: ,d-1) ; firstFix{subjnum}(:,d-1) ; lastFix{subjnum}(:,d-1)];
                CLabel = [ChunklabelW{subjnum}(: , d-1) ; ChunklabelF{subjnum}(: , d-1) ; ChunklabelL{subjnum}(: , d-1)];
                out.PB(subjnum,d-1) = anovan(PB , CLabel , 'display' , 'off');
                out.FD(subjnum,d-1) = anovan(FD , CLabel , 'display' , 'off');
%                 out.PB(subjnum,d-1) = anovaMixed(PB , repmat(subs,3,1) ,'within',[CLabel],{'Chunk place'},'subset' , ~isnan(PB) , 'intercept',1);
%                 out.FD(subjnum,d-1) = anovaMixed(FD , repmat(subs,3,1),'within',[CLabel],{'Chunk place'},'subset' , ~isnan(FD) ,'intercept',1);
                
                
                wf(subjnum , d-1,1)   = nanmean(withinFix{subjnum}(:,d-1));
                wf_s(subjnum , d-1,1) = nanstd(withinFix{subjnum}(:,d-1));
                
                
                ff(subjnum,d-1,1)   = nanmean(firstFix{subjnum}(:,d-1));
                ff_s(subjnum,d-1,1) = nanstd(firstFix{subjnum}(:,d-1));
                
                lf(subjnum,d-1,1)   = nanmean(lastFix{subjnum}(:,d-1));
                lf_s(subjnum,d-1,1) = nanstd(lastFix{subjnum}(:,d-1));
                
                
                pw(subjnum,d-1,1)   = nanmean(PBwithin{subjnum}(:,d-1));
                pw_s(subjnum,d-1,1) = nanstd(PBwithin{subjnum}(:,d-1));
                
                pf(subjnum,d-1,1)   = nanmean(PBfirst{subjnum}(:,d-1));
                pf_s(subjnum,d-1,1) = nanstd(PBfirst{subjnum}(:,d-1));
                
                pl(subjnum,d-1,1)   = nanmean(PBlast{subjnum}(:,d-1));
                pl_s(subjnum,d-1,1) = nanstd(PBlast{subjnum}(:,d-1));
                
                
                isSacc(subjnum,d-1,:)   = [100*nansum(iSfirst{subjnum , d-1})/length(iSfirst{subjnum , d-1}) 100*nansum(iSwithin{subjnum , d-1})/length(iSwithin{subjnum , d-1}) 100*nansum(iSlast{subjnum , d-1})/length(iSlast{subjnum , d-1})];
                
            end
        end
        
        figure('color' , 'white')
        for subjnum = 1:length(SNu)
            subplot(3,3,subjnum)
            errorbar(wf(subjnum , :) , wf_s(subjnum , :), 'LineWidth' , 3);
            hold on
            errorbar(ff(subjnum , :) , ff_s(subjnum , :), 'LineWidth' , 3);
            errorbar(lf(subjnum , :) , lf_s(subjnum , :), 'LineWidth' , 3);
            grid on
            ylabel('mssec')
            title(['Sub ' , num2str(subjnum) , '- Digit Fixation - LDp = ' , num2str(out.FD(subjnum , end))])
            hold on
            ax = gca;
            ax.XTick = [1:3];
            ax.XTickLabel = {'Day 1' , 'Day2' , 'Day 3'};
            legend({'First Chunk Presses' , 'Middle Chunk Presses' , 'Last Chunk Presses'});
         
        end
            
            
        
        
        figure('color' , 'white')
        for subjnum = 1:length(SNu)
            subplot(3,3,subjnum)
            for d = 1:3
                plot(squeeze(isSacc(subjnum , d,  :)) ,'o-', 'LineWidth' , 3);
                hold on
            end
            grid on
            ylabel('percent')
%             title('Percentage of Presses Happend While Fixation')
            title(['Sub ' , num2str(subjnum) , '- Percentage of Presses while Fixated - LDp = ' , num2str(out.FD(subjnum , end))])
            hold on
            ax = gca;
            ax.XTick = [1:3];
            legend({'Day 1' , 'Day2' , 'Day 3'});
            ax.XTickLabel = {'First Chunk Presses' , 'Middle Chunk Presses' , 'Last Chunk Presses'};
            ax.XTickLabelRotation = 45;
        end
        
        figure('color' , 'white')
        for subjnum = 1:length(SNu)
            subplot(3,3,subjnum)
            errorbar(-pf(subjnum , :) , pf_s(subjnum , :), 'LineWidth' , 3);
            hold on
            errorbar(-pw(subjnum , :) , pw_s(subjnum , :), 'LineWidth' , 3);
            errorbar(-pl(subjnum , :) , pl_s(subjnum , :), 'LineWidth' , 3);
            grid on
            ylabel('Digits')
            title(['Sub ' , num2str(subjnum) , '- Preview Benefi - LDp = ' , num2str(out.FD(subjnum , end))])
            hold on
            ax = gca;
            ax.XTick = [1:3];
            ax.XTickLabel = {'Day 1' , 'Day2' , 'Day 3'};
            legend({'First Chunk Presses' , 'Middle Chunk Presses' , 'Last Chunk Presses'});
        end
    
    case 'Random_error'
        for sub = subjnum
            clear C chunks
            ANAChunk = getrow(Dall ,Dall.seqNumb > 7 & ismember(Dall.SN , sub));
            ANAChunk.AllPress(isnan(ANAChunk.AllPress)) = 0;
            C = unique(ANAChunk.AllPress , 'rows');
            DuoCounter  = 1;
            TrpCounter  = 1;
            QudCounter  = 1;
            for c = 1:size(C , 1)
                chunks{c} = C(c, C(c,:)~=0);
                switch length(chunks{c})
                    case 2
                        doubles(DuoCounter , :) = chunks{c};
                        doubleInfo(DuoCounter , 1) = 2; % comes form double chunk
                        doubleInfo(DuoCounter , 2) = 1;
                        DuoCounter = DuoCounter + 1;
                    case 3
                        for pcntr = 1:2
                            doubles(DuoCounter , :) = chunks{c}(pcntr : pcntr+1);
                            doubleInfo(DuoCounter , 1) = 3; % comes form triple chunk
                            doubleInfo(DuoCounter , 2) = pcntr; % to know if it's the first or second duo in the reiplet
                            DuoCounter = DuoCounter + 1;
                        end
                        triples(TrpCounter , :) = chunks{c};
                        tripletInfo(TrpCounter , 1) = 3;
                        tripletInfo(TrpCounter , 2) = 1;
                        TrpCounter = TrpCounter + 1;
                    case 4
                        for pcntr = 1:3
                            doubles(DuoCounter , :) = chunks{c}(pcntr : pcntr+1);
                            doubleInfo(DuoCounter , 1) = 4; % comes form quadruple chunk
                            doubleInfo(DuoCounter , 2) = pcntr; % to know if it's the first, second or third duo in the reiplet
                            DuoCounter = DuoCounter + 1;
                        end
                        for pcntr = 1:2
                            triples(TrpCounter , :) = chunks{c}(pcntr : pcntr+2);
                            tripletInfo(TrpCounter , 1) = 4;
                            tripletInfo(TrpCounter , 2)  = pcntr;
                            TrpCounter = TrpCounter + 1;
                        end
                        quadrpl(QudCounter , :) = chunks{c};
                        quadrplInfo(QudCounter , 1) = 4;
                        quadrplInfo(QudCounter , 2) = 1;
                        QudCounter = QudCounter + 1;
                end
            end
            
            out.DoubInducedError{sub} = [];
            out.TripInducedError{sub} = [];
            out.QuadInducedError{sub} = [];
            out.NotInducedByChunk(sub) = 0;
            out.InducedByChunk(sub) = 0;
            ANARand = getrow(Dall , Dall.seqNumb == 0 & ismember(Dall.SN , sub) & Dall.isgood & Dall.isError & ismember(Dall.Rep , rep) & ismember(Dall.Day , [2:4]));
            out.SumofError(sub) = sum(sum(ANARand.badPress));
            out.WholeChunkErr(sub , :) = [0 0 0];
            for tn = 1:length(ANARand.TN)
                for CL = 2:4
                    BP = find(ANARand.badPress(tn , CL:end))+(CL-1);
                    if ~isempty(BP)
                        for bp = BP
                            switch CL
                                case 2
                                    if ismember(ANARand.AllPress(tn , bp-(CL-1):bp) , doubles , 'rows')
                                        ind = ismember(doubles , ANARand.AllPress(tn , bp-(CL-1):bp) , 'rows');
                                        out.DoubInducedError{sub} = [out.DoubInducedError{sub} ; doubleInfo(ind, :)];
                                        out.InducedByChunk(sub) = out.InducedByChunk(sub)  +1;
                                    else
                                        out.NotInducedByChunk(sub) = out.NotInducedByChunk(sub)  +1;
                                    end
                                    
                                case 3
                                    if ismember(ANARand.AllPress(tn , bp-(CL-1):bp) , triples , 'rows')
                                        ind = ismember(triples , ANARand.AllPress(tn , bp-(CL-1):bp) , 'rows');
                                        out.TripInducedError{sub} = [out.TripInducedError{sub} ; tripletInfo(ind, :)];
                                    end
                                case 4
                                    if ismember(ANARand.AllPress(tn , bp-(CL-1):bp) , quadrpl , 'rows')
                                        ind = ismember(quadrpl , ANARand.AllPress(tn , bp-(CL-1):bp) , 'rows');
                                        out.QuadInducedError{sub} = [out.QuadInducedError{sub} ; quadrplInfo(ind, :)];
                                    end
                            end
                        end
                    end
                end
            end
            
            if length(out.DoubInducedError{sub})>0
                out.WholeChunkErr(sub,1) = out.WholeChunkErr(sub,1) + sum(out.DoubInducedError{sub}(:,1) == 2);
            end
            if length(out.TripInducedError{sub})>0
                out.WholeChunkErr(sub,2) = out.WholeChunkErr(sub,2) + sum(out.TripInducedError{sub}(:,1) == 3);
            end
            if length(out.QuadInducedError{sub})>0
                out.WholeChunkErr(sub,3) =  out.WholeChunkErr(sub,3) + sum(out.QuadInducedError{sub}(:,1) == 4);
            end
        end
        out.InducedByChunk = 100*(out.InducedByChunk./out.SumofError);
        out.NotInducedByChunk = 100*(out.NotInducedByChunk./out.SumofError);
        figure('color' , 'white')
        errorbar([nanmean(out.InducedByChunk) nanmean(out.NotInducedByChunk)] , [nanstd(out.InducedByChunk) nanstd(out.NotInducedByChunk)] , 'LineWidth' , 3)
        grid
        ylabel('percent')
        title('Errors in Random Sequences Induces by Chunks vs not Induced by Chunks')
        hold on
        ax = gca;
        ax.XTick = [1:2];
        ax.XTickLabel = {'Induced By Chunk' , 'Not Induced By Chunk'};
        ax.XTickLabelRotation = 45;
        
    case 'Points'
        ANA1 = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.seqNumb , [1:4])  &  ismember(Dall.Rep , rep) & ismember(Dall.Group , GroupCode));
        ANA0 = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.seqNumb , [0])   &  ismember(Dall.Rep , rep)& ismember(Dall.Group , GroupCode)); 
        
        figure('color' , 'white')
        subplot(1,2,1)
        histogram(ANA1.points , 'Normalization','probability')
        hold on 
        ax = gca;
        ax.XTick = [0 1 2];
        xlabel('Points')
        ylabel('Probability')
        
        title('Chunked')
        subplot(1,2,2)
        histogram(ANA0.points , 'Normalization','probability')
        title('Random')
        ax = gca;
        ax.XTick = [0 1 2];
        xlabel('Points')
        ylabel('Probability')
        out =[];
        
        
end

