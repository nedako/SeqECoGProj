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
subj_name = {'P2', 'P4' , 'XX'};
if subjnum == length(subj_name)
    subjnum = 1:length(subj_name)-1;
end


load([baseDir , '/CMB.mat'])
%load([baseDir , '/se1_all.mat'])
days = {1 ,2 ,3 ,4 ,[2:4] ,[3:4] [1:4]};
% subjnum = 2;

tid = {[1:250] [251 :750] [751:1000] [1:1000]};
switch what
    case 'IPI_MT_seq'
        
        %% within and between IPI ttest with randomization
        %% IPI per day
        %
        ANA = getrow(Dall , Dall.isgood & ~Dall.isError & ismember(Dall.seqNumb , [0:4]));
        ANA.seqNumb(ANA.seqNumb > 0) = 1;
        BN = unique(ANA.BN);
        for bn = 1:length(BN)
            temp = getrow(ANA , ANA.BN == BN(bn));
            boxplot(reshape(temp.IPI , numel(temp.IPI),1) , reshape(temp.IPIarrangement  , numel(temp.IPIarrangement),1))
            hold on
        end
        ANA.seqType = ANA.seqNumb;
        ANA.seqType(ANA.seqNumb>1) = 1;
        
        ipinum = [1:unique(ANA.seqlength) - 1];
        ANA.IPI = ANA.IPI(: , 1:unique(ANA.seqlength) - 1);
        IPI.ipi = reshape(ANA.IPI , numel(ANA.IPI) , 1);
        IPI.ipi = reshape(ANA.IPI , numel(ANA.IPI) , 1);
        IPI.day = repmat(ANA.Day , unique(ANA.seqlength)-1 , 1);
        IPI.BN = repmat(ANA.BN , unique(ANA.seqlength)-1 , 1);
        IPI.ipiNum = reshape(repmat(ipinum , length(ANA.IPI) , 1) , numel(IPI.ipi) , 1);
        IPI.ipiArr = reshape(ANA.IPIarrangement , numel(ANA.IPIarrangement) , 1);
        IPI.seqN = repmat(ANA.seqNumb , unique(ANA.seqlength)-1 , 1);
        IPI.seqT = repmat(ANA.seqType , unique(ANA.seqlength)-1 , 1);
        
        
        
        IPIs = getrow(IPI , IPI.seqN > 0);
        out.IPI = anovan(IPIs.ipi , IPIs.ipiArr,'varnames',{'Within/Between'} , 'display' , 'on' , 'model' , 'full');
        out.IPI_all = anovan(IPI.ipi , IPI.ipiArr,'varnames',{'Random/Within/Between'} , 'display' , 'on' , 'model' , 'full');

        h1 = figure('color' , 'white');
        [xb,ePLOTb,ERRORb] = lineplot(IPI.day , IPI.ipi , 'subset' , IPI.ipiArr == 1);
        hold on
        [xw,ePLOTw,ERRORw] = lineplot(IPI.day , IPI.ipi , 'subset' , IPI.ipiArr == 2);
        [xr,ePLOTr,ERRORr] = lineplot(IPI.day , IPI.ipi , 'subset' , IPI.ipiArr == 0);
        dayss = unique(IPI.day);
        for d = 1:length(dayss)
            IPIs = getrow(IPI , IPI.seqN > 0 & IPI.day == dayss(d));
            out.IPId{d} = anovan(IPIs.ipi , IPIs.ipiArr,'varnames',{'Within/Between'} , 'display' , 'on' , 'model' , 'full');
            [xs{d},PLOTs{d},ERRORs{d}] = lineplot(IPI.ipiNum,IPI.ipi , 'subset' , IPI.day == dayss(d) & IPI.seqN == 1);
            hold on
            [xra{d},PLOTra{d},ERRORra{d}] = lineplot(IPI.ipiNum,IPI.ipi , 'subset' , IPI.day == dayss(d) & IPI.seqN == 0);
        end
        close(h1);
        
        figure('color' , 'white')
        subplot(length(dayss) , 2 , [1:2:2*length(dayss)])
        hb = plotshade(xb',ePLOTb,ERRORb,'transp' , .2 , 'patchcolor' , 'b' , 'linecolor' , 'b' , 'linewidth' , 3 , 'linestyle' , ':')
        hold on
        hw = plotshade(xw',ePLOTw , ERRORw,'transp' , .2 , 'patchcolor' , 'm' , 'linecolor' , 'm' , 'linewidth' , 3 , 'linestyle' , ':')
        hr = plotshade(xr',ePLOTr , ERRORr,'transp' , .2 , 'patchcolor' , 'r' , 'linecolor' , 'r' , 'linewidth' , 3 , 'linestyle' , ':')
        ax = gca;
        ax.FontSize = 20;
        legend([hb,hw,hr] , {'Between chunk IPIs' , 'Within chunk IPIs' , 'Random IPIs'} , 'FontSize' , 20)
        title(['Within chunk vs Between chunk Chunk IPIs , p(W/B) = ' , num2str(out.IPI)],'FontSize' , 16)
        xlabel('Training Blocks', 'FontSize' , 20)
        ylabel('msec'  ,'FontSize' , 20)
        grid on
        figcount = 2;
        for d = 1:length(dayss)
            subplot(length(dayss) , 2 , figcount)
            hold on
            hw = plotshade(xs{d}',PLOTs{d} , ERRORs{d},'transp' , .2 , 'patchcolor' , 'm' , 'linecolor' , 'm' , 'linewidth' , 3 , 'linestyle' , ':')
            hr = plotshade(xra{d}',PLOTra{d} , ERRORra{d},'transp' , .2 , 'patchcolor' , 'r' , 'linecolor' , 'r' , 'linewidth' , 3 , 'linestyle' , ':')
            ax = gca;
            ax.FontSize = 20;
            legend([hb,hw,hr] , {'Between chunk IPIs' , 'Within chunk IPIs' , 'Random IPIs'} , 'FontSize' , 20)
            title(['Within chunk vs Between chunk Chunk IPIs , p(W/B) = ' , num2str(out.IPId{d})],'FontSize' , 16)
            xlabel('IPI Number', 'FontSize' , 20)
            ylabel('msec'  ,'FontSize' , 20)
            grid on
            figcount = figcount+2;
        end
        
        %         out.MT = anovaMixed([ANA1.MT ; ANA0.MT] , [ANA1.SN ;ANA0.SN],'within',[0*ANA1.MT ; 1+0*ANA0.MT],{'Random/Chunked'},'intercept',1)  ;
        out.MT = anovan(ANA.MT, ANA.seqNumb , 'display' , 'off' , 'varnames' , 'Rand/Chunked');
        h1 = figure('color' , 'white');
        [xs,PLOTs,ERRORs] = lineplot(ANA.Day,ANA.MT , 'subset' , ANA.seqType == 1);
        hold on
        [xr,PLOTr,ERRORr] = lineplot(ANA.Day,ANA.MT, 'subset' , ANA.seqType == 0);
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
        xlabel('Training Days', 'FontSize' , 20)
        %         ax.XTick = [1:4];
        %         ax.XTickLabel = {'Day 0 ' ,'Day 1' , 'Day 2' , 'Day 3'};
        legend([hs ,hr] , {'Chunked' , 'Random'})
        grid on
        
  
        
    case 'IPI_MT_chunk'
         %% within and between IPI ttest with randomization
        %% IPI per day
        %
        ANA = getrow(Dall , Dall.isgood & ~Dall.isError & ismember(Dall.seqNumb , [103 104 203 204 303 304]));
        BN = unique(ANA.BN);
      
        ANA3 = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.seqNumb , [103 203 303])  & ~Dall.isError &  ismember(Dall.Rep , rep) );
        ANA4 = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.seqNumb , [104 204 304])  & ~Dall.isError &  ismember(Dall.Rep , rep));
        
        
        IPInums = 1:6;
        BN = unique(ANA3.BN);
        clear IPIwithin IPIbetween
        ipi3num = [1:unique(ANA3.seqlength) - 1];
        ipi4num = [1:unique(ANA4.seqlength) - 1];
        ANA3.IPI = ANA3.IPI(: , 1:unique(ANA3.seqlength) - 1);
        ANA4.IPI = ANA4.IPI(: , 1:unique(ANA4.seqlength) - 1);
        IPI3.ipi = reshape(ANA3.IPI , numel(ANA3.IPI) , 1);
        IPI4.ipi = reshape(ANA4.IPI , numel(ANA4.IPI) , 1);
        IPI3.day = repmat(ANA3.Day , unique(ANA3.seqlength)-1 , 1);
        IPI4.day = repmat(ANA4.Day , unique(ANA4.seqlength)-1 , 1);
        IPI3.BN = repmat(ANA3.BN , unique(ANA3.seqlength)-1 , 1);
        IPI4.BN = repmat(ANA4.BN , unique(ANA4.seqlength)-1 , 1);
        IPI3.ipiNum = reshape(repmat(ipi3num , length(ANA3.IPI) , 1) , numel(IPI3.ipi) , 1);
        IPI4.ipiNum = reshape(repmat(ipi4num , length(ANA4.IPI) , 1) , numel(IPI4.ipi) , 1);

        
        %         out.IPI = anovaMixed([IPIw(:,1) ; IPIb(:,1)] , [IPIw(:,3) ; IPIb(:,3)],'within',[0*IPIw(:,1) ; 1+0*IPIb(:,1)],{'Within/Between'},'intercept',1)  ;

        out.IPI = anovan([IPI3.ipi ; IPI4.ipi] , [[0*IPI3.ipi; 1+0*IPI4.ipi] , [IPI3.day; IPI4.day]],'varnames',{'Quadruple/Triplet' , 'day'} , 'display' , 'on' , 'model' , 'full');

        
        h1 = figure('color' , 'white');
        hold on
        [x3,ePLOT3,ERROR3] = lineplot( IPI3.day, IPI3.ipi);
        [x4,ePLOT4,ERROR4] = lineplot(IPI4.day , IPI4.ipi);
        
        [x3_2,ePLOT3_2,ERROR3_2] = lineplot( IPI3.ipiNum, IPI3.ipi);
        [x4_2,ePLOT4_2,ERROR4_2] = lineplot(IPI4.ipiNum , IPI4.ipi);
        close(h1);
        
        figure('color' , 'white')
        subplot(121)
        hb = plotshade(x3',ePLOT3,ERROR3,'transp' , .2 , 'patchcolor' , 'b' , 'linecolor' , 'b' , 'linewidth' , 3 , 'linestyle' , ':')
        hold on
        hw = plotshade(x4',ePLOT4 , ERROR4,'transp' , .2 , 'patchcolor' , 'm' , 'linecolor' , 'm' , 'linewidth' , 3 , 'linestyle' , ':')
        ax = gca;
        ax.FontSize = 20;
        legend([hb,hw] , {'Triplet IPIs' , 'Quadruple IPIs'} , 'FontSize' , 20)
        title(['Triplet vs Quadruple Chunk IPIs , p(Q/T) = ' , num2str(out.IPI(2))],'FontSize' , 16)
        xlabel('Training Days', 'FontSize' , 20)
        ylabel('msec'  ,'FontSize' , 20)
        grid on
        
        subplot(122)
        hb = plotshade(x3_2',ePLOT3_2,ERROR3_2,'transp' , .2 , 'patchcolor' , 'b' , 'linecolor' , 'b' , 'linewidth' , 3 , 'linestyle' , ':')
        hold on
        hw = plotshade(x4_2',ePLOT4_2,ERROR4_2,'transp' , .2 , 'patchcolor' , 'm' , 'linecolor' , 'm' , 'linewidth' , 3 , 'linestyle' , ':')
        ax = gca;
        ax.FontSize = 20;
        legend([hb,hw] , {'Triplet IPIs' , 'Quadruple IPIs'} , 'FontSize' , 20)
        title(['Triplet vs Quadruple Chunk IPIs , p(Q/T) = ' , num2str(out.IPI(2))],'FontSize' , 16)
        xlabel('IPI number', 'FontSize' , 20)
        ylabel('msec'  ,'FontSize' , 20)
        grid on
        
        
        %         out.MT = anovaMixed([ANA1.MT ; ANA0.MT] , [ANA1.SN ;ANA0.SN],'within',[0*ANA1.MT ; 1+0*ANA0.MT],{'Random/Chunked'},'intercept',1)  ;
        out.MT = anovan([ANA3.MT ; ANA4.MT] , [0*ANA3.MT ; 1+0*ANA4.MT] , 'display' , 'off' , 'varnames' , 'Rand/Chunked');
        h1 = figure('color' , 'white');
        [x3,PLOT3,ERROR3] = lineplot(ANA3.Day,ANA3.MT);
        hold on
        [x4,PLOT4,ERROR4] = lineplot(ANA4.Day,ANA4.MT);
        close(h1);
        
        figure('color' , 'white')
        h3 = plotshade(x3',PLOT3,ERROR3,'transp' , .2 , 'patchcolor' , 'black' , 'linecolor' , 'k' , 'linewidth' , 3 , 'linestyle' , ':')
        hold on
        h4 = plotshade(x4',PLOT4,ERROR4,'transp' , .2 , 'patchcolor' , 'r' , 'linecolor' , 'r' , 'linewidth' , 3 , 'linestyle' , ':')
        hold on
        ax = gca;
        ax.FontSize = 20;
        title(['Movement Time, p = ' , num2str(out.MT)] ,'FontSize' , 20)
        ylabel('msec'  ,'FontSize' , 20)
        xlabel('Training Days', 'FontSize' , 20)
        %         ax.XTick = [1:4];
        %         ax.XTickLabel = {'Day 0 ' ,'Day 1' , 'Day 2' , 'Day 3'};
        legend([h3 ,h4] , {'Triplet' , 'Quadruple'})
        grid on
    case 'RT_seq'
        
        %% within and between IPI ttest with randomization
        %% IPI per day
        ANA = getrow(Dall , Dall.isgood & ~Dall.isError & ismember(Dall.seqNumb , [0:4]));
        ANA.seqNumb(ANA.seqNumb > 0) = 1;
        BN = unique(ANA.BN);
        for bn = 1:length(BN)
            temp = getrow(ANA , ANA.BN == BN(bn));
            boxplot(reshape(temp.IPI , numel(temp.IPI),1) , reshape(temp.IPIarrangement  , numel(temp.IPIarrangement),1))
            hold on
        end
        
        
        ANA1 = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.seqNumb , [1:4])  & ~Dall.isError &  ismember(Dall.Rep , rep) );
        ANA0 = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.seqNumb , [0])  & ~Dall.isError &  ismember(Dall.Rep , rep));
        
        

        ANA1.RT = ANA1.AllPressTimes(:,1);
        ANA0.RT = ANA0.AllPressTimes(:,1);
       
        
  
        
        %         out.RT = anovaMixed([ANA1.RT ; ANA0.RT] , [ANA1.SN ;ANA0.SN],'within',[0*ANA1.RT ; 1+0*ANA0.RT],{'Random/Chunked'},'intercept',1)  ;
        out.RT = anovan([ANA1.RT ; ANA0.RT] , [0*ANA1.RT ; 1+0*ANA0.RT] , 'display' , 'off' , 'varnames' , 'Rand/Chunked');
        h1 = figure('color' , 'white');
        [xs,PLOTs,ERRORs] = lineplot(ANA1.BN,ANA1.RT);
        hold on
        [xr,PLOTr,ERRORr] = lineplot(ANA0.BN,ANA0.RT);
        close(h1);
        
        figure('color' , 'white')
        hs = plotshade(xs',PLOTs,ERRORs,'transp' , .2 , 'patchcolor' , 'black' , 'linecolor' , 'k' , 'linewidth' , 3 , 'linestyle' , ':')
        hold on
        hr = plotshade(xr',PLOTr,ERRORr,'transp' , .2 , 'patchcolor' , 'r' , 'linecolor' , 'r' , 'linewidth' , 3 , 'linestyle' , ':')
        hold on
        ax = gca;
        ax.FontSize = 20;
        title(['Reaction Time, p = ' , num2str(out.RT)] ,'FontSize' , 20)
        ylabel('msec'  ,'FontSize' , 20)
        xlabel('Training Blocks', 'FontSize' , 20)
        %         ax.XTick = [1:4];
        %         ax.XTickLabel = {'Day 0 ' ,'Day 1' , 'Day 2' , 'Day 3'};
        legend([hs ,hr] , {'Chunked' , 'Random'})
        grid on        
    case 'RT_chunk'
        
        %% within and between IPI ttest with randomization
        %% IPI per day
        ANA = getrow(Dall , Dall.isgood & ~Dall.isError & ismember(Dall.seqNumb , [103 203 104 204]));
        ANA.seqNumb(ANA.seqNumb > 0) = 1;
        BN = unique(ANA.BN);
        for bn = 1:length(BN)
            temp = getrow(ANA , ANA.BN == BN(bn));
            boxplot(reshape(temp.IPI , numel(temp.IPI),1) , reshape(temp.IPIarrangement  , numel(temp.IPIarrangement),1))
            hold on
        end
        
        
        ANA1 = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.seqNumb , [104 204])  & ~Dall.isError &  ismember(Dall.Rep , rep) );
        ANA0 = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.seqNumb , [103 203])  & ~Dall.isError &  ismember(Dall.Rep , rep));
        
        

        ANA1.RT = ANA1.AllPressTimes(:,1);
        ANA0.RT = ANA0.AllPressTimes(:,1);
       
        
  
        
        %         out.RT = anovaMixed([ANA1.RT ; ANA0.RT] , [ANA1.SN ;ANA0.SN],'within',[0*ANA1.RT ; 1+0*ANA0.RT],{'Random/Chunked'},'intercept',1)  ;
        out.RT = anovan([ANA1.RT ; ANA0.RT] , [0*ANA1.RT ; 1+0*ANA0.RT] , 'display' , 'off' , 'varnames' , 'Rand/Chunked');
        h1 = figure('color' , 'white');
        [xs,PLOTs,ERRORs] = lineplot(ANA1.BN,ANA1.RT);
        hold on
        [xr,PLOTr,ERRORr] = lineplot(ANA0.BN,ANA0.RT);
        close(h1);
        
        figure('color' , 'white')
        hs = plotshade(xs',PLOTs,ERRORs,'transp' , .2 , 'patchcolor' , 'black' , 'linecolor' , 'k' , 'linewidth' , 3 , 'linestyle' , ':')
        hold on
        hr = plotshade(xr',PLOTr,ERRORr,'transp' , .2 , 'patchcolor' , 'r' , 'linecolor' , 'r' , 'linewidth' , 3 , 'linestyle' , ':')
        hold on
        ax = gca;
        ax.FontSize = 20;
        title(['Reaction Time, p = ' , num2str(out.RT)] ,'FontSize' , 20)
        ylabel('msec'  ,'FontSize' , 20)
        xlabel('Training Blocks', 'FontSize' , 20)
        %         ax.XTick = [1:4];
        %         ax.XTickLabel = {'Day 0 ' ,'Day 1' , 'Day 2' , 'Day 3'};
        legend([hs ,hr] , {'Quadruples' , 'Triplets'})
        grid on
    case 'RT_SF'
        
        %% within and between IPI ttest with randomization
        %% IPI per day
        ANA = getrow(Dall , Dall.isgood & ~Dall.isError & ismember(Dall.seqNumb , [11 22 33 44 55]));
        
        

        ANA.RT = ANA.AllPressTimes(:,1);
        
        %         out.RT = anovaMixed([ANA1.RT ; ANA0.RT] , [ANA1.SN ;ANA0.SN],'within',[0*ANA1.RT ; 1+0*ANA0.RT],{'Random/Chunked'},'intercept',1)  ;
        out.RT = anovan([ANA.RT] , ANA.seqNumb , 'display' , 'on');
        h1 = figure('color' , 'white');
        [xs,PLOTs,ERRORs] = lineplot(ANA.BN,ANA.RT);
        hold on
        close(h1);
        
        figure('color' , 'white')
        hs = plotshade(xs',PLOTs,ERRORs,'transp' , .2 , 'patchcolor' , 'black' , 'linecolor' , 'k' , 'linewidth' , 3 , 'linestyle' , ':')
        hold on
        
        ax = gca;
        ax.FontSize = 20;
        title(['Reaction Time, p = ' , num2str(out.RT)] ,'FontSize' , 20)
        ylabel('msec'  ,'FontSize' , 20)
        xlabel('Training Blocks', 'FontSize' , 20)
        %         ax.XTick = [1:4];
        %         ax.XTickLabel = {'Day 0 ' ,'Day 1' , 'Day 2' , 'Day 3'};
        grid on
        
    case 'MT_SF'
        
        %% within and between IPI ttest with randomization
        %% IPI per day
        ANA = getrow(Dall , Dall.isgood & ~Dall.isError & ismember(Dall.seqNumb , [11 22 33 44 55]));
        
        

        ANA.MT = ANA.AllPressTimes(:,7) - ANA.AllPressTimes(:,1);
        
        %         out.MT = anovaMixed([ANA1.MT ; ANA0.MT] , [ANA1.SN ;ANA0.SN],'within',[0*ANA1.MT ; 1+0*ANA0.MT],{'Random/Chunked'},'intercept',1)  ;
        out.MT = anovan([ANA.MT] , ANA.seqNumb , 'display' , 'on');
        h1 = figure('color' , 'white');
        [xs,PLOTs,ERRORs] = lineplot(ANA.BN,ANA.MT);
        hold on
        close(h1);
        
        figure('color' , 'white')
        hs = plotshade(xs',PLOTs,ERRORs,'transp' , .2 , 'patchcolor' , 'black' , 'linecolor' , 'k' , 'linewidth' , 3 , 'linestyle' , ':')
        hold on
        
        ax = gca;
        ax.FontSize = 20;
        title(['Movement Time, p = ' , num2str(out.MT)] ,'FontSize' , 20)
        ylabel('msec'  ,'FontSize' , 20)
        xlabel('Training Blocks', 'FontSize' , 20)
        %         ax.XTick = [1:4];
        %         ax.XTickLabel = {'Day 0 ' ,'Day 1' , 'Day 2' , 'Day 3'};
        grid on
end
     