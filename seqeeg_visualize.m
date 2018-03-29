function out  = secog_visualize(Dall , subjnum, what,rep)

%%  Cases
%     case 'IPI_MT_seq'
%     case 'IPI_MT_chunk'
%     case 'RT_seq'
%     case 'RT_chunk'
%     case 'RT_SF'
%     case 'MT_SF'

baseDir = '/Users/nkordjazi/Documents/SeqECoG/analyze';
% baseDir = '/Users/nkordjazi/Documents/SeqEye/se1/SeqEye1/se1_data/analyze';
% subj_name = {'XW' , 'ML' , 'DS' , 'BM' , 'HK' , 'BW' 'XX'};
subj_name = {'P2', 'P4' , 'P5' , 'XX'};
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
        ANA.seqNumb(ANA.seqNumb==2) = 1; % 1 and 2 have the same chunking strucutre
        ANA.seqNumb(ANA.seqNumb==4) = 3; % 3 and 4 have the same chunking strucutre
        
        ANA.seqType = ANA.seqNumb; % Random vs structure
        
        ANA.seqType(ANA.seqNumb > 0) = 1;
        BN = unique(ANA.BN);
%         for bn = 1:length(BN)
%             temp = getrow(ANA , ANA.BN == BN(bn));
%             boxplot(reshape(temp.IPI , numel(temp.IPI),1) , reshape(temp.IPIarrangement  , numel(temp.IPIarrangement),1))
%             hold on
%         end
        ANA.seqType = ANA.seqNumb;
        ANA.seqType(ANA.seqNumb>1) = 1;
        IPI = [];
        ipinum = [1:unique(ANA.seqlength) - 1];
        ANA.IPI = ANA.IPI(: , 1:unique(ANA.seqlength) - 1);
        IPI.ipi = reshape(ANA.IPI , numel(ANA.IPI) , 1);
        IPI.day = repmat(ANA.Day , unique(ANA.seqlength)-1 , 1);
        IPI.BN = repmat(ANA.BN , unique(ANA.seqlength)-1 , 1);
        IPI.ipiNum = reshape(repmat(ipinum , length(ANA.IPI) , 1) , numel(IPI.ipi) , 1);
        IPI.ipiArr = reshape(ANA.IPIarrangement , numel(ANA.IPIarrangement) , 1);
        IPI.seqN = repmat(ANA.seqNumb , unique(ANA.seqlength)-1 , 1);
        IPI.seqT = repmat(ANA.seqType , unique(ANA.seqlength)-1 , 1);
        
        % take the first and last IPI out
        IPI = getrow(IPI , ~ismember(IPI.ipiNum , [1,6]));
        
        IPIs = getrow(IPI , IPI.seqT > 0);
        out.IPI = anovan(IPIs.ipi , [IPIs.ipiArr , IPIs.day],'varnames',{'Within/Between' , 'day'} , 'display' , 'on' , 'model' , 'full');
        out.IPI_all = anovan(IPI.ipi , [IPI.seqT , IPI.day],'varnames',{'Random/Chunked' , 'day'} , 'display' , 'on' , 'model' , 'full');

        h1 = figure('color' , 'white');
        [xb,ePLOTb,ERRORb] = lineplot(IPI.day , IPI.ipi , 'plotfcn' , 'nanmedian' ,'subset' , IPI.ipiArr == 1);
        hold on
        [xw,ePLOTw,ERRORw] = lineplot(IPI.day , IPI.ipi , 'plotfcn' , 'nanmedian' ,'subset' , IPI.ipiArr == 2);
        [xr,ePLOTr,ERRORr] = lineplot(IPI.day , IPI.ipi , 'plotfcn' , 'nanmedian' ,'subset' , IPI.ipiArr == 0);
        dayss = unique(IPI.day);
        for d = 1:length(dayss)
            IPIs = getrow(IPI , IPI.seqN > 0 & IPI.day == dayss(d));
            out.IPId{d} = anovan(IPIs.ipi , IPIs.ipiArr,'varnames',{'Within/Between'} , 'display' , 'off' , 'model' , 'full');
            % for strucutre one
            [xs1{d},PLOTs1{d},ERRORs1{d}] = lineplot(IPI.ipiNum,IPI.ipi , 'plotfcn' , 'nanmedian' ,'subset' , IPI.day == dayss(d) & IPI.seqN == 1);
            hold on
            
            % for strucutre two
            [xs2{d},PLOTs2{d},ERRORs2{d}] = lineplot(IPI.ipiNum,IPI.ipi , 'plotfcn' , 'nanmedian' ,'subset' , IPI.day == dayss(d) & IPI.seqN == 3);
            
            % for random
            [xra{d},PLOTra{d},ERRORra{d}] = lineplot(IPI.ipiNum,IPI.ipi , 'plotfcn' , 'nanmedian' ,'subset' , IPI.day == dayss(d) & IPI.seqT == 0);
        end
        close(h1);
        
        figure('color' , 'white')
        subplot(length(dayss) , 2 , [1:2:2*length(dayss)])
        hb = plotshade(xb',ePLOTb,ERRORb,'transp' , .2 , 'patchcolor' , 'r' , 'linecolor' , 'r' , 'linewidth' , 3 );
        hold on
        hw = plotshade(xw',ePLOTw , ERRORw,'transp' , .2 , 'patchcolor' , 'g' , 'linecolor' , 'g' , 'linewidth' , 3  );
        hr = plotshade(xr',ePLOTr , ERRORr,'transp' , .2 , 'patchcolor' , 'b' , 'linecolor' , 'b' , 'linewidth' , 3  );
        legend([hb,hw,hr] , {'Between chunk IPIs' , 'Within chunk IPIs' , 'Random IPIs'});
        title('Within/Between Chunk and Random')
        xlabel('Training Days')
        ylabel('msec')
        set(gca , 'XTick' , unique(IPI.day),'FontSize' , 20)
        grid on
        figcount = 2;
        for d = 1:length(dayss)
            subplot(length(dayss) , 2 , figcount)
            hold on
            hs1 = plotshade(xs1{d}',PLOTs1{d} , ERRORs1{d},'transp' , .2 , 'patchcolor' , 'r' , 'linecolor' , 'r' , 'linewidth' , 3  );
            hs2 = plotshade(xs2{d}',PLOTs2{d} , ERRORs2{d},'transp' , .2 , 'patchcolor' , 'm' , 'linecolor' , 'm' , 'linewidth' , 3  );
            hr = plotshade(xra{d}',PLOTra{d} , ERRORra{d},'transp' , .2 , 'patchcolor' , 'b' , 'linecolor' , 'b' , 'linewidth' , 3 );
            legend([hs1,hs2,hr] , {'Structure 1' , 'Structure 2 ' , 'Random IPIs'})
            title(['Day ',num2str(d) ,', chunking effect: p(W/B) = ' , num2str(out.IPId{d})])
            xlabel('IPI Number')
            ylabel('msec')
            set(gca , 'XTick' , unique(IPI.day),'FontSize' , 20)
            grid on
            figcount = figcount+2;
        end
        
        %         out.MT = anovaMixed([ANA1.MT ; ANA0.MT] , [ANA1.SN ;ANA0.SN],'within',[0*ANA1.MT ; 1+0*ANA0.MT],{'Random/Chunked'},'intercept',1)  ;
        out.MT = anovan(ANA.MT, ANA.seqNumb , 'display' , 'off' , 'varnames' , 'Rand/Chunked');
        h1 = figure('color' , 'white');
        [xs,PLOTs,ERRORs] = lineplot(ANA.Day,ANA.MT , 'plotfcn' , 'nanmedian' ,'subset' , ANA.seqType == 1);
        hold on
        [xr,PLOTr,ERRORr] = lineplot(ANA.Day,ANA.MT, 'plotfcn' , 'nanmedian' ,'subset' , ANA.seqType == 0);
        close(h1);
        
        figure('color' , 'white')
        hs = plotshade(xs',PLOTs,ERRORs,'transp' , .2 , 'patchcolor' , 'r' , 'linecolor' , 'r' , 'linewidth' , 3 );
        hold on
        hr = plotshade(xr',PLOTr,ERRORr,'transp' , .2 , 'patchcolor' , 'b' , 'linecolor' , 'b' , 'linewidth' , 3  );
        hold on
        ax = gca;
        ax.FontSize = 20;
        title('Movement Time')
        ylabel('msec')
        xlabel('Training Days')
        set(gca , 'XTick' , unique(IPI.day),'FontSize' , 20)
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

        out.IPI = anovan([IPI3.ipi ; IPI4.ipi] , [[0*IPI3.ipi; 1+0*IPI4.ipi] , [IPI3.day; IPI4.day]],'varnames',{'Quadruple IPI/Triplet IPI' , 'day'} , 'display' , 'on' , 'model' , 'full');
        h1 = figure('color' , 'white');
        hold on
        [x3,ePLOT3,ERROR3] = lineplot( IPI3.day, IPI3.ipi,'plotfcn' , 'nanmedian');
        [x4,ePLOT4,ERROR4] = lineplot(IPI4.day , IPI4.ipi,'plotfcn' , 'nanmedian');
        
        [x3_2,ePLOT3_2,ERROR3_2] = lineplot( IPI3.ipiNum, IPI3.ipi,'plotfcn' , 'nanmedian');
        [x4_2,ePLOT4_2,ERROR4_2] = lineplot(IPI4.ipiNum , IPI4.ipi,'plotfcn' , 'nanmedian');
        close(h1);
        
        figure('color' , 'white')
        subplot(121)
        hb = plotshade(x3',ePLOT3,ERROR3,'transp' , .2 , 'patchcolor' , 'b' , 'linecolor' , 'b' , 'linewidth' , 3  );
        hold on
        hw = plotshade(x4',ePLOT4 , ERROR4,'transp' , .2 , 'patchcolor' , 'm' , 'linecolor' , 'm' , 'linewidth' , 3  );
        legend([hb,hw] , {'Triplet IPIs' , 'Quadruple IPIs'})
        title(['Triplet vs Quadruple Chunk IPIs , p(Q/T) = ' , num2str(out.IPI(2))])
        xlabel('Training Days')
        ylabel('msec')
        set(gca , 'XTick' , unique(IPI4.day),'FontSize' , 20)
        grid on
        
        subplot(122)
        hb = plotshade(x3_2',ePLOT3_2,ERROR3_2,'transp' , .2 , 'patchcolor' , 'b' , 'linecolor' , 'b' , 'linewidth' , 3  )
        hold on
        hw = plotshade(x4_2',ePLOT4_2,ERROR4_2,'transp' , .2 , 'patchcolor' , 'm' , 'linecolor' , 'm' , 'linewidth' , 3  )
        legend([hb,hw] , {'Triplet IPIs' , 'Quadruple IPIs'})
        title(['Triplet vs Quadruple Chunk IPIs , p(Q/T) = ' , num2str(out.IPI(2))])
        xlabel('IPI number')
        ylabel('msec' )
        set(gca , 'XTick' , unique(IPI4.day),'FontSize' , 20)
        grid on
        
        
        %         out.MT = anovaMixed([ANA1.MT ; ANA0.MT] , [ANA1.SN ;ANA0.SN],'within',[0*ANA1.MT ; 1+0*ANA0.MT],{'Random/Chunked'},'intercept',1)  ;
        out.MT = anovan([ANA3.MT ; ANA4.MT] , [0*ANA3.MT ; 1+0*ANA4.MT] , 'display' , 'off' , 'varnames' , 'Rand/Chunked');
        h1 = figure('color' , 'white');
        [x3,PLOT3,ERROR3] = lineplot(ANA3.Day,ANA3.MT);
        hold on
        [x4,PLOT4,ERROR4] = lineplot(ANA4.Day,ANA4.MT);
        close(h1);
        
        figure('color' , 'white')
        h3 = plotshade(x3',PLOT3,ERROR3,'transp' , .2 , 'patchcolor' , 'b' , 'linecolor' , 'b' , 'linewidth' , 3  )
        hold on
        h4 = plotshade(x4',PLOT4,ERROR4,'transp' , .2 , 'patchcolor' , 'm' , 'linecolor' , 'm' , 'linewidth' , 3  )
        hold on
        title(['Movement Time, p = ' , num2str(out.MT)])
        ylabel('msec')
        xlabel('Training Days')
        set(gca , 'XTick' , unique(IPI4.day),'FontSize' , 20)
        legend([h3 ,h4] , {'Triplet' , 'Quadruple'})
        grid on
    case 'RT_seq'
        
        %% within and between IPI ttest with randomization
        %% IPI per day
        ANA = getrow(Dall , Dall.isgood & ~Dall.isError & ismember(Dall.seqNumb , [0:4]));
        ANA.seqNumb(ANA.seqNumb > 0) = 1;
        BN = unique(ANA.BN);

        ANA1 = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.seqNumb , [1:4])  & ~Dall.isError &  ismember(Dall.Rep , rep) );
        ANA0 = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.seqNumb , [0])  & ~Dall.isError &  ismember(Dall.Rep , rep));

        ANA1.RT = ANA1.AllPressTimes(:,1);
        ANA0.RT = ANA0.AllPressTimes(:,1);

        %         out.RT = anovaMixed([ANA1.RT ; ANA0.RT] , [ANA1.SN ;ANA0.SN],'within',[0*ANA1.RT ; 1+0*ANA0.RT],{'Random/Chunked'},'intercept',1)  ;
        out.RT = anovan([ANA1.RT ; ANA0.RT] , [[0*ANA1.RT ; 1+0*ANA0.RT] [ANA1.Day ; ANA0.Day]] , 'display' , 'on' , 'varnames' , {'Rand/Chunked' , 'Day'} , 'model' , 'full');
        h1 = figure('color' , 'white');
        [xs,PLOTs,ERRORs] = lineplot(ANA1.Day,ANA1.RT,'plotfcn' , 'nanmedian');
        hold on
        [xr,PLOTr,ERRORr] = lineplot(ANA0.Day,ANA0.RT,'plotfcn' , 'nanmedian');
        close(h1);
        
        figure('color' , 'white')
        hs = plotshade(xs',PLOTs,ERRORs,'transp' , .2 , 'patchcolor' , 'r' , 'linecolor' , 'r' , 'linewidth' , 3  );
        hold on
        hr = plotshade(xr',PLOTr,ERRORr,'transp' , .2 , 'patchcolor' , 'b' , 'linecolor' , 'b' , 'linewidth' , 3  );
        hold on
        ax = gca;
        ax.FontSize = 20;
        title(['Reaction Time'])
        ylabel('msec')
        xlabel('Training Days')
        legend([hs ,hr] , {'Chunked' , 'Random'})
        set(gca , 'XTick' , unique(ANA1.Day),'FontSize' , 20)
        grid on        
    case 'RT_chunk'
        
        %% within and between IPI ttest with randomization
        %% IPI per day
        ANA = getrow(Dall , Dall.isgood & ~Dall.isError & ismember(Dall.seqNumb , [103 203 104 204]));
        ANA.seqNumb(ANA.seqNumb > 0) = 1;
        BN = unique(ANA.BN);
 
        ANA1 = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.seqNumb , [104 204])  & ~Dall.isError &  ismember(Dall.Rep , rep) );
        ANA0 = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.seqNumb , [103 203])  & ~Dall.isError &  ismember(Dall.Rep , rep));

        ANA1.RT = ANA1.AllPressTimes(:,1);
        ANA0.RT = ANA0.AllPressTimes(:,1);

        out.RT = anovan([ANA1.RT ; ANA0.RT] , [[0*ANA1.RT ; 1+0*ANA0.RT] [ANA1.Day ; ANA0.Day]]  , 'display' , 'on' , 'varnames' , {'Quadruple/Trilet' , 'Day'} , 'model' , 'full');
        h1 = figure('color' , 'white');
        [xs,PLOTs,ERRORs] = lineplot(ANA1.Day,ANA1.RT,'plotfcn' , 'nanmedian');
        hold on
        [xr,PLOTr,ERRORr] = lineplot(ANA0.Day,ANA0.RT,'plotfcn' , 'nanmedian');
        close(h1);
        
        figure('color' , 'white')
        hs = plotshade(xs',PLOTs,ERRORs,'transp' , .2 , 'patchcolor' , 'black' , 'linecolor' , 'k' , 'linewidth' , 3  )
        hold on
        hr = plotshade(xr',PLOTr,ERRORr,'transp' , .2 , 'patchcolor' , 'r' , 'linecolor' , 'r' , 'linewidth' , 3  )
        hold on
        ax = gca;
        ax.FontSize = 20;
        title(['Reaction Time'])
        ylabel('msec')
        xlabel('Training Days')
        legend([hs ,hr] , {'Quadruples' , 'Triplets'})
        set(gca , 'XTick' , unique(ANA1.Day),'FontSize' , 20)
        grid on
    case 'RT_SF'
        
        %% within and between IPI ttest with randomization
        %% IPI per day
        ANA = getrow(Dall , Dall.isgood & ~Dall.isError & ismember(Dall.seqNumb , [11 22 33 44 55]));
        ANA.RT = ANA.AllPressTimes(:,1);
        out.RT = anovan([ANA.RT] ,[ANA.Day ANA.seqNumb] , 'display' , 'on' , 'model' , 'full','varnames' , {'Day', 'Single Finger Sequences'});
        figure('color' , 'white')
        [xs,PLOTs,ERRORs] = lineplot([ANA.Day ANA.seqNumb],ANA.RT , 'style_shade','plotfcn' , 'nanmedian');
        title(['Reaction Time'])
        ylabel('msec')
        xlabel('Training Days')
        set(gca ,'FontSize' , 20)
        grid on       
    case 'MT_SF'
        
        %% within and between IPI ttest with randomization
        %% IPI per day
        ANA = getrow(Dall , Dall.isgood & ~Dall.isError & ismember(Dall.seqNumb , [11 22 33 44 55]));
        
        

        ANA.MT = ANA.AllPressTimes(:,7) - ANA.AllPressTimes(:,1);
        
        %         out.MT = anovaMixed([ANA1.MT ; ANA0.MT] , [ANA1.SN ;ANA0.SN],'within',[0*ANA1.MT ; 1+0*ANA0.MT],{'Random/Chunked'},'intercept',1)  ;
        out.MT = anovan([ANA.MT] , ANA.seqNumb , 'display' , 'on');
        h1 = figure('color' , 'white');
        [xs,PLOTs,ERRORs] = lineplot(ANA.BN,ANA.MT,'plotfcn' , 'nanmedian');
        hold on
        close(h1);
        
        figure('color' , 'white')
        hs = plotshade(xs',PLOTs,ERRORs,'transp' , .2 , 'patchcolor' , 'black' , 'linecolor' , 'k' , 'linewidth' , 3  )
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
     