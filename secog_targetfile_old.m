function secog_targetfile(SubCode,genchunks, CMB)


% genchunks = 0;
% WantTheStars = 0;
% SubCode = 'SZ';
% GroupCode= 1;
% provode CMB from the subject's target file folder only if you already have the chunks and need to  produce more sequences
% baseDir = '/Users/nkordjazi/Documents/SeqECoG/TargetFiles';
baseDir = '/Users/nedakordjazi/Documents/SeqECoG/TargetFiles';
cd(baseDir)
Fname  = [SubCode , '_tgtFiles'];
mkdir(Fname)


if genchunks
    CMB = secog_getChunks(4,0,2,2, 0);
end

CMB.DesiredCnhkargmnt = [3 4;4 3];
CMB.Seq2Chunk = [1 2 3 1 2 3 4;1 2 3 4 1 2 3];
MaxPress = size(CMB.Seq2Chunk ,2);
save([Fname , '/' , SubCode , '_CMB'] , 'CMB')


% make target file
% switch what
%     case 'targetfile'

rng('shuffle');



SequenceLength = size(CMB.Seq2Chunk , 2);
Chunks = cell2mat(CMB.Chunks);
NumofChunks = size(Chunks , 1);
for cc = 1:NumofChunks
    L(cc) = length(find(Chunks(cc,:)));
end
ccounter = 1;
for cL = 2:4
    j = length(find(L == cL));
    for i = 1:j
        ChunkNumb(ccounter) = (i * 100) + cL;
        ccounter = ccounter + 1;
    end
end




repeating  = 1; % determins that every sequence in the CLAT and Intermixed blocks happens twice
FT = 1;

OrderFields = {'seqNumb','FT','press1','press2','press3','press4','press5','press6','press7','hand','cueS','cueC','cueP','iti','sounds','pulse'};



%% Chunk Training Block
%WantTheStars = 0;
Trials = 1:32+2;  % 32 are the number of trials and 2 is the number of baseline collection trials (the stars)

repeating = 1;
%             NumStarTrials = [1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 9 9 9];
%             OrderFields = {'seqNumb','FT','press1','press2','press3','press4','press5','press6','press7','press8','press9','press10','press11','press12','press13','press14','hand','cueS','cueC','cueP','iti','sounds'};

clear ChunkTrain
% StimTimeLim = zeros(length(Trials) , 1);
if ~ repeating
    for e = 1:15 % number of Blocks
        ChunkTrain.cueS = cellstr(repmat('£',length(Trials),1));
        ChunkTrain.FT = 2 * ones(length(Trials),1);
        ChunkTrain.iti(1:length(Trials),:) = 500;
        ChunkTrain.hand(1:length(Trials),:) = 2;
        ChunkTrain.sounds(1:length(Trials),:) = 1;
        ChunkTrain.pulse(1:length(Trials),:) = 1;
        %     ChunkTrain.StimTimeLim = StimTimeLim;
        
        % making sure that all the chunks are ocurring equal times in each chunk trianing block
        X = [];
        for rep = 1: length(Trials)/NumofChunks
            X = [X  randperm(NumofChunks)];
        end
        X = X(randperm(length(X))); % so that the the cycle of going through all the possible chunks in one run is not detected
        %                 X = kron(X, ones(NumStarTrials(e) + 1 , 1));
        %                 X = X(1:length(Trials));
        % Indeices in where the stras should apear as 0
        %                 starInd = repmat([1 ; zeros(NumStarTrials(e) ,1)] , length(Trials)/(NumStarTrials(e)+1) ,1);
        for i = 1:length(Trials)
            Numdigs = sum((Chunks(X(i),:) ~= 0 ));
            ChunkTrain.cueP{i,1} = char(regexprep(cellstr(num2str(Chunks(X(i),1:Numdigs))),'\s',''));
            
            for press= 1:Numdigs
                comnd  = [' ChunkTrain.press' , num2str(press) , '(i,1) = Chunks(X(i),press);'];
                eval(comnd);
            end
            for press= Numdigs + 1 : SequenceLength
                comnd  = [' ChunkTrain.press' , num2str(press) , '(i,1) = 0;'];
                eval(comnd);
            end
            % Chunk length arrangement number
            ChunkTrain.seqNumb(i,1) = ChunkNumb(X(i));
            ChunkTrain.cueC(i,:) = {'£'} ;
        end
        name = [Fname ,'/'  , SubCode , '_CT_B' , num2str(e) , '.tgt'];
        dsave(name,orderfields(ChunkTrain,OrderFields));
        clear x name
        
        clear ChunkTrain
        
    end
else
    Trials = 1: 32/(repeating+1);
    for e = 1:15 % number of Blocks
        ChunkTrain.cueS = cellstr(repmat('£',length(Trials),1));
        ChunkTrain.FT = 2 * ones(length(Trials),1);
        ChunkTrain.iti(1:length(Trials),:) = 500;
        ChunkTrain.hand(1:length(Trials),:) = 2;
        ChunkTrain.sounds(1:length(Trials),:) = 1;
        ChunkTrain.pulse(1:length(Trials),:) = 1;
        %     ChunkTrain.StimTimeLim = StimTimeLim;
        
        % making sure that all the chunks are ocurring equal times in each chunk trianing block
        X = [];
        for rep = 1: length(Trials)/NumofChunks
            X = [X  randperm(NumofChunks)];
        end
        X = X(randperm(length(X))); % so that the the cycle of going through all the possible chunks in one run is not detected
        %                 X = kron(X, ones(NumStarTrials(e) + 1 , 1));
        %                 X = X(1:length(Trials));
        % Indeices in where the stras should apear as 0
        %                 starInd = repmat([1 ; zeros(NumStarTrials(e) ,1)] , length(Trials)/(NumStarTrials(e)+1) ,1);
        for i = 1:length(Trials)
            Numdigs = sum((Chunks(X(i),:) ~= 0 ));
            ChunkTrain.cueP{i,1} = char(regexprep(cellstr(num2str(Chunks(X(i),1:Numdigs))),'\s',''));
            
            for press= 1:Numdigs
                comnd  = [' ChunkTrain.press' , num2str(press) , '(i,1) = Chunks(X(i),press);'];
                eval(comnd);
            end
            for press= Numdigs + 1 : SequenceLength
                comnd  = [' ChunkTrain.press' , num2str(press) , '(i,1) = 0;'];
                eval(comnd);
            end
            % Chunk length arrangement number
            ChunkTrain.seqNumb(i,1) = ChunkNumb(X(i));
            ChunkTrain.cueC(i,:) = {'£'} ;
        end
        repChunkTrain = getrow(ChunkTrain ,  ChunkTrain.seqNumb == 0);
        for i = Trials
            tempChunkTrain = getrow(ChunkTrain , i);
            for repp = 1:repeating+1
                repChunkTrain  = addstruct(repChunkTrain , tempChunkTrain);
            end
        end
        
        name = [Fname ,'/'  , SubCode , '_CT_B' , num2str(e) , '.tgt'];
        dsave(name,orderfields(repChunkTrain,OrderFields));
        clear x name
        
        clear ChunkTrain repChunkTrain
        
    end
end


%% Test blocks - Random Sequences
clear RandomSeq
if ~repeating
    ElimChunkStart  = 1; % 1 means that random seqs that start with a known chunk will be elimminated
    Trials = 1:24;
    % StimTimeLim = zeros(length(Trials) , 1);
    % h = [1:12]';
    % h =  kron(h, ones(4 , 1));
    % Generate 30 Full Horizons, and then 4 of each for 3 : 12
    % Horizon = [(SequenceLength - 1) * ones(30 , 1) ; h];
    
    % just generate full horizons for now
    % Horizon = [(SequenceLength - 1) * ones(10 , 1)];
    
    
    
    for e= 1:5
        %     RandomSeq.StimTimeLim = StimTimeLim;
        RandomSeq.seqNumb = 0;
        RandomSeq.FT(1:length(Trials),1) = 2;
        %     RandomSeq.Horizon(1:length(Trials),:) = Horizon(e);
        for i = 1:length(Trials)
            %                 ChunkArrange = CMB.DesiredCnhkargmnt(CLA(X(i)) , :);
            if ElimChunkStart
                seq = Chunks(1,1:2); % makes sure to enter the first while loop
                seq = [seq  1 1 1 1];
                tempSeq  = diff(seq);
                tempSeq2 = diff(tempSeq);
                while length(unique([seq(1:2) ; Chunks(:,1:2)] , 'rows')) <= length(Chunks) | sum(tempSeq == 0) > 1 | sum(tempSeq == 1) > 1 | sum(tempSeq == -1) > 1
                    seq = sample_wor([1:5],1,MaxPress);
                    tempSeq = diff(seq);
                    tempSeq2 = diff(tempSeq);
                end
            else
                seq = sample_wor([1:5],1,MaxPress);
                tempSeq = diff(seq);
                tempSeq2 = diff(tempSeq);
                while (sum(tempSeq == 0) > 1)| (sum(tempSeq == 1) > 1) | (sum(tempSeq == -1) > 1 )
                    seq = sample_wor([1:5],1,MaxPress);
                    tempSeq = diff(seq);
                    tempSeq2 = diff(tempSeq);
                end
            end
            RandomSeq.cueP{i,:} = char(regexprep(cellstr(num2str(seq)),'\s',''));
            RandomSeq.cueC(i,:) ={'£'};
            RandomSeq.cueS(i,:) ={'£'};
            RandomSeq.iti(1:i,:) = 500;
            RandomSeq.hand(i,:) = 2;
            RandomSeq.sounds(i,:) = 1;
            RandomSeq.pulse(i,:) = 1;
            
            for press= 1:MaxPress
                comnd  = [' RandomSeq.press' , num2str(press) , '(i,1) = seq(press);'];
                eval(comnd);
            end
        end
        
        name = [Fname ,'/'  , SubCode , '_RAND_B' , num2str(e) , '.tgt'];
        dsave(name,orderfields(RandomSeq,OrderFields));
        
        clear x
        clear RandomSeq
    end
else
    ElimChunkStart  = 1;
    Trials = 1:24;
    uniqueTrials = 1:max(Trials)/(repeating+1);    
    for e= 1:10
        %     RandomSeq.StimTimeLim = StimTimeLim;
        RandomSeq.seqNumb(1:length(uniqueTrials) , 1) = 0;
        RandomSeq.FT(1:length(uniqueTrials) , 1)  = 2;
        %     RandomSeq.Horizon(1:length(Trials),:) = Horizon(e);
        for i = 1:length(uniqueTrials)
            %                 ChunkArrange = CMB.DesiredCnhkargmnt(CLA(X(i)) , :);
            if ElimChunkStart
                seq = Chunks(1,1:2); % makes sure to enter the first while loop
                seq = [seq  1 1 1 1];
                tempSeq  = diff(seq);
                tempSeq2 = diff(tempSeq);
                while length(unique([seq(1:2) ; Chunks(:,1:2)] , 'rows')) <= length(Chunks) | sum(tempSeq == 0) > 1 | sum(tempSeq == 1) > 1 | sum(tempSeq == -1) > 1
                    seq = sample_wor([1:5],1,MaxPress);
                    tempSeq = diff(seq);
                    tempSeq2 = diff(tempSeq);
                end
            else
                seq = sample_wor([1:5],1,MaxPress);
                tempSeq = diff(seq);
                tempSeq2 = diff(tempSeq);
                while (sum(tempSeq == 0) > 1)| (sum(tempSeq == 1) > 1) | (sum(tempSeq == -1) > 1 )
                    seq = sample_wor([1:5],1,MaxPress);
                    tempSeq = diff(seq);
                    tempSeq2 = diff(tempSeq);
                end
            end
            RandomSeq.cueP{i,:} = char(regexprep(cellstr(num2str(seq)),'\s',''));
            RandomSeq.cueC(i,:) ={'£'};
            RandomSeq.cueS(i,:) ={'£'};
            RandomSeq.iti(i,:) = 500;
            RandomSeq.hand(i,:) = 2;
            RandomSeq.sounds(i,:) = 1;
            RandomSeq.pulse(i,:) = 1;
            
            for press= 1:MaxPress
                comnd  = [' RandomSeq.press' , num2str(press) , '(i,1) = seq(press);'];
                eval(comnd);
            end
        end
        
        RS = getrow(RandomSeq , RandomSeq.seqNumb == 10);
        for i = uniqueTrials
            temprs = getrow(RandomSeq , i);
            for repp = 1:repeating+1
                RS  = addstruct(RS , temprs);
            end
        end
    
        name = [Fname ,'/'  , SubCode , '_RAND_B' , num2str(e) , '.tgt'];
        dsave(name,orderfields(RS,OrderFields));
        
        clear x
        clear RandomSeq RS temprs
    end
    
end

%% Chunk Arrangemnet Training Block
tempCL = rem(ChunkNumb , 100);
if ~repeating
    CLA = [1 ;2];
    Trials = 1:24;
%     StimTimeLim = zeros(length(Trials) , 1);
    
    clear ChunkArrangeLearn
    
    for e= 1:15
%         ChunkArrangeLearn.StimTimeLim = StimTimeLim;
        X = [];
        for rep = 1: length(Trials) / length(CLA)
            X = [X ; sample_wor([1:length(CLA)],length(CLA))];
        end
        X = X(randperm(length(X))); % so that the the cycle of going through all the possible CLAs in one run is not detected
        ChunkArrangeLearn.seqNumb = X;
        ChunkArrangeLearn.FT(1:length(Trials),1) = 2;
        for i = 1:length(X)
            ChunkArrange = CMB.DesiredCnhkargmnt(CLA(X(i)) , :);
            seq = [];
            for j = 1: length(find(ChunkArrange))
                numCL(j) = sum(tempCL == ChunkArrange(j));
                chnkind = randi(numCL(j));
                seq = [seq CMB.Chunks{ChunkArrange(j)}(chnkind , 1:ChunkArrange(j))];
            end
            ChunkArrangeLearn.cueP{i,:} = char(regexprep(cellstr(num2str(seq)),'\s',''));
            ChunkArrangeLearn.cueC(i,:) ={'£'};
            ChunkArrangeLearn.cueS(i,:) ={'£'};
            ChunkArrangeLearn.iti(1:i,:) = 500;
            ChunkArrangeLearn.hand(i,:) = 2;
            ChunkArrangeLearn.sounds(i,:) = 1;
            ChunkArrangeLearn.pulse(i,:) = 1;
            
            
            for press= 1:MaxPress
                comnd  = [' ChunkArrangeLearn.press' , num2str(press) , '(i,1) = seq(press);'];
                eval(comnd);
            end
        end
        
        name = [Fname ,'/'  , SubCode , '_CAT_B' , num2str(e) , '.tgt'];
        dsave(name,orderfields(ChunkArrangeLearn,OrderFields));
        
        clear x
        clear ChunkArrangeLearn
    end
else
    CLA = [1;2];
    
    % just generate full horizons for now
    
    Trials = 1:24;
%     StimTimeLim = zeros(length(Trials) , 1);
    
    clear ChunkArrangeLearn
    
    for e= 1:15
%         ChunkArrangeLearn.StimTimeLim = StimTimeLim;
        X = [];
        for rep = 1: .5*length(Trials) / length(CLA)
            X = [X ; sample_wor([1:length(CLA)],length(CLA))];
        end
        X = X(randperm(length(X))); % so that the the cycle of going through all the possible CLAs in one run is not detected
        
        ChunkArrangeLearn.FT(1:length(Trials),1) = 2;
        cn = 1;
        for i = 1:length(X)
            ChunkArrange = CMB.DesiredCnhkargmnt(CLA(X(i)) , :);
            seq = [];
            for j = 1: length(find(ChunkArrange))
                numCL(j) = sum(tempCL == ChunkArrange(j));
                chnkind = randi(numCL(j));
                seq = [seq CMB.Chunks{ChunkArrange(j)}(chnkind , 1:ChunkArrange(j))];
            end
            ChunkArrangeLearn.cueP{cn,:} = char(regexprep(cellstr(num2str(seq)),'\s',''));
            ChunkArrangeLearn.cueC(cn,:) ={'£'};
            ChunkArrangeLearn.cueS(cn,:) ={'£'};
            ChunkArrangeLearn.iti(1:cn,:) = 500;
            ChunkArrangeLearn.hand(cn,:) = 2;
            ChunkArrangeLearn.sounds(cn,:) = 1;
            ChunkArrangeLearn.pulse(cn,:) = 1;
            ChunkArrangeLearn.seqNumb(cn , :) = X(i);
            for press= 1:MaxPress
                comnd  = [' ChunkArrangeLearn.press' , num2str(press) , '(cn,1) = seq(press);'];
                eval(comnd);
            end
            
            cn = cn + 1;
            
            ChunkArrangeLearn.cueP{cn,:} = char(regexprep(cellstr(num2str(seq)),'\s',''));
            ChunkArrangeLearn.cueC(cn,:) ={'£'};
            ChunkArrangeLearn.cueS(cn,:) ={'£'};
            ChunkArrangeLearn.iti(1:cn,:) = 500;
            ChunkArrangeLearn.hand(cn,:) = 2;
            ChunkArrangeLearn.sounds(cn,:) = 1;
            ChunkArrangeLearn.pulse(cn,:) = 1;
            ChunkArrangeLearn.seqNumb(cn , :) = X(i);
            
            for press= 1:MaxPress
                comnd  = [' ChunkArrangeLearn.press' , num2str(press) , '(cn,1) = seq(press);'];
                eval(comnd);
            end
            
            cn  = cn + 1;
        end
        
        name = [Fname ,'/'  , SubCode , '_CAT2_B' , num2str(e) , '.tgt'];
        dsave(name,orderfields(ChunkArrangeLearn,OrderFields));
        
        clear x
        clear ChunkArrangeLearn
    end
    
    
end

%% The inrtermixed CLAT and Random blocks
if ~ repeating    
    CLA = [1;2];
    Trials = 1:36;
%     StimTimeLim = zeros(length(Trials) , 1);
    RandProportion  = 1/3; % Can be 1/3 , 1/2 or 2/3 - portion of the sequences that are going to be random
    ElimChunkStart  = 1; % 1 means that random seqs that start with a known chunk will be elimminated
    % Generate 30 Full Horizons, and then 4 of each for 3 : 12
    % Horizon = [(SequenceLength - 1) * ones(20 , 1) ; h];
    
    % just generate full horizons for now
    clear Intermixed
    for e= 1:15
        RandIndex = randperm(length(Trials));
        RandIndex = RandIndex(1: RandProportion * length(Trials));
        
        CLATIndex = Trials(~ismember(Trials , RandIndex));
        % first create the CLAT sequecnces
        X = [];
        for rep = 1: (length(Trials) - length(RandIndex)) / length(CLA)
            X = [X ; sample_wor([1:length(CLA)],length(CLA))];
        end
        X = X(randperm(length(X))); % so that the the cycle of going through all the possible CLAs in one run is not detected
        Intermixed.seqNumb  = zeros(length(Trials) , 1);
        Intermixed.seqNumb(CLATIndex) = X;
        Intermixed.FT(1:length(Trials),1) = 2;
        Intermixed.iti(1:length(Trials),:) = 500;
        Intermixed.hand(1:length(Trials),:) = 2;
        Intermixed.sounds(1:length(Trials),:) = 1;
        Intermixed.pulse(1:length(Trials),:) = 1;
        
        for i = 1:length(CLATIndex)
            
            ChunkArrange = CMB.DesiredCnhkargmnt(CLA(X(i)) , :);
            seq = [];
            for j = 1: length(find(ChunkArrange))
                numCL(j) = sum(tempCL == ChunkArrange(j));
                chnkind = randi(numCL(j));
                seq = [seq CMB.Chunks{ChunkArrange(j)}(chnkind , 1:ChunkArrange(j))];
            end
            Intermixed.cueP{CLATIndex(i),:} = char(regexprep(cellstr(num2str(seq)),'\s',''));
            Intermixed.cueC(CLATIndex(i),:) ={'£'};
            Intermixed.cueS(CLATIndex(i),:) ={'£'};
            
            for press= 1:MaxPress
                comnd  = [' Intermixed.press' , num2str(press) , '(CLATIndex(i),1) = seq(press);'];
                eval(comnd);
            end
        end
        
        for i = 1:length(RandIndex)
            %                 ChunkArrange = CMB.DesiredCnhkargmnt(CLA(X(i)) , :);
            
            if ElimChunkStart
                seq = Chunks(1,1:2); % makes sure to enter the first while loop
                seq = [seq  1 1 1 1];
                tempSeq  = diff(seq);
                tempSeq2 = diff(tempSeq);
                while length(unique([seq(1:2) ; Chunks(:,1:2)] , 'rows')) <= length(Chunks) | sum(tempSeq == 0) > 1 | sum(tempSeq == 1) > 1 | sum(tempSeq == -1) > 1
                    seq = sample_wor([1:5],1,MaxPress);
                    tempSeq = diff(seq);
                    tempSeq2 = diff(tempSeq);
                end
            else
                seq = sample_wor([1:5],1,MaxPress);
                tempSeq = diff(seq);
                tempSeq2 = diff(tempSeq);
                while (sum(tempSeq == 0) > 1)| (sum(tempSeq == 1) > 1) | (sum(tempSeq == -1) > 1 )
                    seq = sample_wor([1:5],1,MaxPress);
                    tempSeq = diff(seq);
                    tempSeq2 = diff(tempSeq);
                end
            end
            
            
            Intermixed.cueP{RandIndex(i),:} = char(regexprep(cellstr(num2str(seq)),'\s',''));
            Intermixed.cueC(RandIndex(i),:) ={'£'};
            Intermixed.cueS(RandIndex(i),:) ={'£'};
            for press= 1:MaxPress
                comnd  = [' Intermixed.press' , num2str(press) , '(RandIndex(i),1) = seq(press);'];
                eval(comnd);
            end
        end
        switch RandProportion
            case 1/3
                name = [Fname ,'/'  , SubCode , '_TM_B' , num2str(e) , '.tgt'];
            case 1/2
                name = [Fname ,'/'  , SubCode , '_HM_B' , num2str(e) , '.tgt'];
            case 2/3
                name = [Fname ,'/'  , SubCode , '_2TM_B', num2str(e) , '.tgt'];
        end
        
        dsave(name,orderfields(Intermixed,OrderFields));
        
        clear x Intermixed
    end
else

    RandProportion  = 1/3; % Can be 1/3 , 1/2 or 2/3 - portion of the sequences that are going to be random
    ElimChunkStart  = 1;
    CLA = [1;2];
    
    Trials = 1:36;
    
%     StimTimeLim = zeros(length(Trials) , 1);
    
    clear ChunkArrangeLearn
    clear RandomSeq
    
    
    % Test blocks - Random Sequences
    
    ElimChunkStart  = 1; % 1 means that random seqs that start with a known chunk will be elimminated
%     StimTimeLim = zeros(length(Trials) , 1);
       
    for e= 1:15       
%         ChunkArrangeLearn.StimTimeLim = StimTimeLim;
        X = [];
        
        for rep = 1: (1-RandProportion)*length(Trials) / (length(CLA)*(repeating+1))
            X = [X ; sample_wor([1:length(CLA)],length(CLA))];
        end
        X = X(randperm(length(X))); % so that the the cycle of going through all the possible CLAs in one run is not detected
        
        
        cn = 1;
        uniqidx = 1;
        for i = 1:length(X)
            ChunkArrangeLearn.FT(i,1) = 2;
            ChunkArrange = CMB.DesiredCnhkargmnt(CLA(X(i)) , :);
            seq = [];
            for j = 1: length(find(ChunkArrange))
                numCL(j) = sum(tempCL == ChunkArrange(j));
                chnkind = randi(numCL(j));
                seq = [seq CMB.Chunks{ChunkArrange(j)}(chnkind , 1:ChunkArrange(j))];
            end
            ChunkArrangeLearn.cueP{i,:} = char(regexprep(cellstr(num2str(seq)),'\s',''));
            ChunkArrangeLearn.cueC(i,:) ={'£'};
            ChunkArrangeLearn.cueS(i,:) ={'£'};
            ChunkArrangeLearn.iti(1:i,:) = 500;
            ChunkArrangeLearn.hand(i,:) = 2;
            ChunkArrangeLearn.sounds(i,:) = 1;
            ChunkArrangeLearn.pulse(i,:) = 1;
            ChunkArrangeLearn.seqNumb(i , :) = X(i);
            for press= 1:MaxPress
                comnd  = [' ChunkArrangeLearn.press' , num2str(press) , '(i,1) = seq(press);'];
                eval(comnd);
            end
            
           
        end
        
        %name = [Fname ,'/'  , SubCode , num2str(GroupCode) , '_CLAT2_h' ,num2str(Horizon(e)), '_B' , num2str(e) , '.tgt'];
        %dsave(name,orderfields(ChunkArrangeLearn,OrderFields));
        
        clear x
        
        
        
        
        

        rTrials = 1:length(Trials)*RandProportion/(repeating+1);
%         RandomSeq.StimTimeLim = StimTimeLim
        RandomSeq.seqNumb = zeros(length(rTrials) , 1);
        RandomSeq.FT(1:length(rTrials),1) = 2;
        cn  = 1;
        for i = 1:length(rTrials)
            %                 ChunkArrange = CMB.DesiredCnhkargmnt(CLA(X(i)) , :);
            
            if ElimChunkStart
                seq = Chunks(1,1:2); % makes sure to enter the first while loop
                seq = [seq  1 1 1 1];
                tempSeq  = diff(seq);
                tempSeq2 = diff(tempSeq);
                while length(unique([seq(1:2) ; Chunks(:,1:2)] , 'rows')) <= length(Chunks) | sum(tempSeq == 0) > 1 | sum(tempSeq == 1) > 1 | sum(tempSeq == -1) > 1
                    seq = sample_wor([1:5],1,MaxPress);
                    tempSeq = diff(seq);
                    tempSeq2 = diff(tempSeq);
                end
            else
                seq = sample_wor([1:5],1,MaxPress);
                tempSeq = diff(seq);
                tempSeq2 = diff(tempSeq);
                while (sum(tempSeq == 0) > 1)| (sum(tempSeq == 1) > 1) | (sum(tempSeq == -1) > 1 )
                    seq = sample_wor([1:5],1,MaxPress);
                    tempSeq = diff(seq);
                    tempSeq2 = diff(tempSeq);
                end
            end
            
            
            RandomSeq.cueP{i,:} = char(regexprep(cellstr(num2str(seq)),'\s',''));
            RandomSeq.cueC(i,:) ={'£'};
            RandomSeq.cueS(i,:) ={'£'};
            RandomSeq.iti(1:i,:) = 500;
            RandomSeq.hand(i,:) = 2;
            RandomSeq.sounds(i,:) = 1;
            RandomSeq.pulse(i,:) = 1;
            for press= 1:MaxPress
                comnd  = [' RandomSeq.press' , num2str(press) , '(i,1) = seq(press);'];
                eval(comnd);
            end
            
        end
        
        ALLSEQ = addstruct(RandomSeq , ChunkArrangeLearn);
        ids = randperm(length(ALLSEQ.FT));
        ALLSEQ_new = getrow(ALLSEQ , ALLSEQ.seqNumb == 10);
        for j = 1:length(ids)
            temp = getrow(ALLSEQ , ids(j));
            ALLSEQ_new = addstruct(ALLSEQ_new , temp);
        end
        Intermixed = getrow(ALLSEQ_new , ALLSEQ_new.seqNumb == 10);
        for j = 1:length(ALLSEQ_new.FT)
            temp = getrow(ALLSEQ_new , j);
            for repp = 1:repeating+1
            Intermixed = addstruct(Intermixed , temp);
            end
        end
         
        
        
        clear x
        clear RandomSeq
        clear ChunkArrangeLearn
        
        
        switch RandProportion
            case 1/3
                name = [Fname ,'/'  , SubCode , '_TM_B' , num2str(e) , '.tgt'];
            case 1/2
                name = [Fname ,'/'  , SubCode , '_HM_B' , num2str(e) , '.tgt'];
            case 2/3
                name = [Fname ,'/'  , SubCode , '_2TM_B', num2str(e) , '.tgt'];
        end
        
        
        dsave(name,orderfields(Intermixed,OrderFields));
        
        
        clear Intermixed
        
        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Merge






