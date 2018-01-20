function ANA=secog_subj(subjname,fig,block,trial)



% ------------------------- Directories -----------------------------------
baseDir         ='/Volumes/MotorControl/data/SeqECoG/ecog1/data';
% -------------------------------------------------------------------------
cd (baseDir)
mkdir('analyze');

if nargin<2
    fig=0;
end;
datafilename = ['ecog1_' subjname '.dat'];

ANA=[];

D=dload(datafilename);
if (nargin<3)
    % ------ if not specified read all blocks and trials
    block = unique(D.BN);
    outfilename  = ['analyze/ecog1_' subjname '.mat'];
    trial  = [];
elseif (nargin<4)
    trial = [];
    outfilename  = ['analyze/ecog1_' subjname '_B',num2str(block),'.mat'];
else
    outfilename  = ['analyze/ecog1_' subjname '_B',num2str(block),'_T',num2str(trial),'.mat'];
    
end;
trcount = 1;
%define  number of trials
for b  = 1: length(block)
    disp(['Reading Block ' , num2str(block(b))])
    MOV   = movload(['ecog1_' subjname '_' num2str(block(b),'%02d') '.mov']); % all trials of a block  
    if isempty(trial)
        trial = D.TN(find(D.BN==block(b)));
    end
    
    for i=1:length(trial) % loop over all trials
        [C]=secog_trial(MOV{D.TN(trial(i))},getrow(D,(D.BN == block(b) & D.TN == trial(i))),fig);
        fprintf('%d %d\n',block(b),D.TN(trial(i)));
        ANA=addstruct(ANA,C);
    end

    trial  = [];
end

