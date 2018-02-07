tn  = 20;
bn = 1;

% eval(['fname = /Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/iEEG_P1_Sep14_BlockParsed/AVG_Norm_PSD_b' , num2str(bn) , '_P1_Sep14.mat'])
% load(fname)
Fs = 1024;
alltrials = 1;
% definitions, selections...
BN = unique(Dall.BN);
min_freq =  2;
max_freq = 150;
num_frex = 90;

% define wavelet parameters
time = -1:1/Fs:1;
frex = logspace(log10(min_freq),log10(max_freq),num_frex);
Fs = 1024;
fName = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/iEEG_P1_Sep14_BlockParsed/EEG_BN' , num2str(bn) , '.mat'];
load(fName);
GoodElec = [1:30,32:60,63:74,77:80];
eval(['A = EEG_BN' , num2str(bn) , '(GoodElec , :);']);

for i = 1:size(A , 1)
    V(i) = var(A(i , :));
end
[~,Vidx] = sort(V , 'descend');

D = getrow(Dall , ismember(Dall.BN ,bn) & ismember(Dall.TN ,tn));
NEM = find(D.NormEventMarker{1});
a = Pall.PEEG_stim{tn};
b = Pall.PEEG_none{tn};


figure
for i = 1:size(a , 1)
    subplot(1,2,1)
    contourf([1:500],frex,squeeze(a(Vidx(i),:,:)),60,'linecolor','none');
    title(['nonNorm  Block ' , num2str(bn) , ', Trial ' , num2str(tn) , ', SeqNumb = ' , num2str(D.seqNumb)])
    for lin = 1:length(NEM)
        line([NEM(lin) NEM(lin)] , [2 150] , 'color' , 'red' , 'LineStyle' , ':' , 'LineWidth' , 3)
    end
    subplot(1,2,2)
    contourf([1:500],frex,squeeze(b(Vidx(i),:,:)),60,'linecolor','none');
    title(['StimNorm Block ' , num2str(bn) , ', Trial ' , num2str(tn) , ', SeqNumb = ' , num2str(D.seqNumb)])
    for lin = 1:length(NEM)
        line([NEM(lin) NEM(lin)] , [2 150] , 'color' , 'red' , 'LineStyle' , ':' , 'LineWidth' , 3)
    end
    pause()
    drawnow
end