


baseDir = '/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/';
filename = [baseDir , 'p1_Sep-06.dat'];
d= fopen(filename,'r','n','UTF16-LE');
LCounter = 1;

% Skip 15 headerlines
for n=1:15
  tline = fgetl(d)
end
LCounter = 1;

tic
while ischar(tline)
    tline = fgetl(d);
    delimiter = '\t';
    C = strsplit(tline);
    for i = 1 : length(C)
        if isempty(str2num(C{i})) | length(str2num(C{i}))>1
            Data(LCounter, i) = NaN;
        else
            Data(LCounter, i) = str2num(C{i});
        end
     end
    LCounter  = LCounter +1
end
toc
fclose(d)