close all;
ipath = 'K:\20180816 chl1 647 b2 680\DIV18\split\';
pixelSize = 167; 
dbin = 10/1000;  % micr
bb=0;
average_c=[];
txt_file='autocorrelation_amplitude.txt';
% fileID = fopen([ipath txt_file],'w');
% fprintf(fileID,'%8s\t%8s\t%26s\t%8s\t%8s\r\n', 'foldername', 'filename','autocorrelation_amplitude','diameter','length');
d = dir(ipath);
isub = [d(:).isdir]; %# returns logical vector
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'..'})) = [];
aamplitude=[];
for jj=1%:length(nameFolds)
filename = dir(fullfile(strcat([ipath nameFolds{jj,1} '\'],'*_1.txt')));
for ii=1:length(filename)

  data = importdata([ipath nameFolds{jj,1} '\' filename(ii).name]);
  % get the y_corr value to do the FFT
  xcdata = data.data(:,4); %pixel -> nm
  x = (xcdata-min(xcdata))*pixelSize/1000;% micron ' dual objective 141 '148
  ycdata = data.data(:,5);
  y = (ycdata-min(ycdata))*pixelSize/1000;
  % x= X(:,7);
  Mlist=save_molecule_list_all(x,y,data,[ipath nameFolds{jj,1} '\' filename(ii).name(1:end-4) '.bin']);
%   [In, imaxes] = list2img(Mlist,'Zsteps',5);
end
end