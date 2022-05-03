clear all;
close all;
ipath = 'D:\turnkey\fasciculation\synapse density\wt\';
pixelSize = 133; 
dbin = 10/1000;  % micr
bb=0;
average_c=[];
txt_file='autocorrelation_amplitude.txt';
fileID = fopen([ipath txt_file],'w');
fprintf(fileID,'%8s\t%8s\t%26s\t%8s\t%8s\r\n', 'foldername', 'filename','autocorrelation_amplitude','diameter','length');
d = dir(ipath);
isub = [d(:).isdir]; %# returns logical vector
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'..'})) = [];
aamplitude=[];
per=[];
sensitivity=0.6;

filename = dir(fullfile(strcat([ipath '\'],'*.bin')));

for ii=1:2:length(filename)
  
  [MList, memoryMap] = ReadMasterMoleculeList([ipath '\' filename(ii).name]);
  % get the y_corr value to do the FFT
  xcdata_1 = MList.xc(find(MList.c==1)); %pixel -> nm
  ycdata_1 = MList.yc(find(MList.c==1));
  
  [MList, memoryMap] = ReadMasterMoleculeList([ipath '\' filename(ii+1).name]);
  % get the y_corr value to do the FFT
  xcdata_2 = MList.xc(find(MList.c==1)); %pixel -> nm
  ycdata_2 = MList.yc(find(MList.c==1));
  
  D=pdist2([xcdata_1,ycdata_1],[xcdata_2,ycdata_2]);
  [col,row]=find(D<1);
  Im = readSPE([ipath '\'], [filename(ii).name(1:end-9) '.spe']);
  
  Im=double(Im)-median(median(double(Im)));
  Im_contrast = imadjust(Im/max(max(Im)),[0,1],[]);
  image_binary = imbinarize(Im_contrast, 'adaptive','Sensitivity',sensitivity);
  image_binary =(bwareaopen(image_binary, 50));
%   image_binary = imdilate(image_binary, strel('disk', 3))-image_binary;
  imshow(image_binary);hold on
  plot(xcdata_1,ycdata_1,'.');
  plot(xcdata_2,ycdata_2,'.');hold off;
  
  total_area_soma=length(find(image_binary==1))*(pixelSize/1000)^2;
  colo_num(ceil(ii/2))=min([length(unique(col)),length(unique(row))])/total_area_soma;
end


