clear all;
close all;
ipath = 'D:\turnkey\20190807 b2 expression level\b2kd\';
pixelSize = 133; 
calculate_soma=0;
sensitivity=0.5;
clims=[1000,5000];
intensity_axon=[];

filename = dir(fullfile(strcat([ipath '\'],'*.spe')));
for ii=1:length(filename)
        Im = readSPE([ipath '\'], filename(ii).name);
        Im=mean(Im,3);
%           imagesc(Im);
        image_binary = imbinarize(Im/max(Im(:)), 'adaptive', 'sensitivity');
        figure;
         imagesc(image_binary);
        image_intensity=reshape(Im(image_binary==1),[],1)-min(Im(:));
        intensity_axon=[intensity_axon,mean(image_intensity)];
               
end
intensity_axon=intensity_axon;
mean(intensity_axon)
std(intensity_axon)