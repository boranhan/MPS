clear all;
close all;
ipath = '\axon_regions\';
pixelSize = 167; 
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
localization_count=[];
for jj=1:length(nameFolds)
filename = dir(fullfile(strcat([ipath nameFolds{jj,1} '\'],'*storm*.bin')));

for ii=[1:length(filename)]

  [MList, memoryMap] = ReadMasterMoleculeList([ipath nameFolds{jj,1} '\' filename(ii).name]);
  % get the y_corr value to do the FFT
  xcdata = MList.xc(MList.c==1); %pixel -> nm
  ycdata = MList.yc(MList.c==1);
  scatter(xcdata,ycdata,'.')
  x = (xcdata-min(xcdata))*pixelSize/1000;    % micron ' dual objective 141 '148
  % x= X(:,7);
  localization_count=[localization_count,length(x)];
  bin=min(x):dbin:max(x);
  L=size(bin,2);
  nx=hist(x,bin);
  nx = (nx-bb); % important to keep the zero in the middle of account
  
  figure(2);
  hist(x,bin); % do FFT from -800 to 800?
  xlabel('in micron');
  
  ind = find(bin>-3 & bin<3); % to plot and calculat which part
  Nind = size(ind,2);
  
  NFFT = 2^nextpow2(L);
  nx1 = zeros(NFFT,1);
  nx1(1:Nind) = nx(ind);
  fx = fft (nx1,NFFT)/L;
  f = 1/dbin/2*linspace(0,1,NFFT/2+1);
  
  [c,i] = max(abs(fx));
 
    
  figure(3)
  zz=2*abs(fx(1:NFFT/2+1));
  plot(f,zz)

  
%   figure(3)
%   plot(bin(ind),nx(ind))
%   xlabel('in micron');
  
  tlag = 80;
  tt=nx (ind)';
  c=autocorr(tt,tlag);% time lag 65 
  % cc=xcorr(tt+bb,tt+bb);
  xxx=0:(tlag+1)*10/length(c):(tlag+1)*10-(tlag+1)*10/length(c);
  plot(xxx,c,'k','LineWidth',1.25);
  x2 = 600;
  y2 = 0.9;
  txt2 = ['amplitude=' num2str(-(c(11)+c(30))/2+c(20))];
  text(x2,y2,txt2)
  %%%%% distance between adjacent rings %%%%
  gf = fit(xxx(10:end)',c(10:end),'gauss3');
  per(ii)=gf.b1; 

  xlabel('Lag (nm)','FontSize', 25)
  ylabel('Autocorrelation','FontSize', 25);
  set(gca,'XTick',0:200:900);
  set(gca,'YTick',-0.5:0.5:1);
  xlim([0 900])
  set(gca,'FontSize',20)
  average_c=[average_c,c];
  
  %%% diameter calculation. 
  ycdata = MList.yc; %pixel -> nm
  y = (ycdata-min(ycdata))*pixelSize/1000;    % micron ' dual objective 141 '148
  bin=min(y):dbin:max(y);
  L=size(bin,2);
  ny=hist(y,bin);
  threshold=0.005*sum(ny);
  consecutive_index=find(ny>threshold);
  consecutive_index=[0,consecutive_index,consecutive_index(end)+2];
    diameter(ii)=(max(diff(find(diff(consecutive_index)~=1)))-1)*dbin;
    aamplitude(ii)=-(c(11)+c(30))/2+c(20);
%   fprintf(fileID,'%10s\t%20s\t%4.2f\t%4.2f\t%4.2f\r\n',nameFolds{jj,1},filename(ii).name,-(c(11)+c(30))/2+c(20),diameter(ii), max(x));
%   
  
end
end
  aamplitude(aamplitude==0)=[];
  mean(localization_count)
  figure;
  h=hist(aamplitude,-0.2:0.2:0.2);
  average_total=mean(average_c,2);
  
  figure;
  plot(xxx, mean(average_c,2),'k','LineWidth',3);
  xlabel('Lag (nm)','FontSize', 40)
  ylabel('A. C.','FontSize', 40);
  set(gca,'XTick',0:200:650);
  set(gca,'YTick',-0.2:0.2:0.2);
  xlim([0 650])
  ylim([-.2 .2])
  box off
  set(gca,'FontSize',40)
  set(gca,'linewidth',3)
  saveas(gcf,[ipath nameFolds{jj,1} '\folder_average_autoCorrelation.eps']);
