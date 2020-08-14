clear all;
close all;
ipath = 'M:\new molecule validation\20180513 single color CB1 gnb1 ncam1\msDIV21 NCAM-647\';
pixelSize = 167; 
dbin = 10/1000;  % micr
bb=0;
average_c=[];
txt_file='autocorrelation_amplitude.txt';
fileID = fopen([ipath txt_file],'w');
% fprintf(fileID,'%8s\t%8s\t%26s\t%8s\t%8s\r\n', 'foldername', 'filename','autocorrelation_amplitude','diameter','length');
d = dir(ipath);
isub = [d(:).isdir]; %# returns logical vector
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'..'})) = [];
aamplitude=[];
for jj=1%:length(nameFolds)
filename = dir(fullfile(strcat([ipath nameFolds{jj,1} '\'],'*t_1.txt')));

for ii=1:length(filename)

  data = importdata([ipath nameFolds{jj,1} '\' filename(ii).name]);
  % get the y_corr value to do the FFT
  xcdata = data.data(:,4); %pixel -> nm
  ycdata = data.data(:,5);
  scatter(xcdata,ycdata,'.')
  x = (xcdata-min(xcdata))*pixelSize/1000;    % micron ' dual objective 141 '148
  % x= X(:,7);

  bin=min(x):dbin:max(x);
  L=size(bin,2);
  nx=hist(x,bin);
  nx = (nx-bb); % important to keep the zero in the middle of account
  
  figure(2);
  hist(x,bin); % do FFT from -800 to 800?
  xlabel('in micron');
  
  ind = find(bin>-5 & bin<5); % to plot and calculat which part
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
  title(sprintf('period = % .2f nm',1/f(i)*1000));
  xlabel('in 1/micron');
  
  figure(3)
  plot(bin(ind),nx(ind))
  xlabel('in micron');
  
  tlag = 70;
  tt=nx (ind)';
  c=autocorr(tt,tlag);% time lag 65 
  % cc=xcorr(tt+bb,tt+bb);
  xxx=0:(tlag+1)*10/length(c):(tlag+1)*10-(tlag+1)*10/length(c);
  plot(xxx,c,'k','LineWidth',3);
  x2 = 600;
  y2 = 0.9;
%   txt2 = ['amplitude=' num2str(-(c(11)+c(30))/2+c(20))];
%   text(x2,y2,txt2)
  
  
%   figHandle=figure;
% idx = find(hist{i,j}.BinEdges<900);
% plot(hist{i,j}.BinEdges(idx), acfValues{i, j}(idx),'k','LineWidth',1.25);
  xlabel('Lag (nm)','FontSize', 40)
  ylabel('A. C.','FontSize', 40);
  set(gca,'XTick',0:200:650);
  set(gca,'YTick',-0.5:0.5:1);
  xlim([0 650])
  ylim([-0.5 1])
  box off
  set(gca,'FontSize',40)
  set(gca,'linewidth',3)
  saveas(gcf,[ipath nameFolds{jj,1} '\' filename(ii).name(1:end-4) '_autoCorrelation.eps']);
  average_c=[average_c,c];

end
end

figure;
h=hist(aamplitude,-0.2:0.1:1.2);
  average_total=mean(average_c,2);
  -(average_total(11)+average_total(30))/2+average_total(20)
  std(aamplitude)/sqrt(length(aamplitude))
  length(find(aamplitude>=0.27))/length(aamplitude)
  figure;
  plot(xxx, mean(average_c,2),'k','LineWidth',3);
  xlabel('Lag (nm)','FontSize', 40)
  ylabel('A. C.','FontSize', 40);
  set(gca,'XTick',0:200:650);
  set(gca,'YTick',-0.5:0.5:1);
  xlim([0 650])
  ylim([-0.5 1])
  box off
  set(gca,'FontSize',40)
  set(gca,'linewidth',3)
  saveas(gcf,[ipath nameFolds{jj,1} '\folder_average_autoCorrelation.eps']);