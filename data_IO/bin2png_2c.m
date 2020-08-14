clear all;
close all;
ipath = 'R:\20180513 myh10 genetex rab 647 b2 680\DIV20\split\axon_regions\';
filename = dir(fullfile(strcat([ipath '\'],'*.bin')));
for ii=1:length(filename)
    [mList, memoryMap] = ReadMasterMoleculeList([ipath '\' filename(ii).name],'fieldsToLoad',{'xc','yc','z', 'c'},'ZScale',167);
    ROI=[ min(mList.yc (mList.c==2)) max(mList.yc(mList.c==2));min(mList.xc(mList.c==2)) max(mList.xc(mList.c==2))];
    range=[-10:1:0];
    renderedStack_1 = RenderMList([mList.xc(mList.c==2) mList.yc(mList.c==2)],  ...
    'gaussianWidth', 0.1, ...
    'ROI', ROI, ...
    'imageScale', 10);
    renderedStack_2 = RenderMList([mList.xc(mList.c==1) mList.yc(mList.c==1)],  ...
    'gaussianWidth', 0.1, ...
    'ROI', ROI, ...
    'imageScale', 10);
    imbw_1=imbinarize(renderedStack_1);
    imbw_2=imbinarize(renderedStack_2);
    ratio_1(ii)=length(find(imbw_1==1&imbw_2==1))/length(find(imbw_1==1));
    ratio_2(ii)=length(find(imbw_1==1&imbw_2==1))/length(find(imbw_2==1));
    renderedStack=cat(3,renderedStack_1,renderedStack_2);
    renderedStack=cat(3,renderedStack,zeros(size(renderedStack_1,1),size(renderedStack_1,2)));
    imshow(renderedStack, []);
    imwrite(coloredSTORM,[ipath '\' filename(ii).name(1:end-4) '.png'])
end

 'type1 ratio_1'
  aamplitude_random=ratio_1;%(randperm(length(aamplitude)));
  replicate=[mean(aamplitude_random(1:round(length(aamplitude_random)/3))),mean(aamplitude_random(round(length(aamplitude_random)/3)+1:round(2*length(aamplitude_random)/3))),mean(aamplitude_random(round(2*length(aamplitude_random)/3)+1:end))];
  ste=std(replicate)/sqrt(3)
  average=mean(replicate)
  replicate
  

  
  'type2 ratio_1'
  replicate=[mean(aamplitude_random(1:3:end)),mean(aamplitude_random(2:3:end)),mean(aamplitude_random(3:3:end))];
  ste=std(replicate)/sqrt(3)
  replicate
  average=mean(replicate)
  
     'type1 ratio_2'
  aamplitude_random=ratio_2;%(randperm(length(aamplitude)));
  replicate=[mean(aamplitude_random(1:round(length(aamplitude_random)/3))),mean(aamplitude_random(round(length(aamplitude_random)/3)+1:round(2*length(aamplitude_random)/3))),mean(aamplitude_random(round(2*length(aamplitude_random)/3)+1:end))];
  ste=std(replicate)/sqrt(3)
  average=mean(replicate)
  replicate
  
    'type2 ratio_2'
  replicate=[mean(aamplitude_random(1:3:end)),mean(aamplitude_random(2:3:end)),mean(aamplitude_random(3:3:end))];
  ste=std(replicate)/sqrt(3)
  replicate
  average=mean(replicate)
% [renderedImage, edges, parameters] = RenderMList(MList);
% [In, imaxes] = list2img(Mlist,'Zsteps',10);
% Io = STORMcell2img(I,varargin)
% % imshow(renderedImage);
