clear all;
close all;
ipath = 'K:\20180818 worm\exc-2\';
filename = dir(fullfile(strcat([ipath '\'],'*list.bin')));
for ii=1:length(filename)
    [mList, memoryMap] = ReadMasterMoleculeList([ipath '\' filename(ii).name],'fieldsToLoad',{'xc','yc','zc'},'ZScale',167);
    ROI=[ min(mList.xc) max(mList.xc);min(mList.yc) max(mList.yc)];
    range=[0:100:1000];
    renderedStack = RenderStack([mList.xc mList.yc mList.zc], range, ...
    'gaussianWidth', 0.2, ...
    'ROI', ROI, ...
    'imageScale', 5, ...
    'verbose', true);
    coloredSTORM = Gray2RGB(renderedStack, []);
    imshow(coloredSTORM, []);
    imwrite(coloredSTORM,[ipath '\' filename(ii).name(1:end-4) '.png'])
end
% [renderedImage, edges, parameters] = RenderMList(MList);
% [In, imaxes] = list2img(Mlist,'Zsteps',10);
% Io = STORMcell2img(I,varargin)
% % imshow(renderedImage);
