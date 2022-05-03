function FuncMasterRunFile_LR(FolderPath, FileName)
%%Cat 2 (680), Cat 1 (647) 
%Save beads in a folder in the same folder of sample .dax movie.
%xo, yo is original x, y without drift correction to refer to original
%pixel.
close all
% mainDirMoviePath = 'N:\STORM1_Data\20171118 b2-680-1-1 add-647 storm1\';
% %Dir of original dex file.
% folder = 'Splitted\150\';
% FileName = 'storm  add-647 b2-680_0004';
alist_r = [FolderPath FileName '-R_alist.bin'];
alist_l = [FolderPath FileName '-L_alist.bin'];
mlist_r = [FolderPath FileName '-R_mlist.bin'];
mlist_l = [FolderPath FileName '-L_mlist.bin'];
filehead = [FolderPath FileName];
WarpParameters = [FolderPath 'tform.mat'];
RevWarpParameters = [FolderPath 'tformrev.mat'];
%OutputFileName_ra = [mainDirPath 'Right\' FileName '-R_alist_warp.csv'];
%OutputFileName_rm = [mainDirPath 'Right\' FileName '-R_mlist_warp.csv'];
%OutputFileName_la = [mainDirPath 'Left\' FileName '-L_alist.csv'];
%OutputFileName_lm = [mainDirPath 'Left\' FileName '-L_mlist.csv'];
OutputFileName_l = [FolderPath FileName '-L_alist_assigned.txt'];
OutputFileName_r = [FolderPath FileName '-R_alist_assigned.txt'];
IdenticalDist = 1;
Ratiometriccut = 0.5; 
load(WarpParameters);
load(RevWarpParameters);
rmol_alist = readbinfile_n(alist_r); % for autosaved bin file
lmol_alist = readbinfile_n(alist_l);
rmol_mlist = readbinfile_n(mlist_r);
lmol_mlist = readbinfile_n(mlist_l);
% hold on;
[txa,tya] = tforminv(tform,rmol_alist.x,rmol_alist.y);
rmol_alist.x = txa;
rmol_alist.y = tya;
[txm,tym] = tforminv(tform,rmol_mlist.x,rmol_mlist.y);
rmol_mlist.x = txm;
rmol_mlist.y = tym;

delInd = [];
Ratiototal1 = [];
mL = [];
nR = [];

%% Further calibrate using middle 5000-8000 frames.

for i = 5000:8000 % For total of i frames. Build distance matrix for all molecules start in frame i in lmol_alist with all molecules in frame i in rmol (mlist).
     distMatrix = sqrt((repmat(lmol_alist.x(lmol_alist.frame == i), [1,length(rmol_alist.x(rmol_alist.frame == i))]) ...
        -repmat((rmol_alist.x(rmol_alist.frame == i))', [length(lmol_alist.x(lmol_alist.frame == i)),1])).^2 ...
        +(repmat(lmol_alist.y(lmol_alist.frame == i), [1,length(rmol_alist.y(rmol_alist.frame == i))]) ...
        -repmat((rmol_alist.y(rmol_alist.frame == i))', [length(lmol_alist.y(lmol_alist.frame == i)),1])).^2);
    CompMatrix = distMatrix<2.5; %Initial distance allawance 
    [m,n] = find(CompMatrix==1); % find index n of molecule in rmol that is correspondent to the index m in lmol.
    if size(CompMatrix,1)==1
        m=m';
        n=n';
    end
    lmi = sum(lmol_alist.frame<i);
    rmi = sum(rmol_alist.frame<i);
    if isempty(m)==0
    r = lmol_alist.h(m+lmi)./rmol_alist.h(n+rmi);
    A = sortrows([m,r,n]); % sort the array, so that lower ratio appear first, discard higher ratio (take spots with higher rmol h.)
    [C,ia] = unique(A(:,1));
    %if size(m,1) > size(overlapInd,1)
    % when lmol find 2 molecules in rmol, save the first one with larger h.
    %[C, ia] = unique(m); % when lmol find 2 molecules in rmol, save the first one. 
    M = A(ia,1) + lmi;
    N = A(ia,3) + rmi; %to calibrate using matched pairs.
    mL = [mL; M];
    nR = [nR; N]; %to calibrate using matched pairs.
    end
 end
lmol_match = [lmol_alist.x(mL), lmol_alist.y(mL)];
rmol_match = [rmol_alist.x(nR), rmol_alist.y(nR)];
tform2 = cp2tform(lmol_match,rmol_match,'polynomial', 4);
tform2rev = cp2tform(rmol_match,lmol_match,'polynomial', 4);
[txa,tya] = tforminv(tform2,rmol_alist.x,rmol_alist.y);
rmol_alist.x = txa;
rmol_alist.y = tya;
[txm,tym] = tforminv(tform2,rmol_mlist.x,rmol_mlist.y);
rmol_mlist.x = txm;
rmol_mlist.y = tym;

%% Further calibrate using middle 500-1000 frames. End


%%Ke

% warps molecule list according to tform generated in step 1
% iterates through all matched molecules and gives color assignment

IntensityRatio3Total = [];
% these are used to compute output intensity for each molecule
ComponentLZ=.4; 
CompomentRZ=1-ComponentLZ;

% intensity ratio for color assignment - this does 3 category
% classification - matched molecules above IntensityRatioTh are Cat 0,
% otherwise if below, Cat 1 (647); non-matched molecules are Cat 2 (680)
IntensityRatioTh=2.38;

IfAvergaeCoordinates=0;

IfWriteResult=1;
IfWriteZdiff=0;
IfWriteZCor=0;
IfWriteZConcise=0;

SkipData=1; %Use 1 for non-skip

RejectStart=0;
RejectEnd=0;

IfInv2ndZ=0;

IfRead2ndImg=0;
IfRead1stImg=0;

FrameTol=1;

%X and Y tolerance in pixels for molecule matching
XTol=1; %pix0
YTol=1; %pix !!!!!!!!!!!!!!!!!!!!!

MoreMatchMode=0; %0- discard; 1-median; 2-use same frame only

ZCorCount=0;

ZShift=0;
DisplayFrame=1000;


ZTolPlus= 1480;
ZTolMinus=-1270;

outfile = sprintf('%s_Split-L_H%g.txt',filehead,XTol);
outfileLR = sprintf('%s_Split-LR_H%g.txt',filehead,XTol);

Right=rmol_mlist;


Left=lmol_alist;

SplitData=Left;

NumLeft=Left.N;
NumRight=Right.N;

NowFrameNum=-1000;

NumofNoMatch=0;
NumofMoreMatch=0;
NumofMatch=0;

NumofZMismach=0;
NumofNotSameFrame=0;

% NumLeft=1000;
% MatchList=[];
IndicesC=[];
LeftIntensity=[];
RightIntensity=[];


% ZDiffList=zeros(min(NumLeft,NumRight),1);
% ZCorList=zeros(min(NumLeft,NumRight),4);

RightIndexLower=0;
RightIndexUpper=0;

for i=1+RejectStart:SkipData:NumLeft-RejectEnd
    %SplitData.cat(i)=1;
    
    if SplitData.cat(i)==9
        continue;
    end
    
    CurrentFrameNum = Left.frame(i);
    if CurrentFrameNum~=NowFrameNum
        NowFrameNum=CurrentFrameNum;
        
        if ~(mod(NowFrameNum,DisplayFrame))
            fprintf(1,'Current Frame:%d\n', NowFrameNum)
        end
        
        for RInd=RightIndexLower+1:NumRight
            if Right.frame(RInd)>=CurrentFrameNum-FrameTol
%                 RInd=RInd-1;
                break;
            end
        end
        RightIndexLower=max(1,RInd);
        
        for RInd=RightIndexUpper+1:NumRight
            if Right.frame(RInd)>CurrentFrameNum+FrameTol
                RInd=RInd-1;
                break;
            end
        end
        RightIndexUpper=max(1,RInd);
        if isempty(RightIndexUpper)
            RightIndexUpper=NumRight;
        end
        
        IndicesForRight=RightIndexLower:RightIndexUpper;
%         IndicesForRight1=IndicesForRight';
%         IndicesForRight=find(abs(Right.frame-CurrentFrameNum)<=FrameTol);
        
        if isempty(IndicesForRight)
            NumofNoMatch=NumofNoMatch+1;
            SplitData.cat(i)=2;
            continue %!!!!!!!!!!!!!!!!!!!!!!!!!!!
        end
        
        RCoX=Right.x(IndicesForRight);
        RCoY=Right.y(IndicesForRight);
    end

    if isempty(IndicesForRight)
        NumofNoMatch=NumofNoMatch+1;
        SplitData.cat(i)=2;
        continue %!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end

    
    LcoX=Left.x(i);
    LcoY=Left.y(i);
    
    IndicesMeetX=find(LcoX-XTol<RCoX & RCoX<LcoX+XTol);
    if isempty(IndicesMeetX)
        NumofNoMatch=NumofNoMatch+1;
        SplitData.cat(i)=2;
        continue %!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end
    
    YForMolesMeetX=RCoY(IndicesMeetX);
    
    IndicesMeetY=find(LcoY-YTol<YForMolesMeetX & YForMolesMeetX<LcoY+YTol);
    if isempty(IndicesMeetY)
        NumofNoMatch=NumofNoMatch+1;
        SplitData.cat(i)=2;
        continue %!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end
    
    RMatchedMolecule = IndicesForRight(IndicesMeetX(IndicesMeetY));
    ArMatchedXs=Right.x(RMatchedMolecule);
    ArMatchedYs=Right.y(RMatchedMolecule);
    ArMatchedZs=Right.z(RMatchedMolecule);
    ArMatchedCats=Right.cat(RMatchedMolecule);
        
    if length(IndicesMeetY)>1
        NumofMoreMatch=NumofMoreMatch+1;
        if MoreMatchMode==0
            NumofNoMatch=NumofNoMatch+1;
            SplitData.cat(i)=4;
            continue %!!!!!!!!!!!!!!!!!!!!!!!!!!!
        end
        
        ArMatchedFrames=Right.frame(RMatchedMolecule);
        
        if MoreMatchMode==1
            SplitData.cat(i)=7;
            MachedXuse=median(ArMatchedXs);
            MachedYuse=median(ArMatchedYs);
            MachedZuse=median(ArMatchedZs);
        end
        
        if MoreMatchMode==2
            SplitData.cat(i)=7;
            
            ArSameFrame=find(ArMatchedFrames==CurrentFrameNum);
            if length(ArSameFrame)~=1
                SplitData.cat(i)=8;
                NumofNoMatch=NumofNoMatch+1;
                continue %!!!!!!!!!!!!!!!!!!!!!!!!!!!
            end
            
            MachedXuse=ArMatchedXs(ArSameFrame);
            MachedYuse=ArMatchedYs(ArSameFrame);
            MachedZuse=ArMatchedZs(ArSameFrame);
            MatchCatuse=ArMatchedCats(ArSameFrame);
       
        end
        
    
    else
        %SplitData.cat(i)=1;
              
        MachedXuse=ArMatchedXs;
        MachedYuse=ArMatchedYs;
        MachedZuse=ArMatchedZs;
        MatchCatuse=ArMatchedCats;
    end
    
    if MatchCatuse==9
        SplitData.cat(i)=8;
        continue;
    end
    
    %     MatchList=[MatchList RMatchedMolecule];
    NumofMatch=NumofMatch+1;
    
%     Z1st=Left.z(i);
%     
%     if Z1st>0
%         Red1st=Z1st*CoZMatch1stPos;
%     else
%         Red1st=Z1st*CoZMatch1stNeg;
%     end
%     
%     if MachedZuse>0
%         Red2nd=MachedZuse*CoZMatch2ndPos;
%     else
%         Red2nd=MachedZuse*CoZMatch2ndNeg;
%     end
%     
%     ZDiff=Red1st+Red2nd; %"+" here = "--"
%        
%     ZDiffList(NumofMatch)=ZDiff;
%     
%     if ZDiff>ZTolPlus||ZDiff<ZTolMinus
%         NumofZMismach=NumofZMismach+1;
%         if SplitData.cat(i)==7
%             SplitData.cat(i)=6;
%         else
%             SplitData.cat(i)=8;
%         end
%     end
    
   
%     MergedX=(Left.x(i)+MachedXuse)/2;
%     MergedY=(Left.y(i)+MachedYuse)/2;
%     MergedZ=(Left.z(i)-MachedZuse)/2+ZShift;

%     ZCorCount=ZCorCount+1;
%     ZCorList(ZCorCount,:)=[MergedX MergedY Left.z(i) MachedZuse];

    if IfAvergaeCoordinates
        MergedI=Left.I(i)+Right.I(RMatchedMolecule);

        %The weighting below gives additional weight to the width of PSF over the photon count, as the former is perceived to be more important.
        %Alternative weighting methods can also be implemented.

        ErL=(Left.width(i)*Left.width(i))/sqrt(max(Left.I(i),1));
        ErLX=ErL/Left.Ax(i);
        ErLY=ErL*Left.Ax(i);

        ErR=(Right.width(RMatchedMolecule)*Right.width(RMatchedMolecule))/sqrt(Right.I(RMatchedMolecule));
        ErRX=ErR/Right.Ax(RMatchedMolecule);
        ErRY=ErR*Right.Ax(RMatchedMolecule);

        ComponentLX=ErRX*ErRX/(ErRX*ErRX+ErLX*ErLX);
        ComponentRX=1-ComponentLX;

        ComponentLY=ErRY*ErRY/(ErRY*ErRY+ErLY*ErLY);
        ComponentRY=1-ComponentLY;

        MergedX=Left.x(i)*ComponentLX+MachedXuse*ComponentRX;
        MergedY=Left.y(i)*ComponentLY+MachedYuse*ComponentRY;
        MergedZ=Left.z(i)*ComponentLZ+MachedZuse*CompomentRZ; 
        
        SplitData.x(i)=MergedX;
        SplitData.y(i)=MergedY;
        SplitData.z(i)=MergedZ;
    end    
%     MergedList(NumofMatch,:)=[SplitData.cat(i) MergedX MergedY MergedX MergedY Left.h(i) Left.area(i) Left.width(i) ComponentLX ComponentLY ZDiff MergedI Left.frame(i) Left.length(i) Left.link(i) Left.valid(i) MergedZ MergedZ];
%     LeftIntensity(i)=Left.I(i);
%     RightIntensity(i)=Right.I(RMatchedMolecule);
    
    IntensityRatio=Left.I(i)/Right.I(RMatchedMolecule);
    SplitData.Ax(i)=IntensityRatio;
    SplitData.area(i)=Right.I(RMatchedMolecule);
    if IntensityRatio>IntensityRatioTh
        SplitData.cat(i)=0;
%     elseif IntensityRatio<IntensityRatioUpper & IntensityRatio>IntensityRatioLower
%          SplitData.cat(i)=3;
    elseif IntensityRatio<IntensityRatioTh
        SplitData.cat(i)=6;
    else
        SplitData.cat(i)=9;
    end
    
end



%% Further assign non-matched localizations using pixel intensities. L/R ratio above threshold keep in 2, lower goes to 1 (647),  
IntensityRatioTotal = [];
SignalKernel = [0.077847,0.123317,0.077847; 0.123317,0.195346,0.123317; 0.077847,0.123317,0.077847];
BgKernel = ones(11,11)/96;
BgKernel(4:8, 4:8) = 0;
Fragments = floor(max(lmol_alist.frame)/100);
for j = 1:Fragments
    StartFrameNum = 100*j;
    [movieL, infoFileL] = ReadDax([FolderPath, FileName '-L.dax'], 'startFrame', StartFrameNum, 'endFrame', StartFrameNum+99, 'verbose', false);
    [movieR, infoFileR] = ReadDax([FolderPath, FileName '-R.dax'], 'startFrame', StartFrameNum, 'endFrame', StartFrameNum+99, 'verbose', false);
    for k = StartFrameNum:StartFrameNum+99
        IndCurrentFrame = find(SplitData.frame == k);
        if isempty(IndCurrentFrame)
            continue;
        end
        for i = IndCurrentFrame(1):IndCurrentFrame(end)
            if SplitData.cat(i) ~= 2
            continue;
            end
            if SplitData.width(i) > 800
                SplitData.cat(i)=9;
                continue;
            end
            if SplitData.area(i) > 5000
                SplitData.cat(i)=9;
                continue;
            end
            [txo,tyo] = tforminv(tform2rev,SplitData.xo(i),SplitData.yo(i));
            [txoR,tyoR] = tforminv(tformrev, txo, tyo);
            NonValid = txoR < 6  || txoR > 250 || tyoR < 6 || tyoR > 250 || SplitData.xo(i)< 6 || SplitData.yo(i)< 6 || SplitData.xo(i)> 250 || SplitData.yo(i)> 250;
            if NonValid
            SplitData.cat(i)=9;
            continue;
            end
            CurrentFrame = k-StartFrameNum+1;
            Currentxl = round(SplitData.xo(i));
            Currentyl = round(SplitData.yo(i));
            Currentxr = round(txoR);
            Currentyr = round(tyoR);
            % Note: x,y is switched in movie dimensions. (y,x,frame)
        	backgroundL = sum(sum(BgKernel.*double(movieL((Currentyl-5):(Currentyl+5),(Currentxl-5):(Currentxl+5),CurrentFrame))));
            backgroundR = sum(sum(BgKernel.*double(movieR((Currentyr-5):(Currentyr+5),(Currentxr-5):(Currentxr+5),CurrentFrame))));
            IntensityL = sum(sum(SignalKernel.*double(movieL((Currentyr-1):(Currentyr+1),(Currentxr-1):(Currentxr+1),CurrentFrame))))-backgroundL;
            IntensityR = sum(sum(SignalKernel.*double(movieR((Currentyr-1):(Currentyr+1),(Currentxr-1):(Currentxr+1),CurrentFrame))))-backgroundR;
            IntensityRatio = (IntensityL-IntensityR)/(IntensityL+IntensityR);
            IntensityRatioTotal = [IntensityRatioTotal; IntensityRatio];
            if IntensityRatio > Ratiometriccut
            SplitData.cat(i)=2;
            elseif IntensityRatio < Ratiometriccut
            SplitData.cat(i)=5;
            else
            SplitData.cat(i)=9;
            end
        end
    end
    if ~(mod(StartFrameNum,DisplayFrame))
      display(['Final Step Current Frame:' num2str(StartFrameNum)])
    end
end
fprintf(1,'Writing left mol to file...\n');
SplitDataL = SplitData;
writemoltxt(SplitDataL,outfile,0);
Edge = (-3:0.05:3);
fighandle1=histogram(IntensityRatioTotal, Edge);
saveas(fighandle1, [filehead '_Hist.png']);

%% Further assign non-matched localizations using pixel intensities. L/R ratio above threshold keep in 2, lower goes to 1 (647), End


% a=25;
% figure(1)
% scatter(RightIntensity,LeftIntensity,a,'filled')

%% Process rmol, put good to Cat 3 (647), otherwise to 4.


outfiler = sprintf('%s_Split-R_H%g.txt',filehead,XTol);

FrameTol=0;

Right=lmol_mlist; %switch input data, don't change the rest of the code.


Left=rmol_alist;

SplitData=Left;

NumLeft=Left.N;
NumRight=Right.N;

NowFrameNum=-1000;

NumofNoMatch=0;
NumofMoreMatch=0;
NumofMatch=0;

NumofZMismach=0;
NumofNotSameFrame=0;

% NumLeft=1000;
% MatchList=[];
IndicesC=[];
LeftIntensity=[];
RightIntensity=[];


% ZDiffList=zeros(min(NumLeft,NumRight),1);
% ZCorList=zeros(min(NumLeft,NumRight),4);

RightIndexLower=0;
RightIndexUpper=0;

for i=1+RejectStart:SkipData:NumLeft-RejectEnd
    %SplitData.cat(i)=1;
    
    if SplitData.cat(i)==9
        continue;
    end
    
    CurrentFrameNum = Left.frame(i);
    if CurrentFrameNum~=NowFrameNum
        NowFrameNum=CurrentFrameNum;
        
        if ~(mod(NowFrameNum,DisplayFrame))
            fprintf(1,'Current Frame:%d\n', NowFrameNum)
        end
        
        for RInd=RightIndexLower+1:NumRight
            if Right.frame(RInd)>=CurrentFrameNum-FrameTol
%                 RInd=RInd-1;
                break;
            end
        end
        RightIndexLower=max(1,RInd);
        
        for RInd=RightIndexUpper+1:NumRight
            if Right.frame(RInd)>CurrentFrameNum+FrameTol
                RInd=RInd-1;
                break;
            end
        end
        RightIndexUpper=max(1,RInd);
        if isempty(RightIndexUpper)
            RightIndexUpper=NumRight;
        end
        
        IndicesForRight=RightIndexLower:RightIndexUpper;
%         IndicesForRight1=IndicesForRight';
%         IndicesForRight=find(abs(Right.frame-CurrentFrameNum)<=FrameTol);
        
        if isempty(IndicesForRight)
            NumofNoMatch=NumofNoMatch+1;
            SplitData.cat(i)=3;
            continue %!!!!!!!!!!!!!!!!!!!!!!!!!!!
        end
        
        RCoX=Right.x(IndicesForRight);
        RCoY=Right.y(IndicesForRight);
    end

    if isempty(IndicesForRight)
        NumofNoMatch=NumofNoMatch+1;
        SplitData.cat(i)=3;
        continue %!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end

    
    LcoX=Left.x(i);
    LcoY=Left.y(i);
    
    IndicesMeetX=find(LcoX-XTol<RCoX & RCoX<LcoX+XTol);
    if isempty(IndicesMeetX)
        NumofNoMatch=NumofNoMatch+1;
        SplitData.cat(i)=3;
        continue %!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end
    
    YForMolesMeetX=RCoY(IndicesMeetX);
    
    IndicesMeetY=find(LcoY-YTol<YForMolesMeetX & YForMolesMeetX<LcoY+YTol);
    if isempty(IndicesMeetY)
        NumofNoMatch=NumofNoMatch+1;
        SplitData.cat(i)=3;
        continue %!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end
    
    RMatchedMolecule = IndicesForRight(IndicesMeetX(IndicesMeetY));
    ArMatchedXs=Right.x(RMatchedMolecule);
    ArMatchedYs=Right.y(RMatchedMolecule);
    ArMatchedZs=Right.z(RMatchedMolecule);
    ArMatchedCats=Right.cat(RMatchedMolecule);
        
    if length(IndicesMeetY)>1
        NumofMoreMatch=NumofMoreMatch+1;
        if MoreMatchMode==0
            NumofNoMatch=NumofNoMatch+1;
            SplitData.cat(i)=4;
            continue %!!!!!!!!!!!!!!!!!!!!!!!!!!!
        end
        
        ArMatchedFrames=Right.frame(RMatchedMolecule);
        
        if MoreMatchMode==1
            SplitData.cat(i)=7;
            MachedXuse=median(ArMatchedXs);
            MachedYuse=median(ArMatchedYs);
            MachedZuse=median(ArMatchedZs);
        end
        
        if MoreMatchMode==2
            SplitData.cat(i)=7;
            
            ArSameFrame=find(ArMatchedFrames==CurrentFrameNum);
            if length(ArSameFrame)~=1
                SplitData.cat(i)=8;
                NumofNoMatch=NumofNoMatch+1;
                continue %!!!!!!!!!!!!!!!!!!!!!!!!!!!
            end
            
            MachedXuse=ArMatchedXs(ArSameFrame);
            MachedYuse=ArMatchedYs(ArSameFrame);
            MachedZuse=ArMatchedZs(ArSameFrame);
            MatchCatuse=ArMatchedCats(ArSameFrame);
       
        end
        
    
    else
        %SplitData.cat(i)=1;
              
        MachedXuse=ArMatchedXs;
        MachedYuse=ArMatchedYs;
        MachedZuse=ArMatchedZs;
        MatchCatuse=ArMatchedCats;
    end
    
    if MatchCatuse==9
        SplitData.cat(i)=8;
        continue;
    end
    
    %     MatchList=[MatchList RMatchedMolecule];
    NumofMatch=NumofMatch+1;
    
%     Z1st=Left.z(i);
%     
%     if Z1st>0
%         Red1st=Z1st*CoZMatch1stPos;
%     else
%         Red1st=Z1st*CoZMatch1stNeg;
%     end
%     
%     if MachedZuse>0
%         Red2nd=MachedZuse*CoZMatch2ndPos;
%     else
%         Red2nd=MachedZuse*CoZMatch2ndNeg;
%     end
%     
%     ZDiff=Red1st+Red2nd; %"+" here = "--"
%        
%     ZDiffList(NumofMatch)=ZDiff;
%     
%     if ZDiff>ZTolPlus||ZDiff<ZTolMinus
%         NumofZMismach=NumofZMismach+1;
%         if SplitData.cat(i)==7
%             SplitData.cat(i)=6;
%         else
%             SplitData.cat(i)=8;
%         end
%     end
    
   
%     MergedX=(Left.x(i)+MachedXuse)/2;
%     MergedY=(Left.y(i)+MachedYuse)/2;
%     MergedZ=(Left.z(i)-MachedZuse)/2+ZShift;

%     ZCorCount=ZCorCount+1;
%     ZCorList(ZCorCount,:)=[MergedX MergedY Left.z(i) MachedZuse];

    if IfAvergaeCoordinates
        MergedI=Left.I(i)+Right.I(RMatchedMolecule);

        %The weighting below gives additional weight to the width of PSF over the photon count, as the former is perceived to be more important.
        %Alternative weighting methods can also be implemented.

        ErL=(Left.width(i)*Left.width(i))/sqrt(max(Left.I(i),1));
        ErLX=ErL/Left.Ax(i);
        ErLY=ErL*Left.Ax(i);

        ErR=(Right.width(RMatchedMolecule)*Right.width(RMatchedMolecule))/sqrt(Right.I(RMatchedMolecule));
        ErRX=ErR/Right.Ax(RMatchedMolecule);
        ErRY=ErR*Right.Ax(RMatchedMolecule);

        ComponentLX=ErRX*ErRX/(ErRX*ErRX+ErLX*ErLX);
        ComponentRX=1-ComponentLX;

        ComponentLY=ErRY*ErRY/(ErRY*ErRY+ErLY*ErLY);
        ComponentRY=1-ComponentLY;

        MergedX=Left.x(i)*ComponentLX+MachedXuse*ComponentRX;
        MergedY=Left.y(i)*ComponentLY+MachedYuse*ComponentRY;
        MergedZ=Left.z(i)*ComponentLZ+MachedZuse*CompomentRZ; 
        
        SplitData.x(i)=MergedX;
        SplitData.y(i)=MergedY;
        SplitData.z(i)=MergedZ;
    end    
%     MergedList(NumofMatch,:)=[SplitData.cat(i) MergedX MergedY MergedX MergedY Left.h(i) Left.area(i) Left.width(i) ComponentLX ComponentLY ZDiff MergedI Left.frame(i) Left.length(i) Left.link(i) Left.valid(i) MergedZ MergedZ];
%     LeftIntensity(i)=Left.I(i);
%     RightIntensity(i)=Right.I(RMatchedMolecule);
    
    IntensityRatio=Left.I(i)/Right.I(RMatchedMolecule);
    SplitData.Ax(i)=IntensityRatio;
    SplitData.area(i)=Right.I(RMatchedMolecule);
     
    if IntensityRatio>IntensityRatioTh
        SplitData.cat(i)=0;
%     elseif IntensityRatio<IntensityRatioUpper & IntensityRatio>IntensityRatioLower
%          SplitData.cat(i)=3;
    elseif IntensityRatio<IntensityRatioTh
        SplitData.cat(i)=1;
    else
        SplitData.cat(i)=9;
    end
    
end
fprintf(1,'Writing right mol to file...\n');
writemoltxt(SplitData,outfiler,0);



