%% Load Images %%%
clear all;clc;close all;
ProbeFiles = dir('C:\Users\che7oz\Desktop\1x2x3x6x 2493 9312 bcd mrna\w1118\*_c4.jpg');
for i=1:length(ProbeFiles)
ProbeI=imread(ProbeFiles(i).name);
ProbeI=ProbeI(:,:,1);
DapiFiles = dir('C:\Users\che7oz\Desktop\1x2x3x6x 2493 9312 bcd mrna\w1118\*_c3.jpg');
DapiI=imread(DapiFiles(i).name);
DapiI=DapiI(:,:,1);
%% Embryo Identification %%
% figure, imshow(ProbeI); %Shows original Probe image
% figure, imshow(DapiI); %Shows original DAPI stain
Ident=graythresh(DapiI); %Determines the threshold for the outline of the embryo
Ident= .1; %Ident-Ident*(0.5); 
Outline=im2bw(DapiI,Ident); %Converts DAPI to a logical array based on a threshold.
% figure, imshow(Outline);%Shows outline of the embryo
[N,M]=bwlabel(Outline,4); %Labels the objects in 'Outline'
stats=regionprops(N,'all'); %Finds region properties of 'Outline'
embryoArea=[stats.Area]; %Finds the area of the embryo defined by DAPI stain
[embryoSize(i),embryoID]=max(embryoArea); %The object of the largest size is identified as the embryo
N(find(N~=embryoID))=0; %Removes all other smaller objects
NN=(N~=0);
bw2=Outline.*NN;
DapiI=DapiI.*uint8(NN);
% figure, imshow(bw2)
embryoLength = [stats.MajorAxisLength];
EL(i) = max(embryoLength);
embryoWidth = [stats.MinorAxisLength];
EW(i) = max(embryoWidth);
%% Rotate and Crop 
theta=-stats(embryoID).Orientation; %Determines the angle between the major axis of the fitted elipse of the embryo and the x-axis
    I=imrotate(ProbeI,theta); %Orients Probe image horizontally
    ProbeI=imrotate(ProbeI,theta); %Orients 
    W=imrotate(double(bw2),theta);
%        figure,imshow(W);
    BW = edge(W,'canny');
    m0=find(sum(BW,1)>0);
    n0=find(sum(BW,2)>0);
    I1=imcrop(I,[min(m0),min(n0),max(m0)-min(m0),max(n0)-min(n0)]);
    ProbeI1=imcrop(ProbeI,[min(m0),min(n0),max(m0)-min(m0),max(n0)-min(n0)]);
    theta2=180*double(max(max(ProbeI1(:,1:50)))<max(max(ProbeI1(:,max(m0)-min(m0)-48:max(m0)-min(m0)+1)))); % Anterior-posterior
    I1=imrotate(I1,theta2,'crop');
    ProbeI1=imrotate(ProbeI1,theta2,'crop');
    
Outline=imrotate(W, theta);
Outline=imcrop(W,[min(m0),min(n0),max(m0)-min(m0),max(n0)-min(n0)]);
ProbeI1(~Outline>0)=0;
% figure, imshow(ProbeI1)%Shows cropped image
ProbeI1Back = ProbeI1;
ProbeI1Store = ProbeI1;
%

%% Specific Signal
[ProbeIH, ProbeIW] = size(ProbeI1)
DivergeA = ProbeI1([1:ProbeIH], [1:floor(ProbeIW*0.5)])
DivergeP = ProbeI1([1:ProbeIH], [floor(ProbeIW*0.5)+1:ProbeIW])
    Diverge = im2bw(DivergeA, graythresh(DivergeA(DivergeA~=0)));
%     figure, imshow(Diverge);
% Diverge = im2bw(ProbeI1, 53/255); %Histogram Fitted
    DivLabel = bwlabel(Diverge);
    DivProps = regionprops(DivLabel, 'all');
    numberOfBlobs = size(DivProps, 1);
    SignalArea=max([DivProps.Area]);
    DivArea = [DivProps.Area];
    [SignalSize, embryoID]=max(DivArea);
    DivLabel(find(DivLabel~=embryoID))=0;
    DivLabel2=(DivLabel~=0);
%   figure, imshow(DivLabel2);0
    ProbeI1(~DivLabel2>0)=0;
%   figure, imshow(ProbeI1);
    AggInt = sum(ProbeI1);
    SignalSizeStore(i) = SignalSize;
    TotalAggInt(i) = sum(AggInt);
    [count, ex] = imhist(ProbeI1(ProbeI1~=0))
    counts(:,i) = count
% Posterior Region Scanning
[DpSize1, DpSize2] = size(DivergeP)
DivLabelFill = zeros(DpSize1, DpSize2)   
DivLabel2b = [DivLabel2 DivLabelFill]
binaryImage2 = DivLabel2b;
aScan=binaryImage2;
binaryImage3=fliplr(binaryImage2);
% figure, imshow(binaryImage3)
pScan=double(binaryImage3); %Defines posterior scanning blob as a double array
[y1 x1]=size(Outline); %Outputs x and y dimensions of the outline of the cropped image in terms of pixels 
d0=pScan-Outline; % Subtracts the Outline image from the scanning blob
d1=find(d0>0); % Finds where there are still values greater than zero.  
[d2 d3]=size(d1); % Returns matrix dimensions of how many pixels are still greater than zero
d4=1;
DpScan=pScan;%Generates daughter pScan element
while d4>0
        if sum(DpScan(y1,:))<1;
%         figure, imshow(DpScan); % Fig check 1
        DpScan=DpScan'; % Transposes the array 
%         figure, imshow(DpScan); % Fig check 2
        DpScan(:,y1)=[];% Deletes last column of the transposed array; bottom row of the original image
        [s1 s2]=size(DpScan); %Diagnostic
        ZeroFilly=zeros(x1, 1); % Generates a column of zeros with the same x dimension of the original image
        [z1 z2]=size(ZeroFilly) %Diagnostic readout Zeroing
        DpScan=[ZeroFilly DpScan]; % Concatenates the two together
        [s3 s4]=size(DpScan); %Diagnostic
        s3;
        s4;
        DpScan=DpScan'; % Retransposes the array, the original pScan has shifted down one row of pixels. 
        [s5 s6]=size(DpScan); %Diagnostic readout post transposition
        s5;
        s6;
        d0=DpScan-Outline; % Subtracts the Outline image from the scanning blob
%         figure, imshow(d0) % Shows the subtracted image 3
        d1=find(d0>0); % Finds where there are still values greater than zero.  
        [d4 d5]=size(d1);% Finds the size of the matrix which has values greater than zero.
        d4; % Readout of the remaining area difference.
        if sum(DpScan(y1,:))>0
            'xshift'
            pScan(:,1)=[];
            ZeroFillx=zeros(y1,1);
            pScan=[pScan ZeroFillx]; %shift in the x direction
            DpScan=pScan;
            end            
     end            
end
% Anterior Signal and Posterior Background Subtraction Region Overlay
figure, imshow(ProbeI1Store);
hold on;
boundaries = bwboundaries(DpScan);
for k = 1 
thisBoundary = boundaries{k};
plot(thisBoundary(:,2), thisBoundary(:,1), 'r', 'LineWidth', 1);
end
boundaries = bwboundaries(aScan);
for k = 1 
thisBoundary = boundaries{k};
plot(thisBoundary(:,2), thisBoundary(:,1), 'b', 'LineWidth', 1);
end
boundaries = bwboundaries(Outline);
for k = 1 
thisBoundary = boundaries{k};
plot(thisBoundary(:,2), thisBoundary(:,1), 'g', 'LineWidth', 1);
end
title(DapiFiles(i).name)
hold off;
% Background Subtraction
DpScan2 = cast(DpScan, 'logical');
% figure, imshow(DpScan2);
ProbeI1Back(~DpScan2>0)=0;
% figure, imshow(ProbeI1Back);
iBackground = sum(ProbeI1Back);
TotaliBack(i) = sum(iBackground);

end

TotalAggInt = TotalAggInt'
TotaliBack = TotaliBack'
EL = EL*.6431;
EL = EL';
EW = EW*.6431;
EW = EW';
%

Intensity = TotalAggInt - TotaliBack
figure, scatter(EL, Intensity, 'o');
title('EL vs Intensity');
[rho pval] = corr(EL, Intensity)
% Volume
for i = 1:length(ProbeFiles);
Vol(i) = pi*EL(i)*EW(i)*EW(i)/6;
end
figure, scatter(Vol, Intensity);
title('Vol vs Intentsity');
Vol = Vol';
[rho pval] = corr(Vol, Intensity)
%

ELnoise = std(EL)/mean(EL)
Volnoise = std(Vol)/mean(Vol)
Intnoise = std(Intensity)/mean(Intensity)

%
SEL = Intnoise/ELnoise
VEL = Intnoise/Volnoise
%
[x y] = imhist(ProbeI1Store(ProbeI1Store~=0)); set(gca, 'YLim', [0 25000])
%
SignalSizeStore = SignalSizeStore'
IntensityAvg = Intensity./SignalSizeStore
figure, scatter(EL, IntensityAvg)
figure, scatter(embryoSize, IntensityAvg)
[r1 p1] = corrcoef(EL,IntensityAvg)
[r2 p2] = corrcoef(embryoSize, IntensityAvg)
