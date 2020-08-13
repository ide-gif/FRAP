%% 20枚のtiffをimagedata(cell型)に格納
tifffiles = dir('C:\Users\nryu-\OneDrive\ドキュメント\MATLAB\HQ\*.tiff');
numfiles = 20;
imagedata = cell(1,numfiles);
dx = 0.225; %um
dt = 1;     %sec
T = zeros(1,numfiles);  %後でグラフをプロットするために時間のベクトルをつくる
for k = 2:numfiles
    T(k) = T(k-1)+dt;
end
 
for k = 1:numfiles
    imagedata{k} = imread(tifffiles(k).name);   %ここに画像データが入る
end
 
%3*3の範囲に対してメディアンフィルタをかける
for k = 1:numfiles
    imagedata{k} = medfilt2(imagedata{k});
end
figure
imshow(imagedata{1});
 
%% 細胞全体ROIの抽出
%エッジ抽出をしてからROIを抽出してみる
[~,threshold] = edge(imagedata{1},'sobel');
fudgefactor = 0.4; 
ROI1 = edge(imagedata{1},'sobel',threshold * fudgefactor);  %エッジ抽出
 
figure
imshow(ROI1);
se90 = strel('line',2,90);
se0 = strel('line',2,0);
ROI1 = imdilate(ROI1,[se90,se0]);       
figure
imshow(ROI1);
 
ROI1 = imfill(ROI1,'holes');    %内部の塗りつぶし
figure
imshow(ROI1);
 
seD = strel('diamond',1);
ROI1 = imerode(ROI1,seD);
ROI1 = imerode(ROI1,seD);       %平滑化
figure
imshow(ROI1);
 
%このままだと細かいノイズもROIとして認識されてしまう。
%面積が一番大きいROIだけを抽出
cc = bwconncomp(ROI1,8);
numPixels = cellfun(@numel,cc.PixelIdxList);
[~,idx] = max(numPixels);
for k = 1:cc.NumObjects
    if k~=idx
        ROI1(cc.PixelIdxList{k}) = 0;   %このROIにより細胞全体の骨格が同定
    end
end
 
ROIoutline = bwperim(ROI1);
imagedata{1}(ROIoutline) = 255; 
figure
imshow(imagedata{1})
 
%% ROIの蛍光強度を計算
allint=zeros(1,numfiles);
 
%元画像とROIの画素値を乗算して値を求めた。
for k=1:numfiles
    mulim=immultiply(imagedata{k},ROI1);
    allint(k) = mean(mulim(mulim>0));
end
 
figure
plot(T,allint);
xlabel('Time (sec)','FontSize',16);
ylabel('fluorescent intensity','FontSize',16);
 
%% どのタイミングで照射したか調べて照射範囲を特定
D = diff(allint);
[~,T_Laser] = min(D);   %微分値が最も小さいときがレーザーを照射したときとみなせる
T_Laser = T_Laser + 1;  %光を照射したフレーム
 
%T_Laserで蛍光輝度が著しく減少した領域を調べる
%照射する前の画像から照射直後の画像を減算すればレーザーによって退色した領域がわかるはず。
A = immultiply(imagedata{T_Laser-1},ROI1);
ROI2 = immultiply(imagedata{T_Laser},ROI1);
ROI2 = imsubtract(A,ROI2);
figure
imshow(ROI2);
 
ROI2 = imbinarize(ROI2);    %二値化処理
figure
imshow(ROI2);
 
seD = strel('disk',4);
ROI2 = imerode(ROI2,seD);   
ROI2 = imerode(ROI2,seD);   %平滑化
figure
imshow(ROI2);
 
%FRAPから拡散を求めるには照射領域を円とする必要があるので先ほど求めた領域の中心と半径からマスクを作成する
Rmax=60;
Rmin=10;
 
[centersBright, radiiBright] = imfindcircles(ROI2,[Rmin,Rmax],...
    'ObjectPolarity','bright');                             %中心と半径を求める

radiiBright = radiiBright +5;
ROI2 = createCirclesMask(ROI2,centersBright,radiiBright);   %中心と半径からマスクを作成する関数を用いた
 
ROIoutline = bwperim(ROI2);
imagedata{T_Laser}(ROIoutline) = 255; 
figure
imshow(imagedata{T_Laser})
 
 
%% FRAPの時系列データを取得して拡散係数を求める。
frapint=zeros(1,numfiles);
 
for k=1:numfiles
    mulim=immultiply(imagedata{k},ROI2);
    frapint(k) = mean(mulim(mulim>0));
end
 
frapint = frapint - frapint(T_Laser);
frapint = frapint/frapint(1);

%近似式からFRAP信号が1/2に回復する時間tauを求める
xdata = T(T_Laser:numfiles)-(T_Laser-1);
ydata = frapint(T_Laser:numfiles)-frapint(T_Laser);
 
x0 =[ydata(end),(numfiles-T_Laser)/2];  %パラメータの初期値
 
f = @(x,xdata) x(1)*(1-exp(-x(2)*xdata)); %モデル式
 
x = lsqcurvefit(f,x0,xdata,ydata);  %x(1)はFmax,x(2)はtauを表す

figure
plot(T,frapint);
xlabel('Time (sec)','FontSize',16);
ylabel('fluorescent intensity','FontSize',16);
yline(1,'LineWidth',2,'Color','r');
yline(0,'LineWidth',2,'Color','r');
yline(x(1),'-',{'Fmax'},...
    'LineWidth',2,'Color','r');

figure
times = linspace(0,100);
plot(xdata,ydata,'ko',times,f(x,times),'b-');
legend('Data','Fitted line','location','southeast',...
    'FontSize',12);
xlabel('Time (sec)','FontSize',16);
ylabel('fluorescent intensity','FontSize',16);
 
 
%拡散係数を求める
R = dx * radiiBright;   %照射半径[um]
tau = -log(0.5)/x(2);             %[sec]
 
D = 0.224*R^2 / tau;    %[um*um / sec]
