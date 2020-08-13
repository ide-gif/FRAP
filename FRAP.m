%% 20����tiff��imagedata(cell�^)�Ɋi�[
tifffiles = dir('C:\Users\nryu-\OneDrive\�h�L�������g\MATLAB\HQ\*.tiff');
numfiles = 20;
imagedata = cell(1,numfiles);
dx = 0.225; %um
dt = 1;     %sec
T = zeros(1,numfiles);  %��ŃO���t���v���b�g���邽�߂Ɏ��Ԃ̃x�N�g��������
for k = 2:numfiles
    T(k) = T(k-1)+dt;
end
 
for k = 1:numfiles
    imagedata{k} = imread(tifffiles(k).name);   %�����ɉ摜�f�[�^������
end
 
%3*3�͈̔͂ɑ΂��ă��f�B�A���t�B���^��������
for k = 1:numfiles
    imagedata{k} = medfilt2(imagedata{k});
end
figure
imshow(imagedata{1});
 
%% �זE�S��ROI�̒��o
%�G�b�W���o�����Ă���ROI�𒊏o���Ă݂�
[~,threshold] = edge(imagedata{1},'sobel');
fudgefactor = 0.4; 
ROI1 = edge(imagedata{1},'sobel',threshold * fudgefactor);  %�G�b�W���o
 
figure
imshow(ROI1);
se90 = strel('line',2,90);
se0 = strel('line',2,0);
ROI1 = imdilate(ROI1,[se90,se0]);       
figure
imshow(ROI1);
 
ROI1 = imfill(ROI1,'holes');    %�����̓h��Ԃ�
figure
imshow(ROI1);
 
seD = strel('diamond',1);
ROI1 = imerode(ROI1,seD);
ROI1 = imerode(ROI1,seD);       %������
figure
imshow(ROI1);
 
%���̂܂܂��ƍׂ����m�C�Y��ROI�Ƃ��ĔF������Ă��܂��B
%�ʐς���ԑ傫��ROI�����𒊏o
cc = bwconncomp(ROI1,8);
numPixels = cellfun(@numel,cc.PixelIdxList);
[~,idx] = max(numPixels);
for k = 1:cc.NumObjects
    if k~=idx
        ROI1(cc.PixelIdxList{k}) = 0;   %����ROI�ɂ��זE�S�̂̍��i������
    end
end
 
ROIoutline = bwperim(ROI1);
imagedata{1}(ROIoutline) = 255; 
figure
imshow(imagedata{1})
 
%% ROI�̌u�����x���v�Z
allint=zeros(1,numfiles);
 
%���摜��ROI�̉�f�l����Z���Ēl�����߂��B
for k=1:numfiles
    mulim=immultiply(imagedata{k},ROI1);
    allint(k) = mean(mulim(mulim>0));
end
 
figure
plot(T,allint);
xlabel('Time (sec)','FontSize',16);
ylabel('fluorescent intensity','FontSize',16);
 
%% �ǂ̃^�C�~���O�ŏƎ˂��������ׂďƎ˔͈͂����
D = diff(allint);
[~,T_Laser] = min(D);   %�����l���ł��������Ƃ������[�U�[���Ǝ˂����Ƃ��Ƃ݂Ȃ���
T_Laser = T_Laser + 1;  %�����Ǝ˂����t���[��
 
%T_Laser�Ōu���P�x�����������������̈�𒲂ׂ�
%�Ǝ˂���O�̉摜����Ǝ˒���̉摜�����Z����΃��[�U�[�ɂ���đސF�����̈悪�킩��͂��B
A = immultiply(imagedata{T_Laser-1},ROI1);
ROI2 = immultiply(imagedata{T_Laser},ROI1);
ROI2 = imsubtract(A,ROI2);
figure
imshow(ROI2);
 
ROI2 = imbinarize(ROI2);    %��l������
figure
imshow(ROI2);
 
seD = strel('disk',4);
ROI2 = imerode(ROI2,seD);   
ROI2 = imerode(ROI2,seD);   %������
figure
imshow(ROI2);
 
%FRAP����g�U�����߂�ɂ͏Ǝ˗̈���~�Ƃ���K�v������̂Ő�قǋ��߂��̈�̒��S�Ɣ��a����}�X�N���쐬����
Rmax=60;
Rmin=10;
 
[centersBright, radiiBright] = imfindcircles(ROI2,[Rmin,Rmax],...
    'ObjectPolarity','bright');                             %���S�Ɣ��a�����߂�

radiiBright = radiiBright +5;
ROI2 = createCirclesMask(ROI2,centersBright,radiiBright);   %���S�Ɣ��a����}�X�N���쐬����֐���p����
 
ROIoutline = bwperim(ROI2);
imagedata{T_Laser}(ROIoutline) = 255; 
figure
imshow(imagedata{T_Laser})
 
 
%% FRAP�̎��n��f�[�^���擾���Ċg�U�W�������߂�B
frapint=zeros(1,numfiles);
 
for k=1:numfiles
    mulim=immultiply(imagedata{k},ROI2);
    frapint(k) = mean(mulim(mulim>0));
end
 
frapint = frapint - frapint(T_Laser);
frapint = frapint/frapint(1);

%�ߎ�������FRAP�M����1/2�ɉ񕜂��鎞��tau�����߂�
xdata = T(T_Laser:numfiles)-(T_Laser-1);
ydata = frapint(T_Laser:numfiles)-frapint(T_Laser);
 
x0 =[ydata(end),(numfiles-T_Laser)/2];  %�p�����[�^�̏����l
 
f = @(x,xdata) x(1)*(1-exp(-x(2)*xdata)); %���f����
 
x = lsqcurvefit(f,x0,xdata,ydata);  %x(1)��Fmax,x(2)��tau��\��

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
 
 
%�g�U�W�������߂�
R = dx * radiiBright;   %�Ǝ˔��a[um]
tau = -log(0.5)/x(2);             %[sec]
 
D = 0.224*R^2 / tau;    %[um*um / sec]
