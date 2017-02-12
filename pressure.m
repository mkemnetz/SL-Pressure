%%
clearvars;
close all;
clc;

%%
Pt = 14.61;

%%
% PressureVar = hdf5read('/Volumes/LaCie/MATLAB/Research/Shear Layer Pressure/Data/RUN 4/2_5.hws','/wfm_group0/axes/axis1/data_vector/data');



fileNames = {'3_3', '3_6', '3_4', '2_7', '2_1', '2_4', '2_3', '2_6', ...
    '2_5', '2_2', '4_0', '4_2', '4_3', '4_1', '3_2', '3_0', '3_5'};
numFiles = length(fileNames);

dataMatrix = zeros(numFiles, 400000);
for n = 1:numFiles
   str = fileNames{n};
   currentFilePath = sprintf('/Volumes/LaCie/MATLAB/Research/Shear Layer Pressure/Data/RUN 4/%s.hws',str);
   dataMatrix(n, :) = hdf5read(currentFilePath,'/wfm_group0/axes/axis1/data_vector/data')';

end

pressureMatrix = 14.61+((dataMatrix.*10.*1000)./10.2);

%%
pressureMatrixStruc = zeros(22, 400000);

pressureMatrixStruc(1, :) = pressureMatrix(1, :);
pressureMatrixStruc(2, :) = NaN; 
pressureMatrixStruc(3, :) = pressureMatrix(2, :);
pressureMatrixStruc(4, :) = NaN;
pressureMatrixStruc(5, :) = pressureMatrix(3, :);
pressureMatrixStruc(6, :) = NaN;

pressureMatrixStruc(7,  :) = pressureMatrix(4, :);
pressureMatrixStruc(8,  :) = pressureMatrix(5, :);
pressureMatrixStruc(9,  :) = pressureMatrix(6, :);
pressureMatrixStruc(10, :) = pressureMatrix(7, :);
pressureMatrixStruc(11, :) = pressureMatrix(8, :);
pressureMatrixStruc(12, :) = pressureMatrix(9, :);
pressureMatrixStruc(13, :) = pressureMatrix(10, :);
pressureMatrixStruc(14, :) = pressureMatrix(11, :);
pressureMatrixStruc(15, :) = pressureMatrix(12, :);
pressureMatrixStruc(16, :) = pressureMatrix(13, :);
pressureMatrixStruc(17, :) = pressureMatrix(14, :);
pressureMatrixStruc(18, :) = pressureMatrix(15, :);

pressureMatrixStruc(19, :) = NaN;
pressureMatrixStruc(20, :) = pressureMatrix(16, :);
pressureMatrixStruc(21, :) = NaN;
pressureMatrixStruc(22, :) = pressureMatrix(17, :);
% insertNanIndex = [1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1];
% insertValue    = (1-insertNanIndex)./0;
% 
% b_tmp = [pressureMatrix; insertValue];
% 
% 
% a=[1 2 3 5 6 7 9 10 13 14];
% insertNanIndex = [0 diff(a)>1];
% insertValue = (1-insertNanIndex)./0;
% b_tmp = [a; insertValue];
% b = b_tmp(:)';
% b(isinf(b)) = [];

Pcourse = fillmissing(pressureMatrixStruc,'linear');

%%
load('/Volumes/LaCie/MATLAB/Research/Shear Layer Pressure/Data/RUN 4/FORCING_4.mat');

%%
x = 1:9000000;
v1 = data.forcingSignal';
v2 = data.trigger';

xq = 1:0.5:9000000;
vq1 = interp1(x,v1,xq);
vq2 = interp1(x,v2,xq);
% figure
% vq1 = interp1(x,v,xq);
% plot(x,v,'o',xq,vq1,':.');

yF = downsample(vq1,15);
yT = downsample(vq2,15);








