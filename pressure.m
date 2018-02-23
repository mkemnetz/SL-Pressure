%%
clearvars;
close all;
clc;

%%
Pt = 14.48;

%%
% PressureVar = hdf5read('/Volumes/LaCie/MATLAB/Research/Shear Layer Pressure/Data/RUN 4/2_5.hws','/wfm_group0/axes/axis1/data_vector/data');



% fileNames = {'3_3', '3_6', '3_4', '2_7', '2_1', '2_4', '2_3', '2_6', ...
%     '2_5', '2_2', '4_0', '4_2', '4_3', '4_1', '3_2', '3_0', '3_5'};

fileNames = {'3_3', '3_4', '2_7', '2_1', '2_4', '2_3', ...
    '2_5', '2_2', '4_1', '3_2', '4_2', '3_0', '2_0', '4_0', '5_6', '5_3'};
fileBase = '/Volumes/LaCie/MATLAB/Research/Shear Layer Pressure/Data/2-27-2017/Force4/';
numFiles = length(fileNames);

dataMatrix = zeros(numFiles, 300000);
for n = 1:numFiles
   str = fileNames{n};
   currentFilePath = sprintf([fileBase '%s.hws'],str);
   dataMatrix(n, :) = hdf5read(currentFilePath,'/wfm_group0/axes/axis1/data_vector/data')';

end

pressureMatrix = Pt+((dataMatrix.*10.*1000)./10.2);

%%
pressureMatrixStruc = zeros(21, 300000);

pressureMatrixStruc(1, :) = pressureMatrix(1, :);
pressureMatrixStruc(2, :) = NaN; 
pressureMatrixStruc(3, :) = pressureMatrix(2, :);
pressureMatrixStruc(4, :) = NaN;
pressureMatrixStruc(5, :) = pressureMatrix(3, :);
pressureMatrixStruc(6, :) = NaN;

pressureMatrixStruc(7,  :) = pressureMatrix(4, :);
pressureMatrixStruc(8,  :) = pressureMatrix(5, :);
pressureMatrixStruc(9,  :) = pressureMatrix(6, :);
% pressureMatrixStruc(10, :) = pressureMatrix(7, :);
pressureMatrixStruc(10, :) = NaN;
pressureMatrixStruc(11, :) = pressureMatrix(8, :);
pressureMatrixStruc(12, :) = pressureMatrix(9, :);
pressureMatrixStruc(13, :) = pressureMatrix(10, :);
pressureMatrixStruc(14, :) = pressureMatrix(11, :);
% pressureMatrixStruc(15, :) = pressureMatrix(12, :);
pressureMatrixStruc(15, :) = NaN;
pressureMatrixStruc(16, :) = pressureMatrix(13, :);
pressureMatrixStruc(17, :) = pressureMatrix(14, :);
pressureMatrixStruc(18, :) = NaN;

pressureMatrixStruc(19, :) = pressureMatrix(15, :);
pressureMatrixStruc(20, :) = NaN;
pressureMatrixStruc(21, :) = pressureMatrix(16, :);

Pcourse = fillmissing(pressureMatrixStruc,'linear');

%%
currentFilePath = [fileBase 'forcing.hws'];
forcingSignal = hdf5read(currentFilePath,'/wfm_group0/axes/axis1/data_vector/data')';

%% Interpolate onto Finer Mesh
% The y value we want to start at.  This was chosen by Piyush Ranade on pg
% 70 of his dissertation.  It was chosen because it "did a good job
% regularizing pressure" but it should generally be arbitrary.  My value is
% different here because this voltage is scaled to not saturate the NI DAQ.
%  my starting point here and Piyush's should be very similar when scaled.
startY = -0.0075;

% The sampled forceing signal doesnt go from [-1, 1].  The NI DAQ didnt
% have that range (would have saturated).  I multiply by this constant to
% put it on the range [-1, 1].  I actually don't ever use this anymore but
% i'll leave it here incase i want it in the future.
reScale    = 44.5300;

% Now determine the time vector that these samples were taken at
fsamp      = 20000;
x          = 1:length(forcingSignal); 
t          =(1/fsamp).*x;

% Start generating query mesh.  This is the mesh of y values that I want to
% query.  I need to pick these to get good resolution in the peaks.  I
% chose them by hand.
% y_query = [0.0225, 0.0224, 0.0223, 0.0222, 0.0221, 0.0220, 0.0219, ...
%     0.0215, 0.0210, 0.0205, 0.01950, 0.01949, 0.01948, 0.01947, ...
%     0.01946, 0.01945, 0.01930, 0.01920, 0.01910, 0.0190, 0.0189, ...
%     0.0188, 0.0187, 0.0186, 0.0185, 0.0184, 0.0183, 0.0182, 0.0181, ...
%     0.0180];
y_query = [0.0222, 0.0221, 0.0220, 0.0219, ...
    0.0215, 0.0210, 0.0205, 0.01920, 0.01910, 0.0190, 0.0189, ...
    0.0188, 0.0187, 0.0186, 0.0185, 0.0184, 0.0183, 0.0182, 0.0181, ...
    0.0180];

% In the region not near the peak I just do a linear spacing
a         = linspace(0, 0.0180, 5);

% Now I can join everything together and sort
y_query   = [y_query a];
y_query   = sort([y_query -1.*y_query startY]);

% here I remove multiple query points in case I added them.
[~,index] = unique(y_query,'first');
y_query   = y_query(sort(index));

% I use this function crossing.m from MATLAB Central to determine the
% interpolation mesh.  The vector added to time{n} determines the points
% I want to interpolate onto.  
time = cell(1, length(y_query));
for n = 1:length(y_query)
    [~, time{n}] = crossing(forcingSignal-y_query(n),t);
end


% Now I take all those vectors added to time{n} and turn them into one
% vector and sort.  this query mesh represents all the points in time
% corresponding to the y values selected in y_query
tq = sort(cell2mat(time));

% Interpolate the pressure data onto the same mesh chosen by the forcing
% signal
Pcourse       = Pcourse(1:17, :);
PcourseInterp = zeros(size(Pcourse, 1), length(tq)); 

for i = 1:size(Pcourse, 1)
    v = Pcourse(i, :);
    PcourseInterp(i, :) = interp1(t,v,tq, 'spline');
end

% Interpolate the forcing signal onto the chosen mesh
v = forcingSignal;
forcingSignalInterp = interp1(t,v,tq, 'linear');

%% Remove excess Data at start
% I want to start in the same place in phase that Piyush did.  This is
% given in the startY above (the y value pertaining to phi = 0).  So I
% first find where in the y_query vector staryY is.
n  = find(y_query == startY);

% Now I use a placeholder variable to give me all the instances in time
% that my forcing signal is equal to startY.
tt = time{n};

% I really want the indices where my query times (tq) are equal to tt(2).
% I know this should be tt(2) as opposed to tt(k) where k is some number
% from plotting the experimental data.  I could start sometime later in time
% if I wished.  
start         = find(tq == tt(2));

% Now I trim off the front of all my data to start where I defined 0 phase.
T             = tq(start:end);
force_shift   = forcingSignalInterp(start:end);
Pcourse_shift = PcourseInterp(:, start:end);

%% Determine wavelength for phase-locking
% Now from looking at the data we know that this function goes from 
% 0 -> 4pi every 8 indices in t0.
lambda = 4;

% Make t0 evenly divisible by lambda by removing data from the end
tt = tt(2:end-mod(length(tt), lambda));

% Pick out from t0 the locations of the starting point for phase averaging.
%  The arbitrarily chosen 0pi.
startingPhase = tt(1:lambda:end);

% Get T indices where you have the starting points for phase averaging
T_indexes = find(ismember(T, startingPhase));

%% phaselock
a       = T_indexes(2:end)-T_indexes(1:end-1);
lambda2 = mean(a)+1;

if (std(a) ~= 0)
    error('ERROR: your phase locking windows are not of equal length');
end

endBuffer = mod(size(force_shift, 2), lambda2);

force     = force_shift(1:end-endBuffer);
pres      = Pcourse_shift(:, 1:end-endBuffer);
time      = T(1:end-endBuffer);

time_f    = zeros(length(T_indexes), lambda2);
force_f   = zeros(length(T_indexes), lambda2);
pres_f    = zeros(size(pres, 1), lambda2, length(T_indexes));


for i = 1:length(T_indexes)-1
    time_f(i, :)    = time(T_indexes(i):T_indexes(i+1));
    force_f(i, :)   = force(T_indexes(i):T_indexes(i+1));
    pres_f(:, :, i) = pres(:,  T_indexes(i):T_indexes(i+1));
end

pi_f =  (time_f(1, :)- time_f(1, 1)).*(4/(time_f(1, end)- time_f(1, 1)));

%% Plot Pressures
% small perterbation for showing plots on top of one another
delta = 0.001;
figure(); 
plot(pi_f, force_f(1, :)+1*delta);
hold on;
plot(pi_f, force_f(2, :)+2*delta);
plot(pi_f, force_f(100, :)+3*delta);
plot(pi_f, force_f(200, :)+4*delta);
plot(pi_f, force_f(300, :)+5*delta);
plot(pi_f, force_f(1000, :)+6*delta);
% force     = mean(force, 1);

figure()
plot(pi_f, pres_f(10, :, 1));
hold on;
plot(pi_f, pres_f(10, :, 100));
plot(pi_f, pres_f(10, :, 1000));

figure()
plot(pi_f, pres_f(17, :, 1));
hold on;
plot(pi_f, pres_f(17, :, 100));
plot(pi_f, pres_f(17, :, 1000));

pressure_f_mean      = mean(pres_f, 3);
pressure_f_rms       = rms(pres_f, 3);
pressure_f_diff_mean = pressure_f_mean-mean(pressure_f_mean, 2);
pressure_f_diff_rms  = pressure_f_rms-mean(pressure_f_rms, 2);

figure();
surf(pi_f, 1:17, pressure_f_diff_mean, 'EdgeColor', 'none');
view(2);
colormap jet;


%% Interpolate onto finer mesh
[Xq,Yq] = meshgrid(linspace(0,pi_f(end), 1000), linspace(1,17, 1000));
[X,Y]   = meshgrid(pi_f, 1:17);

V                    = pressure_f_mean;
pres_final_mean      = interp2(X,Y,V,Xq,Yq, 'spline');

V                    = pressure_f_diff_mean;
pres_final_diff_mean = interp2(X,Y,V,Xq,Yq, 'spline');

V                    = pressure_f_rms;
pres_final_rms      = interp2(X,Y,V,Xq,Yq, 'spline');

V                    = pressure_f_diff_rms;
pres_final_diff_rms = interp2(X,Y,V,Xq,Yq, 'spline');

% figure();
% surf(Xq, Yq, flip(Vq), 'EdgeColor', 'none');
% view(2);
% colormap jet;
% colorbar;

%% Mean
y_mm                      = linspace(0.85, 136, 1000);
pres_final_flip_mean      = flip(pres_final_mean);
pres_final_diff_flip_mean = flip(pres_final_diff_mean);

figure();
set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);
surf(Xq, y_mm, pres_final_diff_flip_mean./Pt, 'EdgeColor', 'none');
view(2);
colormap jet;
colorbar;

xlim([0 4]);
xticks([0 1/2 1 3/2 2 ... 
    5/2 3 7/2 4])
xticklabels({'0', '\pi/2', '\pi',...
    '3\pi/2', '2\pi', '5\pi/2', ...
    '3\pi', '7\pi/2', '4\pi'})
xlabel('Phase (rad.)', 'interpreter', 'latex');
ylabel('z, mm', 'interpreter', 'latex');
title('$(\tilde{P} - \overline{P})/P_t$', 'interpreter', 'latex');

%% RMS
pres_final_flip_rms      = flip(pres_final_rms);
pres_final_diff_flip_rms = flip(pres_final_diff_rms);
figure();
set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);
surf(Xq, y_mm, pres_final_flip_rms./Pt, 'EdgeColor', 'none');
view(2);
colormap jet;
colorbar;

xlim([0 4]);
xticks([0 1/2 1 3/2 2 ... 
    5/2 3 7/2 4])
xticklabels({'0', '\pi/2', '\pi',...
    '3\pi/2', '2\pi', '5\pi/2', ...
    '3\pi', '7\pi/2', '4\pi'})
xlabel('Phase (rad.)', 'interpreter', 'latex');
ylabel('z, mm', 'interpreter', 'latex');
title('$\tilde{P}_{rms}/P_t$', 'interpreter', 'latex');






