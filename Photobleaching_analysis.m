 clear all
%% load a bin file to MATLAB
%parameters
subfolder = '';
basepath = 'D:\data\test\2022 Photo-stability test\Dish8\cell1\40xWideFieldBin2_1V561nm_1sExposure_600s\';
pathname = [basepath subfolder];

% constants
movname = '\movie.bin';
ncol = 648;         % x'
nrow = 648;          % y
bin = 2;             % bin
bkg = 100*power(bin,2);          % background due to camera bias (100 for bin 1x1)
dt_mov = 1;         % exposure time in second
% DAQname = '\movie_DAQ.txt';
% dnsamp = 10;        % downsampling rate = DAQ rate/camera rate
% dt_daq = 1000/9681.4795;       % in millisecond


% load movie
fname = [pathname movname]
[mov, nframe] = readBinMov(fname, nrow, ncol);
img = mean(mov, 3);
len = size(mov,3);
t_mov = [0:(len-1)]*dt_mov;     % time axis in second
%img = mean(mov(:,:,1:150), 3);
% select ROI for analysis
% [mov_new] = extrabin(mov, 2);
mov = single(mov);
% [roi_points,intens] = clicky(mov, img, 'Please select interested regions')
% [roi_points,intens] = clicky_v1(mov, img, 1000, 10000, 'Please select interested regions')
[roi_points,intens] = clicky_v2(mov, dt_mov, bin, img, 100 , 30000, 'Please select interested regions')
intens_rembkg = intens(:,1)-intens(:,2);
% save clicky figure
saveas(gca,[pathname '\clicky analysis.fig']);
saveas(gca,[pathname '\clicky analysis.png']);
%% I-t figure, check the spontaneous AP
figure()
plot(t_mov,intens_rembkg);
box off
axis tight
ylabel('Intensity')
xlabel('Time (s)')
title('Fast imaging: spontaneous AP')
saveas(gca,[pathname '\spontaneous AP.fig']);
saveas(gca,[pathname '\spontaneous AP.png']);
%% For eliminating the AP (peak) in curve
intens_rembkg_rempeak = imerode(intens_rembkg, ones(200,1));
intens_rembkg_rempeak = smooth(intens_rembkg_rempeak, 200);
%% two-exponential fit
% Initialize arrays to store fits and goodness-of-fit.

% Fit model to data.

param = zeros(size(intens_rembkg,2),6);
intens_rembkg_norm = intens_rembkg./intens_rembkg(1,1);
figure()
plot(intens_rembkg_norm)
[a,b] = ginput(1);
intens_rembkg_norm = intens_rembkg_norm./intens_rembkg_norm(round(a(1)));
half_loc = min(find(intens_rembkg_norm<=0.5));
t_half = t_mov(half_loc);


for n = 1:size(intens_rembkg,2)
fitresult = cell( 3, 1 );
gof = struct( 'sse', cell( 3, 1 ), ...;
    'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );
[xData1, yData1] = prepareCurveData(t_mov,intens_rembkg(:,n) );
% [xData2, yData2] = prepareCurveData( fitted_t2,FlareN_tofit);
% [xData3, yData3] = prepareCurveData( fitted_t3,EGFR_tofit );

% Set up fittype and options.
ft = fittype( 'exp2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'LAR';
[fitresult{n}, gof(n)] = fit( xData1, yData1, ft, opts );
temp = fit( xData1, yData1, ft, opts );
param(n,1) = temp.a;
param(n,2) = -1./temp.b;
param(n,3) = temp.c;
param(n,4) = -1./temp.d;
param(n,5) = param(n,1)./(param(n,1)+param(n,3));
param(n,6) = param(n,3)./(param(n,1)+param(n,3));
% [fitresult{2}, gof(2)] = fit( xData2, yData2, ft, opts );
% [fitresult{3}, gof(3)] = fit( xData3, yData3, ft, opts );
figure();
plot(fitresult{n});hold on
plot(t_mov,yData1);
box off
axis tight
ylabel('Intensity')
xlabel('Time (s)')
title('Photobleaching curve of mNeonGreen')
saveas(gca,[pathname '\region ' num2str(n) '.fig']);
saveas(gca,[pathname '\region ' num2str(n) '.png']);

clear temp
end

for n = 1:size(intens_rembkg_norm,2)
fitresult = cell( 3, 1 );
gof = struct( 'sse', cell( 3, 1 ), ...
    'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );
[xData1, yData1] = prepareCurveData(t_mov,intens_rembkg_norm(:,n) );
% [xData2, yData2] = prepareCurveData( fitted_t2,FlareN_tofit);
% [xData3, yData3] = prepareCurveData( fitted_t3,EGFR_tofit );

% Set up fittype and options.
ft = fittype( 'exp2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'LAR';
[fitresult{n}, gof(n)] = fit( xData1, yData1, ft, opts );
temp = fit( xData1, yData1, ft, opts );
param(n,1) = temp.a;
param(n,2) = -1./temp.b;
param(n,3) = temp.c;
param(n,4) = -1./temp.d;
param(n,5) = param(n,1)./(param(n,1)+param(n,3));
param(n,6) = param(n,3)./(param(n,1)+param(n,3));
% [fitresult{2}, gof(2)] = fit( xData2, yData2, ft, opts );
% [fitresult{3}, gof(3)] = fit( xData3, yData3, ft, opts );
figure();
plot(fitresult{n});hold on
plot(t_mov,yData1);
box off
axis tight
ylabel('Intensity')
xlabel('Time (s)')
title('Photobleaching curve of mNeonGreen')
saveas(gca,[pathname '\Normalized region ' num2str(n) '.fig']);
saveas(gca,[pathname '\Normalized region ' num2str(n) '.png']);
clear temp
end
gof1(1) = gof.rsquare;
gof1(2) = gof.adjrsquare;
gof1(3) = gof.sse;
gof1(4) = gof.rmse;
%%%% two-exponential fit
% Initialize arrays to store fits and goodness-of-fit.
%% Fitting after Removing AP
% Fit model to data.
param = zeros(size(intens_rembkg_rempeak,2),6);
intens_rembkg_norm = intens_rembkg_rempeak./intens_rembkg_rempeak(1,1);
for n = 1:size(intens_rembkg_rempeak,2)
fitresult = cell( 3, 1 );
gof = struct( 'sse', cell( 3, 1 ), ...
    'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );
[xData1, yData1] = prepareCurveData(t_mov,intens_rembkg_rempeak(:,n) );
% [xData2, yData2] = prepareCurveData( fitted_t2,FlareN_tofit);
% [xData3, yData3] = prepareCurveData( fitted_t3,EGFR_tofit );

% Set up fittype and options.
ft = fittype( 'exp2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'LAR';
[fitresult{n}, gof(n)] = fit( xData1, yData1, ft, opts );
temp = fit( xData1, yData1, ft, opts );
param(n,1) = temp.a;
param(n,2) = -1./temp.b;
param(n,3) = temp.c;
param(n,4) = -1./temp.d;
param(n,5) = param(n,1)./(param(n,1)+param(n,3));
param(n,6) = param(n,3)./(param(n,1)+param(n,3));
% [fitresult{2}, gof(2)] = fit( xData2, yData2, ft, opts );
% [fitresult{3}, gof(3)] = fit( xData3, yData3, ft, opts );
figure();
plot(fitresult{n});hold on
plot(t_mov,yData1);
box off
axis tight
ylabel('Intensity')
xlabel('Time (s)')
title('Photobleaching curve of Ace-mNeon')
saveas(gca,[pathname '\(Removing AP) region ' num2str(n) '.fig']);
saveas(gca,[pathname '\(Removing AP) region ' num2str(n) '.png']);

clear temp
end

for n = 1:size(intens_rembkg_rempeak,2)
fitresult = cell( 3, 1 );
gof = struct( 'sse', cell( 3, 1 ), ...
    'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );
[xData1, yData1] = prepareCurveData(t_mov,intens_rembkg_rempeak(:,n) );
% [xData2, yData2] = prepareCurveData( fitted_t2,FlareN_tofit);
% [xData3, yData3] = prepareCurveData( fitted_t3,EGFR_tofit );

% Set up fittype and options.
ft = fittype( 'exp2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'LAR';
[fitresult{n}, gof(n)] = fit( xData1, yData1, ft, opts );
temp = fit( xData1, yData1, ft, opts );
param(n,1) = temp.a;
param(n,2) = -1./temp.b;
param(n,3) = temp.c;
param(n,4) = -1./temp.d;
param(n,5) = param(n,1)./(param(n,1)+param(n,3));
param(n,6) = param(n,3)./(param(n,1)+param(n,3));
% [fitresult{2}, gof(2)] = fit( xData2, yData2, ft, opts );
% [fitresult{3}, gof(3)] = fit( xData3, yData3, ft, opts );
figure();
plot(fitresult{n});hold on
plot(t_mov,yData1);
box off
axis tight
ylabel('Intensity')
xlabel('Time (s)')
title('Photobleaching curve of Ace-mNeon')
saveas(gca,[pathname '\(Removing AP) Normalized region ' num2str(n) '.fig']);
saveas(gca,[pathname '\(Removing AP) Normalized region ' num2str(n) '.png']);
clear temp
end
%% Excel
xlswrite([pathname 'analysis.xlsx'],{'Brightness ','Exposure time (s)','Time (s)', 't_half','Fit parameters'},'Average','A1');
xlswrite(  [pathname 'analysis.xlsx'],intens_rembkg(1,1),'Average','A2');
xlswrite([pathname 'analysis.xlsx'],dt_mov','Average','B2');
xlswrite([pathname 'analysis.xlsx'],len*dt_mov,'Average','C2');
xlswrite([pathname 'analysis.xlsx'],t_half,'Average','D2');
xlswrite([pathname 'analysis.xlsx'],param(1),'Average','E2');
xlswrite([pathname 'analysis.xlsx'],param(2),'Average','F2');
xlswrite([pathname 'analysis.xlsx'],param(3),'Average','G2');
xlswrite([pathname 'analysis.xlsx'],param(4),'Average','H2');
xlswrite([pathname 'analysis.xlsx'],param(5),'Average','I2');
xlswrite([pathname 'analysis.xlsx'],param(6),'Average','J2');
xlswrite([pathname 'analysis.xlsx'],{'A1 ','tau1 (s)','A2', 'tau2 (s)','Normalized  A1','Normalized A1'},'Average','E3');
xlswrite([pathname 'analysis.xlsx'],{'Rsquare','adjRsquare','sse','rmse'},'Average','E5');
xlswrite([pathname 'analysis.xlsx'],gof1(1),'Average','E4');
xlswrite([pathname 'analysis.xlsx'],gof1(2),'Average','F4');
xlswrite([pathname 'analysis.xlsx'],gof1(3),'Average','G4');
xlswrite([pathname 'analysis.xlsx'],gof1(4),'Average','H4');
% xlswrite([pathname 'analysis.xlsx'],intens_mean','Average','I2');
