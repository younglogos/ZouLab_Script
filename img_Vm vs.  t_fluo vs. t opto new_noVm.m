%% load and select
clear all; clc;
% Movie loading path
% Firstly, determine patch or optopatch.1→patch，2→optopatch.If patch, patch

patch = 2;
testmode = 2
if testmode == 1
    dt_mov = 1.218045;
    dnsamp = 8210.1816/(1/dt_mov);        % downsampling rate = DAQ rate/camera rate
else
    dt_mov = 2.153503;
    dnsamp = 4643.6034/(1/dt_mov);        % downsampling rate = DAQ rate/camera rate
end
subfolder = '';
basepath = 'D:\data\test\HY20220514 Cepheid comparisons & applications in neuron\Dish1\Cell2\014355_40x 405 1.93 3V\';
pathname = [basepath subfolder];

load ([pathname '\SNAPT-1.93ms blue light - 3V_4643.6034\matlab variables.mat']);%输入激光文件夹名称
wavelength = 405
;%输入波长信息
% Fs = samprate;

% constants
movname = '\movie.bin';
ncol = 312;         % x
nrow = 208;          % y
camera_bias = 400;          % background due to camera bias (100 for bin 1x1)
DAQname = '\movie_DAQ.txt';

% load movie
fname = [pathname movname];
[mov, nframe] = readBinMov(fname, nrow, ncol);
mov = double(mov);
img = mean(mov, 3);
len = size(mov,3);
t_mov = [0:(len-1)]*dt_mov/1000;     % time axis in second


% select ROI for analysis
[~, intens_sig] = clicky(mov, img, 'Select only 1 ROI, right click when done');
background = mean(intens_sig(:,2));
% save clicky figure
saveas(gca,[pathname '\1 clicky analysis_AP.fig']);
saveas(gca,[pathname '\1 clicky analysis_AP.png']);



% load DAQ data
% load DAQ data
tmp = importdata([pathname DAQname]);   % import data 
data = tmp.data;                    % get array
% data = data(:,2)*100;                 % data in millivolt, column vector
dt_daq = dt_mov/dnsamp;             % DAQ dt in millisecond
t_daq = [0:length(data)-1]*dt_daq;       % DAQ time axis in second

%% plot data-t and Intensity-t
    figure
    subplot(5,1,1)
    plot(mat2gray(img));% gray img plot
    imshow(mat2gray(img));
    title('selected cell');
    subplot()
    time=(1:length(data)')'./9681.4793;
    HCinput = 0;
    subplot(5,1,2:3);
    if exist('Lasers_DO') == 0;
        Lasers_DO = AO_DO(2,:);
    end
        if exist('Lasers_AO') == 0;
        Lasers_AO = AO_DO(1,:);
    end
    plot(t_daq,Lasers_DO.*Lasers_AO,'Color',getrgb(wavelength),'LineWidth',2);hold on;
    ylim([0 1])
    box off
    title(['Light pulse (' num2str(wavelength) ') nm ']);
    axis tight;
    xlabel('time(s)');
    ylabel('Laser DO');
    box off
    subplot(5,1,4:5);
    plot(t_mov,intens_sig(:,1)-intens_sig(:,end));
    hold on
    title('Optical signal');
    xlabel('t(s)');
    ylabel('intensity (W/O background)');
    axis tight;
    xL=[0,max(t_mov)];yL=[min(intens_sig(:,1)-background)-2,max(intens_sig(:,1)-background)+2];% range of figure
    plot(xL,[yL(2),yL(2)],'w',[xL(2),xL(2)],[yL(1),yL(2)],'w')
    box off
    set(gcf,'outerposition',get(0,'screensize'));
    saveas(gca,[pathname '\1 integrated analysis.fig'])
    saveas(gca,[pathname '\1 integrated analysis.png'])
    xlswrite([pathname 'analysis.xlsx'],{'time(s)','optical signal'}','Average','A1');
    xlswrite([pathname 'analysis.xlsx'],t_mov','Average','B1');
    xlswrite([pathname 'analysis.xlsx'],intens_sig(:,1)-intens_sig(:,2),'Average','C1');

%% 10 frames in 1 for R-GCEO signal (smooth)
    intens10 = reshape(sum(reshape(intens_sig,10,[])),[],2)
    t10 = reshape(min(reshape(t_mov,10,[])),1,[])
    figure
    subplot(2,1,1)
    plot(mat2gray(img));% gray img plot
    imshow(mat2gray(img));
    title('selected cell');
    subplot(2,1,2)
    plot(t10,smooth((intens10(:,1)-intens10(:,2)),10));
    hold on
    title('Optical signal');
    xlabel('t(s)');
    ylabel('intensity (W/O background)');
    axis tight;
    xL=[0,max(t_mov)];yL=[min(intens_sig(:,1)-background)-2,max(intens_sig(:,1)-background)+2];% range of figure
    plot(xL,[yL(2),yL(2)],'w',[xL(2),xL(2)],[yL(1),yL(2)],'w')
    box off
    set(gcf,'outerposition',get(0,'screensize'));
    saveas(gca,[pathname '\4 integrated analysis.fig'])
    saveas(gca,[pathname '\4 integrated analysis.png'])
    xlswrite([pathname 'analysis4.xlsx'],{'time(s)','optical signal'}','Average','A1');
    xlswrite([pathname 'analysis4.xlsx'],t_mov','Average','B1');
    xlswrite([pathname 'analysis4.xlsx'],intens_sig(:,1)-intens_sig(:,2),'Average','C1');