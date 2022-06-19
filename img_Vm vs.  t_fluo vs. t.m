%% load and select
clear all; clc;
% Movie loading path
% Firstly, determine patch or optopatch.1¡úpatch£¬2¡úoptopatch
testmode = 1
if testmode == 2
    dt_mov = 2;
    dnsamp = 9681.4793/(1/dt_mov);        % downsampling rate = DAQ rate/camera rate
else
    dt_mov = 2.0284;
    dnsamp = 10471.2/(1/dt_mov);        % downsampling rate = DAQ rate/camera rate
end
subfolder = 'D:\data\test\HY20220602 Cepheid with CheRiff & GCaMP6s\Dish1\Cell5\40x WideField bin2 488_5 561_1 500hz 9mM 60s\';
basepath = '';
pathname = [basepath subfolder];
% constants 
movname = '\movie.bin';
ncol = 1024;         % x 
nrow = 208;         % y
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
[~, intens] = clicky(mov, img, 'Select only 1 ROI, right click when done');
background = mean(intens(:,2));
% save clicky figure
saveas(gca,[pathname '\clicky analysis_V.fig']);
saveas(gca,[pathname '\clicky analysis_V.png']);



% load DAQ data   
% load DAQ data
tmp = importdata([pathname DAQname]);   % import data 
data = tmp.data;                    % get array
Vm = data(:,2)*100;                 % Vm in millivolt, column vector
dt_daq = dt_mov/dnsamp;             % DAQ dt in millisecond
t_daq = [0:length(Vm)-1]*dt_daq;       % DAQ time axis in second
%% plot Vm-t and Intensity-t
figure
subplot(3,1,1)
plot(mat2gray(img));% gray img plotcgvbn 
title('selected cell');
%subplot(3,1,2)
%plot(t_daq,Vm);
%title('membrane voltage signal');
%xlabel('t(s)');
%ylabel('Vm (mV)');
%axis tight;
%hold on
%xL=[0,max(t_mov)];yL=[min(Vm)-10,max(Vm)+10];% range of figure
%plot(xL,[yL(2),yL(2)],'w',[xL(2),xL(2)],[yL(1),yL(2)],'w');hold on
%plot([60:10:100], [-52 -52 -52 -52 -52],'.','Color',[0.99 0 0 ]) ;hold on;
%plot([60 max(t_daq)], [-50 -50],'r','LineWidth',3) ;hold on
%box off
%axis tight
%box off 
subplot(3,1,3)
plot(t_mov,intens(:,1)-background);
hold on
title('Optical signal');
xlabel('t(s)');
ylabel('intensity (W/O background)');
axis tight;
xL=[0,max(t_mov)];yL=[min(intens(:,1)-background)-2,max(intens(:,1)-background)+2];% range of figure
plot(xL,[yL(2),yL(2)],'w',[xL(2),xL(2)],[yL(1),yL(2)],'w')
box off
set(gcf,'outerposition',get(0,'screensize'));
saveas(gca,[pathname '\integrated analysis_1.fig'])
saveas(gca,[pathname '\integrated analysis_1.png'])
%% 10 frames in 1 for glutamate signal (unsmooth)
    intens10 = reshape(sum(reshape(intens,10,[])),[],2)
    t10 = reshape(min(reshape(t_mov,10,[])),1,[])
    figure
    subplot(2,1,1)
    plot(mat2gray(img));% gray img plot
    imshow(mat2gray(img));
    title('selected cell');
    subplot(2,1,2)
    plot(t10,((intens10(:,1)-intens10(:,2))));
    hold on
    title('Optical signal');
    xlabel('t(s)');
    ylabel('intensity (W/O background)');
    axis tight;
    xL=[0,max(t_mov)];yL=[min(intens(:,1)-background)-2,max(intens(:,1)-background)+2];% range of figure
    plot(xL,[yL(2),yL(2)],'w',[xL(2),xL(2)],[yL(1),yL(2)],'w')
    box off
    set(gcf,'outerposition',get(0,'screensize'));
    saveas(gca,[pathname '\4 integrated analysis.fig'])
    saveas(gca,[pathname '\4 integrated analysis.png'])
    xlswrite([pathname 'analysis4.xlsx'],{'time(s)','optical signal'}','Average','A1');
    xlswrite([pathname 'analysis4.xlsx'],t_mov','Average','B1');
    xlswrite([pathname 'analysis4.xlsx'],intens(:,1)-intens(:,2),'Average','C1');