%% parameters
subfolder = '';
basepath = 'C:\Users\dell\Desktop\20210924 AceS81-HaloTag-ST2 in acute mouse slice\slice2\field5\532nm(100%power)_40xWideField_bin4_3';
pathname = basepath;
path = [basepath subfolder];
    

% constants
movname = '\movie.bin'; 

ncol = 512;         % x
nrow = 102;         % y
bin = 4;
bkg = 100 * bin * bin;          % background due to camera bias (100 for bin 1x1)
dt_mov = 2;       % exposure time in millisecond
DAQname = '\movie_DAQ.txt';
dnsamp = 20;       % downsampling rate = DAQ rate
dt_daq = dt_mov/dnsamp;     % in millisecond

% load movie
fname = [path movname]
[mov, nframe] = readBinMov(fname, nrow, ncol);
mov = double(mov);
img = mean(mov, 3);
len = size(mov,3);
t_mov = [0:(len-1)]*dt_mov*0.001;     % time axis in second

% select ROI for analysis
[~, intens] = clicky(mov, img, 'Select only 1 ROI, right click when done');
% save clicky figure
saveas(gca,[path '\clicky analysis_voltage indicator.fig']);
saveas(gca,[path '\clicky analysis_voltage indicator.png']);
%saveas(gca,[path '\clicky analysis_calicium indicator.fig']);
%saveas(gca,[path '\clicky analysis_calicium indicator.png']);
% normalize intensity (in %)

%%
bkg = intens(:,end);
intens_reduBkg =  intens(:,1:end-1) - repmat(bkg, 1, size(intens_reduBkg,2));
for i = 1: size(intens,2)-1
    plot(t_mov, smooth(smooth(intens_reduBkg(:, i))));
    hold on;
end
xlabel('Time (s)');
ylabel('Intensity (au.)');
axis tight;
box off
saveas(gca,[path '\Voltage imaging.fig']);
saveas(gca,[path '\Voltage imaging.png']);

%% Voltage imaging
intensN = zeros(size(intens,1),size(intens,2)-1);
pbleach = zeros(size(intens,1),size(intens,2)-1);
bkg=intens(:,end);
for i = 1: size(intens,2)-1
    %[intensN(:,i), pbleach(:,i)] = rem_pbleach((intens(:,i)-bkg), round(size(intens(:,i),1)/10)+1); % positive indicators
    [intensN(:,i), pbleach(:,i)] = rem_pbleach(-1*(intens(:,i)-bkg), round(size(intens(:,i),1)/10)+1); % negative indicators;
    %intensN1(:,i) = intens(:,i) - bkg   + power(pbleach1(:,i) - min(pbleach1(:,i)),0.5);
    %pbleach(:,i) = -1*power(abs(pbleach(:,i)+min(pbleach(:,i))),1);
    %[intensN(:,i), ~] = rem_pbleach(intensN1(:,i), round(size(intensN1(:,i)/50000,1)));
    %intensN(:,i) = intensN(:,i) + power(pbleach(:,i) - min(pbleach(:,i)),1);
    %intensN(:,i) = intensN(:,i)./power(abs(pbleach(:,i)+max(pbleach(:,i))),1);
end
intens_base = sort(-1*intensN, 1); % negative indicators;
%intens_base = sort(intensN, 1); % positive indicators;
intens_base = mean(-1*intens_base(1:500,:)); % negative indicators;
%intens_base = mean(intens_base(1n  :10000,:)); % positive indicators
intens_norm = (intensN./repmat(intens_base, size(intensN,1),1)-1)*100; 
% intens_norm = intens_norm./(intens_norm(1,1));
% load DAQ data
%tmp = importdata([path DAQname]);   % import data 
%data = tmp.data;                    % get array
%Vm = data(:,2)*100;     % Vm in millivolt, column vector
% downsample Vm
%len = length(Vm);
%t_daq = [0:(len-1)]*dt_daq./1000;                 % time axis in second
%Vm_dnsamp = mean(reshape(Vm,dnsamp,len/dnsamp));
%intens_rembkg = intens(:,1)-intens(:,2);
%intens_rembkg = intens;

% plot traces
figure;
% plot F-t      
%subplot(2,1,1);
itens_size = size(intens)
%plot(t_mov, smooth(smooth(intens_rembkg)));
for i = 1: itens_size(2)-1
    plot(t_mov, smooth(smooth(intens_norm(:,i))));
    hold on;
end
xlabel('time (s)');
ylabel('¦¤F/F0 (%)');
axis tight;
box off
saveas(gca,[path '\Voltage imaging.fig']);
saveas(gca,[path '\Voltage imaging.png']);