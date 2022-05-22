% parameters
subfolder = '';
basepath = 'D:\data\test\20220122 msx\Dish1\Field1\532nm_40xWideField_bin2';
pathname = basepath;
path = [basepath subfolder];

% constants'
movname = '\movie.bin';
ncol = 1024;         % x
  
nrow = 1024;         % y
bkg = 400;          % background due to camera bias (100 for bin 1x1)
dt_mov = 10;       % exposure time in millisecond
DAQname = '\movie_DAQ.txt';
dnsamp = 20;       % downsampling rate = DAQ rate
dt_daq = dt_mov/dnsamp;     % in millisecond

% load movie
%filter = ones(1, 1, 4);
fname = [path movname];
[mov, nframe] = readBinMov(fname, nrow, ncol);
%mov_1 = convn(mov, filter, 'valid');
%mov_2 = mov_1(:,:,1:4:end);
%mov_2 = double(mov_2);
%img_2 = mean(mov_2, 3);
mov = double(mov);
img = mean(mov, 3);
len = size(mov,3);
t_mov = [0:(len-1)]*dt_mov*0.001;     % time axis in second

% select ROI for analysis
[~, intens] = clicky(mov, img, 'Select only 1 ROI, right click when done');
% save clicky figure
%saveas(gca,[path '\clicky analysis_voltage indicator.fig']);
%saveas(gca,[path '\clicky analysis_voltage indicator.png']);
saveas(gca,[path '\clicky analysis_calicium indicator.fig']);
saveas(gca,[path '\clicky analysis_calicium indicator.png']);
% normalize intensity (in %)

%% Calcium imaging
intensN = zeros(size(intens,1),size(intens,2)-1);
bkg=mean(intens(:,end));
for i = 1: size(intens,2)-1
    intensN(:,i) = intens(:,i)-bkg;
end
intens_base = sort(intensN, 1);
intens_base = mean(intens_base(1:7000,:));
intens_norm = (intensN-repmat(intens_base, size(intensN,1),1))./repmat(intens_base, size(intensN,1),1)*100;
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
    plot(t_mov, smooth(intens_norm(:,i), 50));
    hold on;
end
xlabel('time (s)');
ylabel('¦¤F/F0 (%)');
axis tight;
box off
saveas(gca,[path '\Calcium imaging.fig']);
saveas(gca,[path '\Calcium imaging.png']);
%%
% load DAQ data
tmp = importdata([path DAQname]);   % import data 
data = tmp.data;                    % get array
Vm = data(:,2)*100;     % Vm in millivolt, column vector
% downsample Vm
len = length(Vm);
t_daq = [0:(len-1)]*dt_daq./1000;                 % time axis in second
Vm_dnsamp = mean(reshape(Vm,dnsamp,len/dnsamp));
% plot Vm-t
subplot(2,1,2);
plot(t_daq,smooth(smooth(smooth(Vm))));
xlabel('time (s)');
ylabel('V_m (mV)');
axis tight;
box off
% save F-V figure
saveas(gca,[path '\F-V analysis.fig']);
saveas(gca,[path '\F-V analysis.png']);