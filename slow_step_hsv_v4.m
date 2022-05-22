% parameters
clear all;
subfolder = '';
basepath = 'H:\ZouLab\YJunqi\Screening method\System test\20220507 Depolarization titration\Dish3\cell1\205956_-100 to +100_500ms\';
pathname = basepath
path = [basepath subfolder];

% constants
movname = '\movie.bin';
ncol = 176;         % x
nrow = 96;         % y
bkg = 400;          % background due to camera bias (100 for bin 1x1)
dt_mov = 0.9452;       % exposure time in millisecond
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
saveas(gca,[path '\clicky analysis_gap junction.fig']);
saveas(gca,[path '\clicky analysis_gap junction.png']);
% normalize intensity wrt -100 mV (in %)
bkg=mean(intens(:,end));
background_reference=sort(reshape(mean(mov,3),ncol*nrow,1));
background_reference=mean(background_reference(10:round(0.05*length(background_reference))))
intens_norm = (intens-background_reference)/(max(intens)-background_reference)*100;
% intens_norm = intens_norm./(intens_norm(1,1));
% load DAQ data

tmp = importdata([path DAQname]);   % import data 
data = tmp.data;                    % get array
Vm = data(:,2)*100;     % Vm in millivolt, column vector
% downsample Vm
len = length(Vm);
t_daq = [0:(len-1)]*dt_daq./1000;                 % time axis in second
Vm_dnsamp = mean(reshape(Vm,dnsamp,len/dnsamp));
intens_rembkg = intens(:,1:end-1)-intens(:,end);
%intens_rembkg = intens;

%%
% plot traces
figure;
% plot F-t
subplot(2,1,1);
itens_size = size(intens)
plot(t_mov, smooth(smooth(intens_rembkg)));
% for i = 1: itens_size(2)
%     plot(t_mov, smooth(smooth(intens_rembkg(:,i))));
%     hold on;
% end
xlabel('time (s)')';
ylabel('Intens (a.u.)');
axis tight;
box off
% plot Vm-t
subplot(2,1,2);
plot(t_daq,Vm);
xlabel('time (s)');
ylabel('V_m (mV)');
axis tight;
box off
% save F-V figure
saveas(gca,[path '\F-V analysis.fig']);
saveas(gca,[path '\F-V analysis.png']);
%%
color_map = hsv(11);
color_map = flipud(color_map)
intens_stack = zeros(1000,6);
for n = 1:11
  intens_stack(:,n) = intens_rembkg((251+1000*(n-1):250+1000*n));
end
figure();
for n = 1:11
plot([0:dt_mov:(1000-1)*dt_mov]', intens_stack(:,n),'color',color_map(n,:));hold on
end
box off
axis tight
xlabel('Time (ms)')
ylabel('Intensity')
saveas(gca,[path '\F-V stack.fig']);
saveas(gca,[path '\F-V stack.png']);

norm_stack = intens_stack./mean(intens_stack(1:ceil(100/dt_mov),:),1);
dy_stack = intens_stack./mean(intens_stack(1:ceil(100/dt_mov),:),1) - 1;

figure();
for n = 1:11
plot(smooth(smooth([0:dt_mov:(1000-1)*dt_mov]', norm_stack(:,n))),'color',color_map(n,:));hold on
end
box off
axis tight
xlabel('Time (ms)')
ylabel('Normalized Intensity (from -70 mV)')
saveas(gca,[path '\Norm.F-V stack.fig']);
saveas(gca,[path '\Norm.F-V stack.png']);

figure();
for n = 1:11
plot(smooth(smooth([0:dt_mov:(1000-1)*dt_mov]',dy_stack(:,n))),'color',color_map(n,:));hold on
end
box off
axis tight
xlabel('Time (ms)')
ylabel('Dynamic range (from -70 mV)')
saveas(gca,[path '\Dynamic_F-V stack.fig']);
saveas(gca,[path '\Dynamic_F-V stack.png']);

%% for the colorbar
figure();
temp = [-100:(200/10):100]';
temp = repmat(temp',3,1);
imshow(temp',[])
colormap(hsv(11))
saveas(gca,[path '\colorbar.fig']);
saveas(gca,[path '\colorbar.png']);

%%
dy_plat = [];
dy_peak = [];
dy_base = [];
for n = 1:11
  dy_stack_single = dy_stack(:,n);
  dy_plat(:,n) = dy_stack_single(737:747); % 697-706 ms
  dy_peak(:,n) = dy_stack_single(269:379); % 254-358 ms
%   dy_peak(:,n) = dy_stack_single(263:363); % 249-343 ms
  dy_base(:,n) = dy_stack_single(105:215); % 99-203 ms
end
% peak1 = -sort(-dy_peak(:,1:2));
% peak2 = sort(dy_peak(:,3:11));

xlswrite([pathname '1analysis.xlsx'],{'Vm', 'peak_delta(254:358 points) F/F0', 'plat(737-747 points)_delta F/F0', 'base_delta(105-215 potints) F/F0'},'A1:D1');
xlswrite([pathname '1analysis.xlsx'],{-100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100}','A2:A12');
% xlswrite([pathname '1analysis.xlsx'],mean(peak1(1:100,:))','B2:B3');
% xlswrite([pathname '1analysis.xlsx'],mean(peak2(1:30,:))','B4:B12');
xlswrite([pathname '1analysis.xlsx'],mean(dy_peak,1)','B2:B12');
xlswrite([pathname '1analysis.xlsx'],mean(dy_plat,1)','C2:C12');
xlswrite([pathname '1analysis.xlsx'],mean(dy_base,1)','D2:D12');

save([pathname 'All variants.mat']);
%%
figure();
for n = 1:11
plot(smooth([0:dt_mov:(1000-1)*dt_mov]',dy_stack(:,n), 30),'color',color_map(n,:));hold on
end
box off
axis tight
xlabel('Time (ms)')   
ylabel('Dynamic range (from -70 mV)')
saveas(gca,[path '\Dynamic_F-V stack_smoothStep30.fig']);
saveas(gca,[path '\Dynamic_F-V stack_smoothStep30.png']);
