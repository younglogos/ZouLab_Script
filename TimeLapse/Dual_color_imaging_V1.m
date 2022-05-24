clear all
%% loading two channels

name1 = 'RhoVR1';
name2 = 'mNeon-CAAX';
wvlt1 = 580;
wvlt2 = 520;
name1color = getrgb(wvlt1);
name2color = getrgb(wvlt2);
bin = 1;
subfolder = '';
basepath = 'X:\91 Data and analysis\YJunqi\Screening method\Stable HEK293T cell line_Kir2.1-P2A-ArcLight\20210305 KCl stimulus pre-experiment\Dish2\test\';
pathname = [basepath subfolder];
FileTif=[pathname name1 '.tif'];
InfoImage=imfinfo(FileTif);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
mov=zeros(nImage,mImage,NumberImages,'uint16');
start_frame = 1;

TifLink = Tiff(FileTif, 'r');
for i=1:NumberImages
    TifLink.setDirectory(i);
    mov(:,:,i)=TifLink.read();
end
TifLink.close();
dt_mov = 10000;         % exposure time in millisecond
mov = double(mov);
mov = mov(:,:,start_frame:end);
img = mean(mov, 3);
len = size(mov,3);
% t_mov = [0:(len-1)]*dt_mov/1000;     % time axis in second

% second channel
FileTif=[pathname name2 '.tif'];
InfoImage=imfinfo(FileTif);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
mov2=zeros(nImage,mImage,NumberImages,'uint16');

TifLink = Tiff(FileTif, 'r');
for i=1:NumberImages
    TifLink.setDirectory(i);
    mov2(:,:,i)=TifLink.read();
end
TifLink.close();
% ncol = 1024./bin;         % x
% nrow = 1024./bin;          % y
ncol = size(mov,2);         % x
nrow = size(mov,1);          % y
bkg = 100*power(bin,2);          % background due to camera bias (100 for bin 1x1)
dt_mov = 10000;         % exposure time in millisecond
mov2 = double(mov2);
mov2 = mov2(:,:,start_frame:end);
img2 = mean(mov2, 3);
len2 = size(mov2,3);
t_mov = [0:(len2-1)]*dt_mov/1000;     % time axis in second
% mov = mov(:,:,3:end);
% mov2 = mov2(:,:,3:end);
mkdir([pathname 'analysis']);

%% save 2 movies
mkdir([pathname 'analysis']);
cmin1 = 100;
cmax1 = 200;
figure(20)
for j = 1:(size(mov,3))
    imshow(mov(:,:,j),[cmin1,cmax1], 'Border', 'tight');hold on
    colormap(pseudocolor(wvlt1));
    text(10,30,[num2str(t_mov(j)) ' s'], 'FontSize', 18, 'color', [0.99 0.99 0.99]);
    %text(550,20,[name1], 'FontSize', 20, 'color', name2color)
    text(10,nrow-20,[name1], 'FontSize', 18, 'color', name1color);
    M(j) = getframe(gca);
    %     saveas(gca, ['SNAPT_movie_' num2str((j)) '.tif']);  % Use this line to save each frame as a separate TIF file.
end
close(gcf)
myObj = VideoWriter([pathname  '\analysis\' name1 ' .avi'],'Uncompressed AVI');%create an AVI file
myObj.FrameRate = 5;
open(myObj);
writeVideo(myObj,M)
close(myObj);

myObj = VideoWriter([pathname  '\analysis\' name1  '_compressed_H.264 encoding.mp4'],'MPEG-4');%create an AVI file
myObj.FrameRate = 5;
open(myObj);
writeVideo(myObj,M)
close(myObj);

cmin2 = 100;
cmax2 = 500;
clear M
figure(20)
for j = 1:(size(mov2,3))
    imshow(mov2(:,:,j),[cmin2,cmax2], 'Border', 'tight');hold on;
    colormap(pseudocolor(wvlt2));
    text(10,30,[num2str(t_mov(j)) ' s'], 'FontSize', 18, 'color', [0.99 0.99 0.99])
    %text(550,20,[name2], 'FontSize', 20, 'color', getrgb(598));
    text(10,nrow-20,[name2], 'FontSize', 18, 'color', name2color);
    M(j) = getframe(gca);
    %     saveas(gca, ['SNAPT_movie_' num2str((j)) '.tif']);  % Use this line to save each frame as a separate TIF file.
end
close(gcf)
myObj = VideoWriter([pathname  '\analysis\' name2 ' .avi'],'Uncompressed AVI');%create an AVI file
myObj.FrameRate = 5;
open(myObj);
writeVideo(myObj,M)
close(myObj);
myObj = VideoWriter([pathname  '\analysis\' name2  '_compressed_H.264 encoding.mp4'],'MPEG-4');%create an AVI file
myObj.FrameRate = 5;
open(myObj);
writeVideo(myObj,M)
close(myObj);
%% background checking
[~, intens_forbkg,intens2_bkg] = clicky_fullscreen(mov, bin, mean(mov(:,:,1:end),3), bkg ,200 ,name1, mov2,wvlt1);
saveas(gca,[pathname 'analysis\selected regions of bkg.fig']);
saveas(gca,[pathname 'analysis\selected regions of bkg.png']);

% figure(); 
% plot(t_mov,intens_forbkg);
% box off
% axis tight
% xlabel('time(s)')
% title('Selected background')
% saveas(gca,[pathname 'analysis' '\selected bkg.fig']);
% saveas(gca,[pathname 'analysis' '\selected bkg.png']);


intens_forbkg = intens_forbkg-mean(intens_forbkg(1:8));
intens_forbkg = intens_forbkg';

figure();
plot(t_mov,intens_forbkg)
axis tight
box off
xlabel('Time (s)')
ylabel('normalized intensity for bkg correction')
saveas(gca,[pathname 'analysis' '\bkg_trace_' num2str(name1) '_.fig']);
saveas(gca,[pathname 'analysis' '\bkg_trace_' num2str(name1) '_.png']);

% figure();
% plot(t_mov,intens_bkg);
% box off
% axis tight
% xlabel('time(s)')
% title('Selected background')
% saveas(gca,[pathname 'analysis' '\selected bkg' num2str(name1) '.fig']);
% saveas(gca,[pathname 'analysis' '\selected bkg' num2str(name1) '.png']);

intens2_bkg = intens2_bkg-mean(intens2_bkg(1:8));
intens2_bkg = intens2_bkg';
figure();
plot(t_mov,intens2_bkg)
axis tight
box off
xlabel('Time (s)')
ylabel('normalized intensity for bkg correction')
saveas(gca,[pathname 'analysis' '\selected bkg' num2str(name2) '.fig']);
saveas(gca,[pathname 'analysis' '\selected bkg' num2str(name2) '.png']);
%% analysis
[roi_points, intens,intens2] = clicky_fullscreen(mov, bin, mean(mov(:,:,1:end),3), bkg ,cmax1 ,name1, mov2,wvlt1);
saveas(gca,[pathname 'analysis\selected regions of ' name1 '.fig']);
saveas(gca,[pathname 'analysis\selected regions of ' name1 '.png']);

figure();
set(gcf,'outerposition',get(0,'screensize'));
imshow(img,[cmin1,cmax1]);hold on;
colormap(pseudocolor(wvlt1));
for region_num=1:max(size(roi_points))
    plot([roi_points{1,region_num}(:,1)],[roi_points{1,region_num}(:,2)],'LineWidth',3,'Color', getrgb(400+round((region_num-1)/size(intens,2)*270)));hold on
    text(mean(roi_points{1,region_num}(:,1)),mean(roi_points{1,region_num}(:,2)),num2str(region_num),'Color',[0.99 0.99 0.99],'FontSize',20);hold on;
    plot([8 69.53],[nrow-20 nrow-20],'LineWidth',5,'Color',[0.99 0.99 0.99]);
end
saveas(gca,[pathname 'analysis\selected regions of ' name1 '.fig']);
saveas(gca,[pathname 'analysis\selected regions of ' name1 '.png']);
close(gcf)


figure();
set(gcf,'outerposition',get(0,'screensize'));
imshow(img2,[cmin2,cmax2]);hold on;
colormap(pseudocolor(wvlt2));
for region_num=1:max(size(roi_points))
    plot([roi_points{1,region_num}(:,1)],[roi_points{1,region_num}(:,2)],'LineWidth',3,'Color', getrgb(400+round((region_num-1)/size(intens,2)*270)));hold on
    text(mean(roi_points{1,region_num}(:,1)),mean(roi_points{1,region_num}(:,2)),num2str(region_num),'Color',[0.99 0.99 0.99],'FontSize',20);hold on;
    plot([8 69.53],[nrow-20 nrow-20],'LineWidth',5,'Color',[0.99 0.99 0.99]);
end
saveas(gca,[pathname 'analysis\selected regions of ' name2 '.fig']);
saveas(gca,[pathname 'analysis\selected regions of ' name2 '.png']);
close(gcf)
for n = 1:size(intens,2)
    figure();
    plot(t_mov,(intens(:,n))./(intens(1,n)),'Color',name1color,'LineWidth',1.5);hold on
    plot(t_mov,(intens2(:,n))./(intens2(1,n)),'Color',name2color,'LineWidth',1.5);
    title(['Region ' num2str(n)]);
    legend(name1, name2 ,'Location','northoutside','Orientation','horizontal');
    xlabel('Time(s)');
    ylabel('Relative intensity')
    box off
    axis tight
    saveas(gca,[pathname 'analysis\Dual-color plot region ' num2str(n) '.fig']);
    saveas(gca,[pathname 'analysis\Dual-color plot region ' num2str(n) '.png']);
    close(gcf)
    figure();
    [hAx,hLine1,hLine2] = plotyy(t_mov,(intens(:,n))./(intens(1,n)),t_mov,(intens2(:,n))./(intens2(1,n)));
    ylabel(hAx(1),name1) % left y-axis
    ylabel(hAx(2),name2) % right y-axis
    title(['Region ' num2str(n)]);
    xlabel('Time(s)')
    box off
    xlim(hAx(1),[0 max(t_mov)]);
    xlim(hAx(2),[0 max(t_mov)]);
    axis tight
    hLine1.Color = name1color;
    hLine2.Color = name2color;
    set(hAx(:),'Ycolor','k') %设定两个Y轴的颜色为黑色
    legend( name1 , name2 ,'Location','northoutside','Orientation','horizontal')
    saveas(gca,[pathname 'analysis\2 y Dual-color plot region ' num2str(n) '.fig']);
    saveas(gca,[pathname 'analysis\2 y Dual-color plot region ' num2str(n) '.png']);
    close(gcf)
end
%% Channel 1 "tight" plot
% ColorMatrix = get(gca,'colororder');
intens_norm = zeros(size(intens));
intens_norm_rela = zeros(size(intens));
for n = 1:size(intens,2)
    colororder = mod(n,7);
    if colororder == 0
        colororder = 7;
    end
    intens_norm(:,n) = intens(:,n)./intens(1,n);
end
figure()
for n = 1:size(intens,2)
    colororder = mod(n,7);
    if colororder == 0
        colororder = 7;
    end
    plot(t_mov,intens_norm(:,n),'Color', getrgb(400+round((n-1)/size(intens,2)*270)));hold on
end
title(['Normalized intensity of ' num2str(name1)]);
xlabel('time (s)')
legend('Cell 1', 'Cell 2','Cell 3', 'Cell 4','Cell 5', 'Cell 6','Cell 7', 'Cell 8','Cell 9', 'Cell 10','Cell 11', 'Cell 12','Cell 13', 'Cell 14','Cell 15', 'Cell 16','Cell 17', 'Cell 18','Cell 19','Location','Best')
legend('boxoff')
ylim([0.2 max(intens_norm(:))+0.1])
box off
% axis tight
saveas(gca,[pathname 'analysis' '\norm_intensity of '  num2str(name1) '.fig']);
saveas(gca,[pathname 'analysis' '\norm_intensity of '  num2str(name1) '.png']);
figure()
for n = 1:size(intens,2)
    colororder = mod(n,7);
    if colororder == 0
        colororder = 7;
    end
    plot(t_mov,intens_norm(:,n),'Color', getrgb(400+round((n-1)/size(intens,2)*270)));hold on
end
legend('Cell 1', 'Cell 2','Cell 3', 'Cell 4','Cell 5', 'Cell 6','Cell 7', 'Cell 8','Cell 9', 'Cell 10','Cell 11', 'Cell 12','Cell 13', 'Cell 14','Cell 15','Location','Best')
legend('boxoff')
title(['Normalized intensity of ' num2str(name1)]);
xlabel('time (s)')
ylim([min(intens_norm(:))-0.05 max(intens_norm(:))+0.1])
box off
% axis tight
saveas(gca,[pathname 'analysis' '\norm_intensity of '  num2str(name1) ' tight.fig']);
saveas(gca,[pathname 'analysis' '\norm_intensity of '  num2str(name1) ' tight.png']);

figure();
plot(t_mov,mean(intens_norm,2),'Color',getrgb(wvlt1),'LineWidth',2);
title(['Mean normalized intensity of ' num2str(name1)]);
xlabel('time (s)')
ylim([min(mean(intens_norm,2))-0.05 max(mean(intens_norm,2))+0.1])
box off
saveas(gca,[pathname 'analysis' '\mean_norm_intensity of '  num2str(name1) ' tight.fig']);
saveas(gca,[pathname 'analysis' '\mean_norm_intensity of '  num2str(name1) ' tight.png']);
%% Channel 2 "tight" plot
intens2_norm = zeros(size(intens2));
intens2_norm_rela = zeros(size(intens2));
for n = 1:size(intens2,2)
    colororder = mod(n,7);
    if colororder == 0
        colororder = 7;
    end
    intens2_norm(:,n) = intens2(:,n)./intens2(1,n);
    
end
figure()
for n = 1:size(intens2,2)
    colororder = mod(n,7);
        if colororder == 0
        colororder = 7;
    end
    plot(t_mov,intens2_norm(:,n),'Color', getrgb(400+round((n-1)/size(intens,2)*270)));hold on
end
legend('Cell 1', 'Cell 2','Cell 3', 'Cell 4','Cell 5', 'Cell 6','Cell 7', 'Cell 8','Cell 9', 'Cell 10','Cell 11', 'Cell 12','Cell 13', 'Cell 14','Cell 15','Location','Best')
legend('boxoff')
title(['Normalized intensity of ' num2str(name2)]);
xlabel('time (s)')
ylim([0.1 max(intens_norm(:))+0.1])
box off
% axis tight
saveas(gca,[pathname 'analysis' '\norm_intensity of '  num2str(name2) ' tight.fig']);
saveas(gca,[pathname 'analysis' '\norm_intensity of '  num2str(name2) ' tight.png']);

figure()
for n = 1:size(intens2,2)
    colororder = mod(n,7);
    if colororder == 0
        colororder = 7;
    end
    plot(t_mov,intens2_norm(:,n),'Color', getrgb(400+round((n-1)/size(intens,2)*270)));hold on
end
legend('Cell 1', 'Cell 2','Cell 3', 'Cell 4','Cell 5', 'Cell 6','Cell 7', 'Cell 8','Cell 9', 'Cell 10','Cell 11', 'Cell 12','Cell 13', 'Cell 14','Cell 15','Location','Best')
legend('boxoff')
title(['Normalized intensity of ' num2str(name2)]);
xlabel('time (s)')
ylim([min(intens2_norm(:))-0.05 max(intens2_norm(:))+0.1])
box off
% axis tight
saveas(gca,[pathname 'analysis' '\norm_intensity of '  num2str(name2) ' tight.fig']);
saveas(gca,[pathname 'analysis' '\norm_intensity of '  num2str(name2) ' tight.png']);


figure();
plot(t_mov,mean(intens2_norm,2),'Color',getrgb(wvlt2),'LineWidth',2);
title(['Mean normalized intensity of ' num2str(name2)]);
xlabel('time (s)')
ylim([min(mean(intens2_norm,2))-0.05 max(mean(intens2_norm,2))+0.1])
box off
saveas(gca,[pathname 'analysis' '\mean_norm_intensity of '  num2str(name2) ' tight.fig']);
saveas(gca,[pathname 'analysis' '\mean_norm_intensity of '  num2str(name2) ' tight.png']);
%% intens1 (without bkg)
mkdir([pathname 'analysis\without background\'])
% ColorMatrix = get(gca,'colororder');
intens_norm_wobkg = zeros(size(intens));
intens_norm_rela_wobkg = zeros(size(intens));
for n = 1:size(intens,2)
    colororder = mod(n,7);
    if colororder == 0
        colororder = 7;
    end
    intens_norm_wobkg(:,n) = (intens(:,n)-mean(intens_forbkg,1)')./(intens(1,n)-intens_forbkg(1,1));
end
figure()
for n = 1:size(intens,2)
    colororder = mod(n,7);
    if colororder == 0
        colororder = 7;
    end
    plot(t_mov,intens_norm_wobkg(:,n),'Color', getrgb(400+round((n-1)/size(intens,2)*270)));hold on
end
legend('Cell 1', 'Cell 2','Cell 3', 'Cell 4','Cell 5', 'Cell 6','Cell 7', 'Cell 8','Cell 9', 'Cell 10','Cell 11', 'Cell 12','Cell 13', 'Cell 14','Cell 15','Location','Best')
legend('boxoff')
title(['Normalized intensity of ' num2str(name1) ' (without background)']);
xlabel('time (s)')
ylim([0.2 max(intens_norm_wobkg(:))+0.1])
box off
% axis tight
saveas(gca,[pathname 'analysis\without background\norm_intensity of '  num2str(name1) '_w_o bkg.fig']);
saveas(gca,[pathname 'analysis\without background\norm_intensity of '  num2str(name1) '_w_o bkg.png']);

figure()
for n = 1:size(intens,2)
    colororder = mod(n,7);
    if colororder == 0
        colororder = 7;
    end
    plot(t_mov,intens_norm_wobkg(:,n),'Color', getrgb(400+round((n-1)/size(intens,2)*270)));hold on
end
legend('Cell 1', 'Cell 2','Cell 3', 'Cell 4','Cell 5', 'Cell 6','Cell 7', 'Cell 8','Cell 9', 'Cell 10','Cell 11', 'Cell 12','Cell 13', 'Cell 14','Cell 15','Location','Best')
legend('boxoff')
title(['Normalized intensity of ' num2str(name1) ' (without background)']);
xlabel('time (s)')
ylim([min(intens_norm_wobkg(:))-0.05 max(intens_norm_wobkg(:))+0.1])
box off
% axis tight
saveas(gca,[pathname 'analysis\without background\norm_intensity of '  num2str(name1) ' tight_w_o bkg.fig']);
saveas(gca,[pathname 'analysis\without background\norm_intensity of '  num2str(name1) ' tight_w_o bkg.png']);

figure();
plot(t_mov,mean(intens_norm_wobkg,2),'Color',getrgb(wvlt1),'LineWidth',2);
title(['Mean normalized intensity of ' num2str(name1) ' (without background)']);
xlabel('time (s)')
ylim([min(mean(intens_norm_wobkg,2))-0.05 max(mean(intens_norm_wobkg,2))+0.1])
box off
saveas(gca,[pathname 'analysis\without background\mean_norm_intensity of '  num2str(name1) ' tight_w_o bkg.fig']);
saveas(gca,[pathname 'analysis\without background\mean_norm_intensity of '  num2str(name1) ' tight_w_o bkg.png']);
%% intens2 (without bkg)
intens2_norm_wobkg = zeros(size(intens2));
% intens2_norm_rela = zeros(size(intens2));
for n = 1:size(intens2,2)
    colororder = mod(n,7);
    if colororder == 0
        colororder = 7;
    end
    intens2_norm_wobkg(:,n) = (intens2(:,n)-mean(intens2_bkg,1)')./(intens2(1,n)-intens2_bkg(1,1));
    
end
figure()
for n = 1:size(intens2,2)
    colororder = mod(n,7);
    if colororder == 0
        colororder = 7;
    end
    plot(t_mov,intens2_norm_wobkg(:,n),'Color', getrgb(400+round((n-1)/size(intens,2)*270)));hold on
end
legend('Cell 1', 'Cell 2','Cell 3', 'Cell 4','Cell 5', 'Cell 6','Cell 7', 'Cell 8','Cell 9', 'Cell 10','Cell 11', 'Cell 12','Cell 13', 'Cell 14','Cell 15','Location','Best')
legend('boxoff')
title(['Normalized intensity of ' num2str(name2) ' without background']);
xlabel('time (s)')
ylim([0.1 max(intens2_norm_wobkg(:))+0.1])
box off
% axis tight
saveas(gca,[pathname 'analysis\without background\norm_intensity of '  num2str(name2) ' tight.fig']);
saveas(gca,[pathname 'analysis\without background\norm_intensity of '  num2str(name2) ' tight.png']);

figure()
for n = 1:size(intens2,2)
    colororder = mod(n,7);
    if colororder == 0
        colororder = 7;
    end
    plot(t_mov,intens2_norm_wobkg(:,n),'Color', getrgb(400+round((n-1)/size(intens,2)*270)));hold on
end
legend('Cell 1', 'Cell 2','Cell 3', 'Cell 4','Cell 5', 'Cell 6','Cell 7', 'Cell 8','Cell 9', 'Cell 10','Cell 11', 'Cell 12','Cell 13', 'Cell 14','Cell 15','Location','Best')
legend('boxoff')
title(['Normalized intensity of ' num2str(name2) ' with out background']);
xlabel('time (s)')
ylim([min(intens2_norm_wobkg(:))-0.05 max(intens2_norm_wobkg(:))+0.1])
box off
% axis tight
saveas(gca,[pathname 'analysis\without background\norm_intensity of '  num2str(name2) ' tight.fig']);
saveas(gca,[pathname 'analysis\without background\norm_intensity of '  num2str(name2) ' tight.png']);


figure();
plot(t_mov,mean(intens2_norm_wobkg,2),'Color',getrgb(wvlt2),'LineWidth',2);
title(['Mean normalized intensity of ' num2str(name2) ' without background']);
xlabel('time (s)')
ylim([min(mean(intens2_norm_wobkg,2))-0.05 max(mean(intens2_norm_wobkg,2))+0.1])
box off
saveas(gca,[pathname 'analysis\without background\mean_norm_intensity of '  num2str(name2) ' tight.fig']);
saveas(gca,[pathname 'analysis\without background\mean_norm_intensity of '  num2str(name2) ' tight.png']);

intens_std = std(intens_norm_wobkg')';
intens2_std = std(intens2_norm_wobkg')';

intens_sem = intens_std./sqrt(size(intens,2));
intens2_sem = intens2_std./sqrt(size(intens2,2));

%% 1st order of difference (intergrated)
figure()
for n = 1:size(intens2,2)
    colororder = mod(n,7);
    if colororder == 0
        colororder = 7;
    end
    plot(t_mov(1:end-1),diff(intens2_norm(:,n),1),'Color', getrgb(400+round((n-1)/size(intens,2)*270)));hold on
end
legend('Cell 1', 'Cell 2','Cell 3', 'Cell 4','Cell 5', 'Cell 6','Cell 7', 'Cell 8','Cell 9', 'Cell 10','Cell 11', 'Cell 12','Cell 13', 'Cell 14','Cell 15','Location','Best')
legend('boxoff')
title(['1st order of difference: ' num2str(name2)]);
xlabel('time (s)')
% ylim([min(intens2_norm(:))-0.05 max(intens2_norm(:))+0.1])
box off
axis tight
saveas(gca,[pathname 'analysis' '\1st order of difference_' num2str(name2) ' tight.fig']);
saveas(gca,[pathname 'analysis' '\1st order of difference_' num2str(name2) ' tight.png']);



%% 1st order of difference
figure()
for n = 1:size(intens,2)
    colororder = mod(n,7);
    if colororder == 0
        colororder = 7;
    end
    plot(t_mov(1:end-1),diff(intens_norm(:,n),1),'Color', getrgb(400+round((n-1)/size(intens,2)*270)));hold on
end
legend('Cell 1', 'Cell 2','Cell 3', 'Cell 4','Cell 5', 'Cell 6','Cell 7', 'Cell 8','Cell 9', 'Cell 10','Cell 11', 'Cell 12','Cell 13', 'Cell 14','Cell 15','Location','Best')
legend('boxoff')
title(['1st order of difference: ' num2str(name1)]);
xlabel('time (s)')
% ylim([min(intens2_norm(:))-0.05 max(intens2_norm(:))+0.1])
box off
axis tight
saveas(gca,[pathname 'analysis' '\1st order of difference_' num2str(name1) ' tight.fig']);
saveas(gca,[pathname 'analysis' '\1st order of difference_' num2str(name1) ' tight.png']);
%% Export to EXCEL
% xlswrite([pathname 'data_analysis.xlsx'],{'Cdensity', 'Current','Threshold Potential', 'Maximum Rising Vm','Spike Of Each Injection', 'Resting Potential','Spike Number','Max Rising Speed','Amplitude','Vpeak','FWHM','APD_threshold'},'AIS','A1');
for n = 1:size(intens,2)
    xlswrite([pathname 'data_analysis.xlsx'],{['Cell ' num2str(n)]},name1,[char(65+n) '1']);
end

for n = 1:size(intens,2)
    xlswrite([pathname 'data_analysis.xlsx'],{['Cell ' num2str(n)]},name2,[char(65+n) '1']);
end
xlswrite([pathname 'data_analysis.xlsx'],{'Standard dev.'},name1,[char(66+n) '1']);
xlswrite([pathname 'data_analysis.xlsx'],{'Standard dev.'},name2,[char(66+n) '1']);

xlswrite([pathname 'data_analysis.xlsx'],intens_std,name1,[char(66+n) '2']);
xlswrite([pathname 'data_analysis.xlsx'],intens2_std,name2,[char(66+n) '2']);

xlswrite([pathname 'data_analysis.xlsx'],{'S.E.M.'},name1,[char(67+n) '1']);
xlswrite([pathname 'data_analysis.xlsx'],{'S.E.M.'},name2,[char(67+n) '1']);

xlswrite([pathname 'data_analysis.xlsx'],{['Average, n = ' num2str(size(intens,2))]},name1,[char(68+n) '1']);
xlswrite([pathname 'data_analysis.xlsx'],{['Average, n = ' num2str(size(intens,2))]},name2,[char(68+n) '1']);


xlswrite([pathname 'data_analysis.xlsx'],intens_sem,name1,[char(67+n) '2']);
xlswrite([pathname 'data_analysis.xlsx'],intens2_sem,name2,[char(67+n) '2']);

xlswrite([pathname 'data_analysis.xlsx'],mean(intens_norm_wobkg,2),name1,[char(68+n) '2']);
xlswrite([pathname 'data_analysis.xlsx'],mean(intens2_norm_wobkg,2),name2,[char(68+n) '2']);


xlswrite([pathname 'data_analysis.xlsx'],{'Time (s)'},name1,'A1');
xlswrite([pathname 'data_analysis.xlsx'],{'Time (s)'},name2,'A1');

xlswrite([pathname 'data_analysis.xlsx'],t_mov',name1,'A2');
xlswrite([pathname 'data_analysis.xlsx'],t_mov',name2,'A2');

xlswrite([pathname 'data_analysis.xlsx'],intens_norm_wobkg,name1,'B2');
xlswrite([pathname 'data_analysis.xlsx'],intens2_norm_wobkg,name2,'B2');

%% dual plot of averge dynamic intensity with shadow
figure()
seShoadow1 = area(t_mov, [(mean(intens_norm_wobkg, 2) - intens_sem), (2 * intens_sem)]);
hold on;
set(seShoadow1(1),'Visible','off');
set(seShoadow1(2),'FaceColor',[0.85, 0.32, 0]' ,'FaceAlpha',0.5, 'EdgeColor', 'none');
seShoadow2 = area(t_mov, [(mean(intens2_norm_wobkg, 2) - intens2_sem), (2 * intens2_sem)]);
hold on;
set(seShoadow2(1),'Visible','off');
set(seShoadow2(2),'FaceColor',[0.46, 0.68, 0.17]' ,'FaceAlpha',0.5, 'EdgeColor', 'none');
plot(t_mov, mean(intens_norm_wobkg,2), 'Color', name1color,'LineWidth',1.5);
plot(t_mov, mean(intens2_norm_wobkg,2), 'Color', name2color,'LineWidth',1.5);
legend(name1, name2, 'Location','northoutside','Orientation','horizontal');
box off;
axis tight;
xlabel('Time(s)');
ylabel('Relative intensity')
xlabel('Time (ms)');
%ylim([0,1]);
ylabel('Normalized Intensity');
%saveas(gca,[pathname '\shadow_mean_norm_intensity of mNeon-CAAX tight.fig']);
%saveas(gca,[pathname '\shadow_mean_norm_intensity of mNeon-CAAX tight.png']);

%% ratio calculation
mkdir([pathname '\Ratio analysis\'])
intens_ratio = (intens)./(intens2);% EM630/EM585; Higher pH, higher ratio;
intens_ratio_wobkg = (intens-intens_forbkg')./(intens2-intens2_bkg');% EM630/EM585; Higher pH, higher ratio;

figure()
for n = 1:size(intens,2)
    colororder = mod(n,7);
    if colororder == 0
        colororder = 7;
    end
    plot(t_mov,intens_ratio(:,n),'Color', getrgb(400+round((n-1)/size(intens,2)*270)));hold on
end
legend('Cell 1', 'Cell 2','Cell 3', 'Cell 4','Cell 5', 'Cell 6','Cell 7', 'Cell 8','Cell 9', 'Cell 10','Cell 11', 'Cell 12','Cell 13', 'Cell 14','Cell 15','Location','Best')
legend('boxoff')
title(['Intensitiy ratio: ' num2str(wvlt1) ' / ' num2str(wvlt2)]);
xlabel('time (s)')
ylabel(['Em. ' num2str(wvlt1) ' / ' num2str(wvlt2)])
box off
axis tight
saveas(gca,[pathname '\Ratio analysis\intensity ratio' num2str(wvlt1) '_' num2str(wvlt2) '_tight.fig']);
saveas(gca,[pathname '\Ratio analysis\intensity ratio' num2str(wvlt1) '_' num2str(wvlt2) '_tight.png']);

figure()
for n = 1:size(intens2,2)
    colororder = mod(n,7);
    if colororder == 0
        colororder = 7;
    end
    plot(t_mov,intens_ratio_wobkg(:,n),'Color', getrgb(400+round((n-1)/size(intens,2)*270)));hold on
end
legend('Cell 1', 'Cell 2','Cell 3', 'Cell 4','Cell 5', 'Cell 6','Cell 7', 'Cell 8','Cell 9', 'Cell 10','Cell 11', 'Cell 12','Cell 13', 'Cell 14','Cell 15','Location','Best')
legend('boxoff')
title(['Intensitiy ratio: ' num2str(wvlt1) ' / ' num2str(wvlt2)]);
xlabel('time (s)')
ylabel(['Em. ' num2str(wvlt1) ' / ' num2str(wvlt2) ', w/o bkg.'])
box off
axis tight
saveas(gca,[pathname '\Ratio analysis\intensity ratio' num2str(wvlt1) '_' num2str(wvlt2) '_tight_without bkg.fig']);
saveas(gca,[pathname '\Ratio analysis\intensity ratio' num2str(wvlt1) '_' num2str(wvlt2) '_tight_without bkg.png']);

%
k = 1;

for n = 438.5:0.5:488
    
    ratio_colormap(k,:) =getrgb(n);% from blue to green
    k=k+1;
end

for n = 488:0.05:557.95
    
    ratio_colormap(k,:) =getrgb(n);% from blue to green
    k=k+1;
end

figure()
set(figure(),'Position',get(0,'screensize'))
subplot(3,1,1)
imshow(mov(:,:,1),[500,20000])
ax1 = subplot(3,1,1);
colormap(ax1,pseudocolor(wvlt1))
colorbar
title(['1st frame: EM. ' num2str(wvlt1)])
subplot(3,1,2)
imshow(mov2(:,:,1),[500,5000])
ax2 = subplot(3,1,2);
colormap(ax2,pseudocolor(wvlt2))
colorbar
title(['1st frame: EM. ' num2str(wvlt2)])
subplot(3,1,3)
imshow((mov(:,:,1)-400)./(mov2(:,:,1)-400),[0.01 15]);
ax3 = subplot(3,1,3);
colormap(ax3,ratio_colormap)
title(['Ratio image: ' num2str(wvlt1) ' / ' num2str(wvlt2)]);
colorbar
title(colorbar,['Em. ' num2str(wvlt1) ' / ' num2str(wvlt2)])
saveas(gca,[pathname '\Ratio analysis\1st frame_Snapshot.fig']);
saveas(gca,[pathname '\Ratio analysis\1st frame_Snapshot.png']);

figure()
set(figure(),'Position',get(0,'screensize'))
subplot(3,1,1)
imshow(mov(:,:,end),[500,20000])
ax1 = subplot(3,1,1);
colormap(ax1,pseudocolor(wvlt1))
colorbar
title(['Last frame: EM. ' num2str(wvlt1)])
subplot(3,1,2)
imshow(mov2(:,:,end),[500,5000])
ax2 = subplot(3,1,2);
colormap(ax2,pseudocolor(wvlt2))
title(['Last frame: EM. ' num2str(wvlt2)])
colorbar
subplot(3,1,3)
imshow((mov(:,:,end)-400)./(mov2(:,:,end)-400),[0.01 15]);
ax3 = subplot(3,1,3);
colormap(ax3,ratio_colormap)
colorbar
title(['Ratio image: ' num2str(wvlt1) ' / ' num2str(wvlt2)]);
colorbar
title(colorbar,['Em. ' num2str(wvlt1) ' / ' num2str(wvlt2)])
saveas(gca,[pathname '\Ratio analysis\last frame_Snapshot.fig']);
saveas(gca,[pathname '\Ratio analysis\last frame_Snapshot.png']);
%%
mov_ratio = (mov-100*power(bin,1))./(mov2-100*power(bin,1));
clear M
figure(20)
for j = 1:(size(mov_ratio,3))
    imshow(mov_ratio(:,:,j),[0.01 15], 'Border', 'tight');hold on;
    colormap(ratio_colormap);
    text(10,30,[num2str(t_mov(j)) ' s'], 'FontSize', 18, 'color', [0.99 0 0])
    %text(550,20,[name2], 'FontSize', 20, 'color', getrgb(598));
    text(10,nrow-20,[ num2str(wvlt1) ' / ' num2str(wvlt2)], 'FontSize', 18, 'color', [0.99 0 0]);
    colorbar
    title(colorbar,['Em. ' num2str(wvlt1) ' / ' num2str(wvlt2)])
    M(j) = getframe(gcf);
    %     saveas(gca, ['SNAPT_movie_' num2str((j)) '.tif']);  % Use this line to save each frame as a separate TIF file.
end
close(gcf)
myObj = VideoWriter([pathname  '\Ratio analysis\'  num2str(wvlt2) '_' num2str(wvlt1) '_colorbar.avi'],'Uncompressed AVI');%create an AVI file
myObj.FrameRate = 5;
open(myObj);
writeVideo(myObj,M)
close(myObj);
myObj = VideoWriter([pathname  '\Ratio analysis\'  num2str(wvlt2) '_' num2str(wvlt1)  '_compressed_H.264 encoding_colorbar.mp4'],'MPEG-4');%create an AVI file
myObj.FrameRate = 5;
open(myObj);
writeVideo(myObj,M)
close(myObj);

% clear M
% figure(20)
% for j = 1:(size(mov_ratio,3))
%     imshow(mov_ratio(:,:,j),[0.4 1], 'Border', 'tight');hold on;
%     colormap jet;
%     text(10,30,[num2str(t_mov(j)) ' s'], 'FontSize', 18, 'color', [0.99 0.99 0.99])
%     %text(550,20,[name2], 'FontSize', 20, 'color', getrgb(598));
%     text(10,nrow-20,[ num2str(wvlt2) ' / ' num2str(wvlt1)], 'FontSize', 18, 'color', [0.99 0.99 0.99]);
% %     colorbar
%     M(j) = getframe(gcf);
%     %     saveas(gca, ['SNAPT_movie_' num2str((j)) '.tif']);  % Use this line to save each frame as a separate TIF file.
% end
% close(gcf)
% myObj = VideoWriter([pathname  '\Ratio analysis\'  num2str(wvlt2) '_' num2str(wvlt1) ' .avi'],'Uncompressed AVI');%create an AVI file
% myObj.FrameRate = 5;
% open(myObj);
% writeVideo(myObj,M)
% close(myObj);
% myObj = VideoWriter([pathname  '\Ratio analysis\'  num2str(wvlt2) '_' num2str(wvlt1)  '_compressed_H.264 encoding.mp4'],'MPEG-4');%create an AVI file
% myObj.FrameRate = 5;
% open(myObj);
% writeVideo(myObj,M)
% close(myObj);
save([pathname name1 '_' name2 ' .mat']);% save the variant



