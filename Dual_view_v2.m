% Analyze response of voltage sensors in two different channel with dual-view imaging
% Created by Luxing Peng on 20210716
% Modified by Junqi Yang on 20210717 and 20210822
%%

clear all;
%parameters
name1 = 'Dual-color';
% name2 = 'mito-pHTomato';
name_LEFT = 'ArcLight';
name_RIGHT = 'HASAP1-JF635';

WL_LEFT = 525;
WL_RIGHT = 655;

% WL_RIGHT = 600;
nameLEFTcolor = getrgb(WL_LEFT);
nameRIGHTcolor = getrgb(WL_RIGHT);


subfolder = '';
basepath = 'D:\data\20x old\20x dry\';
pathname = [basepath subfolder];
[file_list] = RScanDir(pathname,  '\matlab variables.mat')
wavename = file_list{1}(length(basepath):end-length('matlab variables.mat'))
load (file_list{1});
mkdir([pathname 'analysis']);

% constants
movname = '\movie.bin';
fileID = fopen([pathname 'movie_info.txt']);
movie_info = textscan(fileID,'%q');
ncol = str2num(movie_info{1,1}{14,1}) ;         % x
nrow = str2num(movie_info{1,1}{17,1}) ;          % y
bin = str2num(movie_info{1,1}{28,1});
bkg = 100*power(bin,2);          % background due to camera bias (100 for bin 1x1)
dt_mov = str2num(movie_info{1,1}{32,1});         % exposure time in millisecond
DAQname = '\movie_DAQ.txt';
dt_daq = 1000./samprate;       % in millisecond
dnsamp = dt_mov./dt_daq;        % downsampling rate = DAQ rate/camera rate

% load movie
fname = [pathname movname]
[mov, nframe] = readBinMov(fname, nrow, ncol);
%% Extra bin?
% [mov] = extrabin(mov, 1);
% bin = bin + 1;
% ncol = size(mov,2)  ;         % x
% nrow = size(mov,1) ;          % y
% bkg = 100*power(bin,2);          % background due to camera bias (100 for bin 1x1)

%%
%mov = double(mov);
%toc;
img = mean(mov, 3);
len = size(mov,3);
t_mov = [0:(len-1)]*dt_mov/1000;     % time axis in second
%img = mean(mov(:,:,1:1100), 3);
% select ROI for analysis
intens = squeeze(mean(mean(mov)));
intens = intens-bkg;
figure();
plot(t_mov,intens);
box off
axis tight
xlabel('Time (s)');
ylabel('Intensity');
saveas(gca,[pathname 'analysis' '\global intensity.fig']);
saveas(gca,[pathname 'analysis' '\global intensity.png']);

%% Pick up the key frames
threshold_down = 50; % related to the global intens, Need to be auto
%threshold_up = 200;
%intens(intens>=threshold_up)=0
% threshold_pixel = 15; % for clicky: upper
[~,peak] = findpeaks(intens,'MinPeakDistance',2,'MinPeakHeight',threshold_down);
t_keyframe = t_mov(peak); 
% peak = spikefind(intens,threshold_down);
mov2 = double([]);
intens_light = [];
% t_keyframe = t_keyframe(2:end);
% peak = peak(2:end)
for p = 1:length(t_keyframe)
    %     intens_light(p) = intens(peak(p)-1)+intens(peak(p))+intens(peak(p)+1)+intens(peak(p)+2)-4*mean(intens([peak(p)-6:peak(p)-2 peak(p)+3:peak(p)+8]));% without background
    %     mov2(:,:,p)=double(mov(:,:,peak(p))+mov(:,:,peak(p)+1)+mov(:,:,peak(p)-1)+mov(:,:,peak(p)+2))-4*mean(mov(:,:,[peak(p)-6:peak(p)-2 peak(p)+3:peak(p)+7]),3);
    intens_light(p) = intens(peak(p)-1)+intens(peak(p))+intens(peak(p)+1)+intens(peak(p)+2)+intens(peak(p)+3)-5*mean(intens(20:40));% without background
    mov2(:,:,p)=double(mov(:,:,peak(p)-1)+mov(:,:,peak(p))+mov(:,:,peak(p)+1)+mov(:,:,peak(p)+2)+mov(:,:,peak(p)+3))-5*mean(mov(:,:,20:40),3);
    %     intens_light(p) = intens(peak(p)-1)+intens(peak(p))+intens(peak(p)+1)-3*mean(intens(20:40));% without background
    %     mov2(:,:,p)=double(mov(:,:,peak(p)-1)+mov(:,:,peak(p))+mov(:,:,peak(p)+1))-3*mean(mov(:,:,20:40),3);
    
    
    % if Lasers_AO_DO(4,(round(t_keyframe(p)*samprate)))==1;%4 = elimanating points with blue-light illumination
    %         intens_light(p) = 0;
    %         mov2(:,:,p)= 0;
    %     end
    %     if Lasers_AO_DO(4,(ceil(t_keyframe(p)*samprate)))==1;%4 = elimanating points with blue-light illumination
    %         intens_light(p) = 0;
    %         mov2(:,:,p)= 0;
    %     end
end
%elimanating points with blue-light illumination in time and intensity
%sequence
temp3 = find(intens_light==0);
t_keyframe(temp3)=[];
intens_light(temp3)=[];
mov2(:,:,temp3)=[];
clear temp3

img2 = mean(mov2,3);
t_loc = spikefind(intens,threshold_down);
plot(t_keyframe,intens_light,'.');% without background
box off
axis tight
xlabel('time (s)')
ylabel('Intensity (W/O camera bias)')
title('global intensity with excitation')
saveas(gca,[pathname 'analysis' '\global dual color intensity.fig']);
saveas(gca,[pathname 'analysis' '\global dual color intensity.png']);
figure()
plot(t_keyframe,intens_light./(intens_light(1,1)),'.');% without background
box off
axis tight
xlabel('time (s)')
ylabel('Normalized intensity (W/O camera bias)')
title('Normalized global intensity of excitation')
saveas(gca,[pathname 'analysis' '\Normalized dual color intensity.fig']);
saveas(gca,[pathname 'analysis' '\Normalized dual color intensity.png']);
%% Select the region of mito in the cells with expression and non-expression cells
cmin_RIGHT = 50;
cmax_RIGHT = 1000;
[ROI_points_RIGHT,intens_RIGHT] = clicky_v2(mov2, imaging_flash_interval, 0, img2, cmin_RIGHT , cmax_RIGHT, '');
saveas(gca,[pathname 'analysis' '\Expressing Cells_RIGHT.fig']);
saveas(gca,[pathname 'analysis' '\Expressing Cells_RIGHT.png']);
cmin_LEFT = 50;
cmax_LEFT = 2000;
[ROI_points_LEFT,intens_LEFT] = clicky_v2(mov2, imaging_flash_interval, 0, img2, cmin_LEFT, cmax_LEFT, '');
saveas(gca,[pathname 'analysis' '\Expressing Cells_LEFT.fig']);
saveas(gca,[pathname 'analysis' '\Expressing Cells_LEFT.png']);


%% Select the region of cyto in the cells with expression and non-expression cells
%
% cmin = 100;
% cmax = 4000;
% [ROI_points_notmito,intens_region_notmito] = clicky_v2(mov2, imaging_flash_interval, 0, img2, cmin , cmax, '')
% saveas(gca,[pathname 'analysis' '\Expressing Cells_not mito.fig']);
% saveas(gca,[pathname 'analysis' '\Expressing Cells_not mito.png']);
% cmin = 100;
% cmax = 4000;
% [ROI_points_neg_notmito,intens_negative_notmito] = clicky_v2(mov2, imaging_flash_interval, 0, img2, cmin , cmax, '')
% saveas(gca,[pathname 'analysis' '\NO_Expressing Cells_not mito.fig']);
% saveas(gca,[pathname 'analysis' '\NO_Expressing Cells_not mito.png']);
%% Display the ROI of 4-time selection
% merge_mask = imread(['X:\91 Data and analysis\MRR\20210414 MRR\HEK 293T cells\pc3.1-Kir2.1-P2A-CheRiff-EGFP_TMRM\Dish 2\Cell 3\151703_after stimulation\Composite.jpg']);
img2_LEFT = img2(:,1:size(img2,2)/2,:);
img2_RIGHT = img2(:,size(img2,2)/2+1:end,:);
figure();
set(gcf,'outerposition',get(0,'screensize'));
subplot(1,2,1)
imshow(img2_LEFT,[cmin_LEFT cmax_LEFT]);hold on;
h =subplot(1,2,1);
colormap(h,pseudocolor(WL_LEFT))
title(h,name_LEFT)
for region_num=1:max(size(ROI_points_LEFT))
    plot([ROI_points_LEFT{1,region_num}(:,1)],[ROI_points_LEFT{1,region_num}(:,2)],'LineWidth',3,'Color', getrgb(400+round((region_num-1)/size(intens_LEFT,2)*270)));hold on
    text(mean(ROI_points_LEFT{1,region_num}(:,1)),mean(ROI_points_LEFT{1,region_num}(:,2)),num2str(region_num),'Color',[0.99 0.99 0.99],'FontSize',20);hold on;
    plot([8 69.53],[nrow-20 nrow-20],'LineWidth',5,'Color',[0.99 0.99 0.99]);
end
colorbar
subplot(1,2,2)
imshow(img2_RIGHT,[cmin_RIGHT cmax_RIGHT]);hold on;
for region_num=1:max(size(ROI_points_LEFT))
    plot([ROI_points_RIGHT{1,region_num}(:,1)]-size(img,2)/2,[ROI_points_RIGHT{1,region_num}(:,2)],'LineWidth',3,'Color', getrgb(400+round((region_num-1)/size(intens_RIGHT,2)*270)));hold on
    text(mean(ROI_points_RIGHT{1,region_num}(:,1))-size(img,2)/2,mean(ROI_points_RIGHT{1,region_num}(:,2)),num2str(region_num),'Color',[0.99 0.99 0.99],'FontSize',20);hold on;
    plot([8 69.53],[nrow-20 nrow-20],'LineWidth',5,'Color',[0.99 0.99 0.99]);
end
h =subplot(1,2,2);
colormap(h,pseudocolor(WL_RIGHT))
title(h,name_RIGHT)
colorbar

% colormap(pseudocolor(WL_LEFT));

% title(['Selected regions: LEFT: ' name_LEFT ', RIGHT: ' name_LEFT])
saveas(gca,[pathname 'analysis\selected regions of ' name1 '.fig']);
saveas(gca,[pathname 'analysis\selected regions of ' name1 '.png']);
close(gcf)


% figure();
% set(gcf,'outerposition',get(0,'screensize'));
% imshow(merge_mask,[cmin cmax]);hold on;
% colormap(pseudocolor(WL_LEFT));
% for region_num=1:max(size(ROI_points_notmito))
%     plot([ROI_points_notmito{1,region_num}(:,1)],[ROI_points_notmito{1,region_num}(:,2)],'LineWidth',3,'Color', getrgb(400+round((region_num-1)/size(intens_region,2)*270)));hold on
%     text(mean(ROI_points_notmito{1,region_num}(:,1)),mean(ROI_points_notmito{1,region_num}(:,2)),num2str(region_num),'Color',[0.99 0.99 0.99],'FontSize',20);hold on;
%     plot([8 69.53],[nrow-20 nrow-20],'LineWidth',5,'Color',[0.99 0.99 0.99]);
% end
% saveas(gca,[pathname 'analysis\selected regions of ' name1 '_not mito.fig']);
% saveas(gca,[pathname 'analysis\selected regions of ' name1 '_not mito.png']);
% close(gcf)
%
% figure();
% set(gcf,'outerposition',get(0,'screensize'));
% imshow(merge_mask,[cmin cmax]);hold on;
% colormap(pseudocolor(WL_LEFT));
% for region_num=1:max(size(ROI_points_neg_notmito))
%     plot([ROI_points_neg_notmito{1,region_num}(:,1)],[ROI_points_neg_notmito{1,region_num}(:,2)],'LineWidth',3,'Color', getrgb(400+round((region_num-1)/size(intens_region,2)*270)));hold on
%     text(mean(ROI_points_neg_notmito{1,region_num}(:,1)),mean(ROI_points_neg_notmito{1,region_num}(:,2)),num2str(region_num),'Color',[0.99 0.99 0.99],'FontSize',20);hold on;
%     plot([8 69.53],[nrow-20 nrow-20],'LineWidth',5,'Color',[0.99 0.99 0.99]);
% end
% saveas(gca,[pathname 'analysis\selected regions of ' name1 ' negative_not mito.fig']);
% saveas(gca,[pathname 'analysis\selected regions of ' name1 ' negative_not mito.png']);
% close(gcf)


%%
[roi_points_background,background,~,region_square] = clicky_fullscreen(mov2, 0,img2, cmin_LEFT , cmax_LEFT, '', WL_LEFT,WL_LEFT);
saveas(gca,[pathname 'analysis\selected regions of background.fig']);
saveas(gca,[pathname 'analysis\selected regions of background.png']);
for n = 1:size(background,2)/2
    background_LEFT(:,n) = background(:,n)*region_square(n);
end
background_LEFT = sum(background_LEFT,2)./sum(region_square(1:size(background,2)/2));

for n = size(background,2)/2+1:size(background,2)
    background_RIGHT(:,n) = background(:,n)*region_square(n);
end
background_RIGHT = sum(background_RIGHT,2)./sum(region_square(size(background,2)/2+1:size(background,2)));

intens_LEFT_noBKG = intens_LEFT - background_LEFT;
intens_RIGHT_noBKG = intens_RIGHT- background_RIGHT;
%% Normalized intensity
figure()
set(gcf,'outerposition',get(0,'screensize'));
subplot(1,2,1)
for  region_num=1:max(size(ROI_points_LEFT))
    plot(t_keyframe,intens_LEFT_noBKG(:,region_num)./intens_LEFT_noBKG(1,region_num),'.','Color', getrgb(400+round((region_num-1)/size(intens_LEFT_noBKG,2)*270)),'MarkerSize',12);hold on
end
box off
% axis ([0 max(t) 0.90 1.30])
axis tight
title({'Normalized intensity of ' name_LEFT })
xlabel('Time (s)')
ylabel('Normalized intensity')
legend('Cell 1', 'Cell 2','Cell 3', 'Cell 4','Cell 5', 'Cell 6','Cell 7', 'Cell 8','Cell 9', 'Cell 10','Cell 11', 'Cell 12','Cell 13', 'Cell 14','Cell 15', 'Cell 16','Cell 17', 'Cell 18','Cell 19','Location','best')
legend('boxoff')


subplot(1,2,2)
for  region_num=1:max(size(ROI_points_RIGHT))
    plot(t_keyframe,intens_RIGHT_noBKG(:,region_num)./intens_RIGHT_noBKG(1,region_num),'.','Color', getrgb(400+round((region_num-1)/size(intens_RIGHT_noBKG,2)*270)),'MarkerSize',12);hold on
end
box off
% axis ([0 max(t) 0.90 1.30])
axis tight

title({'Normalized intensity of ' name_RIGHT })
xlabel('Time (s)')
ylabel('Normalized intensity')
legend('Cell 1', 'Cell 2','Cell 3', 'Cell 4','Cell 5', 'Cell 6','Cell 7', 'Cell 8','Cell 9', 'Cell 10','Cell 11', 'Cell 12','Cell 13', 'Cell 14','Cell 15', 'Cell 16','Cell 17', 'Cell 18','Cell 19','Location','best')
legend('boxoff')

% plot([61 119],[1.05 1.05],'LineWidth',3,'Color',getrgb(488))
saveas(gca,[pathname 'analysis\Normalized intensity of ' name1 '.fig']);
saveas(gca,[pathname 'analysis\Normalized intensity of ' name1 '.png']);
close(gcf)


%%
%
% figure()
% for  region_num=1:max(size(ROI_points_notmito))
% plot(t_keyframe,intens_region_notmito(:,region_num)./intens_region_notmito(1,region_num),'o','Color', getrgb(400+round((region_num-1)/size(intens_region,2)*270)),'MarkerSize',3);hold on
%
% end
% % plot([61 119],[1.05 1.05],'LineWidth',3,'Color',getrgb(488))
% box off
% axis tight
% title({'Normalized intensity of TMRM from expressing cells (not mito regions)';'Blue light: 0.01 s, cycle length = 0.1s'})
% xlabel('Time (s)')
% ylabel('Normalized intensity')
% legend('Cell 1', 'Cell 2','Cell 3', 'Cell 4','Cell 5', 'Cell 6','Cell 7', 'Cell 8','Cell 9', 'Cell 10','Cell 11', 'Cell 12','Cell 13', 'Cell 14','Cell 15', 'Cell 16','Cell 17', 'Cell 18','Cell 19','Location','Best')
% legend('boxoff')
% saveas(gca,[pathname 'analysis\Normalized intensity of ' name1 '_not mito.fig']);
% saveas(gca,[pathname 'analysis\Normalized intensity of ' name1 '_not mito.png']);
%
%
% close(gcf)
%
% figure()
% for  region_num=1:max(size(ROI_points_neg_notmito))
% plot(t_keyframe,intens_negative_notmito(:,region_num)./intens_negative_notmito(1,region_num),'o','Color', getrgb(400+round((region_num-1)/size(intens_negative,2)*270)),'MarkerSize',3);hold on
% end
% % plot([61 119],[1.05 1.05],'LineWidth',3,'Color',getrgb(488))
% axis tight
% title({'Normalized intensity of TMRM from non-expressing cells (not mito regions)';'Blue light: 0.01 s, cycle length = 0.1s'})
% xlabel('Time (s)')
% ylabel('Normalized intensity')
%
%
% legend('Cell 1', 'Cell 2','Cell 3', 'Cell 4','Cell 5', 'Cell 6','Cell 7', 'Cell 8','Cell 9', 'Cell 10','Cell 11', 'Cell 12','Cell 13', 'Cell 14','Cell 15', 'Cell 16','Cell 17', 'Cell 18','Cell 19','Location','Best')
% legend('boxoff')
% box off
% saveas(gca,[pathname 'analysis\Normalized intensity of ' name1 ' negative_not mito.fig']);
% saveas(gca,[pathname 'analysis\Normalized intensity of ' name1 ' negative_not mito.png']);
% close(gcf)
%% Normalized intensity with standard error shadow
intens_LEFT_norm = intens_LEFT_noBKG./intens_LEFT_noBKG(1,:);
intens_RIGHT_norm = intens_RIGHT_noBKG./intens_RIGHT_noBKG(1,:);
intens_LEFT_sem = std(intens_LEFT_norm')'/sqrt(size(intens_LEFT_norm, 2));
intens_RIGHT_sem = std(intens_RIGHT_norm')'/sqrt(size(intens_RIGHT_norm, 2));

figure()
set(gcf,'outerposition',get(0,'screensize'));
subplot(1,2,1)
seShoadow_LEFT = area(t_keyframe, [(mean(intens_LEFT_norm, 2) - intens_LEFT_sem), (2 * intens_LEFT_sem )]);
hold on;
set(seShoadow_LEFT(1),'Visible','off');
% set(seShoadow2(1),'Visible','off');
%set(seShoadow2(2),'FaceColor',[0.46, 0.68, 0.17]' ,'FaceAlpha',0.5, 'EdgeColor', 'none');
set(seShoadow_LEFT(2),'FaceColor',[0, 1, 1]' ,'FaceAlpha',0.3, 'EdgeColor', 'none');
plot(t_keyframe, mean(intens_LEFT_norm, 2), 'Color', getrgb(WL_LEFT),'LineWidth',2);hold on
box off
% axis ([0 max(t) 0.90 1.30])
axis tight
title({'Normalized intensity of ' name_LEFT })
xlabel('Time (s)')
ylabel('Normalized intensity')
% legend('Cell 1', 'Cell 2','Cell 3', 'Cell 4','Cell 5', 'Cell 6','Cell 7', 'Cell 8','Cell 9', 'Cell 10','Cell 11', 'Cell 12','Cell 13', 'Cell 14','Cell 15', 'Cell 16','Cell 17', 'Cell 18','Cell 19','Location','best')
% legend('boxoff')


subplot(1,2,2)
seShoadow_RIGHT = area(t_keyframe, [(mean(intens_RIGHT_norm, 2) - intens_LEFT_sem), (2 * intens_LEFT_sem )]);
hold on;
set(seShoadow_RIGHT(1),'Visible','off');
% set(seShoadow2(1),'Visible','off');
%set(seShoadow2(2),'FaceColor',[0.46, 0.68, 0.17]' ,'FaceAlpha',0.5, 'EdgeColor', 'none');
set(seShoadow_RIGHT(2),'FaceColor',[1, 0, 1]' ,'FaceAlpha',0.2, 'EdgeColor', 'none');
plot(t_keyframe, mean(intens_RIGHT_norm, 2), 'Color', getrgb(WL_RIGHT),'LineWidth',2);hold on
box off
% axis ([0 max(t) 0.90 1.30])
axis tight
title({'Normalized intensity of ' name_RIGHT })
xlabel('Time (s)')
ylabel('Normalized intensity')
% legend('Cell 1', 'Cell 2','Cell 3', 'Cell 4','Cell 5', 'Cell 6','Cell 7', 'Cell 8','Cell 9', 'Cell 10','Cell 11', 'Cell 12','Cell 13', 'Cell 14','Cell 15', 'Cell 16','Cell 17', 'Cell 18','Cell 19','Location','best')
% legend('boxoff')

% plot([61 119],[1.05 1.05],'LineWidth',3,'Color',getrgb(488))
saveas(gca,[pathname 'analysis\Normalized mean intensity of ' name1 '.fig']);
saveas(gca,[pathname 'analysis\Normalized mean intensity of ' name1 '.png']);
close(gcf)


%% Dynamic range with standard error shadow
intens_LEFT_dy = (intens_LEFT_norm ./ mean(intens_LEFT_norm(1:8, :),1) -1) * 100;
intens_RIGHT_dy = (intens_RIGHT_norm./mean(intens_RIGHT_norm(1:8, :),1) - 1) * 100;
intens_LEFT_dysem = std(intens_LEFT_dy')'/sqrt(size(intens_LEFT_dy, 2));
intens_RIGHT_dysem = std(intens_RIGHT_dy')'/sqrt(size(intens_RIGHT_dy, 2));

figure()
set(gcf,'outerposition',get(0,'screensize'));
subplot(1,2,1)
seShoadow_LEFT = area(t_keyframe, [(mean(intens_LEFT_dy, 2) - intens_LEFT_dysem), (2 * intens_LEFT_dysem)]);
hold on;
set(seShoadow_LEFT(1),'Visible','off');
% set(seShoadow2(1),'Visible','off');
%set(seShoadow2(2),'FaceColor',[0.46, 0.68, 0.17]' ,'FaceAlpha',0.5, 'EdgeColor', 'none');
set(seShoadow_LEFT(2),'FaceColor',[0, 1, 1]' ,'FaceAlpha',0.3, 'EdgeColor', 'none');
plot(t_keyframe, mean(intens_LEFT_dy, 2), 'Color', getrgb(WL_LEFT),'LineWidth',2);hold on
box off
% axis ([0 max(t) 0.90 1.30])
axis tight
title({'Dyanamic range of ' name_LEFT })
xlabel('Time (s)')
ylabel('¦¤F/F0 (%)')
% legend('Cell 1', 'Cell 2','Cell 3', 'Cell 4','Cell 5', 'Cell 6','Cell 7', 'Cell 8','Cell 9', 'Cell 10','Cell 11', 'Cell 12','Cell 13', 'Cell 14','Cell 15', 'Cell 16','Cell 17', 'Cell 18','Cell 19','Location','best')
% legend('boxoff')


subplot(1,2,2)
seShoadow_RIGHT = area(t_keyframe, [(mean(intens_RIGHT_dy, 2) - intens_RIGHT_dysem), (2 * intens_RIGHT_dysem )]);
hold on;
set(seShoadow_RIGHT(1),'Visible','off');
% set(seShoadow2(1),'Visible','off');
%set(seShoadow2(2),'FaceColor',[0.46, 0.68, 0.17]' ,'FaceAlpha',0.5, 'EdgeColor', 'none');
set(seShoadow_RIGHT(2),'FaceColor',[1, 0, 1]' ,'FaceAlpha',0.2, 'EdgeColor', 'none');
plot(t_keyframe, mean(intens_RIGHT_dy, 2), 'Color', getrgb(WL_RIGHT),'LineWidth',2);hold on
box off
% axis ([0 max(t) 0.90 1.30])
axis tight
title({'Dyanamic range of ' name_RIGHT })
xlabel('Time (s)')
ylabel('¦¤F/F0 (%)')
% legend('Cell 1', 'Cell 2','Cell 3', 'Cell 4','Cell 5', 'Cell 6','Cell 7', 'Cell 8','Cell 9', 'Cell 10','Cell 11', 'Cell 12','Cell 13', 'Cell 14','Cell 15', 'Cell 16','Cell 17', 'Cell 18','Cell 19','Location','best')
% legend('boxoff')

% plot([61 119],[1.05 1.05],'LineWidth',3,'Color',getrgb(488))
saveas(gca,[pathname 'analysis\Mean dynamic of ' name1 '.fig']);
saveas(gca,[pathname 'analysis\Mean dynamic of ' name1 '.png']);
close(gcf)

%% Merge dynamic range with standard error shadow
figure()
seShoadow_LEFT = area(t_keyframe, [(mean(intens_LEFT_dy, 2) - intens_LEFT_dysem), (2 * intens_LEFT_dysem)]);
hold on;
set(seShoadow_LEFT(1),'Visible','off');
% set(seShoadow2(1),'Visible','off');
%set(seShoadow2(2),'FaceColor',[0.46, 0.68, 0.17]' ,'FaceAlpha',0.5, 'EdgeColor', 'none');
set(seShoadow_LEFT(2),'FaceColor',[0, 1, 1]' ,'FaceAlpha',0.3, 'EdgeColor', 'none');
line1 = plot(t_keyframe, mean(intens_LEFT_dy, 2), 'Color', getrgb(WL_LEFT),'LineWidth',2);hold on
box off
% axis ([0 max(t) 0.90 1.30])
% legend('Cell 1', 'Cell 2','Cell 3', 'Cell 4','Cell 5', 'Cell 6','Cell 7', 'Cell 8','Cell 9', 'Cell 10','Cell 11', 'Cell 12','Cell 13', 'Cell 14','Cell 15', 'Cell 16','Cell 17', 'Cell 18','Cell 19','Location','best')
% legend('boxoff')

seShoadow_RIGHT = area(t_keyframe, [(mean(intens_RIGHT_dy, 2) - intens_RIGHT_dysem), (2 * intens_RIGHT_dysem )]);
hold on;
set(seShoadow_RIGHT(1),'Visible','off');
% set(seShoadow2(1),'Visible','off');
%set(seShoadow2(2),'FaceColor',[0.46, 0.68, 0.17]' ,'FaceAlpha',0.5, 'EdgeColor', 'none');
set(seShoadow_RIGHT(2),'FaceColor',[1, 0, 1]' ,'FaceAlpha',0.2, 'EdgeColor', 'none');
line2 = plot(t_keyframe, mean(intens_RIGHT_dy, 2), 'Color', getrgb(WL_RIGHT),'LineWidth',2);hold on
legend([line1, line2], {'ArcLightA242', 'HASAP1JF635'})
box off
% axis ([0 max(t) 0.90 1.30])
axis tight
title({'Dyanamic range of indicators'})
xlabel('Time (s)')
ylabel('¦¤F/F0 (%)')
% legend('Cell 1', 'Cell 2','Cell 3', 'Cell 4','Cell 5', 'Cell 6','Cell 7', 'Cell 8','Cell 9', 'Cell 10','Cell 11', 'Cell 12','Cell 13', 'Cell 14','Cell 15', 'Cell 16','Cell 17', 'Cell 18','Cell 19','Location','best')
% legend('boxoff')

% plot([61 119],[1.05 1.05],'LineWidth',3,'Color',getrgb(488))
saveas(gca,[pathname 'analysis\Merge dynamic of ' name1 '.fig']);
saveas(gca,[pathname 'analysis\Merge dynamic of ' name1 '.png']);
close(gcf)

%%
mkdir([pathname  'analysis\movie\tiff\' num2str(name1) '_lapse\'])
clear M j myObj
figure();
for j = 1:size(mov2,3)
    imshow(mov2(:,1:size(mov2,2)/2,j),[cmin_LEFT cmax_LEFT], 'Border', 'tight')
    colormap(pseudocolor(WL_LEFT));
    text(10,20,[num2str(t_keyframe(j)) ' s'], 'FontSize', 20, 'color', [0.99 0.99 0.99])
    %     text(62,20,[' min'], 'FontSize', 20, 'color', [0.99 0.99 0.99])
    text(10,nrow-20,[name_LEFT], 'FontSize', 18, 'color', getrgb(WL_LEFT))
    %     text(10,480,['180 uM GM1'], 'FontSize', 16, 'color', [0.99 0.99 0.99])
    M(j) = getframe(gca);
    %saveas(gca, [pathname '\Blink_' num2str((j)) '.tif']);  % Use this line to save each frame as a separate TIF file.
end
speedup = 1;
myObj = VideoWriter([pathname '\analysis\movie\'  '\' name_LEFT ' movie_' num2str(speedup) 'x.avi'],'Uncompressed AVI');%create an AVI file
myObj.FrameRate = speedup*(1000/dt_mov);
open(myObj);
writeVideo(myObj,M)
close(myObj);
myObj = VideoWriter([pathname '\analysis\movie\' '\' name_LEFT ' movie_compressed_H.264 encoding_' num2str(speedup) 'x.mp4'],'MPEG-4');%create an AVI file
myObj.FrameRate = speedup*(1000/dt_mov);
myObj.Quality = 100;
open(myObj);
writeVideo(myObj,M)
close(myObj);


for j = 1:size(M,2)
    im = frame2im(M(j));
    %     im = (im(41:860,646:1305,:));
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if j == 1
        imwrite(imind,cm,[pathname '\analysis\movie\Movie_5fps_LEFT.gif'],'gif', 'Loopcount',inf,'DelayTime',0.2);
        imwrite(imind,cm,[pathname '\analysis\movie\Movie_10fps_LEFT.gif'],'gif', 'Loopcount',inf,'DelayTime',0.1);
        
    else
        imwrite(imind,cm,[pathname '\analysis\movie\Movie_5fps_LEFT.gif'],'gif','WriteMode','append','DelayTime',0.2);
        imwrite(imind,cm,[pathname '\analysis\movie\Movie_10fps_LEFT.gif'],'gif','WriteMode','append','DelayTime',0.1);
        
    end
    
end


clear M j myObj
figure();
for j = 1:size(mov2,3)
    imshow(mov2(:,size(mov2,2)/2+1:end,j),[cmin_RIGHT cmax_RIGHT], 'Border', 'tight')
    colormap(pseudocolor(WL_RIGHT));
    text(10,20,[num2str(t_keyframe(j)) ' s'], 'FontSize', 20, 'color', [0.99 0.99 0.99])
    %     text(62,20,[' min'], 'FontSize', 20, 'color', [0.99 0.99 0.99])
    text(nrow/2-5,nrow-20,[name_RIGHT], 'FontSize', 18, 'color', getrgb(WL_RIGHT))
    %     text(10,480,['180 uM GM1'], 'FontSize', 16, 'color', [0.99 0.99 0.99])
    M(j) = getframe(gca);
    %saveas(gca, [pathname '\Blink_' num2str((j)) '.tif']);  % Use this line to save each frame as a separate TIF file.
end
speedup = 1;
myObj = VideoWriter([pathname '\analysis\movie\'  '\' name_RIGHT ' movie_' num2str(speedup) 'x.avi'],'Uncompressed AVI');%create an AVI file
myObj.FrameRate = speedup*(1000/dt_mov);
open(myObj);
writeVideo(myObj,M)
close(myObj);
myObj = VideoWriter([pathname '\analysis\movie\' '\' name_RIGHT ' movie_compressed_H.264 encoding_' num2str(speedup) 'x.mp4'],'MPEG-4');%create an AVI file
myObj.FrameRate = speedup*(1000/dt_mov);
myObj.Quality = 100;
open(myObj);
writeVideo(myObj,M)
close(myObj);


for j = 1:size(M,2)
    im = frame2im(M(j));
    %     im = (im(41:860,646:1305,:));
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if j == 1
        imwrite(imind,cm,[pathname '\analysis\movie\Movie_5fps_RIGHT.gif'],'gif', 'Loopcount',inf,'DelayTime',0.2);
        imwrite(imind,cm,[pathname '\analysis\movie\Movie_10fps_RIGHT.gif'],'gif', 'Loopcount',inf,'DelayTime',0.1);
        
    else
        imwrite(imind,cm,[pathname '\analysis\movie\Movie_5fps_RIGHT.gif'],'gif','WriteMode','append','DelayTime',0.2);
        imwrite(imind,cm,[pathname '\analysis\movie\Movie_10fps_RIGHT.gif'],'gif','WriteMode','append','DelayTime',0.1);
        
    end
    
end

% mkdir([pathname  '\movie\tiff\' num2str(name1) '_lapse\'])
% Tiff([pathname  '\movie\tiff\' num2str(name1) '.tif'],'w')
for num_tiff=1:size(mov2,3)
    imwrite(uint16(mov2(:,:,num_tiff)),[pathname  'analysis\movie\tiff\' num2str(name1) '.tif'],'WriteMode', 'append',  'Compression','none')
    imwrite(uint16(mov2(:,:,num_tiff)),[pathname  'analysis\movie\tiff\' num2str(name1) '_lapse\' num2str(name1)  '_' num2str(num_tiff,2) '.tif'],'WriteMode', 'append',  'Compression','none')
end
%% Export the data to excel: mito regions
intens_LEFT_std = std(intens_LEFT_noBKG')';
intens_LEFT_sem = intens_LEFT_std./sqrt(size(intens_LEFT_noBKG,2));

intens_RIGHT_std = std(intens_RIGHT_noBKG')';
intens_RIGHT_sem = intens_RIGHT_std./sqrt(size(intens_RIGHT_noBKG,2));

% intens2_sem = intens2_std./sqrt(size(intens2,2));

% xlswrite([pathname 'data_analysis.xlsx'],{'Cdensity', 'Current','Threshold Potential', 'Maximum Rising Vm','Spike Of Each Injection', 'Resting Potential','Spike Number','Max Rising Speed','Amplitude','Vpeak','FWHM','APD_threshold'},'AIS','A1');
for n = 1:size(intens_LEFT_noBKG,2)
    xlswrite([pathname 'data_analysis.xlsx'],{['Cell ' num2str(n)]},name_LEFT,[char(65+n) '1']);
end

% for n = 1:size(intens,2)
%     xlswrite([pathname 'data_analysis.xlsx'],{['Cell ' num2str(n)]},name2,[char(65+n) '1']);
% end
xlswrite([pathname 'data_analysis.xlsx'],{'Standard dev.'},name_LEFT,[char(66+n) '1']);
% xlswrite([pathname 'data_analysis.xlsx'],{'Standard dev.'},name2,[char(66+n) '1']);

xlswrite([pathname 'data_analysis.xlsx'],intens_LEFT_std,name_LEFT,[char(66+n) '2']);
% xlswrite([pathname 'data_analysis.xlsx'],intens2_std,name2,[char(66+n) '2']);

xlswrite([pathname 'data_analysis.xlsx'],{'S.E.M.'},name_LEFT,[char(67+n) '1']);
% xlswrite([pathname 'data_analysis.xlsx'],{'S.E.M.'},name2,[char(67+n) '1']);

xlswrite([pathname 'data_analysis.xlsx'],{['Average, n = ' num2str(size(intens_LEFT_noBKG,2))]},name_LEFT,[char(68+n) '1']);
% xlswrite([pathname 'data_analysis.xlsx'],{['Average, n = ' num2str(size(intens,2))]},name2,[char(68+n) '1']);


xlswrite([pathname 'data_analysis.xlsx'],intens_LEFT_sem, name_LEFT,[char(67+n) '2']);
% xlswrite([pathname 'data_analysis.xlsx'],intens2_sem,name2,[char(67+n) '2']);

xlswrite([pathname 'data_analysis.xlsx'],mean(intens_LEFT_noBKG,2),name_LEFT,[char(68+n) '2']);
% xlswrite([pathname 'data_analysis.xlsx'],mean(intens2_norm_wobkg,2),name2,[char(68+n) '2']);


xlswrite([pathname 'data_analysis.xlsx'],{'Time (s)'},name_LEFT,'A1');
% xlswrite([pathname 'data_analysis.xlsx'],{'Time (s)'},name2,'A1');

xlswrite([pathname 'data_analysis.xlsx'],t_keyframe',name_LEFT,'A2');
% xlswrite([pathname 'data_analysis.xlsx'],t_mov',name2,'A2');

xlswrite([pathname 'data_analysis.xlsx'],intens_LEFT_noBKG,name_LEFT,'B2');
% xlswrite([pathname 'data_analysis.xlsx'],intens2_norm_wobkg,name2,'B2');
for n = 1:size(intens_RIGHT_noBKG,2)
    xlswrite([pathname 'data_analysis.xlsx'],{['Cell ' num2str(n)]},name_RIGHT,[char(65+n) '1']);
end

% for n = 1:size(intens,2)
%     xlswrite([pathname 'data_analysis.xlsx'],{['Cell ' num2str(n)]},name2,[char(65+n) '1']);
% end
xlswrite([pathname 'data_analysis.xlsx'],{'Standard dev.'},name_RIGHT,[char(66+n) '1']);
% xlswrite([pathname 'data_analysis.xlsx'],{'Standard dev.'},name2,[char(66+n) '1']);

xlswrite([pathname 'data_analysis.xlsx'],intens_RIGHT_std,name_RIGHT,[char(66+n) '2']);
% xlswrite([pathname 'data_analysis.xlsx'],intens2_std,name2,[char(66+n) '2']);

xlswrite([pathname 'data_analysis.xlsx'],{'S.E.M.'},name_RIGHT,[char(67+n) '1']);
% xlswrite([pathname 'data_analysis.xlsx'],{'S.E.M.'},name2,[char(67+n) '1']);

xlswrite([pathname 'data_analysis.xlsx'],{['Average, n = ' num2str(size(intens_RIGHT_noBKG,2))]},name_RIGHT,[char(68+n) '1']);
% xlswrite([pathname 'data_analysis.xlsx'],{['Average, n = ' num2str(size(intens,2))]},name2,[char(68+n) '1']);


xlswrite([pathname 'data_analysis.xlsx'],intens_RIGHT_sem,name_RIGHT,[char(67+n) '2']);
% xlswrite([pathname 'data_analysis.xlsx'],intens2_sem,name2,[char(67+n) '2']);

xlswrite([pathname 'data_analysis.xlsx'],mean(intens_RIGHT_noBKG,2),name_RIGHT,[char(68+n) '2']);
% xlswrite([pathname 'data_analysis.xlsx'],mean(intens2_norm_wobkg,2),name2,[char(68+n) '2']);


xlswrite([pathname 'data_analysis.xlsx'],{'Time (s)'},name_RIGHT,'A1');
% xlswrite([pathname 'data_analysis.xlsx'],{'Time (s)'},name2,'A1');

xlswrite([pathname 'data_analysis.xlsx'],t_keyframe',name_RIGHT,'A2');
% xlswrite([pathname 'data_analysis.xlsx'],t_mov',name2,'A2');

xlswrite([pathname 'data_analysis.xlsx'],intens_RIGHT_noBKG,name_RIGHT,'B2');
% xlswrite([pathname 'data_analysis.xlsx'],intens2_norm_wobkg,name2,'B2')
%%
clear mov
% intens_LEFT_dy(7) is the fluorescent dynamic range at the 140th second
% intens_LEFT_dy(22) is the fluorescent dynamic range at the 300th second
window_LEFT_dy = intens_LEFT_dy(7:15,:);
window_RIGHT_dy = intens_RIGHT_dy(7:15, :); 
peak_LEFT_dy = max(abs(window_LEFT_dy));
peak_RIGHT_dy = max(abs(window_RIGHT_dy));
peak_ratio = (peak_RIGHT_dy./peak_LEFT_dy)';
save([pathname 'All variants.mat']);
fclose('all');


%%
mean_ArcLightA242 = mean(ArcLightA242, 2);
mean_HASAP1JF635 = mean(HASAP1JF635, 2);
ArcLightA242_sem = std(ArcLightA242')'/sqrt(size(ArcLightA242, 2));
HASAP1JF635_sem = std(HASAP1JF635')'/sqrt(size(HASAP1JF635, 2));
traceRatio = abs(HASAP1JF635./ArcLightA242);
meanRatio = mean(traceRatio, 2);
semRatio = std(traceRatio')'/sqrt(size(traceRatio, 2));

figure()
seShoadow_ArcLightA242 = area(t_keyframe, [mean_ArcLightA242 - ArcLightA242_sem, (2 * ArcLightA242_sem)]);
hold on;
set(seShoadow_ArcLightA242(1),'Visible','off');
% set(seShoadow2(1),'Visible','off');
%set(seShoadow2(2),'FaceColor',[0.46, 0.68, 0.17]' ,'FaceAlpha',0.5, 'EdgeColor', 'none');
set(seShoadow_ArcLightA242(2),'FaceColor',[0, 1, 1]' ,'FaceAlpha',0.3, 'EdgeColor', 'none');
line1 = plot(t_keyframe, mean_ArcLightA242, 'Color', getrgb(WL_LEFT),'LineWidth',2);hold on


seShoadow_HASAP1JF635 = area(t_keyframe, [mean_HASAP1JF635 - HASAP1JF635_sem, (2 * HASAP1JF635_sem)]);
hold on;
set(seShoadow_HASAP1JF635(1),'Visible','off');
% set(seShoadow2(1),'Visible','off');
%set(seShoadow2(2),'FaceColor',[0.46, 0.68, 0.17]' ,'FaceAlpha',0.5, 'EdgeColor', 'none');
set(seShoadow_HASAP1JF635(2),'FaceColor',[1, 0, 1]' ,'FaceAlpha',0.2, 'EdgeColor', 'none');
line2 = plot(t_keyframe, mean_HASAP1JF635, 'Color', getrgb(WL_RIGHT),'LineWidth',2);hold on
legend([line1, line2], {'ArcLightA242', 'HASAP1JF635'})

box off
% axis ([0 max(t) 0.90 1.30])
axis tight
title({'Dynamic range of voltage indicators'})
xlabel('Time (s)')
ylabel('¦¤F/F0 (%)')

saveas(gca,['X:\91 Data and analysis\YJunqi\Screening method\Stimulus pre-experiments\Statistics\20210822_35mMsummary' '.fig']);
saveas(gca,['X:\91 Data and analysis\YJunqi\Screening method\Stimulus pre-experiments\Statistics\20210822_35mMsummary' '.png']);


figure()
seShoadow_ratio = area(t_keyframe, [meanRatio - semRatio, (2 * semRatio)]);
hold on;
set(seShoadow_ratio(1),'Visible','off');
% set(seShoadow2(1),'Visible','off');
%set(seShoadow2(2),'FaceColor',[0.46, 0.68, 0.17]' ,'FaceAlpha',0.5, 'EdgeColor', 'none');
set(seShoadow_ratio(2),'FaceColor',[0.67, 0, 1]' ,'FaceAlpha',0.3, 'EdgeColor', 'none');
line1 = plot(t_keyframe, meanRatio,'Color', [0, 0, 1], 'LineWidth',2);hold on

box off
% axis ([0 max(t) 0.90 1.30])
axis tight
title({'Trace ratio of voltage indicators'})
xlabel('Time (s)')
ylabel('Absolute ratio')

saveas(gca,['X:\91 Data and analysis\YJunqi\Screening method\Stimulus pre-experiments\Statistics\20210822_35mMsummary_traceRatio' '.fig']);
saveas(gca,['X:\91 Data and analysis\YJunqi\Screening method\Stimulus pre-experiments\Statistics\20210822_35mMsummary_traceRatio' '.png']);
