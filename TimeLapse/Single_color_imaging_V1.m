clear all
%% loading two channels
    
name1 = 'mNeonGreen-CaaX';
% name2 = 'HASAP1JF635';
wvlt1 = 520;
% wvlt2 = 655;
name1color = getrgb(wvlt1);
% name2color = getrgb(wvlt2);
bin = 1;
basepath = 'H:\ZouLab\YJunqi\Screening method\System test\20220307 Gramicidin test\Dish4\test';
subfolder = '\';

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
dt_mov = 30000;         % exposure time in millisecond
mov = double(mov);
mov = mov(:,:,start_frame:end);
img = mean(mov, 3);
len = size(mov,3);
% t_mov = [0:(len-1)]*dt_mov/1000;     % time axis in second

% second channel
% FileTif=[pathname name2 '.tif'];
% InfoImage=imfinfo(FileTif);
% mImage=InfoImage(1).Width;
% nImage=InfoImage(1).Height;
% NumberImages=length(InfoImage);
% mov2=zeros(nImage,mImage,NumberImages,'uint16');

% TifLink = Tiff(FileTif, 'r');
% for i=1:NumberImages
%     TifLink.setDirectory(i);
%     mov2(:,:,i)=TifLink.read();
% end
% TifLink.close();
% ncol = 1024./bin;         % x
% nrow = 1024./bin;          % y
ncol = size(mov,2);         % x
nrow = size(mov,1);          % y
bkg = 100*power(bin,2);          % background due to camera bias (100 for bin 1x1)
dt_mov = 30000;         % exposure time in millisecond
% mov2 = double(mov2);
% mov2 = mov2(:,:,start_frame:end);
% img2 = mean(mov2, 3);
% len2 = size(mov2,3);
t_mov = [0:(len-1)]*dt_mov/1000;     % time axis in second
t_mov(7:end) = t_mov(7:end)+1200;
% mov = mov(:,:,3:end);
% mov2 = mov2(:,:,3:end);
% mkdir([pathname 'analysis']);

%% save movies
expGroup = 'analysis\mNeonGreen-CaaX\'
mkdir([pathname expGroup]);
cmin1 = 100;
cmax1 = 500;
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
myObj = VideoWriter([pathname expGroup name1 '.avi'],'Uncompressed AVI');%create an AVI file
myObj.FrameRate = 5;
open(myObj);
writeVideo(myObj,M)
close(myObj);

myObj = VideoWriter([pathname expGroup name1  '_compressed_H.264 encoding.mp4'],'MPEG-4');%create an AVI file
myObj.FrameRate = 5;
open(myObj);
writeVideo(myObj,M)
close(myObj);

% cmin2 = 100;
% cmax2 = 300;
% clear M
% figure(20)
% for j = 1:(size(mov2,3))
%     imshow(mov2(:,:,j),[cmin2,cmax2], 'Border', 'tight');hold on;
%     colormap(pseudocolor(wvlt2));
%     text(10,30,[num2str(t_mov(j)) ' s'], 'FontSize', 18, 'color', [0.99 0.99 0.99])
%     %text(550,20,[name2], 'FontSize', 20, 'color', getrgb(598));
%     text(10,nrow-20,[name2], 'FontSize', 18, 'color', name2color);
%     M(j) = getframe(gca);
%     %     saveas(gca, ['SNAPT_movie_' num2str((j)) '.tif']);  % Use this line to save each frame as a separate TIF file.
% end
% close(gcf)
% myObj = VideoWriter([pathname expGroup name2 '.avi'],'Uncompressed AVI');%create an AVI file
% myObj.FrameRate = 5;
% open(myObj);
% writeVideo(myObj,M)
% close(myObj);
% myObj = VideoWriter([pathname expGroup name2  '_compressed_H.264 encoding.mp4'],'MPEG-4');%create an AVI file
% myObj.FrameRate = 5;
% open(myObj);
% writeVideo(myObj,M)
% close(myObj);
%% background checking
colordef white
[~, intens_forbkg] = clicky_v2(mov, dt_mov, bin, mean(mov(:,:,1:end),3), bkg ,200 ,[name1, "backgroud"]);
saveas(gca,[pathname expGroup 'selected regions of bkg.fig']);
saveas(gca,[pathname expGroup 'selected regions of bkg.png']);

% figure(); 
% plot(t_mov,intens_forbkg);
% box off
% axis tight
% xlabel('time(s)')
% title('Selected background')
% saveas(gca,[pathname 'analysis' '\selected bkg.fig']);
% saveas(gca,[pathname 'analysis' '\selected bkg.png']);



% intens_forbkg = intens_forbkg-mean(intens_forbkg(1:8));
% intens_forbkg = intens_forbkg';

figure();
plot(t_mov,intens_forbkg)
axis tight
box off
xlabel('Time (s)')
ylabel('normalized intensity for bkg correction')
saveas(gca,[pathname expGroup 'selected bkg_' num2str(name1) '.fig']);
saveas(gca,[pathname expGroup 'selected bkg_' num2str(name1) '.png']);


% figure();
% plot(t_mov,intens_bkg);
% box off
% axis tight
% xlabel('time(s)')
% title('Selected background')
% saveas(gca,[pathname 'analysis' '\selected bkg' num2str(name1) '.fig']);
% saveas(gca,[pathname 'analysis' '\selected bkg' num2str(name1) '.png']);

% intens2_bkg = intens2_bkg-mean(intens2_bkg(1:8));
% intens2_bkg = intens2_bkg';
% figure();
% plot(t_mov,intens2_bkg)
% axis tight
% box off
% xlabel('Time (s)')
% ylabel('normalized intensity for bkg correction')
% saveas(gca,[pathname expGroup 'selected bkg_' num2str(name2) '.fig']);
% saveas(gca,[pathname expGroup 'selected bkg_' num2str(name2) '.png']);
%% analysis
[roi_points, intens] = clicky_v2(mov, dt_mov, bin, mean(mov(:,:,1:end),3), bkg, cmax1, name1);
saveas(gca,[pathname expGroup 'selected regions of ' name1 '.fig']);
saveas(gca,[pathname expGroup 'selected regions of ' name1 '.png']);

figure();
set(gcf,'outerposition',get(0,'screensize'));
imshow(img,[cmin1,cmax1]);hold on;
colormap(pseudocolor(wvlt1));
for region_num=1:max(size(roi_points))
    plot([roi_points{1,region_num}(:,1)],[roi_points{1,region_num}(:,2)],'LineWidth',3,'Color', getrgb(400+round((region_num-1)/size(intens,2)*270)));hold on
    text(mean(roi_points{1,region_num}(:,1)),mean(roi_points{1,region_num}(:,2)),num2str(region_num),'Color',[0.99 0.99 0.99],'FontSize',20);hold on;
    plot([8 69.53],[nrow-20 nrow-20],'LineWidth',5,'Color',[0.99 0.99 0.99]);
end
saveas(gca,[pathname expGroup 'selected regions of ' name1 '.fig']);
saveas(gca,[pathname expGroup 'selected regions of ' name1 '.png']);
close(gcf)


% figure();
% set(gcf,'outerposition',get(0,'screensize'));
% imshow(img2,[cmin2,cmax2]);hold on;
% colormap(pseudocolor(wvlt2));
% for region_num=1:max(size(roi_points))
%     plot([roi_points{1,region_num}(:,1)],[roi_points{1,region_num}(:,2)],'LineWidth',3,'Color', getrgb(400+round((region_num-1)/size(intens,2)*270)));hold on
%     text(mean(roi_points{1,region_num}(:,1)),mean(roi_points{1,region_num}(:,2)),num2str(region_num),'Color',[0.99 0.99 0.99],'FontSize',20);hold on;
%     plot([8 69.53],[nrow-20 nrow-20],'LineWidth',5,'Color',[0.99 0.99 0.99]);
% end
% saveas(gca,[pathname expGroup 'selected regions of ' name2 '.fig']);
% saveas(gca,[pathname expGroup 'selected regions of ' name2 '.png']);
% close(gcf)
for n = 1:size(intens,2)
    figure();
    plot(t_mov,(intens(:,n))./(intens(1,n)),'Color',name1color,'LineWidth',1.5);hold on
%     plot(t_mov,(intens2(:,n))./(intens2(1,n)),'Color',name2color,'LineWidth',1.5);
    title(['Region ' num2str(n)]);
    legend(name1, 'Location','northoutside','Orientation','horizontal');
    xlabel('Time(s)');
    ylabel('Relative intensity')
    box off
    axis tight
    saveas(gca,[pathname expGroup 'plot region ' num2str(n) '.fig']);
    saveas(gca,[pathname expGroup 'plot region ' num2str(n) '.png']);
    close(gcf)
%     figure();
%     [hAx,hLine1,hLine2] = plot(t_mov,(intens(:,n))./(intens(1,n)),t_mov));
%     ylabel(hAx(1),name1) % left y-axis
%     ylabel(hAx(2),name2) % right y-axis
%     title(['Region ' num2str(n)]);
%     xlabel('Time(s)')
%     box off
%     xlim(hAx(1),[0 max(t_mov)]);
%     xlim(hAx(2),[0 max(t_mov)]);
%     axis tight
%     hLine1.Color = name1color;
%     hLine2.Color = name2color;
%     set(hAx(:),'Ycolor','k') %设定两个Y轴的颜色为黑色
%     legend(name1, 'Location','northoutside','Orientation','horizontal')
%     saveas(gca,[pathname expGroup 'plot region ' num2str(n) '.fig']);
%     saveas(gca,[pathname expGroup '2 y Dual-color plot region ' num2str(n) '.png']);
%     close(gcf)
end
%% Channel 1 "tight" plot
% ColorMatrix = get(gca,'colororder');
intens_norm = zeros(size(intens));
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
legend('Cell 1', 'Cell 2','Cell 3', 'Cell 4','Cell 5', 'Cell 6','Cell 7', 'Cell 8','Cell 9', 'Cell 10','Cell 11', 'Cell 12','Cell 13', 'Cell 14','Cell 15', 'Cell 16','Cell 17', 'Cell 18','Cell 19','Cell 20','Cell 21','Cell 22','Cell 23','Cell 24', 'Cell 25','Cell 26','Cell 27', 'Cell 28','Cell 29','Cell 30', 'Location','Best')
legend('boxoff')
ylim([0.2 max(intens_norm(:))+0.1])
box off
% axis tight
saveas(gca,[pathname expGroup 'norm_intensity of '  num2str(name1) '.fig']);
saveas(gca,[pathname expGroup 'norm_intensity of '  num2str(name1) '.png']);
figure()
for n = 1:size(intens,2)
    colororder = mod(n,7);
    if colororder == 0
        colororder = 7;
    end
    plot(t_mov,intens_norm(:,n),'Color', getrgb(400+round((n-1)/size(intens,2)*270)));hold on
end
legend('Cell 1', 'Cell 2','Cell 3', 'Cell 4','Cell 5', 'Cell 6','Cell 7', 'Cell 8','Cell 9', 'Cell 10','Cell 11', 'Cell 12','Cell 13', 'Cell 14','Cell 15','Cell 16','Cell 17', 'Cell 18','Cell 19','Cell 20','Cell 21','Cell 22','Cell 23','Cell 24', 'Cell 25','Cell 26','Cell 27', 'Cell 28','Cell 29','Cell 30', 'Location','Best')
legend('boxoff')
title(['Normalized intensity of ' num2str(name1)]);
xlabel('time (s)')
ylim([min(intens_norm(:))-0.05 max(intens_norm(:))+0.1])
box off
% axis tight
saveas(gca,[pathname expGroup 'norm_intensity of '  num2str(name1) ' tight.fig']);
saveas(gca,[pathname expGroup 'norm_intensity of '  num2str(name1) ' tight.png']);

figure();
plot(t_mov,mean(intens_norm,2),'Color',getrgb(wvlt1),'LineWidth',2);
title(['Mean normalized intensity of ' num2str(name1)]);
xlabel('time (s)')
ylim([min(mean(intens_norm,2))-0.05 max(mean(intens_norm,2))+0.1])
box off
saveas(gca,[pathname expGroup 'mean_norm_intensity of '  num2str(name1) ' tight.fig']);
saveas(gca,[pathname expGroup 'mean_norm_intensity of '  num2str(name1) ' tight.png']);

%% intens (without bkg)
mkdir([pathname expGroup 'without background\'])
% ColorMatrix = get(gca,'colororder');
intens_wobkg = zeros(size(intens));
intens_norm_wobkg = zeros(size(intens));
for n = 1:size(intens,2)
    colororder = mod(n,7);
    if colororder == 0
        colororder = 7;
    end
    intens_wobkg(:,n) = (intens(:,n)-intens_forbkg(n));
    
end
% pbleach1 = [imerode(intens_wobkg(1:6,:), ones(3,size(intens, 2))); imerode(intens_wobkg(7:end,:), ones(6,size(intens, 2)))];
% pbleach1 = [smoothdata(intens_wobkg(1:6,:)./mean(intens_wobkg(1:3,:)), 1, 'movmedian',3); smoothdata(intens_wobkg(7:end,:)./mean(intens_wobkg(7:9,:)), 1, 'movmedian',3)];
% pbleach1 = [intens_wobkg(1:6,:)./intens_wobkg(1,:); intens_wobkg(7:end,:)./intens_wobkg(7,:)];
% intens_wobkg = intens_wobkg./pbleach1;
intens_norm_wobkg = intens_wobkg./intens_wobkg(1,:);

figure()
for n = 1:size(intens,2)
    colororder = mod(n,7);
    if colororder == 0
        colororder = 7;
    end
    plot(t_mov,intens_norm_wobkg(:,n),'Color', getrgb(400+round((n-1)/size(intens,2)*270)));hold on
end
legend('Cell 1', 'Cell 2','Cell 3', 'Cell 4','Cell 5', 'Cell 6','Cell 7', 'Cell 8','Cell 9', 'Cell 10','Cell 11', 'Cell 12','Cell 13', 'Cell 14','Cell 15','Cell 16','Cell 17','Cell 18','Cell 19','Cell 20','Cell 20','Cell 21','Cell 22','Cell 23','Cell 24', 'Cell 25','Cell 26','Cell 27', 'Cell 28','Cell 29','Cell 30', 'Location','Best')
legend('boxoff')
title(['Normalized intensity of ' num2str(name1) ' (without background)']);
xlabel('time (s)')
ylim([0.2 max(intens_norm_wobkg(:))+0.1])
box off
% axis tight
saveas(gca,[pathname expGroup 'without background\norm_intensity of '  num2str(name1) '_w_o bkg.fig']);
saveas(gca,[pathname expGroup 'without background\norm_intensity of '  num2str(name1) '_w_o bkg.png']);

figure()
for n = 1:size(intens,2)
    colororder = mod(n,7);
    if colororder == 0
        colororder = 7;
    end
    plot(t_mov,intens_norm_wobkg(:,n),'Color', getrgb(400+round((n-1)/size(intens,2)*270)));hold on
end
legend('Cell 1', 'Cell 2','Cell 3', 'Cell 4','Cell 5', 'Cell 6','Cell 7', 'Cell 8','Cell 9', 'Cell 10','Cell 11', 'Cell 12','Cell 13', 'Cell 14','Cell 15','Cell 16','Cell 17','Cell 18','Cell 19','Cell 20','Cell 20','Cell 21','Cell 22','Cell 23','Cell 24', 'Cell 25','Cell 26','Cell 27', 'Cell 28','Cell 29','Cell 30', 'Location','Best')
legend('boxoff')
title(['Normalized intensity of ' num2str(name1) ' (without background)']);
xlabel('time (s)')
ylim([min(intens_norm_wobkg(:))-0.05 max(intens_norm_wobkg(:))+0.1])
box off
% axis tight
saveas(gca,[pathname expGroup 'without background\norm_intensity of '  num2str(name1) ' tight_w_o bkg.fig']);
saveas(gca,[pathname expGroup 'without background\norm_intensity of '  num2str(name1) ' tight_w_o bkg.png']);

figure();
plot(t_mov,mean(intens_norm_wobkg,2),'Color',getrgb(wvlt1),'LineWidth',2);
title(['Mean normalized intensity of ' num2str(name1) ' (without background)']);
xlabel('time (s)')
ylim([min(mean(intens_norm_wobkg,2))-0.05 max(mean(intens_norm_wobkg,2))+0.1])
box off
saveas(gca,[pathname expGroup 'without background\mean_norm_intensity of '  num2str(name1) ' tight_w_o bkg.fig']);
saveas(gca,[pathname expGroup 'without background\mean_norm_intensity of '  num2str(name1) ' tight_w_o bkg.png']);

intens_std = std(intens_norm_wobkg')';
intens_sem = intens_std./sqrt(size(intens,2));
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
saveas(gca,[pathname expGroup '1st order of difference_' num2str(name1) ' tight.fig']);
saveas(gca,[pathname expGroup '1st order of difference_' num2str(name1) ' tight.png']);
%% Export to EXCEL
% xlswrite([pathname 'data_analysis.xlsx'],{'Cdensity', 'Current','Threshold Potential', 'Maximum Rising Vm','Spike Of Each Injection', 'Resting Potential','Spike Number','Max Rising Speed','Amplitude','Vpeak','FWHM','APD_threshold'},'AIS','A1');
for n = 1:size(intens,2)
    xlswrite([pathname expGroup 'data_analysis.xlsx'],{['Cell ' num2str(n)]},name1,[char(65+n) '1']);
end


xlswrite([pathname expGroup 'data_analysis.xlsx'],{'Standard dev.'},name1,[char(66+n) '1']);

xlswrite([pathname expGroup 'data_analysis.xlsx'],intens_std,name1,[char(66+n) '2']);

xlswrite([pathname expGroup 'data_analysis.xlsx'],{'S.E.M.'},name1,[char(67+n) '1']);

xlswrite([pathname expGroup 'data_analysis.xlsx'],{['Average, n = ' num2str(size(intens,2))]},name1,[char(68+n) '1']);

xlswrite([pathname expGroup 'data_analysis.xlsx'],intens_sem,name1,[char(67+n) '2']);

xlswrite([pathname expGroup 'data_analysis.xlsx'],mean(intens_norm_wobkg,2),name1,[char(68+n) '2']);

xlswrite([pathname expGroup 'data_analysis.xlsx'],{'Time (s)'},name1,'A1');

xlswrite([pathname expGroup 'data_analysis.xlsx'],t_mov',name1,'A2');

xlswrite([pathname expGroup 'data_analysis.xlsx'],intens_norm_wobkg,name1,'B2');

%% dual plot of averge dynamic normalized intensity with shadow
figure()
seShoadow1_1 = area(t_mov(1:6), [(mean(intens_norm_wobkg(1:6,:), 2) - intens_sem(1:6)), (2 * intens_sem(1:6))]);
set(seShoadow1_1(1),'Visible','off');
%set(seShoadow1(2),'FaceColor',[0.85, 0.32, 0]' ,'FaceAlpha',0.4, 'EdgeColor', 'none');
set(seShoadow1_1(2),'FaceColor',[0.46, 0.68, 0.17]' ,'FaceAlpha',0.5, 'EdgeColor', 'none');
hold on;
seShoadow1_2 = area(t_mov(7:end), [(mean(intens_norm_wobkg(7:end,:), 2) - intens_sem(7:end)), (2 * intens_sem(7:end))]);
set(seShoadow1_2(1),'Visible','off');
%set(seShoadow1(2),'FaceColor',[0.85, 0.32, -2]' ,'FaceAlpha',0.4, 'EdgeColor', 'none');
set(seShoadow1_2(2),'FaceColor',[0.46, 0.68, 0.17]' ,'FaceAlpha',0.5, 'EdgeColor', 'none');
hold on;
plot1_1 = plot(t_mov(1:6), mean(intens_norm_wobkg(1:6,:),2), 'Color', name1color,'LineWidth',1.5);
plot1_2 = plot(t_mov(7:end), mean(intens_norm_wobkg(7:end,:),2), 'Color', name1color,'LineWidth',1.5);
legend(plot1_1, name1, 'Location','northoutside','Orientation','horizontal');
box off;
axis tight;
xlabel('Time(s)');
ylabel('Relative intensity')
xlabel('Time (s)');
% ylim([0.8 ,1.15]);
ylabel('Normalized Intensity');
saveas(gca,[pathname expGroup 'shadow_mean_norm_intensity.fig']);
saveas(gca,[pathname expGroup 'shadow_mean_norm_intensity.png']);

delta_intens = (intens_norm_wobkg - mean(intens_norm_wobkg(1:5, :)))./mean(intens_norm_wobkg(1:5,:));
delta_sem = std(delta_intens')'./sqrt(size(delta_intens, 2));

figure()
seShoadow2_1 = area(t_mov(1:6), [(mean(delta_intens(1:6,:), 2) - delta_sem(1:6)), (2 * delta_sem(1:6))]);
set(seShoadow2_1(1),'Visible','off');
%set(seShoadow1_2(2),'FaceColor',[0.85, 0.32, 0]' ,'FaceAlpha',0.4, 'EdgeColor', 'none');
set(seShoadow2_1(2),'FaceColor',[0, 1, 1]' ,'FaceAlpha',0.3, 'EdgeColor', 'none');
hold on;
seShoadow2_2 = area(t_mov(7:end), [(mean(delta_intens(7:end,:), 2) - delta_sem(7:end)), (2 * delta_sem(7:end))]);
set(seShoadow2_2(1),'Visible','off');
%set(seShoadow2_2(2),'FaceColor',[0.46, 0.68, 0.17]' ,'FaceAlpha',0.4, 'EdgeColor', 'none');
set(seShoadow2_2(2), 'FaceColor', [0, 1, 1]' , 'FaceAlpha',0.3, 'EdgeColor', 'none');
hold on;
plot2_1 = plot(t_mov(1:6), mean(delta_intens(1:6,:),2), 'Color', name1color,'LineWidth',1.5);
plot2_2 = plot(t_mov(7:end), mean(delta_intens(7:end,:),2), 'Color', name1color,'LineWidth',1.5);
legend(plot2_1, name1, 'Location','northoutside','Orientation','horizontal');
box off;
axis tight;
xlabel('Time(s)');
ylabel('Relative intensity')
xlabel('Time (s)');
% ylim([-0.2,0.15]);
ylabel('ΔF/F0');
saveas(gca,[pathname expGroup 'shadow_deltaF_F0.fig']);
saveas(gca,[pathname expGroup 'shadow_deltaF_F0.png']);

%%
save([pathname 'All variants.mat']);
%% Reatio calculation
dyArcLightA242 = ArcLightA242 ./ mean(ArcLightA242(1:6, :)) -1;
dyHASAP1JF635 = HASAP1JF635 ./ mean(HASAP1JF635(1:6, :)) -1;
mean_dyArcLightA242 = mean(dyArcLightA242, 2);
mean_dyHASAP1JF635 = mean(dyHASAP1JF635, 2);
dyArcLightA242_sem = std(dyArcLightA242')'/sqrt(size(dyArcLightA242, 2));
dyHASAP1JF635_sem = std(dyHASAP1JF635')'/sqrt(size(dyHASAP1JF635, 2));
traceRatio = abs(dyHASAP1JF635./dyArcLightA242);
meanRatio = mean(traceRatio, 2);
semRatio = std(traceRatio')'/sqrt(size(traceRatio, 2));

figure()
seShoadowArcLightA242_1 = area(t_mov(1:6), [(mean(dyArcLightA242(1:6,:), 2) - dyArcLightA242_sem(1:6)), (2 * dyArcLightA242_sem(1:6))]);
set(seShoadowArcLightA242_1(1),'Visible','off');
%set(seShoadow1_2(2),'FaceColor',[0.85, 0.32, 0]' ,'FaceAlpha',0.4, 'EdgeColor', 'none');
set(seShoadowArcLightA242_1(2),'FaceColor',[0, 1, 1]' ,'FaceAlpha',0.3, 'EdgeColor', 'none');
hold on;
seShoadowArcLightA242_2 = area(t_mov(7:end), [(mean(dyArcLightA242(7:end,:), 2) - dyArcLightA242_sem(7:end)), (2 * dyArcLightA242_sem(7:end))]);
set(seShoadowArcLightA242_2(1),'Visible','off');
%set(seShoadow2_2(2),'FaceColor',[0.46, 0.68, 0.17]' ,'FaceAlpha',0.4, 'EdgeColor', 'none');
set(seShoadowArcLightA242_2(2), 'FaceColor', [0, 1, 1]' , 'FaceAlpha',0.3, 'EdgeColor', 'none');
hold on;
seShoadowHASAP1JF635_1 = area(t_mov(1:6), [(mean(dyHASAP1JF635(1:6,:), 2) - dyHASAP1JF635_sem(1:6)), (2 * dyHASAP1JF635_sem(1:6))]);
set(seShoadowHASAP1JF635_1(1),'Visible','off');
%set(seShoadow1_2(2),'FaceColor',[0.85, 0.32, 0]' ,'FaceAlpha',0.4, 'EdgeColor', 'none');
set(seShoadowHASAP1JF635_1(2),'FaceColor',[1, 0, 1]' ,'FaceAlpha',0.2, 'EdgeColor', 'none');
hold on;
seShoadowHASAP1JF635_2 = area(t_mov(7:end), [(mean(dyHASAP1JF635(7:end,:), 2) - dyHASAP1JF635_sem(7:end)), (2 * dyHASAP1JF635_sem(7:end))]);
set(seShoadowHASAP1JF635_2(1),'Visible','off');
%set(seShoadow2_2(2),'FaceColor',[0.46, 0.68, 0.17]' ,'FaceAlpha',0.4, 'EdgeColor', 'none');
set(seShoadowHASAP1JF635_2(2), 'FaceColor', [1, 0, 1]' , 'FaceAlpha',0.2, 'EdgeColor', 'none');
plot1_1 = plot(t_mov(1:6), mean(dyArcLightA242(1:6,:),2), 'Color', name1color,'LineWidth',1.5);
plot1_2 = plot(t_mov(7:end), mean(dyArcLightA242(7:end,:),2), 'Color', name1color,'LineWidth',1.5);
plot2_1 = plot(t_mov(1:6), mean(dyHASAP1JF635(1:6,:),2), 'Color', name2color,'LineWidth',1.5);
plot2_2 = plot(t_mov(7:end), mean(dyHASAP1JF635(7:end,:),2), 'Color', name2color,'LineWidth',1.5);
legend([plot1_1, plot2_1], {'ArcLightA242', 'HASAP1JF635'}, 'Location','northoutside','Orientation','horizontal');

box off
% axis ([0 max(t) 0.90 1.30])
axis tight
title({'Dynamic range of voltage indicators'})
xlabel('Time (s)')
ylabel('ΔF/F0');
set(gcf,'color','w');
breakxaxis([150 1380]);
saveas(gca,['X:\91 Data and analysis\YJunqi\Screening method\System test\Statistics\20220404_gramicidinSummary' '.fig']);
saveas(gca,['X:\91 Data and analysis\YJunqi\Screening method\System test\Statistics\20220404_gramicidinSummary' '.png']);


figure()
seShoadow_ratio_1 = area(t_mov(1:6), [(mean(traceRatio(1:6,:), 2) - semRatio(1:6)), (2 * semRatio(1:6))]);
set(seShoadow_ratio_1(1),'Visible','off');
%set(seShoadow1_2(2),'FaceColor',[0.85, 0.32, 0]' ,'FaceAlpha',0.4, 'EdgeColor', 'none');
set(seShoadow_ratio_1(2),'FaceColor', [0.67, 0, 1]' , 'FaceAlpha',0.3, 'EdgeColor', 'none');
hold on;
seShoadow_ratio_2 = area(t_mov(7:end), [(mean(traceRatio(7:end,:), 2) - semRatio(7:end)), (2 * semRatio(7:end))]);
set(seShoadow_ratio_2(1),'Visible','off');
%set(seShoadow1_2(2),'FaceColor',[0.85, 0.32, 0]' ,'FaceAlpha',0.4, 'EdgeColor', 'none');
set(seShoadow_ratio_2(2),'FaceColor', [0.67, 0, 1]' , 'FaceAlpha',0.3, 'EdgeColor', 'none');
hold on;
% set(seShoadow2(1),'Visible','off');
%set(seShoadow2(2),'FaceColor',[0.46, 0.68, 0.17]' ,'FaceAlpha',0.5, 'EdgeColor', 'none');
line1 = plot(t_mov(1:6), meanRatio(1:6),'Color', [0, 0, 1], 'LineWidth',2);hold on
line2 = plot(t_mov(7:end), meanRatio(7:end),'Color', [0, 0, 1], 'LineWidth',2);

box off
% axis ([0 max(t) 0.90 1.30])
axis tight
title({'Trace ratio of voltage indicators'})
xlabel('Time (s)')
ylabel('Absolute ratio')
set(gcf,'color','w');
breakxaxis([150 1380]);
saveas(gca,['X:\91 Data and analysis\YJunqi\Screening method\System test\Statistics\20220404_gramicidinSummary_deltaTraceRatio' '.fig']);
saveas(gca,['X:\91 Data and analysis\YJunqi\Screening method\System test\Statistics\20220404_gramicidinSummary_deltaTraceRatio' '.png']);

%%
mean_ArcLightA242 = mean(ArcLightA242, 2);
mean_HASAP1JF635 = mean(HASAP1JF635, 2);
ArcLightA242_sem = std(ArcLightA242')'/sqrt(size(ArcLightA242, 2));
HASAP1JF635_sem = std(HASAP1JF635')'/sqrt(size(HASAP1JF635, 2));
traceRatio = HASAP1JF635./ArcLightA242;
meanRatio = mean(traceRatio, 2);
semRatio = std(traceRatio')'/sqrt(size(traceRatio, 2));
deltaRatio = traceRatio./mean(traceRatio(1:6,:))-1;
meanDeltaRatio = mean(deltaRatio, 2);
semDeltaRatio = std(deltaRatio')'/sqrt(size(deltaRatio, 2));

figure()
seShoadowRatio_1 = area(t_mov(1:6), [(mean(deltaRatio(1:6,:), 2) - semDeltaRatio(1:6)), (2 * semDeltaRatio(1:6))]);
set(seShoadowRatio_1(1),'Visible','off');
%set(seShoadow1_2(2),'FaceColor',[0.85, 0.32, 0]' ,'FaceAlpha',0.4, 'EdgeColor', 'none');
set(seShoadowRatio_1(2),'FaceColor',[0, 1, 1]' ,'FaceAlpha',0.3, 'EdgeColor', 'none');
hold on;
seShoadowRatio_2 = area(t_mov(7:end), [(mean(deltaRatio(7:end,:), 2) - semDeltaRatio(7:end)), (2 * semDeltaRatio(7:end))]);
set(seShoadowRatio_2(1),'Visible','off');
%set(seShoadow2_2(2),'FaceColor',[0.46, 0.68, 0.17]' ,'FaceAlpha',0.4, 'EdgeColor', 'none');
set(seShoadowRatio_2(2), 'FaceColor', [0, 1, 1]' , 'FaceAlpha',0.3, 'EdgeColor', 'none');
hold on;
plot1 = plot(t_mov(1:6), mean(deltaRatio(1:6,:),2),'Color', [0, 0, 1], 'LineWidth',2);
plot2 = plot(t_mov(7:end), mean(deltaRatio(7:end,:),2), 'Color', [0, 0, 1], 'LineWidth',2);
box off
% axis ([0 max(t) 0.90 1.30])
axis tight
title({'Dynamic range of voltage indicators'})
xlabel('Time (s)')
ylabel('ΔR/R0');
set(gcf,'color','w');
breakxaxis([150 1380]);
saveas(gca,['X:\91 Data and analysis\YJunqi\Screening method\Stimulus pre-experiments\System test\Statistics\20220404_gramicidinSummary_DeltaR2R0' '.fig']);
saveas(gca,['X:\91 Data and analysis\YJunqi\Screening method\Stimulus pre-experiments\System test\Statistics\20220404_gramicidinSummary_DeltaR2R0' '.png']);

