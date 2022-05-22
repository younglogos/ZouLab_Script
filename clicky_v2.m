function [roi_points,intens] = clicky_v2(movie_in, dtmov, bin, refimg, bottom, top, titleofplot)
% function [roi_points, intens] = clicky(movie_in, message, refimg);
% movie_in: input movie.
% refimg: optional image (black and white or color) on which to click to
% select ROIs.
% bin was the binning settings, it should be 1,2, or 4 for removing
% background and setting contrust of refimg
% bottom and top are for contrust setting;
% titleofplot and titleofselection are for title of region selected by
% users and the plot of intensity of ROIs
% Click to define an ROI, right click to complete.  Plots the intensity
% trace averaged over the ROI.  Keeps generating ROIs until you generate an
% ROI with no points (i.e. just a right click on the image).
% Returns coordinates of all ROIs and the intensity traces.
% AEC and JMK, 12 Oct. 2011.
% modified by PZ, 2014-08-20
% modified by Luxin Peng, 2017-11-21

if nargin == 1;
    dtmov = 1;
    refimg = mean(movie_in, 3);
    bottom = min(refimg);
    top = max(refimg);
    titleofselection = 'Selected cells'
elseif nargin == 2;
    bin = 2;
    bottom = min(refimg);
    top = max(refimg);
    titleofselection = 'Selected cells'
elseif nargin == 3;
    bottom = min(refimg);
    top = max(refimg);
    titleofselection = 'Selected cells'
end;

nframes = size(movie_in, 3);

figure; clf;
subplot(3,1,1:2)
imshow(refimg, [bottom,top], 'InitialMagnification', 'fit')
title('Select ROIs, right click when done');
hold on;


[ysize, xsize] = size(refimg(:,:,1));
npts = 1;
colorindex = 0;
order = get(gca,'ColorOrder');
nroi = 1;
intens = [];
[x, y] = meshgrid(1:xsize, 1:ysize);

while(npts > 0)
       subplot(3,1,1:2)
    [xv, yv] = (getline(gca, 'closed'));
    if size(xv,1) < 3  % exit loop if only a line is drawn
       subplot(3,1,1:2)
       title('Selected cells');
       hold on;
    break
    end
    inpoly = inpolygon(x,y,xv,yv);

    %draw the bounding polygons and label them
    currcolor = order(1+mod(colorindex,size(order,1)),:);
    plot(xv, yv, 'Linewidth',1.5,'Color',currcolor);
    text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor,'FontSize',14);

    itrace = squeeze(sum(sum(movie_in.*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:))-100*power(bin,2);

    subplot(3,1,3) % plot the trace
    hold on;
    plot([0:dtmov:dtmov*(length(itrace)-1)]',itrace,'Color',currcolor);
    xlabel('times (s)');
    ylabel(strcat('Intensity W/O background, bin = ',num2str(bin)))
    title(titleofplot);
    colorindex = colorindex+1;

    intens = [intens; itrace'];
    roi_points{nroi} = [xv, yv];
    nroi = nroi + 1;
end
intens = intens';