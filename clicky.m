function [roi_points, intens] = clicky(movie_in, refimg, message)
% function [roi_points, intens] = clicky(movie_in, refimg, message);
% movie_in: input movie.
% refimg: optional image (black and white or color) on which to click to
% select ROIs.
% Click to define an ROI, right click to complete.  Plots the intensity
% trace averaged over the ROI.  Keeps generating ROIs until you generate an
% ROI with no points (i.e. just a right click on the image).
% Returns coordinates of all ROIs and the intensity traces.
% AEC and JMK, 12 Oct. 2011.
% modified by PZ, 2014-08-20

if nargin == 1;
    refimg = mean(movie_in, 3);
    message = [];
elseif nargin == 2;
    message = [];
end;

nframes = size(movie_in, 3);

figure; clf;
subplot(2,1,1)
imshow(refimg, [], 'InitialMagnification', 'fit')
title(message);
hold on;

[ysize, xsize] = size(refimg(:,:,1));
npts = 1;
colorindex = 0;
order = get(gca,'ColorOrder');
nroi = 1;
intens = [];
[x, y] = meshgrid(1:xsize, 1:ysize);

while(npts > 0)
    subplot(2,1,1)
    [xv, yv] = (getline(gca, 'closed'));
    if size(xv,1) < 3  % exit loop if only a line is drawn
    break
    end
    inpoly = inpolygon(x,y,xv,yv);

    %draw the bounding polygons and label them
    currcolor = order(1+mod(colorindex,size(order,1)),:);
    plot(xv, yv, 'Linewidth', 1,'Color',currcolor);
    text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor,'FontSize',12);

    itrace = squeeze(sum(sum(movie_in.*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));

    subplot(2,1,2) % plot the trace
    hold on;
    plot(itrace,'Color',currcolor);
    colorindex = colorindex+1;

    intens = [intens; itrace'];
    roi_points{nroi} = [xv, yv];
    nroi = nroi + 1;
end
intens = intens';