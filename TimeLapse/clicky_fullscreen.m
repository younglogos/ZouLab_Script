function [roi_points,intens,intens2,region_square] = clicky_fullscreen(movie_in, bin, refimg, bottom , top, titleofplot, movie_in2, wvlt)
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



nframes = size(movie_in, 3);

figure; clf;
set(gcf,'outerposition',get(0,'screensize'));
imshow(refimg, [bottom,top], 'InitialMagnification', 'fit')
colormap(pseudocolor(wvlt));
title('Select ROIs, right click when done');
hold on;


[ysize, xsize] = size(refimg(:,:,1));
npts = 1;
colorindex = 0;
order = get(gca,'ColorOrder');
nroi = 1;
intens = [];
intens2 = [];
itrace = [];
itrace2 = [];
region_square = [];
[x, y] = meshgrid(1:xsize, 1:ysize);

while(npts > 0)
    [xv, yv] = (getline(gca, 'closed'));
    if size(xv,1) < 3  % exit loop if only a line is drawn
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
    itrace2 = squeeze(sum(sum(movie_in2.*repmat(inpoly, [1, 1, size(movie_in2, 3)]))))/sum(inpoly(:))-100*power(bin,2);
    colorindex = colorindex+1;

    intens = [intens; itrace'];
    intens2 = [intens2; itrace2'];
    roi_points{nroi} = [xv, yv];
    region_square = [region_square sum(inpoly(:))];
    nroi = nroi + 1;
end
intens = intens';
intens2 = intens2';
