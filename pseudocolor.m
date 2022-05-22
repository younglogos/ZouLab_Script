function [pseudo_cmap]=pseudocolor(WL);
base = getrgb(WL);
pseudo_cmap = zeros(256,3);
if base(1) ~= 0
    pseudo_cmap(:,1) = [0:base(1)/255:base(1)]';
else
    pseudo_cmap(:,1) = 0;
end
if base(2) ~= 0
    pseudo_cmap(:,2) = [0:base(2)/255:base(2)]';
else
    pseudo_cmap(:,2) = 0;
end
if base(3) ~= 0
    pseudo_cmap(:,3) = [0:base(3)/255:base(3)]';
else
    pseudo_cmap(:,3) = 0;
end