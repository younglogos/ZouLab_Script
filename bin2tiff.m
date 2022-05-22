% parameters
subfolder = '';
basepath = 'X:\91 Data and analysis\YJunqi\Sensitivity\20200104 INS-1 cells\Dish2\cell1\213434_OD 2(488 nm)_long term';
pathname = basepath
path = [basepath subfolder];

% constants
movname = '\movie.bin';
ncol = 512;         % x
nrow = 104;         % y
bkg = 1600;          % background due to camera bias (100 for bin 1x1)
dt_mov = 2.0658;       % exposure time in millisecond

fname = [path movname];
[mov, nframe] = readBinMov(fname, nrow, ncol);

for i = 1:size(mov,3)
    imwrite(uint16(mov(:,:,i)), ['image', num2str(i), '.tif'], 'tif');
end

clear all
