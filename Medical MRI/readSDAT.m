% The program uses the function freadVAXD-function from
% http://www.mathworks.com/matlabcentral/fileexchange/22675-read-vaxd-and-vaxg-files-in-r2008b-and-later
% to read a Philips SDAT spectroscopy file (an FID). The SPAR file
% describes the acquisition parameters. Adapt filename and number of points.
function [complexFID] = readSDAT(filename,want_to_see_plot)

fileid = fopen(filename, 'r', 'ieee-le');
Npts = 1024;

FID = freadVAXD(fileid, 2*Npts, 'single');
realFID = FID(1:2:Npts*2-1);
imagFID = FID(2:2:Npts*2);

complexFID = realFID + 1i*imagFID;


if want_to_see_plot
    figure
    plot(1:Npts,abs(complexFID))
end
% The figure shows that first point has already been divided by 2 by scanner to
% correct baseline after FFT.
% It also shows that there are rapidly decaying signal components that only
% affect the first fews samples. 
end