function average = geg_calc_constants(Vdata,Vmask,wmcsfthresh)
% This computes everything on a slice-basis intead of image-basis, which is
% supposed to be more correct. Also, works when images have different
% dimensions, but are in the same space.
%
% Based on Ashburners script to measure a ROI volume
% Modified by M. Grothe to include only positive values
% Adapted & Integrated by G. Gonzalez-Escamilla
%_______________________________________________________________________
% $Id: geg_calc_constants.m 001 2015-10-01 14:55:42Z $
% 
% 
% rev = '$Rev: 002 $'; % 01-October-2015

tot = 0;
vox = 0;
for i=1:Vdata(1).dim(3),% this takes slices above z dimension (up->down)
    M    = spm_matrix([0 0 i]);
    img  = spm_slice_vol(Vdata(1),M,Vdata(1).dim(1:2),0);%counts the dimentions x and y on the brain image
    Mmsk = Vmask.mat\Vdata(1).mat*M;
    mask = spm_slice_vol(Vmask,Mmsk,Vdata(1).dim(1:2),0);%counts the dimentions x and y on the mask
    mask(mask<wmcsfthresh) = 0;
    mask(mask>wmcsfthresh) = 1;
    mask = round(mask); %found to be necessary for Uint8 image data    
    img  = img.*mask; % restrings the image to the mask
    img(~isfinite(img)) = 0;% give the total slices to obtain volume
    tot = tot + sum(img(img>0)); %to obtain total sum of each slice; sums the anterior to the current slice
                             % problematic if negative voxel values exist,
                             % then the (:) was replaced by (img>0)
    img(img>0)=1; %sets positive voxels to 1
    vox = vox + sum(img(img>0)); %counts positive voxels by summing 1s
end
signal = tot;
% divide sum of voxel-values by total number of contributing voxels
% (i.e. positive voxels)
average = signal/vox;


