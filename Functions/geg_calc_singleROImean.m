function geg_calc_singleROImean()
% Returns on screen the mean SUV of a binary mask for all subjects
% This is the scrip of Ashburner to measure volume, but modified by M
% Grothe and converted to function by me
% 
%_______________________________________________________________________
% Copyright (C) 2015
%
% Gabriel Gonzalez-Escamilla
% $Id: geg_calc_singleROImean.m 001 2015-10-03 10:13:53Z $
% 
% 
% rev = '$Rev: 010 $'; % 25-October-2015

Vmask  = spm_vol(spm_select(Inf,'image','Select the binary masks'));
Vdata  = spm_vol(spm_select(Inf,'image','Select warped functionals'));
signal = zeros(numel(Vdata),numel(Vmask));
average = zeros(numel(Vdata),numel(Vmask));
for k=1:numel(Vmask),
    % sum up voxel-values within ROI defined by mask (slice-wise)
    disp(Vmask(k).fname)
    for j=1:numel(Vdata),
        fprintf('%d,',j)
        tot = 0;
        vox = 0;
        for i=1:Vdata(j).dim(3),
            M    = spm_matrix([0 0 i]);
            img  = spm_slice_vol(Vdata(j),M,Vdata(j).dim(1:2),0);
            Mmsk = Vmask(k).mat\Vdata(j).mat*M;
            mask = spm_slice_vol(Vmask(k),Mmsk,Vdata(j).dim(1:2),0);
            img  = img.*mask;
            img(~isfinite(img)) = 0;
            tot = tot + sum(img(img>0)); %sums only positive voxels
            img(img>0)=1; %sets positive voxels to 1
            vox = vox + sum(img(img>0)); %counts only positive voxels by summing 1s
        end
        signal(j,k) = tot;
        % divide sum of voxel-values by total number of contributing voxels
        % (i.e. voxels with positive value in (PET) image)
        average(j,k) = signal(j,k)/vox;
    end
    fprintf('\n')
    disp (average);
    disp(' ')
end
