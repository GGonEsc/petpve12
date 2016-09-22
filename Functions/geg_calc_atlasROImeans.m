function geg_calc_atlasROImeans()
% This computes everything on a slice-basis intead of image-basis, which is
% supposed to be more correct. It assumes that the parcelled-atlas and the
% functional volumes are in the same space or at least coregistered.
% It outputs an xls file with the values of each ROI (columns) and each
% subject (rows). 
% 
% Saves computed values to a xls file.
% 
% Based on Ashburners script to measure a ROI volume
% Modified by M. Grothe to include only positive values
% Adapted & Integrated by G. Gonzalez-Escamilla
%_______________________________________________________________________
% $Id: geg_calc_atlasROImeans.m 001 2015-10-05 22:44:53Z $
% 
% 
% rev = '$Rev: 002 $'; % 07-October-2015

Vmask  = spm_vol(spm_select([1,1],'image','Select the template atlas in avg space'));
AtlasDesc = spm_select([1,1],'^*\.txt','Select Atlas descriptor file');
Vdata  = spm_vol(spm_select(Inf,'image','Select warped functionals'));
opth = spm_select([1,1],'Dir','Select output directory');

AtlDesc=importdata(char(AtlasDesc));
ParcAtlasNames=AtlDesc.textdata;
ParcAtlasIds=AtlDesc.data;

signal = zeros(numel(Vdata),numel(ParcAtlasIds));
average = zeros(numel(Vdata),numel(ParcAtlasIds));
snames = cell(numel(Vdata),1);

[~,nam,~] = spm_fileparts(Vmask.fname);
tots_xls = fullfile(opth,['w',nam, '_ROItots.xls']);
xlswrite(tots_xls,['SubjID', ParcAtlasNames'],'ROItotals','A1')

for k=1:numel(ParcAtlasIds)
    actp = ParcAtlasIds(k);
    fprintf('%d,',actp) % progression indicator
    for j=1:numel(Vdata),
        tot = 0;
        vox = 0;
        [~,nam,~] = spm_fileparts(Vdata(j).fname);
        snames{j,1} = nam;
        % sum up voxel-values within ROI defined by mask (slice-wise)
        for i=1:Vdata(j).dim(3),
            M    = spm_matrix([0 0 i]);
            img  = spm_slice_vol(Vdata(j),M,Vdata(j).dim(1:2),0);
            Mmsk = Vmask(1).mat\Vdata(j).mat*M;
            atlas = spm_slice_vol(Vmask,Mmsk,Vdata(j).dim(1:2),0);
            atlas = round(atlas); %found to be necessary for Uint8 image data
            mask = atlas==actp;
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
    
end
fprintf('\n')
% print to a xls file
xlswrite(tots_xls,snames,'ROItotals','A2')
xlswrite(tots_xls,average,'ROItotals','B2')
disp('Done')
