function geg_composite_from_atlas()
% This function creates a binary mask containing one or series of selected
% regions from a parcellation atlas. This is usefull to create masks
% excluding some regions. i.e. cortical composite used for amyloid-PET.
% 
%_______________________________________________________________________
% Copyright (C) 2015
%
% Gabriel Gonzalez-Escamilla
% $Id: geg_composite_from_atlas.m 001 2015-10-03 10:13:53Z $
% 
% 
% rev = '$Rev: 010 $'; % 25-October-2015

Vatlas  = spm_vol(spm_select(Inf,'image','Select the template atlas'));
AtlasDesc = spm_select([1,1],'^*\.txt','Select Atlas descriptor file (txt)');
saveName = spm_input('output name: ','+1','s','region_');
if strncmpi(saveName(end-4:end),'.nii',4)
    [pth,~,~] = spm_fileparts(Vatlas(1).fname);
    saveName = saveName(1:end-4); ext = saveName(end-4:end);
else
    [pth,~,ext] = spm_fileparts(Vatlas(1).fname);
end

AtlDesc=importdata(char(AtlasDesc));
ParcAtlasIds=AtlDesc.data;

%-Start progress plot
spm_progress_bar('Init',numel(ParcAtlasIds),'constructing composite','ROIs completed');

atlasi = spm_read_vols(Vatlas);
natlasi = zeros(size(atlasi));
for k=1:numel(ParcAtlasIds)
    actp = ParcAtlasIds(k);
    natlasi(atlasi==actp) = 1;    
    spm_progress_bar('Set',k); % progression indicator
end
spm_progress_bar('Clear')
newAtlas = Vatlas;
newAtlas.fname = fullfile(pth,[saveName,ext]);
newAtlas.descrip='PETPVE12';
spm_write_vol(newAtlas,natlasi);
disp('Done')