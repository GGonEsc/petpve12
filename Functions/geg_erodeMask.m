function eroi = geg_erodeMask(ROI,erothresh,fwhm,save_eroded)
% Erodes a mask, using a Gaussian weighting function of FWHM and a threshold
% (erothresh). This is slightly different from only identify the edges of
% the area and erase them N times (N-voxels edge erosion). Even when both
% methods erode the edges of the image.
% As higher is the fwhm and the threshold more erosion will be the
% resultant image. 
% It has been design and programmed thinking on probability tissue
% segments, but it also works with binary images.
% 
%_______________________________________________________________________
% Copyright (C) 2015
%
% Gabriel Gonzalez-Escamilla
% $Id: geg_Segment.m 001 2015-03-30 09:28:03Z $
% 
% 
% rev = $Rev: 005 $'; 04-April-2015

fprintf('Eroding %s \n',ROI.fname)
volhandle = 1;
Vroi = ROI;%spm_vol(ROI);%spm_vol(varargin{3});
roi = spm_read_vols(Vroi)>0;
[x,y,z] = ndgrid(1:Vroi.dim(1),1:Vroi.dim(2),1:Vroi.dim(3));
xyz = [x(roi(:))'; y(roi(:))'; z(roi(:))'];
% reset data type to save disk space
Vroi.dt(1) = spm_type('uint8');
% reset scaling factor of Vroi handle
Vroi.pinfo(1:2) = Inf;
clear x y z
vols{volhandle}.roi = struct('Vroi',Vroi, 'xyz',xyz, 'roi',roi,...
    'hroi',1, 'mode','set', 'tool', 'box', ...
    'thresh',[60 140], 'box',[4 4 4],...
    'cb',[], 'polyslices',1, 'csize',5,...
    'erothresh',erothresh);%'erothresh',.5 %.99
V = zeros(size(vols{volhandle}.roi.roi));
spm_smooth(double(vols{volhandle}.roi.roi), V, fwhm);% smooths X mm a roi without saving, but filling V with the new matrix
[ero(1,:),ero(2,:),ero(3,:)] = ind2sub(vols{volhandle}.roi.Vroi.dim(1:3),...
    find(V(:)>vols{volhandle}.roi.erothresh));
if strcmp(vols{volhandle}.roi.mode,'set')
    toset = ero;
    toclear = vols{volhandle}.roi.xyz;
else
    toclear = ero;
end;
itoclear = sub2ind(vols{volhandle}.roi.Vroi.dim(1:3), ...
    toclear(1,:), toclear(2,:), toclear(3,:));
vols{volhandle}.roi.roi(itoclear) = false;
vols{volhandle}.roi.xyz = setdiff(vols{volhandle}.roi.xyz',toclear','rows')'; % check that is empty after this
itoset = round(sub2ind(vols{volhandle}.roi.Vroi.dim(1:3), ...
    toset(1,:), toset(2,:), toset(3,:)));
vols{volhandle}.roi.roi(itoset) = true;
eroi=vols{1, 1}.roi.roi;

if save_eroded    
    newm=Vroi;
    [pth,nam,ext] = spm_fileparts(newm.fname);
    newm.fname = fullfile(pth,['e' nam, ext]);% '_eroded'
    spm_write_vol(newm,eroi);
end