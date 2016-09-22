function [Vres, Vo] = geg_reslice(Vref,Vdata,interp)
% This function applies a reslicing when two images are not identical.
% It assumes that images are at least coregistered.
% 
%_______________________________________________________________________
% Copyright (C) 2015
%
% Gabriel Gonzalez-Escamilla
% $Id: geg_reslice.m 001 2015-30-10 18:43:36Z $
% 

m1 = struct2cell(Vref); m2 = struct2cell(Vdata);
Vi = [m1, m2];% new order of inputs to always match 1st image space
Headings = {'fname'; 'dim'; 'dt'; 'pinfo'; 'mat'; 'n'; 'descrip'; 'private'};
Vi = cell2struct(Vi, Headings, 1);
Vo = Vref.fname; cvt = Vdata.fname;
% interp = 0;%1=trilinear; 0=nearest;
dtype = spm_type('float32');
[sts, str] = spm_check_orientations(Vi, false);
if isstruct(Vo)
    Vi   = spm_cat_struct(Vo,Vi);
    refstr = 'output';
else
    Vi   = Vi(:);
    refstr = '1st';
end
Vo = struct('fname', cvt, 'dim', Vi(1).dim(1:3), 'dt', [dtype spm_platform('bigend')],...
    'pinfo', [Inf Inf Inf]', 'mat', Vi(1).mat, 'n', 1,...
    'descrip', 'PETPVE12');

if ~sts
    for i=1:size(str,1)
        fprintf('Warning: %s - using %s image.\n',strtrim(str(i,:)),refstr);
    end
    
    %-Computation
    Vres = zeros(Vo.dim(1:3));
    
    %-Start progress plot
    spm_progress_bar('Init',Vo.dim(3),'reslicing','planes completed');
    
    %-Loop over planes computing result Y
    for p = 1:Vo.dim(3)
        B = spm_matrix([0 0 -p 0 0 0 1 1 1]);%B    = spm_matrix([0 0 i]);%
        M = inv(B * inv(Vo.mat) * Vi(2).mat);
        d = spm_slice_vol(Vi(2), M, Vo.dim(1:2), [interp,NaN]);
        if prod(Vo.dim(1:2)) ~= numel(d)
            error(['"',f,'" produced incompatible image.']); end
        Vres(:,:,p) = reshape(d,Vo.dim(1:2));
        spm_progress_bar('Set',p);
    end
    spm_progress_bar('Clear')
else
    fprintf('Warning: The images do all have same dimensions/orientation. - not doing anything. \n');
    Vres = spm_read_vols(Vdata);
end