%==========================================================================
% Inverse Deformation utilities
%==========================================================================
function out = geg_iydeformation(job)
% applies the inversed deformation to the image, called 'pull'
% job - a struct with fields
%       job.out{cell}.pull = structure with pull parameters (fnames, savedir, interp, mask, fwhm)
%       job.comp{cell}.def = Name of the deformation field (iy_*.nii)
% out - a struct with fields
%       .def        - file name of created deformation field
%       .warpedImg  - the warped Image (3D matrix)
%       .warpedName - file name of warped images
%_______________________________________________________________________
% Based on John Ashburner's: spm_deformations.m 6137
% 
% Copyright (C) 2015
%
% Gabriel Gonzalez-Escamilla
% $Id: geg_iydeformation.m 001 2015-03-02 10:16:23Z $
% 
% 
% rev = '$Rev: 004 $'; % 07-April-2015

[Def,mat] = get_comp(job.comp);
out = struct('warped',{{}});
fn = fieldnames(job.out{1});
fn = fn{1};
[out.warpedImg,out.warpedName] = pull_def(Def,mat,job.out{1}.(fn));

%==========================================================================
% Sub-functions
%==========================================================================
%==========================================================================
% function [Def,mat] = get_comp(job)
%==========================================================================
function [Def,mat] = get_comp(job)
% Return the composition of a number of deformation fields.
if isempty(job)
    error('Empty list of jobs in composition');
end
[Def,mat] = get_job(job{1});
for i=2:numel(job)
    Def1         = Def;
    mat1         = mat;
    [Def,mat]    = get_job(job{i});
    M            = inv(mat1);
    tmp          = zeros(size(Def),'single');
    tmp(:,:,:,1) = M(1,1)*Def(:,:,:,1)+M(1,2)*Def(:,:,:,2)+M(1,3)*Def(:,:,:,3)+M(1,4);
    tmp(:,:,:,2) = M(2,1)*Def(:,:,:,1)+M(2,2)*Def(:,:,:,2)+M(2,3)*Def(:,:,:,3)+M(2,4);
    tmp(:,:,:,3) = M(3,1)*Def(:,:,:,1)+M(3,2)*Def(:,:,:,2)+M(3,3)*Def(:,:,:,3)+M(3,4);
    Def(:,:,:,1) = single(spm_diffeo('bsplins',Def1(:,:,:,1),tmp,[1,1,1,0,0,0]));
    Def(:,:,:,2) = single(spm_diffeo('bsplins',Def1(:,:,:,2),tmp,[1,1,1,0,0,0]));
    Def(:,:,:,3) = single(spm_diffeo('bsplins',Def1(:,:,:,3),tmp,[1,1,1,0,0,0]));
    clear tmp
end
%==========================================================================
% function [Def,mat] = get_job(job)
%==========================================================================
function [Def,mat] = get_job(job)
% Determine what is required, and pass the relevant bit of the
% job out to the appropriate function.
fn = fieldnames(job);
fn = fn{1};
switch fn
    case {'comp'}
        [Def,mat] = get_comp(job.(fn));
    case {'def'}
        [Def,mat] = get_def(job.(fn));
    case {'dartel'}
        [Def,mat] = get_dartel(job.(fn));
    case {'id'}
        [Def,mat] = get_id(job.(fn));
    otherwise
        error('Unrecognised job type');
end
%==========================================================================
% function [Def,mat] = get_def(job)
%==========================================================================
function [Def,mat] = get_def(job)
% Load a deformation field saved as an image
Nii = nifti(job);%job{1}
Def = single(Nii.dat(:,:,:,1,:));
d   = size(Def);
if d(4)~=1 || d(5)~=3, error('Deformation field is wrong!'); end
Def = reshape(Def,[d(1:3) d(5)]);
mat = Nii.mat;
%==========================================================================
% function [Def,mat] = get_dartel(job)
%==========================================================================
function [Def,mat] = get_dartel(job)
Nii    = nifti(job.flowfield{1});
if ~isempty(job.template{1})
    [pth,nam] = fileparts(job.template{1});
    if exist(fullfile(pth,[nam '_2mni.mat']),'file')
        load(fullfile(pth,[nam '_2mni.mat']),'mni');
    else
        % Affine registration of Dartel Template with MNI space.
        fprintf('** Affine registering "%s" with MNI space **\n', nam);
        tpm        = fullfile(spm('Dir'),'tpm','TPM.nii');
        Mmni       = spm_get_space(tpm);
        Nt         = nifti(job.template{1});
        mni.affine = Mmni/spm_klaff(Nt,tpm);
        mni.code   = 'MNI152';
        save(fullfile(pth,[nam '_2mni.mat']),'mni', spm_get_defaults('mat.format'));
    end
    Mat    = mni.affine;
else
    Mat    = Nii.mat;
end
% Integrate a Dartel flow field
y0      = spm_dartel_integrate(Nii.dat,job.times,job.K);
if all(job.times == [0 1]),
    mat = Nii.mat0;
    Def = affine(y0,single(Mat));
else
    mat = Mat;
    Def = affine(y0,single(Nii.mat0));
end
%==========================================================================
% function out = pull_def(Def,mat,job)
%==========================================================================
function [outWarp,outNam] = pull_def(Def,mat,job)

PI      = job.fnames;
intrp   = job.interp;
intrp   = [intrp*[1 1 1], 0 0 0];
outNam     = cell(1,numel(PI));
if numel(PI)==0, return; end

if job.mask
    oM  = zeros(4,4);
    odm = zeros(1,3);
    dim = size(Def);
    msk = true(dim);
    for m=1:numel(PI)
        [pth,nam,ext,num] = spm_fileparts(PI{m});
        NI = nifti(fullfile(pth,[nam ext]));
        dm = NI.dat.dim(1:3);
        if isempty(num)
            j_range = 1:size(NI.dat,4);
        else
            num     = sscanf(num,',%d');
            j_range = num(1);
        end
        for j=j_range,
            M0 = NI.mat;
            if ~isempty(NI.extras) && isstruct(NI.extras) && isfield(NI.extras,'mat')
                M1 = NI.extras.mat;
                if size(M1,3) >= j && sum(sum(M1(:,:,j).^2)) ~=0
                    M0 = M1(:,:,j);
                end
            end
            M   = inv(M0);
            if ~all(M(:)==oM(:)) || ~all(dm==odm)
                tmp = affine(Def,M);
                msk = tmp(:,:,:,1)>=1 & tmp(:,:,:,1)<=size(NI.dat,1) ...
                    & tmp(:,:,:,2)>=1 & tmp(:,:,:,2)<=size(NI.dat,2) ...
                    & tmp(:,:,:,3)>=1 & tmp(:,:,:,3)<=size(NI.dat,3);
            end
            oM  = M;
            odm = dm;
        end
    end
end

oM = zeros(4,4);
spm_progress_bar('Init',numel(PI),'Resampling','volumes completed');
for m=1:numel(PI)

    % Generate headers etc for output images
    [pth,nam,ext,num] = spm_fileparts(PI{m});
    NI = nifti(fullfile(pth,[nam ext]));
    j_range = 1:size(NI.dat,4);
    k_range = 1:size(NI.dat,5);
    l_range = 1:size(NI.dat,6);
    if ~isempty(num)
        num = sscanf(num,',%d');
        if numel(num)>=1, j_range = num(1); end
        if numel(num)>=2, k_range = num(2); end
        if numel(num)>=3, l_range = num(3); end
    end
    NO = NI;
    if isfield(job.savedir,'savepwd'), wd = pwd;
    elseif isfield(job.savedir,'saveusr'), wd = job.savedir.saveusr{1};
    elseif isfield(job.savedir,'savesrc'), wd = pth;
    else wd = pwd;
    end
    if sum(job.fwhm.^2)==0
        NO.dat.fname = fullfile(wd,['w' nam ext]);
        NO.descrip   = sprintf('Warped');
    else
        NO.dat.fname = fullfile(wd,['sw' nam ext]);
        NO.descrip   = sprintf('Smoothed (%gx%gx%g subopt) warped',job.fwhm);
    end
    dim            = size(Def);
    dim            = dim(1:3);
    NO.dat.dim     = [dim NI.dat.dim(4:end)];
    NO.dat.offset  = 0; % For situations where input .nii images have an extension.
    NO.mat         = mat;
    NO.mat0        = mat;
    NO.mat_intent  = 'Aligned';
    NO.mat0_intent = 'Aligned';
    if isempty(num), outNam{m} = NO.dat.fname;
    else  outNam{m} = [NO.dat.fname, ',', num2str(num(1))]; 
    end
    NO.extras      = [];
    create(NO);

    % Smoothing settings
    vx  = sqrt(sum(mat(1:3,1:3).^2));
    krn = max(job.fwhm./vx,0.25);

    % Loop over volumes within the file
    for j=j_range
        M0 = NI.mat;
        if ~isempty(NI.extras) && isstruct(NI.extras) && isfield(NI.extras,'mat')
            M1 = NI.extras.mat;
            if size(M1,3) >= j && sum(sum(M1(:,:,j).^2)) ~=0
                M0 = M1(:,:,j);
            end
        end
        M  = inv(M0);
        if ~all(M(:)==oM(:))
            % Generate new deformation (if needed)
            Y = affine(Def,M);
        end
        oM = M;
        % Write the warped data for this time point
        for k=k_range
            for l=l_range
                C   = spm_diffeo('bsplinc',single(NI.dat(:,:,:,j,k,l)),intrp);
                dat = spm_diffeo('bsplins',C,Y,intrp);
                if job.mask
                    dat(~msk) = NaN;
                end
                if sum(job.fwhm.^2)~=0
                    spm_smooth(dat,dat,krn); % Side effects
                end
                NO.dat(:,:,:,j,k,l) = dat;
            end
        end
    end
    spm_progress_bar('Set',m);
    outWarp=NO.dat;
end
spm_progress_bar('Clear');
function Def = affine(y,M)
Def          = zeros(size(y),'single');
Def(:,:,:,1) = y(:,:,:,1)*M(1,1) + y(:,:,:,2)*M(1,2) + y(:,:,:,3)*M(1,3) + M(1,4);
Def(:,:,:,2) = y(:,:,:,1)*M(2,1) + y(:,:,:,2)*M(2,2) + y(:,:,:,3)*M(2,3) + M(2,4);
Def(:,:,:,3) = y(:,:,:,1)*M(3,1) + y(:,:,:,2)*M(3,2) + y(:,:,:,3)*M(3,3) + M(3,4);
%==========================================================================
% function [Def,mat] = get_id(job)
%==========================================================================
function [Def,mat] = get_id(job)
% Get an identity transform based on an image volume
[mat, dim] = spm_get_matdim(job.space{1});
Def = identity(dim, mat);
%==========================================================================
% function Def = identity(d,M)
%==========================================================================
function Def = identity(d,M)
[y1,y2]   = ndgrid(single(1:d(1)),single(1:d(2)));
Def       = zeros([d 3],'single');
for y3=1:d(3)
    Def(:,:,y3,1) = y1*M(1,1) + y2*M(1,2) + (y3*M(1,3) + M(1,4));
    Def(:,:,y3,2) = y1*M(2,1) + y2*M(2,2) + (y3*M(2,3) + M(2,4));
    Def(:,:,y3,3) = y1*M(3,1) + y2*M(3,2) + (y3*M(3,3) + M(3,4));
end
