function varargout = geg_segment(job,arg)
% Segment a bunch of images
% FORMAT geg_segment(job)
%
% job is a structure containing the following fields:
% job.channel(n).vols{m}
% job.channel(n).biasreg
% job.channel(n).biasfwhm
% job.channel(n).write
% job.tissue(k).tpm
% job.tissue(k).ngaus
% job.tissue(k).native
% job.tissue(k).warped
% job.warp.affreg
% job.warp.reg
% job.warp.samp
% job.warp.write
% job.warp.SS
% job.warp.print
%
% See the user interface for a description of the fields.
%_______________________________________________________________________
% Copyright (C) 2015
%
% Gabriel Gonzalez-Escamilla
% $Id: geg_Segment.m 012 2015-03-16 20:13:03Z $

% rev = '$Rev: 012 $'; 24-March-2015

% check whether estimation & write of tissue segments should be performed
estwrite = isfield(job,'opts');

% set some defaults if segmentations are not estimated
if ~estwrite
    job.opts = struct('biasreg',0.001,'biasfwhm',60,'affreg','mni',...
        'warpreg',4,'samp',3,'ngaus',[2 2 2 3 4 2]);
end

channel = struct('vols',{job.data});

warp = struct('affreg', job.opts.affreg,...
    'samp', job.opts.samp,...
    'reg', job.opts.warpreg,...
    'write', job.output.warps,...
    'sanlm', job.extopts.sanlm,...
    'mrf', job.extopts.mrf,...
    'cleanup', job.extopts.cleanup,...
    'dartelwarp', isfield(job.extopts.dartelwarp,'normhigh'));
if isfield(job.extopts.dartelwarp,'normhigh')
    warp.darteltpm = job.extopts.dartelwarp.normhigh.darteltpm{1};
end

% prepare tissue priors and number of gaussians for all 6 classes
if estwrite
    tissue=struct;
    [pth,nam,ext,~] = spm_fileparts(job.opts.tpm{1});
    for i=1:6
        tissue(i).ngaus = job.opts.ngaus(i);
        tissue(i).tpm = [fullfile(pth,[nam ext]) ',' num2str(i)];
    end
end

tissue(1).warped = [0  0  0 ];
tissue(1).native = [job.output.GM.native  0     0    ];
tissue(2).warped = [0  0  0 ];
tissue(2).native = [job.output.WM.native  0     0    ];
tissue(3).warped = [0 0 0];
tissue(3).native = [job.output.CSF.native 0    0   ];

% never write class 4-6
for i=4:6
    tissue(i).warped = [0 0 0];
    tissue(i).native = [0 0 0];
end

job.bias     = job.output.bias.native;
job.label    = job.output.label.native;
job.jacobian = 0;
job.SS       = job.output.SS;
job.biasreg  = job.opts.biasreg;
job.biasfwhm = job.opts.biasfwhm;
job.channel  = channel;
job.warp     = warp;
job.warps    = job.output.warps;
job.tissue   = tissue;

if nargin == 1, arg = 'run'; end

switch lower(arg)
    case 'run'
        varargout{1} = run_job(job,estwrite);
    case 'check'
        varargout{1} = check_job(job);
    case 'vfiles'
        varargout{1} = vfiles_job(job);
    case 'vout'
        varargout{1} = vout_job(job);
    otherwise
        error('Unknown argument ("%s").', arg);
end
return
%_______________________________________________________________________
%=======================================================================
%_______________________________________________________________________
function vout = run_job(job, estwrite)
% This function computes the segmentation parameters from SPM12, and calls
% for the improved bias, skull-stripping and cleanup processes
%
vout   = vout_job(job);

% load tpm priors only for estimate and write
if estwrite % estimate and write segmentations
    tpm    =  char(cat(1,job.tissue(:).tpm));
    tpm    = spm_load_priors8(tpm);
end

for subj=1:numel(job.channel(1).vols), % Subject by subject
    if estwrite % estimate and write segmentations
        images = '';
        for n=1:numel(job.channel),
            images = strvcat(images,job.channel(n).vols{subj});
        end
        obj.image    = spm_vol(images);
        spm_check_orientations(obj.image);
        
        obj.fwhm    = 1;
        obj.biasreg  = cat(1,job.biasreg);
        obj.biasfwhm = cat(1,job.biasfwhm);
        obj.tpm      = tpm;
        obj.lkp      = [];
        if all(isfinite(cat(1,job.tissue.ngaus))),
            for k=1:numel(job.tissue),
                obj.lkp = [obj.lkp ones(1,job.tissue(k).ngaus)*k];
            end;
        end
        obj.reg      = job.warp.reg;
        obj.samp     = job.warp.samp;
        
        [pth,nam] = fileparts(job.channel(1).vols{subj});
        if ~exist(fullfile(pth,[nam '_seg8.mat']),'file');%iter==1,
            % Initial affine registration.
            Affine  = eye(4); % Start creating an affine matrix such that the mapping
            % from voxels in the individual to those in the template
            % is by tpm.M\Affine*obj.image(1).mat. in the Combined
            % Segmentation and Spatial Normalisation --> spm_preproc8.m
            % function
            if ~isempty(job.warp.affreg),
                fprintf('Estimating segmentation parameters..\n');
                % Use old SPM Template (more blurred)
                % VG = spm_vol(fullfile(spm('Dir'),'toolbox','OldNorm','T1.nii'));
                % Use SPM12 template (More defined)
                VG = spm_vol(fullfile(spm('Dir'),'canonical','avg152T1.nii'));
                VF = spm_vol(obj.image(1));
                
                % smooth source with 8mm
                VF1 = spm_smoothto8bit(VF,8);
                
                % Rescale images so that globals are better conditioned
                VF1.pinfo(1:2,:) = VF1.pinfo(1:2,:)/spm_global(VF1);
                VG.pinfo(1:2,:)  = VG.pinfo(1:2,:)/spm_global(VG);
                
                fprintf('Initial Coarse Affine Registration..\n');
                aflags    = struct('sep',8, 'regtype',job.warp.affreg,...
                    'WG',[],'WF',[],'globnorm',0);
                aflags.sep = max(aflags.sep,max(sqrt(sum(VG(1).mat(1:3,1:3).^2))));
                aflags.sep = max(aflags.sep,max(sqrt(sum(VF(1).mat(1:3,1:3).^2))));
                
                M = eye(4);
                try
                    spm_plot_convergence('Init','Coarse Affine Registration','Mean squared difference','Iteration');
                catch
                    spm_chi2_plot('Init','Coarse Affine Registration','Mean squared difference','Iteration');
                end
                [Affine, scale]  = spm_affreg(VG, VF1, aflags, M);
                
                aflags.WG  = spm_vol(fullfile(spm('Dir'),'toolbox','FieldMap','brainmask.nii'));
                aflags.sep = aflags.sep/2;
                try
                    spm_plot_convergence('Init','Fine Affine Registration','Mean squared difference','Iteration');
                catch
                    spm_chi2_plot('Init','Fine Affine Registration','Mean squared difference','Iteration');
                end
                Affine  = spm_affreg(VG, VF1, aflags, Affine, scale);
                
                fprintf('Fine Affine Registration..\n');
                Affine  = spm_maff8(obj.image(1),job.warp.samp,obj.fwhm,  tpm,Affine,job.warp.affreg);
            end;
            obj.Affine = Affine; 
            res = spm_preproc8(obj);
            savefields(fullfile(pth,[nam '_seg8.mat']),res);% writes the segmentation parameters matrix
        else
            % Load results from previous iteration of the new group-specific tissue probability map.
            fprintf('Previous Segmentation parameters for image %s founded ...\n',nam)
            res       = load(fullfile(pth,[nam '_seg8.mat']));
            obj.Affine = res.Affine;
            obj.Twarp  = res.Twarp;
            obj.Tbias  = res.Tbias;
            if ~isempty(obj.lkp),
                obj.mg     = res.mg;
                obj.mn     = res.mn;
                obj.vr     = res.vr;
            end
        end
    end % This end was introduced to eliminate some code
    % Write out the required data.
    tc = [cat(1,job.tissue(:).native) cat(1,job.tissue(:).warped)];
    bf = job.bias;
    df = job.warp.write;
    lb = job.label;
    geg_seg_write(res, tc, bf, df, lb, job.warp, tpm, job) % segment, SS and PVE-labels
    disp('Done')
    fprintf('\n')
end
return
%_______________________________________________________________________
%_______________________________________________________________________
function msg = check_job(job)
msg = {};
if numel(job.channel) >1,
    k = numel(job.channel(1).vols);
    for i=2:numel(job.channel),
        if numel(job.channel(i).vols)~=k,
            msg = {['Incompatible number of images in channel ' num2str(i)]};
            break
        end
    end
elseif numel(job.channel)==0,
    msg = {'No data'};
end
return
%_______________________________________________________________________
%_______________________________________________________________________
function savefields(fnam,p)
if length(p)>1, error('Can''t save fields.'); end;
fn = fieldnames(p);
if numel(fn)==0, return; end;
for i=1:length(fn),
    eval([fn{i} '= p.' fn{i} ';']);
end;
try
    mat_ver = spm_check_version('matlab','7');
catch
    mat_ver = spm_matlab_version_chk('7');
end
if mat_ver >= 0
    save(fnam,'-V6',fn{:});
else
    save(fnam,fn{:});
end;

return;
%_______________________________________________________________________
%_______________________________________________________________________
function vout = vout_job(job)
% setting-up output file names definitions
n     = numel(job.channel(1).vols);
parts = cell(n,4);
biascorr  = {};
label  = {};

for j=1:n,
    [parts{j,:}] = spm_fileparts(job.channel(1).vols{j});
end

if job.bias(1),
    biascorr = cell(n,1);
    for j=1:n
        biascorr{j} = fullfile(parts{j,1},['m',parts{j,2},'.nii']);
    end
end

if job.label(1),
    label = cell(n,1);
    for j=1:n
        label{j} = fullfile(parts{j,1},['c0',parts{j,2},'.nii']);
    end
end

for j=1:n,
    [parts{j,:}] = spm_fileparts(job.channel(1).vols{j});
end
param = cell(n,1);
for j=1:n
    param{j} = fullfile(parts{j,1},[parts{j,2},'_seg8.mat']);
end

tiss = struct('c',{});
for i=1:numel(job.tissue),
    if job.tissue(i).native(1),
        tiss(i).c = cell(n,1);
        for j=1:n
            tiss(i).c{j} = fullfile(parts{j,1},['c',num2str(i),parts{j,2},'.nii']);
        end
    end
end

if job.warp.write(1),
    fordef = cell(n,1);
    for j=1:n
        fordef{j} = fullfile(parts{j,1},['y_',parts{j,2},'.nii']);
    end
else
    fordef = {};
end

if job.warp.write(2),
    invdef = cell(n,1);
    for j=1:n
        invdef{j} = fullfile(parts{j,1},['iy_',parts{j,2},'.nii']);
    end
else
    invdef = {};
end

if job.SS.nativeSS(1),
    T1ss = cell(n,1);
    for j=1:n
        T1ss{j} = fullfile(parts{j,1},['ss',parts{j,2},'.nii']);
    end
else
    T1ss = {};
end

if job.SS.nativeSS(2),
    BIASss = cell(n,1);
    for j=1:n
        BIASss{j} = fullfile(parts{j,1},['ss',parts{j,2},'.nii']);
    end
else
    BIASss = {};
end

vout  = struct('tiss',tiss,'label',{label},...
    'biascorr',{biascorr},'param',{param},...
    'invdef',{invdef},'fordef',{fordef},'T1ss',{T1ss},'BIASss',{BIASss});
%_______________________________________________________________________
%_______________________________________________________________________
function vf = vfiles_job(job)
vout = vout_job(job);
vf   = vout.param;
if ~isempty(vout.invdef),     vf = {vf{:}, vout.invdef{:}}; end
if ~isempty(vout.fordef),     vf = {vf{:}, vout.fordef{:}}; end
if ~isempty(vout.jacobian),   vf = {vf{:}, vout.jacobian{:}}; end

if ~isempty(vout.biascorr),   vf = {vf{:}, vout.biascorr{:}}; end
if ~isempty(vout.wbiascorr),  vf = {vf{:}, vout.wbiascorr{:}}; end
if ~isempty(vout.label),      vf = {vf{:}, vout.label{:}}; end
if ~isempty(vout.wlabel),     vf = {vf{:}, vout.wlabel{:}}; end
if ~isempty(vout.rlabel),     vf = {vf{:}, vout.rlabel{:}}; end
if ~isempty(vout.alabel),     vf = {vf{:}, vout.alabel{:}}; end

for i=1:numel(vout.tiss)
    if ~isempty(vout.tiss(i).c),   vf = {vf{:}, vout.tiss(i).c{:}};   end
    if ~isempty(vout.tiss(i).rc),  vf = {vf{:}, vout.tiss(i).rc{:}};  end
    if ~isempty(vout.tiss(i).rca), vf = {vf{:}, vout.tiss(i).rca{:}}; end
    if ~isempty(vout.tiss(i).wc),  vf = {vf{:}, vout.tiss(i).wc{:}};  end
    if ~isempty(vout.tiss(i).mwc), vf = {vf{:}, vout.tiss(i).mwc{:}}; end
    if ~isempty(vout.tiss(i).m0wc),vf = {vf{:}, vout.tiss(i).m0wc{:}};end
end
vf = reshape(vf,numel(vf),1);
%_______________________________________________________________________
%_______________________________________________________________________
%=======================================================================
%=======================================================================
%_______________________________________________________________________
%_______________________________________________________________________
function cls = geg_seg_write(res,tc,bf,df,lb,warp,tpm,job)
% Write out preprocessed data
% FORMAT cls = geg_seg_write(res,tc,bf,df)
%__________________________________________________________________________
% Copyright (C) 2015

% based on Christian Gaser's version of
% cg_vbm8_write.m 439 2012-03-23 15:46:43Z gaser $
% $Id: Gabriel Gonzalez-Escamilla 001 2015-03-16 $

% tc - tissue classes: native, dartel-rigid, dartel-affine, warped, warped-mod, warped-mod0
% bf - bias field: corrected, warp corrected, affine corrected
% df - deformations: forward, inverse
% lb - label: native, warped label, rigid label, affine label
% jc - jacobian: no, normalized

% get current release number
A = ver;
for i=1:length(A)
    if strcmp(A(i).Name,'PET Partial Volume Effects correction Toolbox')
        r = str2double(A(i).Version);
    end
end

if exist('r','var')
    fprintf('PETPVE12 r%d: %s\n',r,res.image.fname);
end

if ~isstruct(tpm) || (~isfield(tpm, 'bg1') && ~isfield(tpm, 'bg')),
    tpm = spm_load_priors8(tpm);
end

% Defining what to do (default & chosen parameters):
SS_opt=job.SS.nativeSS;
bias_fwhm   = geg_petpve12_get_defaults('extopts.bias_fwhm');
init_kmeans = geg_petpve12_get_defaults('extopts.kmeans');
finalmask   = job.extopts.finalmask;
gcut        = job.extopts.gcut;
cleanup     = job.extopts.cleanup;
write_bmask = job.SS.BMoutput;
warp.open_th = 0.25; % initial threshold for skull-stripping
warp.dilate = 1; % number of final dilations for skull-stripping
apply_dartel = warp.dartelwarp;   % apply dartel normalization
save_warped = job.extopts.dartelwarp.normhigh.svwarped;
svdartel_exports = job.extopts.dartelwarp.normhigh.svexports;
svjc = job.extopts.dartelwarp.normhigh.svjacobian;

do_cls   = any(tc(:)) || any(lb) || any(df) || nargout>1;
do_defs = any(df) ||  any(tc(:,2));
do_defs = do_cls || do_defs;

histeq_deep = 0;
try
    histeq_deep = geg_petpve12_get_defaults('extopts.histeq_deep');
catch
end

d1        = size(tpm.dat{1});
d1        = d1(1:3);
M1        = tpm.M;

if isfield(res,'mg'),
    lkp = res.lkp;
    Kb  = max(lkp);
else
    Kb  = size(res.intensity(1).lik,2);
end

N   = numel(res.image);
if N > 1
    warning('PETPVE12 does not support multiple channels. Only the first channel will be used.');
end

[pth,nam] = fileparts(res.image(1).fname);
ind  = res.image(1).n;
d    = res.image(1).dim(1:3);

[x1,x2,o] = ndgrid(1:d(1),1:d(2),1);
x3  = 1:d(3);

% chan is going to be the bias-corrected image (m*.nii)
chan(N) = struct('B1',[],'B2',[],'B3',[],'T',[],'Nc',[],'Nf',[],'ind',[]);
for n=1:N,
    d3         = [size(res.Tbias{n}) 1];
    chan(n).B3 = spm_dctmtx(d(3),d3(3),x3);
    chan(n).B2 = spm_dctmtx(d(2),d3(2),x2(1,:)');
    chan(n).B1 = spm_dctmtx(d(1),d3(1),x1(:,1));
    chan(n).T  = res.Tbias{n};
    
    [pth1,nam1,~] = fileparts(res.image(n).fname);%[pth1,nam1,ext1]
    chan(n).ind      = res.image(n).n;
    
    if bf(n,1),
        chan(n).Nc      = nifti;
        chan(n).Nc.dat  = file_array(fullfile(pth1,['m1', nam1, '.nii']),...
            res.image(n).dim(1:3),...
            [spm_type('float32') spm_platform('bigend')],...
            0,1,0);
        chan(n).Nc.mat  = res.image(n).mat;
        chan(n).Nc.mat0 = res.image(n).mat;
        chan(n).Nc.descrip = 'Bias corrected';
        create(chan(n).Nc);
    end
end

tiss(Kb) = struct('Nt',[]);
cls      = cell(1,Kb);
for k1=1:Kb, % kb is 6 because there are 6 compartments in the TPM
    cls{k1} = zeros(d(1:3),'uint8');
    if tc(k1,1),
        tiss(k1).Nt      = nifti;
        tiss(k1).Nt.dat  = file_array(fullfile(pth,['c', num2str(k1), nam, '.nii']),...
            res.image(1).dim(1:3),...
            [spm_type('int16') spm_platform('bigend')],...
            0,1/255,0);
        tiss(k1).Nt.mat  = res.image(n).mat;
        tiss(k1).Nt.mat0 = res.image(n).mat;
        tiss(k1).Nt.descrip = ['Tissue class ' num2str(k1)];
        create(tiss(k1).Nt);
        do_cls = true;
    end;
end

prm     = [3 3 3 0 0 0];
Coef    = cell(1,3);
Coef{1} = spm_bsplinc(res.Twarp(:,:,:,1),prm);
Coef{2} = spm_bsplinc(res.Twarp(:,:,:,2),prm);
Coef{3} = spm_bsplinc(res.Twarp(:,:,:,3),prm);

if do_defs, % create the deformation fields
    if df(2),
        [pth,nam,~]=fileparts(res.image(1).fname);%[pth,nam,ext1]
        Ndef      = nifti;
        Ndef.dat  = file_array(fullfile(pth,['iy_', nam1, '.nii']),...
            [res.image(1).dim(1:3),1,3],...
            [spm_type('float32') spm_platform('bigend')],...
            0,1,0);
        if apply_dartel
            Ndef.dat.fname = fullfile(pth,['iy_r', nam1, '.nii']);
        end        
        Ndef.mat  = res.image(1).mat;
        Ndef.mat0 = res.image(1).mat;
        Ndef.descrip = 'Inverse Deformation';
        create(Ndef);
    end
end

spm_progress_bar('init',length(x3),['Working on ' nam],'Planes completed');
M = tpm.M\res.Affine*res.image(1).mat;
if histeq_deep
    tmp_histeq_mask = spm_vol(char(geg_petpve12_get_defaults('extopts.histeq_mask')));
    histeq_mask = zeros(d(1:3),'uint8');
    M2 = tmp_histeq_mask.mat\res.Affine*res.image(1).mat;
end
for z=1:length(x3),
    
    % Bias corrected image
    cr = cell(1,N);
    for n=1:N,% Write plane by plane
        f = spm_sample_vol(res.image(n),x1,x2,o*x3(z),0);
        bf1 = exp(transf1(chan(n).B1,chan(n).B2,chan(n).B3(z,:),chan(n).T));
        bf1(bf1>100) = 100;
        cr{n} = bf1.*f;
        
        % Write a plane of bias corrected data
        if bf(n,1),
            chan(n).Nc.dat(:,:,z,chan(n).ind(1),chan(n).ind(2)) = cr{n};
        end
        if ~isempty(chan(n).Nf),
            % Write a plane of bias field
            chan(n).Nf.dat(:,:,z,chan(n).ind(1),chan(n).ind(2)) = bf1;
        end;
    end
    
    if do_defs,
        [t1,t2,t3] = defs(Coef,z,res.MT,prm,x1,x2,x3,M);
        if exist('Ndef','var'), % creates the inverse flow_field for every TPM
            tmp = tpm.M(1,1)*t1 + tpm.M(1,2)*t2 + tpm.M(1,3)*t3 + tpm.M(1,4);
            Ndef.dat(:,:,z,1,1) = tmp;
            tmp = tpm.M(2,1)*t1 + tpm.M(2,2)*t2 + tpm.M(2,3)*t3 + tpm.M(2,4);
            Ndef.dat(:,:,z,1,2) = tmp;
            tmp = tpm.M(3,1)*t1 + tpm.M(3,2)*t2 + tpm.M(3,3)*t3 + tpm.M(3,4);
            Ndef.dat(:,:,z,1,3) = tmp;
        end
        
        if do_cls,
            msk = (f==0) | ~isfinite(f);
            
            if isfield(res,'mg'),
                % Bayes theorem allows one to formally incorporate prior
                % knowledge into computing statistical probabilities.
                %
                % The “posterior” probability of the parameters given the
                % data is an optimal combination of prior knowledge and new
                % data, weighted by their relative precision.
                %
                %       posterior ~= likelihood * prior_distribution
                %   so, posterior ~= new_Data * existing_Data
                % then, joint_probability = new_Data * existing_Data
                %
                q   = zeros([d(1:2) Kb]);
                q1  = likelihoods(cr,[],res.mg,res.mn,res.vr);
                q1  = reshape(q1,[d(1:2),numel(res.mg)]);
                b   = spm_sample_priors8(tpm,t1,t2,t3);
                for k1=1:Kb,
                    q(:,:,k1) = sum(q1(:,:,lkp==k1),3).*b{k1};
                end
            else
                q   = spm_sample_priors8(tpm,t1,t2,t3);
                q   = cat(3,q{:});
                for n=1:N,
                    tmp = round(cr{n}*res.intensity(n).interscal(2) + res.intensity(n).interscal(1));
                    tmp = min(max(tmp,1),size(res.intensity(n).lik,1));
                    for k1=1:Kb,
                        likelihood = res.intensity(n).lik(:,k1);
                        q(:,:,k1)  = q(:,:,k1).*likelihood(tmp);
                    end
                end
            end
            
            if histeq_deep % To apply a histogram equalization to every slice
                [t01,t02,t03] = defs(Coef,z,res.MT,prm,x1,x2,x3,M2);
                histeq_mask(:,:,z) = uint8(round(spm_sample_vol(tmp_histeq_mask,t01,t02,t03,0)));
            end
            
            sq = sum(q,3) + eps^2;
            for k1=1:Kb,
                tmp            = q(:,:,k1);
                tmp(msk)       = 0;
                tmp            = tmp./sq;
                if ~isempty(cls{k1}),
                    cls{k1}(:,:,z) = uint8(round(255 * tmp));
                end
            end
        end
        if save_warped || svdartel_exports
            % initialize y only at first slice
            if z==1
                y = zeros([res.image(1).dim(1:3),3],'single');
            end
            % Here it creates the field
            y(:,:,z,1) = t1;
            y(:,:,z,2) = t2;
            y(:,:,z,3) = t3;
        end
        
    end
    spm_progress_bar('set',z);
end
spm_progress_bar('clear');

clear q q1 Coef b cr

% Here makes a copy of the bias corrected image (m0), the next bias image
% will be histogram equalized (m1), and the user might not want it like
% that. This keeps the bias image before the clean upn (m1 is always SS).
if exist(fullfile(pth,['m1', nam, '.nii']),'file')
    copyfile(fullfile(pth,['m1', nam, '.nii']),fullfile(pth,['m0', nam, '.nii']))
end

% % Preparation for skull-stripping
% load bias corrected image
src = zeros(res.image(1).dim(1:3),'single');
for z=1:length(x3),
    f = spm_sample_vol(res.image(1),x1,x2,o*x3(z),0);
    bf1 = exp(transf1(chan(1).B1,chan(1).B2,chan(1).B3(z,:),chan(1).T));
    % restrict bias field to maximum of 100
    % (sometimes artefacts at the borders can cause huge values in bias field)
    bf1(bf1>100) = 100;
    src(:,:,z) = single(bf1.*f);
end

% Optionally apply non local means denoising filter (NLMdf)
% This is the noise reduction function needed for skull-stripping
%
% prevent NaNs
src(isnan(src)) = 0;
% for windows disable multi-threading
if strcmp(mexext,'mexw32') || strcmp(mexext,'mexw64')
    warp.sanlm = min(1,warp.sanlm);
end
switch warp.sanlm
    case 0
    case 1 % use single-threaded version
        fprintf('Applying Spatial Adaptive Non-Local Means (SANLM) Filter\n')
        SANLMdf(src,'noMT')
    otherwise % use multi-threaded version
        fprintf('Applying Spatial Adaptive Non-Local Means (SANLM) Filter with multi-threading\n')
        SANLMdf(src,'MT')
end

% The skull-stripping
if do_cls && do_defs,
    vx_vol = sqrt(sum(res.image(1).mat(1:3,1:3).^2));
    scale_morph = 1/mean(vx_vol);
    if gcut
        % skull-stripping using graph-cuts
        opt.verb = 0; % display process (0=nothing, 1=only points, 2=times)
        fprintf('Skull-stripping using graph-cuts\n');
        cls_old = cls;
        try
            % src = bias corrected T1 image after noise reduction (single)
            % cls = 6 uint8 tissue class images (GM,WM,CSF,.,.,head) (cell)
            % res = .mat containing the parameters of the segmentation
            % opt =     main parameter options (optional input)(struct)
            %  .RSS           remove sinus sagittalis
            %                 controls the detection (has to be greater or equal one)
            %                 (RSS=1 more tissue, RSS=2.0 medium tissue (default), RSS>2 less tissue)
            %  .BVH           remove high intensity blood vessels
            %  .BVN           remove low  intensity blood vessels
            %  .FILL          detect ventricle and deep GM structures
            %                 create segment map with filled ventrile and deep GM structures
            [src,cls,mask] = SSGC(src,cls,res,opt); % Skull-Stripping with Graph-Cuts
        catch
            fprintf('Graph-cut failed\n');
            gcut = 0;
        end
        % check whether graph-cut failed
        if (sum(mask(:))/sum((single(cls_old{1}(:))+single(cls_old{2}(:))+single(cls_old{3}(:)))/255)<0.8)
            fprintf('Graph-cut failed\n');
            gcut = 0;
            cls = cls_old;
        end
        clear cls_old
    end
    if ~gcut
        fprintf('Skull-stripping using morphological operations\n');
        % use mask of GM and WM
        mask = single(cls{1});
        mask = mask + single(cls{2});
        
        % keep largest connected component after at least 1 iteration of opening
        n_initial_openings = max(1,round(scale_morph*warp.cleanup));
        mask = cg_morph_vol(mask,'open',n_initial_openings,warp.open_th);
        mask = mask_largest_cluster(mask,0.5);
        
        % dilate and close to fill ventricles
        mask = cg_morph_vol(mask,'dilate',warp.dilate,0.5);
        mask = cg_morph_vol(mask,'close',round(scale_morph*10),0.5);
        
        % remove sinus
        mask = mask & ((single(cls{5})<single(cls{1})) | ...
            (single(cls{5})<single(cls{2})) | ...
            (single(cls{5})<single(cls{3})));
        
        % fill holes that may remain
        mask = cg_morph_vol(mask,'close',round(scale_morph*2),0.5);
    end
    
    % calculate label image for all classes
    cls2 = zeros([d(1:2) Kb]);
    label2 = zeros(d,'uint8');
    for i=1:d(3)
        for k1 = 1:Kb
            cls2(:,:,k1)  = cls{k1}(:,:,i);
        end
        % find maximum for reordered segmentations
        [maxi,maxind] = max(cls2(:,:,[3,1,2,4:Kb]),[],3);
        for k1 = 1:Kb
            label2(:,:,i) = label2(:,:,i) + uint8((maxind == k1).*(maxi~=0)*k1);
        end
    end
    
    % set all non-brain tissue outside mask to 0
    label2(mask == 0)  = 0;
    
    % and for skull/bkg tissue classes to 0
    label2(label2 > 3) = 0;
    
    % fill remaining holes in label with 1
    mask = cg_morph_vol(label2,'close',round(scale_morph*2),0);
    label2((label2 == 0) & (mask > 0)) = 1;
    
    % use index to speed up and save memory
    sz = size(mask);
    [indx, indy, indz] = ind2sub(sz,find(mask>0));
    indx = max((min(indx) - 1),1):min((max(indx) + 1),sz(1));
    indy = max((min(indy) - 1),1):min((max(indy) + 1),sz(2));
    indz = max((min(indz) - 1),1):min((max(indz) + 1),sz(3));
    
    label = label2(indx,indy,indz);
    
    clear cls2 label2
    
    vol = double(src(indx,indy,indz));
    
    % mask source image because Amap needs a skull stripped image
    % set label and source inside outside mask to 0
    vol(mask(indx,indy,indz)==0) = 0;
    
    % use local histogram equalization for bias corrected image
    if histeq_deep
        
        clear tmp_histeq_mask t01 t02 t03
        
        histeq_mask_ind = histeq_mask(indx,indy,indz);
        clear histeq_mask;
        
        outermask_gt_0 = (vol > 0 ) & (histeq_mask_ind == 0);
        max_vol = max(vol(outermask_gt_0));
        vol = vol/max_vol;
        h1 = hist(vol(outermask_gt_0),512);
        
        histeq_mask_gt_0 = histeq_mask_ind >0;
        vol(histeq_mask_gt_0) = histeq_deep*histeq(vol(histeq_mask_gt_0),h1) + (1-histeq_deep)*vol(histeq_mask_gt_0);
        
        vol = vol*max_vol;
        
        src(:) = 0;
        src(indx, indy, indz) = vol;
        
        % Here applies the histogram equalization and the NLMdf to the bias corrected image
        for z=1:length(x3),
            % Bias corrected image
            % Write a plane of bias corrected data
            if bf(1,1),
                chan(1).Nc.dat(:,:,z,chan(1).ind(1),chan(1).ind(2)) = src(:,:,z);
            end
        end
        
    end    
    clear chan
    
    % Adaptive maximum a posteriori estimations (Amap)(Gaser, 2009)
    % Amap parameters
    n_iters = 200; sub = 16; n_classes = 3; pve = 5;
    iters_icm = 20;    
    if init_kmeans, fprintf('Adaptive Maximum A Posteriori (AMAP) Approach with Kmeans\n');
    else            fprintf('Adaptive Maximum A Posteriori (AMAP) without Kmeans\n');
    end    
    [prob, ~] = AmapMex(vol, label, n_classes, n_iters, sub, pve, init_kmeans, warp.mrf, vx_vol, iters_icm, bias_fwhm);%[prob, means]
    
    % reorder probability maps according to spm order
    % For Tohka Partial volume estimation in brain MRI the segments are:
    % csf = label 1; gm = label 2; wm = label 3
    prob = prob(:,:,:,[2 3 1]);
    clear vol
    
    % use cleanup
    if cleanup
        % get sure that all regions outside mask are zero
        for i=1:3
            cls{i}(:) = 0;
        end
        fprintf('Cleanning up... \n');
        [cls{1}(indx,indy,indz), cls{2}(indx,indy,indz), cls{3}(indx,indy,indz)] = cg_cleanup_gwc(prob(:,:,:,1), ...
            prob(:,:,:,2), prob(:,:,:,3), warp.cleanup);
        sum_cls = cls{1}(indx,indy,indz)+cls{2}(indx,indy,indz)+cls{3}(indx,indy,indz);
        label(sum_cls<0.15*255) = 0;
    else
        for i=1:3
            cls{i}(:) = 0;
            cls{i}(indx,indy,indz) = prob(:,:,:,i);
        end
    end;
    clear prob
    
    if finalmask
        fprintf('Applying final masking\n');
        % create final mask
        mask = single(cls{1});
        mask = mask + single(cls{2});
        
        % keep largest connected component after at least 1 iteration of opening
        n_initial_openings = max(1,round(scale_morph*2));
        mask = cg_morph_vol(mask,'open',n_initial_openings,0.5);
        mask = mask_largest_cluster(mask,0.5);
        
        % dilate and close to fill ventricles
        mask = cg_morph_vol(mask,'dilate',2,0.5);
        mask = cg_morph_vol(mask,'close',20,0.5);
        
        ind_mask = find(mask == 0);
        for i=1:3
            cls{i}(ind_mask) = 0;
        end
        
        % mask label
        label2 = zeros(d,'uint8');
        label2(indx,indy,indz) = label;
        label2(ind_mask) = 0;
        
        label = label2(indx,indy,indz);
        clear label2
    end
    % clear last 3 tissue classes to save memory
    for i=4:6
        cls{i} = [];
    end
end

M0 = res.image(1).mat;

% Here goes the DARTEL & preparation for affine and rigid transforms 
% run dartel registration to GM/WM dartel template
if apply_dartel
    darteltpm = warp.darteltpm;
    % find position of '_1_'
    numpos = strfind(darteltpm,'Template_1.nii');
    numpos = numpos+8;
    if isempty(numpos)
        numpos = strfind(darteltpm,'_1_');
    end
    if isempty(numpos)
        error('Could not find _1_ that indicates the first Dartel template in geg_vbm8_defaults.');
    end
    if strcmp(darteltpm(1,end-1:end),',1') >0
        darteltpm = darteltpm(1,1:end-2);
    end
    for j=1:6
        for i=1:2
            run2(i).tpm =  [darteltpm(1:numpos) num2str(j) darteltpm(numpos+2:end) ',' num2str(i)];
        end
        tpm2{j} = spm_vol(strvcat(cat(1,run2(:).tpm)));
    end
end

% prepare transformations for rigidly or affine aligned images
vx = NaN;
bb = nan(2,3);

% Sort out bounding box etc
[bb1,vx1] = bbvox_from_V(tpm.V(1));
bb(~isfinite(bb)) = bb1(~isfinite(bb));
if ~isfinite(vx), vx = abs(prod(vx1))^(1/3); end;
bb(1,:) = vx*round(bb(1,:)/vx);
bb(2,:) = vx*round(bb(2,:)/vx);
if any(tc(:,2)) || any(tc(:,3)) || apply_dartel
    
    % figure out the mapping from the volumes to create to the original
    mm = [[
        bb(1,1) bb(1,2) bb(1,3)
        bb(2,1) bb(1,2) bb(1,3)
        bb(1,1) bb(2,2) bb(1,3)
        bb(2,1) bb(2,2) bb(1,3)
        bb(1,1) bb(1,2) bb(2,3)
        bb(2,1) bb(1,2) bb(2,3)
        bb(1,1) bb(2,2) bb(2,3)
        bb(2,1) bb(2,2) bb(2,3)]'; ones(1,8)];    
    vx2  = M1\mm;
    odim = abs(round((bb(2,1:3)-bb(1,1:3))/vx))+1;
    vx3  = [[
        1       1       1
        odim(1) 1       1
        1       odim(2) 1
        odim(1) odim(2) 1
        1       1       odim(3)
        odim(1) 1       odim(3)
        1       odim(2) odim(3)
        odim(1) odim(2) odim(3)]'; ones(1,8)];
    
    % rigid transformation
    if svdartel_exports
        x      = affind(rgrid(d),M0);
        y1     = affind(y,M1);        
        [~,R]  = spm_get_closest_affine(x,y1,single(cls{1})/255);
        clear x y1
        
        % rigid parameters
        Mr      = M0\inv(R)*M1*vx2/vx3;
        mat0r   =    R\M1*vx2/vx3;
        matr    = mm/vx3;
    end
    
    % affine parameters
    Ma      = M0\inv(res.Affine)*M1*vx2/vx3;
    mat0a   = res.Affine\M1*vx2/vx3;
    mata    = mm/vx3;
end

% dartel spatial normalization to given template
if apply_dartel
    disp('Dartel normalization...');
    % use GM/WM for dartel
    n1 = 2;    
    f = zeros([odim(1:3) 2],'single');
    g = zeros([odim(1:3) 2],'single');
    u = zeros([odim(1:3) 3],'single');
    for k1=1:n1
        for i=1:odim(3),
            f(:,:,i,k1) = single(spm_slice_vol(single(cls{k1}),Ma*spm_matrix([0 0 i]),odim(1:2),[1,NaN])/255);
        end
    end
    
    rform = 0;    % regularization form: 0 - Linear Elastic Energy
    code  = 2;    % multinomial
    lmreg = 0.01; % LM regularization
    cyc = 3;      % cycles
    its = 3;      % relaxation iterations
    for i=1:6, param(i).its = 3; end % inner iterations    
    param(1).rparam = [4 2 1e-6]; % regularization parameters: mu, lambda, id
    param(1).K = 0;               % time steps
    param(2).rparam = [2 1 1e-6];
    param(2).K = 0;
    param(3).rparam = [1 0.5 1e-6];
    param(3).K = 1;
    param(4).rparam = [0.5 0.25 1e-6];
    param(4).K = 2;
    param(5).rparam = [0.25 0.125 1e-6];
    param(5).K = 4;
    param(6).rparam = [0.25 0.125 1e-6];
    param(6).K = 6;   
    
    it0 = 0;
    for it = 1:numel(param)
        prm   = [rform, param(it).rparam, lmreg, cyc, its, param(it).K, code];
        % load new template for this iteration
        for k1=1:n1
            for i=1:odim(3),
                g(:,:,i,k1) = single(spm_slice_vol(tpm2{it}(k1),spm_matrix([0 0 i]),odim(1:2),[1,NaN]));
            end
        end
        for j = 1:param(it).its,
            it0 = it0 + 1;
            [u,ll] = dartel3(u,f,g,prm);
            fprintf('%d \t%g\t%g\t%g\t%g\n',...
                it0,ll(1),ll(2),ll(1)+ll(2),ll(3));
        end
    end    
    [pth,nam,~]=fileparts(res.image(1).fname);    
    y0 = spm_dartel_integrate(reshape(u,[odim(1:3) 1 3]),[0 1], 6);
    clear f g
    [t1,t2] = ndgrid(1:d(1),1:d(2),1);
    t3 = 1:d(3);    
    prm     = [3 3 3 0 0 0];
    Coef    = cell(1,3);
    Coef{1} = spm_bsplinc(y0(:,:,:,1),prm);
    Coef{2} = spm_bsplinc(y0(:,:,:,2),prm);
    Coef{3} = spm_bsplinc(y0(:,:,:,3),prm);
    for z=1:d(3)
        % inverse deformation fields (iy_*)
        [t11,t22,t33] = defs2(Coef,z,Ma,prm,t1,t2,t3);
        y(:,:,z,1) = t11;
        y(:,:,z,2) = t22;
        y(:,:,z,3) = t33;
    end
    clear Coef y0 t1 t2 t3 y1 y2 y3 t11 t22 t33 x1a y1a z1a
end
% % Here ends the DARTEL


% Starting with data ouputs writing
fprintf('Writing data outputs \n');

% write raw segmented images
for k1=1:3,
    if ~isempty(tiss(k1).Nt),
        tiss(k1).Nt.dat(:,:,:,ind(1),ind(2)) = double(cls{k1})/255;
    end
end
clear tiss

% Writting brainmask (binary)
mask2=zeros(size(src)); mask2(src>0)=1;
if write_bmask
    BM.Nt      = nifti;
    BM.Nt.dat  = file_array(fullfile(pth,['bm_', nam, '.nii']),...
        res.image(1).dim(1:3),...
        [spm_type('int16') spm_platform('bigend')],...
        0,1/255,0);
    BM.Nt.mat  = res.image(n).mat;
    BM.Nt.mat0 = res.image(n).mat;
    BM.Nt.descrip = 'Brain mask';
    create(BM.Nt);
    BM.Nt.dat(:,:,:,ind(1),ind(2))=double(mask2);
end

% write skull-stripped images
if any(SS_opt)
    if SS_opt(1)
        % after this I have to open the mask and the T1, then multiply and save the
        % skull-stripped image, and erase the mask
        BSSname=['ss0', nam, '.nii'];
        images = fullfile(pth,[nam, '.nii']);
        Oimage    = spm_vol(images);
        Oima=spm_read_vols(Oimage);
        biasSS=Oima.*mask2;
    elseif SS_opt(2)
        BSSname = ['ss1', nam, '.nii'];
        Bname = fullfile(pth,['m0', nam, '.nii']);
        Bimage = spm_vol(Bname);
        Bima = spm_read_vols(Bimage);
        biasSS=Bima.*mask2;
    end
    BSS.Nt      = nifti;
    BSS.Nt.dat  = file_array(fullfile(pth,BSSname),...
        res.image(1).dim(1:3),...
        [spm_type('int32') spm_platform('bigend')],...
        0,1/255,0);
    BSS.Nt.mat  = res.image(n).mat;
    BSS.Nt.mat0 = res.image(n).mat;
    BSS.Nt.descrip = 'Skull-stripped image';
    create(BSS.Nt);
    BSS.Nt.dat(:,:,:,ind(1),ind(2))=double(biasSS);
end

% save raw tissue class volumes in ml in log-file
if do_cls
    volfactor = abs(det(M0(1:3,1:3)))/1000;
    vol_txt = fullfile(pth,['c', nam1, '_seg8.txt']);
    fid = fopen(vol_txt, 'w');
    for i=1:3
        vol = volfactor*sum(cls{i}(:))/255;
        fprintf(fid,'%5.3f\t',vol);
    end;
    fclose(fid);
end

% native label
if lb(1),
    N      = nifti;
    N.dat  = file_array(fullfile(pth1,['c0', nam, '.nii']),...
        res.image(1).dim(1:3),...
        'float32',0,1,0);
    N.mat  = res.image(1).mat;
    N.mat0 = res.image(1).mat;
    N.descrip = 'tissue-labeled';
    create(N);
    N.dat(:,:,:) = 0;
    N.dat(indx,indy,indz) = double(label)*3/255;
end

% Write DARTEL ...
% write dartel exports
if svdartel_exports,
    % write rigid aligned label
    tmp1 = zeros(res.image(1).dim(1:3),'single');
    tmp1(indx,indy,indz) = double(label)*3;
    VT      = struct('fname',fullfile(pth,['rp0', nam, '.nii']),...
        'dim',  odim,...
        'dt',   [spm_type('int16') spm_platform('bigend')],...
        'pinfo',[1/255 0]',...
        'mat',matr);
    VT = spm_create_vol(VT);
    
    Ni             = nifti(VT.fname);
    Ni.mat0        = mat0r;
    Ni.mat_intent  = 'Aligned';
    Ni.mat0_intent = 'Aligned';
    create(Ni);
    
    for i=1:odim(3),
        tmp = spm_slice_vol(tmp1,Mr*spm_matrix([0 0 i]),odim(1:2),[1,NaN])/255;
        VT  = spm_write_plane(VT,tmp,i);
    end
    clear tmp1

    % write rigid aligned tissues    
    for k1=1:size(tc,1),        
        if svdartel_exports && tc(k1,1)~=0
            VT      = struct('fname',fullfile(pth,['rp', num2str(k1), nam, '.nii']),...
                'dim',  odim,...
                'dt',   [spm_type('int16') spm_platform('bigend')],...
                'pinfo',[1/255 0]',...
                'mat',matr);
            VT = spm_create_vol(VT);
            
            Ni             = nifti(VT.fname);
            Ni.mat0        = mat0r;
            Ni.mat_intent  = 'Aligned';
            Ni.mat0_intent = 'Aligned';
            create(Ni);
            
            for i=1:odim(3),
                tmp = spm_slice_vol(single(cls{k1}),Mr*spm_matrix([0 0 i]),odim(1:2),[1,NaN])/255;
                VT  = spm_write_plane(VT,tmp,i);
            end
        end
        
        % write affine normalized tissues 
        % only usefull when creating templates
        if tc(k1,3),
            VT      = struct('fname',fullfile(pth,['rp', num2str(k1), nam, '_affine.nii']),...
                'dim',  odim,...
                'dt',   [spm_type('int16') spm_platform('bigend')],...
                'pinfo',[1/255 0]',...
                'mat',mata);
            VT = spm_create_vol(VT);
            
            Ni             = nifti(VT.fname);
            % get rid of the QFORM0 rounding warning
            warning('off','')
            Ni.mat0        = mat0a;
            warning('on','')
            Ni.mat_intent  = 'Aligned';
            Ni.mat0_intent = 'Aligned';
            create(Ni);
            
            for i=1:odim(3),
                tmp = spm_slice_vol(single(cls{k1}),Ma*spm_matrix([0 0 i]),odim(1:2),[1,NaN])/255;
                VT  = spm_write_plane(VT,tmp,i);
            end
        end
    end
end

% write jacobian determinant
if svjc
    if ~apply_dartel
        warning('Jacobian can only be saved if dartel normalization was used.');
    else
        [~, dt] = spm_dartel_integrate(reshape(u,[odim(1:3) 1 3]),[1 0], 6);
        N      = nifti;
        N.dat  = file_array(fullfile(pth,['jac_wrp1', nam, '.nii']),...
            d1,...
            [spm_type('float32') spm_platform('bigend')],...
            0,1,0);
        N.mat  = M1;
        N.mat0 = M1;
        N.descrip = 'Jacobian';
        create(N);
        
        N.dat(:,:,:) = dt;
    end
end

if any(tc(:,4)),
    C = zeros([d1,3],'single');
end

clear u

% warped tissue classes
if save_warped,    
    spm_progress_bar('init',3,'Warped Tissue Classes','Classes completed');
    
    % estimate a more accurate jacobian determinant
    y2 = spm_invert_def(y,M1,d1,M0,[1 0]);
    dt = spm_def2det(y2(:,:,:,1),y2(:,:,:,2),y2(:,:,:,3),M1);
    dt = dt./abs(det(M0(1:3,1:3))/det(M1(1:3,1:3)));
    clear y2
    
    M2 = M1\res.Affine*M0;
    
    for k1 = 1:3,
        if ~isempty(cls{k1}),
            c = single(cls{k1})/255;
            [c,w]  = dartel3('push',c,y,d1(1:3));
            C(:,:,:,k1) = c./(w+eps);
            clear w
            c = C(:,:,:,k1).*dt;
            if nargout>=1,
                cls{k1} = c;
            end
            
            if save_warped == 3,
                N      = nifti;
                N.dat  = file_array(fullfile(pth,['mwp', num2str(k1), nam, '.nii']),...
                    d1,...
                    [spm_type('float32') spm_platform('bigend')],...
                    0,1,0);
                if apply_dartel
                    N.dat.fname = fullfile(pth,['mwrp', num2str(k1), nam, '.nii']);
                end
                N.mat  = M1;
                N.mat0 = M1;
                N.descrip = ['Jac. sc. warped tissue class ' num2str(k1)];
                create(N);
                N.dat(:,:,:) = c*abs(det(M0(1:3,1:3))/det(M1(1:3,1:3)));
            end
            if save_warped == 2,
                N      = nifti;
                N.dat  = file_array(fullfile(pth,['m0wp', num2str(k1), nam, '.nii']),...
                    d1,...
                    [spm_type('float32') spm_platform('bigend')],...
                    0,1,0);
                if apply_dartel
                    N.dat.fname = fullfile(pth,['m0wrp', num2str(k1), nam, '.nii']);
                end
                N.mat  = M1;
                N.mat0 = M1;
                N.descrip = ['Jac. sc. warped tissue class non-lin only' num2str(k1)];
                create(N);
                
                N.dat(:,:,:) = c*abs(det(M2(1:3,1:3)));
            end
            spm_progress_bar('set',k1);
        end
    end    
    spm_progress_bar('Clear');
end

% warped label
if save_warped == 2,
    c = zeros(res.image(n).dim(1:3),'single');
    c(indx,indy,indz) = single(label);
    [c,w]  = dartel3('push',c,y,d1(1:3));
    C = optimNn(w,c,[1  vx vx vx 1e-4 1e-6 0  3 2]);
    clear w
    N      = nifti;
    N.dat  = file_array(fullfile(pth,['wp0', nam, '.nii']),...
                            d1,'float32',0,1,0);
    if apply_dartel
        N.dat.fname = fullfile(pth,['wrp0', nam, '.nii']);
    end
    N.mat  = M1;
    N.mat0 = M1;
    N.descrip = 'Warped bias corrected image ';
    create(N);
    N.dat(:,:,:) = double(C)*3/255;
end
% End of DARTEL outputs

if nargout == 0
    clear cls
end

% deformations
if df(1),
    y         = spm_invert_def(y,M1,d1,M0,[1 0]);
    N         = nifti;
    N.dat     = file_array(fullfile(pth,['y_', nam1, '.nii']),...
        [d1,1,3],'float32',0,1,0);
     if do_dartel
        N.dat.fname = fullfile(pth,['y_r', nam1, '.nii']);
     end
    N.mat     = M1;
    N.mat0    = M1;
    N.descrip = 'Deformation';
    create(N);
    N.dat(:,:,:,:,:) = reshape(y,[d1,1,3]);
end


return;
%=======================================================================
%=======================================================================
function [x1,y1,z1] = defs(sol,z,MT,prm,x0,y0,z0,M)
iMT = inv(MT);
x1  = x0*iMT(1,1)+iMT(1,4);
y1  = y0*iMT(2,2)+iMT(2,4);
z1  = (z0(z)*iMT(3,3)+iMT(3,4))*ones(size(x1));
x1a = x0    + spm_bsplins(sol{1},x1,y1,z1,prm);
y1a = y0    + spm_bsplins(sol{2},x1,y1,z1,prm);
z1a = z0(z) + spm_bsplins(sol{3},x1,y1,z1,prm);
x1  = M(1,1)*x1a + M(1,2)*y1a + M(1,3)*z1a + M(1,4);
y1  = M(2,1)*x1a + M(2,2)*y1a + M(2,3)*z1a + M(2,4);
z1  = M(3,1)*x1a + M(3,2)*y1a + M(3,3)*z1a + M(3,4);
return;
%=======================================================================
%=======================================================================
function t = transf1(B1,B2,B3,T)
if ~isempty(T)
    d2 = [size(T) 1];
    t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
    t  = B1*t1*B2';
else
    t = zeros(size(B1,1),size(B2,1),size(B3,1));
end;
return;
%=======================================================================
%=======================================================================
function p = likelihoods(f,bf,mg,mn,vr)
K  = numel(mg);
N  = numel(f);
M  = numel(f{1});
cr = zeros(M,N);
for n=1:N,
    if isempty(bf),
        cr(:,n) = double(f{n}(:));
    else
        cr(:,n) = double(f{n}(:).*bf{n}(:));
    end
end
p  = ones(numel(f{1}),K);
for k=1:K,
    amp    = mg(k)/sqrt((2*pi)^N * det(vr(:,:,k)));
    d      = cr - repmat(mn(:,k)',M,1);
    p(:,k) = amp * exp(-0.5* sum(d.*(d/vr(:,:,k)),2));
end
%=======================================================================
%=======================================================================
function y = mask_largest_cluster(y, th)

if nargin < 2
    th = 0.5;
end

sz = size(y);

th = th*max(single(y(:)));

mask = y > th;
Q = find(mask);

Qth = find(y <= th & y>0);
yth = y(Qth);

% save memory by using bounding box where y > th
[indx, indy, indz] = ind2sub(size(mask),Q);
indx = max((min(indx) - 1),1):min((max(indx) + 1),sz(1));
indy = max((min(indy) - 1),1):min((max(indy) + 1),sz(2));
indz = max((min(indz) - 1),1):min((max(indz) + 1),sz(3));

[A,~] = spm_bwlabel(double(mask(indx,indy,indz)),26);

clear mask

if isempty(A)
    error('No cluster found!');
end

% interrupt if cluster was > 7.5% of whole image to save time
max_A = max(A(:));
sz_cluster = zeros(max_A,1);
for i=1:max_A
    QA = find(A == i);
    ind = i;
    if length(QA)/numel(A) > 0.075 
        break
    end
    sz_cluster(i) = length(QA);
end

if length(QA)/numel(A) <= 0.075 
    [~, ind] = max(sz_cluster);
    QA = find(A == ind);
end

indices=(1:numel(A))';
QA0 = indices(A ~= ind);
A = y(indx,indy,indz);
A(QA0) = 0;

y(indx,indy,indz) = A;
y(Qth) = yth;

return;
%=======================================================================
%=======================================================================
function [bb,vx] = bbvox_from_V(V)
vx = sqrt(sum(V(1).mat(1:3,1:3).^2));
if det(V(1).mat(1:3,1:3))<0, vx(1) = -vx(1); end;

o  = V(1).mat\[0 0 0 1]';
o  = o(1:3)';
bb = [-vx.*(o-1) ; vx.*(V(1).dim(1:3)-o)];
return;
%=======================================================================
%=======================================================================
function [x1,y1,z1] = defs2(sol,z,M,prm,x0,y0,z0)
iM = inv(M);
z01 = z0(z)*ones(size(x0));

x1a  = iM(1,1)*x0 + iM(1,2)*y0 + iM(1,3)*z01 + iM(1,4);
y1a  = iM(2,1)*x0 + iM(2,2)*y0 + iM(2,3)*z01 + iM(2,4);
z1a  = iM(3,1)*x0 + iM(3,2)*y0 + iM(3,3)*z01 + iM(3,4);

x1 = spm_bsplins(sol{1},x1a,y1a,z1a,prm);
y1 = spm_bsplins(sol{2},x1a,y1a,z1a,prm);
z1 = spm_bsplins(sol{3},x1a,y1a,z1a,prm);

return;
%=======================================================================
%=======================================================================
function y1 = affind(y0,M)
y1 = zeros(size(y0),'single');
for d=1:3,
    y1(:,:,:,d) = y0(:,:,:,1)*M(d,1);
    y1(:,:,:,d) = y1(:,:,:,d) + y0(:,:,:,2)*M(d,2);
    y1(:,:,:,d) = y1(:,:,:,d) + y0(:,:,:,3)*M(d,3) + M(d,4);
end
%=======================================================================
%=======================================================================
function x = rgrid(d)
x = zeros([d(1:3) 3],'single');
[x1,x2] = ndgrid(single(1:d(1)),single(1:d(2)));
for i=1:d(3),
    x(:,:,i,1) = x1;
    x(:,:,i,2) = x2;
    x(:,:,i,3) = single(i);
end
%=======================================================================