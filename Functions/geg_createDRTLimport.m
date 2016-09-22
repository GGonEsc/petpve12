function geg_createDRTLimport(job)
%
%_______________________________________________________________________
% Gabriel Gonzalez-Escamilla
% $Id: geg_createDRTLimport.m 001 2015-03-28 16:20:04Z $

%rev = '$Rev: 003 $';

S=job.resF;

for subj=1:length(S)
    actres=S{subj};
    [pth,nam,ext] = fileparts(actres);
    if exist(fullfile(pth,[nam,ext]),'file')
        res       = load(fullfile(pth,[nam ext]));
        geg_DRTLimp_write(res,pth)        
    end
    fprintf('\n');
end

function geg_DRTLimp_write(res,newpth)

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

do_defs=1;
do_cls=1;
tc=[0 1 0; 0 1 0; 0 1 0];
lb = [0 0 1 0];
bf = [0 0 0];

% Read essentials from tpm (it will be cleared later)
tpm = res.tpm;
if ~isstruct(tpm) || (~isfield(tpm, 'bg1') && ~isfield(tpm, 'bg')),
    tpm = spm_load_priors8(tpm);
end

if isfield(res,'mg'),
    lkp = res.lkp;
    Kb  = max(lkp);
else
    Kb  = size(res.intensity(1).lik,2);
end
N   = numel(res.image);

M1        = tpm.M;
[bb1,vx1] = bbvox_from_V(tpm.V(1));

vx = NaN;
bb = nan(2,3);

% Sort out bounding box etc
bb(~isfinite(bb)) = bb1(~isfinite(bb));
if ~isfinite(vx), vx = abs(prod(vx1))^(1/3); end;
bb(1,:) = vx*round(bb(1,:)/vx);
bb(2,:) = vx*round(bb(2,:)/vx);

d    = res.image(1).dim(1:3);

[x1,x2,o] = ndgrid(1:d(1),1:d(2),1);
x3  = 1:d(3);

chan(N) = struct('B1',[],'B2',[],'B3',[],'T',[],'Nc',[],'Nf',[],'ind',[]);
for n=1:N,
    d3         = [size(res.Tbias{n}) 1];
    chan(n).B3 = spm_dctmtx(d(3),d3(3),x3);
    chan(n).B2 = spm_dctmtx(d(2),d3(2),x2(1,:)');
    chan(n).B1 = spm_dctmtx(d(1),d3(1),x1(:,1));
    chan(n).T  = res.Tbias{n};    
    chan(n).ind      = res.image(n).n;
end

cls      = cell(1,Kb);
for k1=1:Kb,
    cls{k1} = zeros(d(1:3),'uint8');
end

prm     = [3 3 3 0 0 0];
Coef    = cell(1,3);
Coef{1} = spm_bsplinc(res.Twarp(:,:,:,1),prm);
Coef{2} = spm_bsplinc(res.Twarp(:,:,:,2),prm);
Coef{3} = spm_bsplinc(res.Twarp(:,:,:,3),prm);

M = tpm.M\res.Affine*res.image(1).mat;

warp.open_th = 0.25; % initial threshold for skull-stripping
warp.dilate = 1; % number of final dilations for skull-stripping
warp.cleanup = geg_petpve12_get_defaults('extopts.cleanup');

histeq_deep = 0;
try
    histeq_deep = geg_petpve12_get_defaults('extopts.histeq_deep');
catch
end
if histeq_deep
    tmp_histeq_mask = spm_vol(char(geg_petpve12_get_defaults('extopts.histeq_mask')));
    histeq_mask = zeros(d(1:3),'uint8');
    M2 = tmp_histeq_mask.mat\res.Affine*res.image(1).mat;
end

for z=1:length(x3),
    
    for n=1:N,
        f = spm_sample_vol(res.image(n),x1,x2,o*x3(z),0);
        bf1 = exp(transf(chan(n).B1,chan(n).B2,chan(n).B3(z,:),chan(n).T));
        bf1(bf1>100) = 100;
        cr{n} = bf1.*f;
    end
    
    if do_defs,
        [t1,t2,t3] = defs(Coef,z,res.MT,prm,x1,x2,x3,M);
        if exist('Ndef','var'), % This the iy_ image
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

            if histeq_deep
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
        if bf(1,2) || bf(1,3) || any(lb([2,3,4])) || df(1) || any(any(tc(:,[2,3,4,5,6]))) || nargout>=1,
            % initialize y only at first slice
            if z==1
                y = zeros([res.image(1).dim(1:3),3],'single');
            end
            y(:,:,z,1) = t1;
            y(:,:,z,2) = t2;
            y(:,:,z,3) = t3;
        end
    end
end

% load bias corrected image
src = zeros(res.image(1).dim(1:3),'single');
for z=1:length(x3),
    f = spm_sample_vol(res.image(1),x1,x2,o*x3(z),0);
    bf1 = exp(transf(chan(1).B1,chan(1).B2,chan(1).B3(z,:),chan(1).T));
    % restrict bias field to maximum of 100 
    % (sometimes artefacts at the borders can cause huge values in bias field)
    bf1(bf1>100) = 100;
    src(:,:,z) = single(bf1.*f);
end
% prevent NaN
src(isnan(src)) = 0;
if do_cls && do_defs,
    
    % default parameters
    bias_fwhm   = geg_petpve12_get_defaults('extopts.bias_fwhm');
    init_kmeans = geg_petpve12_get_defaults('extopts.kmeans');
    finalmask   = geg_petpve12_get_defaults('extopts.finalmask');
    gcut        = geg_petpve12_get_defaults('extopts.gcut');
    
    vx_vol = sqrt(sum(res.image(1).mat(1:3,1:3).^2));
    scale_morph = 1/mean(vx_vol);
    
    if gcut
        % skull-stripping using graph-cut
        opt.verb = 0; % display process (0=nothing, 1=only points, 2=times)
        fprintf('Skull-stripping using graph-cut\n');
        cls_old = cls;
        try
          [src,cls,mask] = SSGC(src,cls,res,opt);
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
    
    % % New
    vol = double(src(indx,indy,indz));    
    
    % mask source image because Amap needs a skull stripped image
    % set label and source inside outside mask to 0
    vol(mask(indx,indy,indz)==0) = 0;
    
    % use local histogram equalization
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

        for z=1:length(x3),
            % Bias corrected image
            % Write a plane of bias corrected data
            if bf(1,1),
                chan(1).Nc.dat(:,:,z,chan(1).ind(1),chan(1).ind(2)) = src(:,:,z);
            end
        end

    end

    clear chan
    % Amap parameters
    n_iters = 200; sub = 16; n_classes = 3; pve = 5; iters_icm = 20;
    warp.mrf = geg_petpve12_get_defaults('extopts.mrf');
    
    if init_kmeans, fprintf('Amap with Kmeans\n');   
    else            fprintf('Amap without Kmeans\n');   
    end

    [prob, ~] = AmapMex(vol, label, n_classes, n_iters, sub, pve, init_kmeans, warp.mrf, vx_vol, iters_icm, bias_fwhm);%[prob, means]
    
    % reorder probability maps according to spm order
    prob = prob(:,:,:,[2 3 1]);
    clear vol 
    
    % use cleanup
    if warp.cleanup
        % get sure that all regions outside mask are zero
        for i=1:3
            cls{i}(:) = 0;
        end
       % disp('Clean up...');        
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
        fprintf('Final masking\n');
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

% prepare transformations for rigidly or affine aligned images
M0 = res.image(1).mat;
if any(tc(:,2)) || any(tc(:,3)) || lb(1,3) || lb(1,4) || bf(1,3)

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
    if (any(tc(:,2)) || lb(1,3))
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

% % Write outputs
% write affine label
if lb(1,4),
    tmp1 = zeros(res.image(1).dim(1:3),'single');
    tmp1(indx,indy,indz) = single(label)*3;
    [~,nam,~]=fileparts(res.image(1).fname);
    VT      = struct('fname',fullfile(newpth,['rp0', nam, '_affine.nii']),...
            'dim',  odim,...
            'dt',   [spm_type('int16') spm_platform('bigend')],...
            'pinfo',[1/255 0]',...
            'mat',mata);
    VT = spm_create_vol(VT);

    N             = nifti(VT.fname);
    % get rid of the QFORM0 rounding warning
    warning('off',' ')
    N.mat0        = mat0a;
    warning('on',' ')
    N.mat_intent  = 'Aligned';
    N.mat0_intent = 'Aligned';
    create(N);

    for i=1:odim(3),
        tmp = spm_slice_vol(tmp1,Ma*spm_matrix([0 0 i]),odim(1:2),[1,NaN])/255;
        VT  = spm_write_plane(VT,tmp,i);
    end
    clear tmp1
end

% write rigid aligned label
if lb(1,3),
    tmp1 = zeros(res.image(1).dim(1:3),'single');
    tmp1(indx,indy,indz) = double(label)*3;
    [~,nam,~]=fileparts(res.image(1).fname);
    VT      = struct('fname',fullfile(newpth,['rc0', nam, '.nii']),...
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
end


% write dartel exports
for k1=1:size(tc,1),

    % write rigid aligned tissues
    if tc(k1,2),
        [~,nam,~]=fileparts(res.image(1).fname);
        VT      = struct('fname',fullfile(newpth,['rc', num2str(k1), nam, '.nii']),...
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
    if tc(k1,3),
        [~,nam,~]=fileparts(res.image(1).fname);
        VT      = struct('fname',fullfile(newpth,['rc', num2str(k1), nam, '_affine.nii']),...
                'dim',  odim,...
                'dt',   [spm_type('int16') spm_platform('bigend')],...
                'pinfo',[1/255 0]',...
                'mat',mata);
        VT = spm_create_vol(VT);

        Ni             = nifti(VT.fname);
        % get rid of the QFORM0 rounding warning
        warning('off',' ')
        Ni.mat0        = mat0a;
        warning('on',' ')
        Ni.mat_intent  = 'Aligned';
        Ni.mat0_intent = 'Aligned';
        create(Ni);

        for i=1:odim(3),
            tmp = spm_slice_vol(single(cls{k1}),Ma*spm_matrix([0 0 i]),odim(1:2),[1,NaN])/255;
            VT  = spm_write_plane(VT,tmp,i);
        end
    end
end

%=======================================================================
%=======================================================================
% % SUB-FUNCITONS
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
function t = transf(B1,B2,B3,T)
if ~isempty(T)
    d2 = [size(T) 1];
    t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
    t  = B1*t1*B2';
else
    t = zeros(size(B1,1),size(B2,1),size(B3,1));
end;
return;
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
[A,~] = spm_bwlabel(double(mask(indx,indy,indz)),26);%[A,num]
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
	if length(QA)/numel(A) > 0.075 %prod(size(A))
		break
	end
	sz_cluster(i) = length(QA);
end
if length(QA)/numel(A) <= 0.075 %prod(size(A))
	[~, ind] = max(sz_cluster);%[mx, ind]
	QA = find(A == ind);
end
QA0 = find(A ~= ind);% not changed: A ~= ind;
A = y(indx,indy,indz);
A(QA0) = 0;
y(indx,indy,indz) = A;
y(Qth) = yth;
return;
%=======================================================================
function y1 = affind(y0,M)
y1 = zeros(size(y0),'single');
for d=1:3,
    y1(:,:,:,d) = y0(:,:,:,1)*M(d,1);
    y1(:,:,:,d) = y1(:,:,:,d) + y0(:,:,:,2)*M(d,2);
    y1(:,:,:,d) = y1(:,:,:,d) + y0(:,:,:,3)*M(d,3) + M(d,4);
end
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