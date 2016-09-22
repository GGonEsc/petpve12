function geg_SSusingSegments(job)
%
% Inputs:
% res matrix = fullfile(pth,[nam '_seg8.mat']),'file')
%
%
% Ask for the T1 or bias image, the tissue images (c1, c2,
% c3) and the res file
% If no res file as input, do michel procedure:
% skull-strip automatico: el MRI lo multiplico con una mascara (c1 + c2) > 0.5
%
% Is added to the batch in the spatial tools as skull-stripping
%_______________________________________________________________________
% Copyright (C) 2015
%
% Gabriel Gonzalez-Escamilla
% $Id: geg_SSusingSegments.m 015 2015-03-20 20:13:03Z $

do_morph=job.SStype;

% Check that for every MRI has a GM and WM segment
if numel(job.data) ~= (numel(job.tiss1)) || numel(job.data) ~= numel(job.tiss2)
    error('Not equal number of inputs: %s MRI, %s GM and %s WM images', num2str(numel(job.data)),num2str(numel(job.tiss1)),num2str(length(job.tiss2)))
end
fprintf('Preparing data for Skull-Stripping\n');
for subj=1:length(job.data)
    
    % selecting subject tissue data
    tiss1n = char(job.tiss1{subj});
    %tiss1=fullfile(pth,['c1', nam, '.nii']); %tiss1n = char(tiss1);
    tiss1image    = spm_vol(tiss1n);
    tiss1ima=spm_read_vols(tiss1image);
    
    tiss2n = char(job.tiss2{subj});
    %tiss2=fullfile(pth,['c2', nam, '.nii']); %tiss2n = char(tiss2);
    tiss2image    = spm_vol(tiss2n);
    tiss2ima=spm_read_vols(tiss2image);
    
    % construct cls with the three tissues
    cls={tiss1ima tiss2ima};
    geg_makeBmask(cls,do_morph,subj,job)
    fprintf('\n');
end
fprintf('Skull-Stripping finished\n')

function geg_makeBmask(cls,do_morph,subj,job)

[pth,nam,~,~] = spm_fileparts(job.data{subj});
fprintf('Working with subject: %s \n', nam);
image = spm_vol(job.data{subj});
ind  = image(1).n;
src = spm_read_vols(image);

cleanup=job.SScleanup;% 0=no, 1= light, 2=thorough
write_bmask=job.SSoutput;

if do_morph
    lkp = [ones(1,2) ones(1,2)*2];
    Kb  = max(lkp);
    d    = image(1).dim(1:3);
    open_th = 0.25; % initial threshold for skull-stripping 
    dilate = 1; % number of final dilations for skull-stripping
    %[pth,nam] = fileparts(image(1).fname);
    
    % Defining what to write    
    do_finalmask = job.SStypeOpt.SStypeO;%geg_petpve12_get_defaults('extopts.finalmask');
    SS_opt=1;
    
    % The skull-stripping    
    vx_vol = sqrt(sum(image(1).mat(1:3,1:3).^2));
    scale_morph = 1/mean(vx_vol);
    
    fprintf('Skull-stripping using morphological operations\n');
    % use mask of GM and WM
    mask = single(cls{1});
    mask = mask + single(cls{2});
    
    % keep largest connected component after at least 1 iteration of opening
    n_initial_openings = max(1,round(scale_morph*cleanup));
    mask = cg_morph_vol(mask,'open',n_initial_openings,open_th);
    mask = mask_largest_cluster(mask,0.5);
    
    % dilate and close to fill ventricles
    mask = cg_morph_vol(mask,'dilate',dilate,0.5);
    mask = cg_morph_vol(mask,'close',round(scale_morph*10),0.5);
    
    %{
    % remove sinus (only possible if have the CSF and skull tissues)
    mask = mask & ((single(cls{5})<single(cls{1})) | ...
    (single(cls{5})<single(cls{2})) | ...
    (single(cls{5})<single(cls{3})));
    %}
    
    % fill holes that may remain
    mask = cg_morph_vol(mask,'close',round(scale_morph*2),0.5);
        
    % calculate label image for all classes
    cls2 = zeros([d(1:2) Kb]);
    label2 = zeros(d,'uint8');
    for i=1:d(3)
        for k1 = 1:Kb
            cls2(:,:,k1)  = cls{k1}(:,:,i);
        end
        % find maximum for reordered segmentations
        [maxi,maxind] = max(cls2(:,:,[1,2,4:Kb]),[],3);
        for k1 = 1:Kb
            label2(:,:,i) = label2(:,:,i) + uint8((maxind == k1).*(maxi~=0)*k1);
        end
    end
    
    % set all non-brain tissue outside mask to 0
    label2(mask == 0)  = 0;
    
    % and for skull/bkg tissue classes to 0
    label2(label2 > 3) = 0;
        
    % fill remaining holes in label with 1
    mask = cg_morph_vol(label2,'close',round(scale_morph*6),0);% 2 lefts the ventricles unfilled
    label2((label2 == 0) & (mask > 0)) = 1;
    
    % use index to speed up and save memory (This was not in here)
    sz = size(mask);
    [indx, indy, indz] = ind2sub(sz,find(mask>0));
    indx = max((min(indx) - 1),1):min((max(indx) + 1),sz(1));
    indy = max((min(indy) - 1),1):min((max(indy) + 1),sz(2));
    indz = max((min(indz) - 1),1):min((max(indz) + 1),sz(3));
    label = label2(indx,indy,indz);
    % I normally pass from here to apply the final mask and then save
    
    %{
    % % This if for aplying the Amap and the cleanup 
    % (It does not apply local histogram equalization)
    % It is actually not suitable, because I have only two tissue segments
    if do_amap && length(cls)==3
        bias_fwhm   = geg_petpve12_get_defaults('extopts.bias_fwhm');    
        mrf = geg_petpve12_get_defaults('extopts.mrf');
        
        vol = double(src(indx,indy,indz));
        
        % mask source image because Amap needs a skull stripped image
        % set label and source inside outside mask to 0
        vol(mask(indx,indy,indz)==0) = 0;
        
        % Adaptive maximum a posteriori estimations (Amap)(Gaser, 2009)
        % Amap parameters
        n_iters = 200; sub = 16; n_classes = 3; pve = 5; iters_icm = 20;
        if init_kmeans, fprintf('Adaptive Maximum A Posteriori (AMAP) Approach with Kmeans\n');
        else            fprintf('Adaptive Maximum A Posteriori (AMAP) without Kmeans\n');
        end
        [prob, ~] = AmapMex(vol, label, n_classes, n_iters, sub, pve, init_kmeans, mrf, vx_vol, iters_icm, bias_fwhm);%[prob, means]
        
        % reorder probability maps according to spm order
        % For Tohka Partial volume estimation in brain MRI the segments are:
        % csf = label 1; gm = label 2; wm = label 3
        prob = prob(:,:,:,[2 3 1]);
        clear vol
        
        % use cleanup
        if do_cleanup
            % get sure that all regions outside mask are zero
            for i=1:2
                cls{i}(:) = 0;
            end
            fprintf('Cleanning up... \n');
            [cls{1}(indx,indy,indz), cls{2}(indx,indy,indz), cls{3}(indx,indy,indz)] = cg_cleanup_gwc(prob(:,:,:,1), ...
                prob(:,:,:,2), prob(:,:,:,3), cleanup);
            sum_cls = cls{1}(indx,indy,indz)+cls{2}(indx,indy,indz)+cls{3}(indx,indy,indz);
            label(sum_cls<0.15*255) = 0;
        else
            for i=1:3
                cls{i}(:) = 0;
                cls{i}(indx,indy,indz) = prob(:,:,:,i);
            end
        end;
        clear prob        
    end
    %}
    
    if do_finalmask
        fprintf('Applying final masking\n');
        % create final mask
        mask2 = single(cls{1});
        mask2 = mask2 + single(cls{2});
        
        % keep largest connected component after at least 1 iteration of opening
        n_initial_openings = max(1,round(scale_morph*2));
        mask2 = cg_morph_vol(mask2,'open',n_initial_openings,0.5);
        mask2 = mask_largest_cluster(mask2,0.5);
        
        % dilate and close to fill ventricles
        mask2 = cg_morph_vol(mask2,'dilate',2,0.5);
        mask2 = cg_morph_vol(mask2,'close',20,0.5);
        
        ind_mask = find(mask2 == 0);
        for i=1:2
            cls{i}(ind_mask) = 0;
        end
        
        % mask label
        label2 = zeros(d,'uint8');
        label2(indx,indy,indz) = label;
        label2(ind_mask) = 0;
        label = label2(indx,indy,indz);
        clear label2        
    end
    
    % Starting with data ouputs writing
    fprintf('Writing data outputs \n');
    
    % Writting brainmask (binary)
    if write_bmask
        BM.Nt      = nifti;
        BM.Nt.dat  = file_array(fullfile(pth,['bm_', nam, '.nii']),...
            image(1).dim(1:3),...
            [spm_type('int16') spm_platform('bigend')],...
            0,1/255,0);
        BM.Nt.mat  = image(1).mat;
        BM.Nt.mat0 = image(1).mat;
        BM.Nt.descrip = 'Brain mask';
        create(BM.Nt);
        if do_finalmask
            BM.Nt.dat(:,:,:,ind(1),ind(2))=double(mask2);
        else
            BM.Nt.dat(:,:,:,ind(1),ind(2))=double(mask);%double(label2)/255
        end
    end
    
    % write skull-stripped images
    if any(SS_opt)
        SSname = ['ss', nam, '.nii'];
        %check that mask is binary isbin=find(mask~=0 & mask~=1);
        if do_finalmask
            ISS=src.*mask2;
        else
            ISS=src.*double(mask);
        end
        BSS.Nt      = nifti;
        BSS.Nt.dat  = file_array(fullfile(pth,SSname),...
            image(1).dim(1:3),...
            [spm_type('int32') spm_platform('bigend')],...
            0,1/255,0);
        BSS.Nt.mat  = image(1).mat;
        BSS.Nt.mat0 = image(1).mat;
        BSS.Nt.descrip = 'Skull-stripped image';
        create(BSS.Nt);
        BSS.Nt.dat(:,:,:,ind(1),ind(2))=double(ISS);%double(label2)/255
    end
    
else
    
    fprintf('Skull-stripping using thresholding operations\n');    
    
    fprintf('Aplying threshold %s \n',num2str(job.SS_Thresh));
    Bmaski=cls{1}+cls{2};
    ThrMask=(Bmaski) > job.SS_Thresh;
    
    fprintf('Filling ventricles\n');
    % These are some operations similar to the ones above, but simpler
    vx_vol = sqrt(sum(image(1).mat(1:3,1:3).^2));
    scale_morph = 1/mean(vx_vol);
    n_initial_openings = max(1,round(scale_morph*cleanup));    
    Bmaskb=zeros(size(ThrMask));
    
    % Binarize the mask before filling it
    Bmaskb(ThrMask>0)=1; 
    
    % Fill the ventricles
    mask = cg_morph_vol(Bmaskb,'open',n_initial_openings,0.5);
    mask = mask_largest_cluster(mask,0.5);
    mask = cg_morph_vol(mask,'close',20,0.5);        
        
    fprintf('Writing data outputs \n');
    % Writing outputs
    if write_bmask
        SSname = ['bm', nam, '.nii'];
        BM.Nt      = nifti;
        BM.Nt.dat  = file_array(fullfile(pth,SSname),...
            image(1).dim(1:3),...
            [spm_type('int16') spm_platform('bigend')],...
            0,1/255,0);
        BM.Nt.mat  = image(1).mat;
        BM.Nt.mat0 = image(1).mat;
        BM.Nt.descrip = 'Skull-stripped image';
        create(BM.Nt);
        BM.Nt.dat(:,:,:,ind(1),ind(2))=double(mask);
    end
    
    SSname = ['ss', nam, '.nii'];
    %check that mask is binary isbin=find(mask~=0 & mask~=1);
    ISS=src.*mask;
    BSS1.Nt      = nifti;
    BSS1.Nt.dat  = file_array(fullfile(pth,SSname),...
        image(1).dim(1:3),...
        [spm_type('int32') spm_platform('bigend')],...
        0,1/255,0);
    BSS1.Nt.mat  = image(1).mat;
    BSS1.Nt.mat0 = image(1).mat;
    BSS1.Nt.descrip = 'Skull-stripped image';
    create(BSS1.Nt);
    BSS1.Nt.dat(:,:,:,ind(1),ind(2))=double(ISS);
    
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
    %return;
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
indices=(1:numel(A))';
QA0 = indices(A ~= ind);%find(A ~= ind); %
A = y(indx,indy,indz);
A(QA0) = 0;
y(indx,indy,indz) = A;
y(Qth) = yth;

return;
