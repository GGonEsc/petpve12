function geg_calc_singleROImean_native()
% This function converts a standard mask into native space and computes the
% corresponding mean signal.
% This usefull to extract values on the uncorrected PET data, so intensity
% normalization can be performed on both, uncorrected and PVEc PETs, for
% comparison.

% data
mask  = spm_select([1,1],'image','Select binary mask');
S  = spm_select(Inf,'image','Select native PET images');

% 'Atlas in standard space'
djob.out{1}.pull = struct('fnames',{{}},'savedir',{{}},'interp',{{}},'mask',{{}},'fwhm',{{}});
djob.out{1}.pull.mask   = 0;
djob.out{1}.pull.interp = 0;% it should be always zero
djob.out{1}.pull.fwhm   = [0,0,0];
deftype = spm_input('Input Fields type','1','b',{'Deformations','FlowFields'},{'defs','dartel'});
deftype = deftype{1};
switch deftype
    case 'defs'
        InDefs  = spm_select(size(S,1),'^iy_.*\.nii','Select deformation-fields for every subject');
        if isempty(InDefs)
            fprintf(2,'Operation terminated. No deformation matrices selected \n');
            return
        end
    case 'dartel'
        Inflowfields = spm_select(size(S,1),'^u_.*\.nii','Select flow-fields for every subject');
        if isempty(Inflowfields)
            fprintf(2,'Operation terminated. No deformation matrices selected \n');
            return
        end
        DRTLtemplate = {''};
        DRTLsteps    = {6};
end

disp(mask)
signal = zeros(size(S,1),1);
average = zeros(size(S,1),1);

%-Start progress plot
spm_progress_bar('Init',size(S,1),'Deform + Computing mean','Subjects completed');

% Running for every subject
for subj=1:size(S,1)
    % Loading data
    actPT = char(S(subj,:));% PET name
    [odir,~,~] = spm_fileparts(actPT);
    PTm = spm_vol(actPT);% PET header
    fprintf('%d,',subj)
    
    % 'Parcellation atlas in standard space'
    djob.out{1}.pull.fnames = {mask};
    djob.out{1}.pull.savedir.saveusr = {odir};
    switch deftype
        case 'defs'
            djob.comp{1}.def = InDefs(subj,:); % Name of the deformation field (iy_*.nii)
        case 'dartel'
            djob.comp{1}.dartel = struct('flowfield',{{}},'template',{{}},'times',{{}},'K',{{}});
            djob.comp{1}.dartel.flowfield = {Inflowfields(subj,:)};% Name of the flow_field (u_*.nii)
            djob.comp{1}.dartel.template  = DRTLtemplate;
            djob.comp{1}.dartel.times     = [0 1];
            djob.comp{1}.dartel.K         = DRTLsteps{1};
            djob.comp{1}.dartel = orderfields(djob.comp{1}.dartel, {'flowfield','template','times','K'});
            djob.comp{2}.id = struct('space',{{}});% identity matrix
            djob.comp{2}.id.space = {deblank(actPT)};
    end
    dmask = geg_iydeformation(djob);
    dmaskm = spm_vol(dmask.warpedName{1});
    
    % make sure dimensions are equal (fix it if not)
    Vi = [PTm,dmaskm];% new order of inputs to always match first image space
    [sts,~] = spm_check_orientations(Vi, false);
    tot = 0; vox = 0;
    if ~sts        
        % Computing mean values (images with diff dim/vox)
        S1 = spm_vol(S);
        for i=1:S1(subj).dim(3),
            M    = spm_matrix([0 0 i]);
            img  = spm_slice_vol(S1(subj),M,S1(subj).dim(1:2),0);
            Mmsk = dmaskm(1).mat\S1(subj).mat*M;
            maskm = spm_slice_vol(dmaskm(1),Mmsk,S1(subj).dim(1:2),0);
            img  = img.*maskm;
            img(~isfinite(img)) = 0;
            tot = tot + sum(img(img>0)); %sums only positive voxels
            img(img>0)=1; %sets positive voxels to 1
            vox = vox + sum(img(img>0)); %counts only positive voxels by summing 1s
        end        
    else
        % Computing mean values (images with equal dim/vox)
        img = spm_read_vols(PTm);
        dmaski = spm_read_vols(dmaskm);
        dmaski(dmaski>0) = 1;% binarize the mask just in case is not
        img  = img.*dmaski;
        img(~isfinite(img)) = 0;
        tot = tot + sum(img(img>0)); %sums only positive voxels
        img(img>0)=1; %sets positive voxels to 1
        vox = vox + sum(img(img>0)); %counts only positive voxels by summing 1s        
    end
    signal(subj,1) = tot;
    % divide sum of voxel-values by total number of contributing voxels
    % (i.e. voxels with positive value in (PET) image)
    average(subj,1) = signal(subj,1)/vox;
    
    spm_progress_bar('Set',subj);
    delete(dmaskm.fname)
end
spm_progress_bar('Clear')
fprintf('\n')
disp (average);
disp('Done')

