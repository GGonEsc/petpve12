function geg_calc_atlasROImeans_native()
% This function converts a standard mask into native space and computes the
% corresponding mean signal.
% This usefull to extract values on the uncorrected PET data, so intensity
% normalization can be performed on both, uncorrected and PVEc PETs, for
% comparison.

% data
atlas  = spm_select([1,1],'image','Select template atlas');
AtlasDesc = spm_select([1,1],'^*\.txt','Select Atlas descriptor file');
S  = spm_select(Inf,'image','Select native PET images');
opth = spm_select([1,1],'Dir','Select output directory');

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

AtlDesc=importdata(char(AtlasDesc));
ParcAtlasNames=AtlDesc.textdata;
ParcAtlasIds=AtlDesc.data;
disp(atlas)
signal = zeros(size(S,1),numel(ParcAtlasIds));
average = zeros(size(S,1),numel(ParcAtlasIds));
snames = cell(size(S,1),1);

[~,nam,~] = spm_fileparts(atlas);
tots_xls = fullfile(opth,[nam, '_ROItots.xls']);
xlswrite(tots_xls,['SubjID', ParcAtlasNames'],'ROItotals','A1')

% Running for every subject
for subj=1:size(S,1)
    % Loading data
    actPT = char(S(subj,:));% PET name
    [odir,snam,~] = spm_fileparts(actPT);
    PTm = spm_vol(actPT);% PET header    
    snames{subj,1} = snam;
    
    % 'Template atlas in standard space'
    djob.out{1}.pull.fnames = {atlas};
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
    fprintf('Deforming %d of %d, ',subj,size(S,1));
    dmask = geg_iydeformation(djob);
    dmaskm = spm_vol(dmask.warpedName{1});
    
    % make sure dimensions are equal (fix it if not)
    Vi = [PTm,dmaskm];% new order of inputs to always match first image space
    [sts,~] = spm_check_orientations(Vi, false);
    if ~sts
        % Adjust the images
        [dmaski,~] = geg_reslice(PTm,dmaskm,0);
    else
        dmaski = spm_read_vols(dmaskm);
    end
    PTi = spm_read_vols(PTm);
    
    fprintf('Estimating Regions: ');
    for k=1:numel(ParcAtlasIds)
        actp = ParcAtlasIds(k);
        fprintf('%d,',actp) % progression indicator
        img = PTi;
        
        % Identify actual parcel
        actparc = round(dmaski)==(actp);% logical
        actparc = double(actparc);% numerical
        
        % Computing mean values
        img  = img.*actparc;
        img(~isfinite(img)) = 0;
        tot = sum(img(img>0)); %sums only positive voxels
        img(img>0)=1; %sets positive voxels to 1
        vox = sum(img(img>0)); %counts only positive voxels by summing 1s
        signal(subj,k) = tot;
        % divide sum of voxel-values by total number of contributing voxels
        % (i.e. voxels with positive value in (PET) image)
        average(subj,k) = signal(subj,k)/vox;
    end
    delete(dmask.warpedName{1}(1:end-2));  
    fprintf('\n')
end
fprintf('\n')
% print to a xls file
xlswrite(tots_xls,snames,'ROItotals','A2')
xlswrite(tots_xls,average,'ROItotals','B2')
disp('Done')

