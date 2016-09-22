function geg_create_atlas_WMrois(job)
% This function automatically parcells the brain WM into various regions. 
% It uses previously computed transformations from MRI to MNI space
% to resample the selected atlas to native space, thus, allowing to perform
% the WM parcellation in MRI native space. 
% It also defines the new ROIs according to the segmentation process. i.e.
% voxels declared as grey/white matter or CSF during the segmentation.
% Note that the atlas must be deep enough to be able to declare WM regions,
% e.g. AAL or Hammers atlases.
% It also creates a the txt description file, necessary to perform the GTM.
% This function does not add the background as a label, that is always
% performed during the GTM process.
% 
% Methodology:
% 1) The atlas is deformed to native space
% 2) The tissue compartments are thresholded to avoid overlapping.
% 3) The atlas is multiplied by the GM tissue map, obtaining segmented GM
% ROIs. This operation therefore defines voxels that were assigned to grey
% matter regions.
% 4) The atlas is multiplied by the WM tissue map, obtaining segmented WM
% ROIs.
% 5) The atlas is multiplied by the CSF tissue map, this creates
% coimplementary regions from the atlas definition.
% 6) Finally, complementary segmented regions are created by joining the
% three segmented tissues toghether.
%
% As atlas is the same for all subjects it assumes a single txt
% file containing 2 columns (names and assigned-ID-value)
%_______________________________________________________________________
% Copyright (C) 2015
%
% Gabriel Gonzalez-Escamilla
% $Id: geg_PVEcorrection2.m 001 2016-03-29 22:27:36Z $
% 
% 
% rev = '$Rev: 002 $'; % 04-May-2016


if nargin==1
    tim = fieldnames(job.SegImgs);
    tim = tim{1};
    switch tim
        case 'Tsegs'
            N1     = job.SegImgs.Tsegs.tiss1;
            N2     = job.SegImgs.Tsegs.tiss2;
            N3     = job.SegImgs.Tsegs.tiss3;
            labels = 0; S = size(N1,1);
        case 'lblsegs'
            N      = job.SegImgs.lblsegs.pvelbl;
            labels = 1; S = size(N,1);
    end
    gmthresh     = geg_petpve12_get_defaults('PVEopts.GTMtissthr');
    ParcAtlasdes = job.ParcDescF;
    erothr       = fieldnames(job.wmOpts);
    erothr       = erothr{1};
    switch erothr
        case 'wmthresh'
            whthr  = geg_petpve12_get_defaults('PVEopts.GTMtissthr');
            thresh = true;
        case 'wmero'
            erothresh = geg_petpve12_get_defaults('PVEopts.EroThresh');
            erofwhm   = geg_petpve12_get_defaults('PVEopts.Erofwhm');
            thresh = false;
    end
    cin = fieldnames(job.csfOpts);
    cin = cin{1};
    switch cin
        case 'oneCSF'
            singleCSF = true;
        case 'CSFrois'
            singleCSF = false;
    end
    dfn          = fieldnames(job.PVE_AParc_opts);
    dfn          = dfn{1};
    switch dfn
        case 'type2'
            % 'Atlas in standard space'
            ParcAtlas               = job.PVE_AParc_opts.type2.ParcAtlas;
            djob.out{1}.pull        = struct('fnames',{{}},'savedir',{{}},'interp',{{}},'mask',{{}},'fwhm',{{}});
            djob.out{1}.pull.mask   = 0;
            djob.out{1}.pull.interp = geg_petpve12_get_defaults('PVEopts.RoiAtlasInterp');
            djob.out{1}.pull.fwhm   = [0 0 0];
            deftype = fieldnames(job.PVE_AParc_opts.type2.InDefs);
            deftype = deftype{1};
            switch deftype
                case 'defs'
                    InDefs  = job.PVE_AParc_opts.type2.InDefs.defs;
                case 'dartel'
                    Inflowfields = job.PVE_AParc_opts.type2.InDefs.dartel.flowfield;
                    DRTLtemplate = {''}; 
                    DRTLsteps    = job.PVE_AParc_opts.type2.InDefs.dartel.K;
            end            
        otherwise
            fprintf('Not a recognized option!\n')
            return;
    end
elseif nargin < 1
    labels = spm_input('Type of segmented image(s)','1','b',{'Tissue_Segments','Labeled_pveMRI'},[0,1]);
    if labels
        N  = spm_select(Inf,{'^c0.*\.nii','^p0.*\.nii'},'Select native MRI labeled image');
        S  = size(N,1);
    else
        N1 = spm_select(Inf,{'^c1.*\.nii','^p1.*\.nii'},'Select native GM maps');
        N2 = spm_select(size(N1,1),{'^c2.*\.nii','^p2.*\.nii'},'Select native WM maps');
        N3 = spm_select(size(N1,1),{'^c3.*\.nii','^p3.*\.nii'},'Select native CSF maps');
        S  = size(N1,1);
    end
    gmthresh   = spm_input('GM threshold','1','e',0,1);
    thresh  = spm_input('Option to include WM/CSF in the atlas?','+1','b',{'threhsold','erode'},[1,0]);
    if ~thresh
        erothresh = geg_petpve12_get_defaults('PVEopts.EroThresh');
        erofwhm   = geg_petpve12_get_defaults('PVEopts.Erofwhm');
        thresh = false;
    else
        whthr  = geg_petpve12_get_defaults('PVEopts.GTMtissthr');
        thresh = true;
    end
    dfn = spm_input('Pacellation atlas','+1','b',{'AtlasStandard'},{'type2'});
    dfn = dfn{1};
    switch dfn        
        case 'type2'
            % 'Atlas in standard space'
            ParcAtlas        = {spm_select([1,1],'^*\.nii','Select Template Atlas in standard space ')};
            ParcAtlasdes     = spm_select([1,1],'^*\.txt','Select Atlas descriptor file');
            djob.out{1}.pull = struct('fnames',{{}},'savedir',{{}},'interp',{{}},'mask',{{}},'fwhm',{{}});
            djob.out{1}.pull.mask   = 0;
            djob.out{1}.pull.interp = @(val)geg_petpve12_get_defaults('PVEopts.RoiAtlasInterp', val{:});
            djob.out{1}.pull.fwhm   = [0,0,0];            
            deftype = spm_input('Input Fields type','+1','b',{'Deformations','FlowFields'},{'defs','dartel'});
            deftype = deftype{1};
            switch deftype
                case 'defs'
                    InDefs  = spm_select([S,1],'^iy_.*\.nii','Select deformation-fields for every subject');
                case 'dartel'
                    Inflowfields = spm_select([S,1],'^u_.*\.nii','Select flow-fields for every subject');
                    DRTLtemplate = {''}; 
                    DRTLsteps    = 6;
            end 
    end
end
save_segs = geg_petpve12_get_defaults('WMrois.savesegs'); % option to save the extracted tissue compartments from PVC-labeled MRimage
save_extras = geg_petpve12_get_defaults('WMrois.saveextras'); % option to save the comprobatory images at different stages

% Preparing atlas description (same for all subjects)
% Must acknowledge how many ROIs are and the names.
AtlDesc=importdata(char(ParcAtlasdes));
SubjParcAtlasNames=AtlDesc.textdata;
SubjParcAtlasIds=AtlDesc.data;

% get current release number
A = ver;
for i=1:length(A)
    if strcmp(A(i).Name,'PET Partial Volume Effects correction Toolbox')
        r = str2double(A(i).Version);
    end
end

% Running for every subject
for subj=1:S
    
    if exist('r','var')
        fprintf('PETPVE12 r%d: %s\n',r,char(subj));
    end

    % Set starting time
    fprintf('Starting time: %s \n',datestr(clock,'local'))
    
    % Loading data
    if ~labels
            actGM = char(N1(subj,:));% unsmoothed GM name
            actWM = char(N2(subj,:));% unsmoothed WM name
            actCSF = char(N3(subj,:));% unsmoothed CSF name
            GMm = spm_vol(actGM);% GM id matrix
            GMi = spm_read_vols(GMm);% GM image
            CSFm = spm_vol(actCSF);% CSF id matrix
            WMm = spm_vol(actWM);% WM id matrix
            WMi = spm_read_vols(WMm);% WM image
            CSFi = spm_read_vols(CSFm);% CSF image
    elseif labels
            actlbl = char(N(subj,:));% unsmoothed labeled image name
            lblm = spm_vol(actlbl);% labeled image id matrix
            lbli = spm_read_vols(lblm);% labeled image
            GMi = single(round(lbli)==2);
            WMi = single(round(lbli)==3); 
            CSFi = single(round(lbli)==1); 
            GMm = lblm;
            GMm.fname = fullfile(pth,['c1' nam(3:end),ext]);
            WMm = lblm;
            WMm.fname = fullfile(pth,['c2' nam(3:end),ext]);
            CSFm = lblm;
            CSFm.fname = fullfile(pth,['c3' nam(3:end),ext]);
            if save_segs
                spm_write_vol(GMm,GMi);
                spm_write_vol(WMm,WMi);
            end
    end
    [pth,nam,ext] = spm_fileparts(GMm.fname);
    odir = pth; nam = nam(3:end);
    if strncmpi(ext,'.nii',4)
        cvt = [nam ext];
    else
        cvt = [nam '.nii'];
    end 
    
    if ~labels, [~,namw,~,~] = spm_fileparts(WMm.fname); [~,namc,~,~] = spm_fileparts(CSFm.fname); end
    switch dfn
        case 'type2'
            % 'Parcellation atlas in standard space'
            fprintf('Deforming ROI atlas to %s,. \n',actGM)
            djob.out{1}.pull.fnames = ParcAtlas;
            djob.out{1}.pull.savedir.saveusr = {odir};
            switch deftype
                case 'defs'
                    djob.comp{1}.def = InDefs(subj); % Name of the deformation field (iy_*.nii)
                case 'dartel'
                    djob.comp{1}.dartel = struct('flowfield',{{}},'template',{{}},'times',{{}},'K',{{}});
                    djob.comp{1}.dartel.flowfield = Inflowfields(subj);% Name of the flow_field (u_*.nii)
                    djob.comp{1}.dartel.template  = DRTLtemplate;
                    djob.comp{1}.dartel.times     = [0 1];% always perform template to subject projection
                    djob.comp{1}.dartel.K         = DRTLsteps;
                    djob.comp{1}.dartel = orderfields(djob.comp{1}.dartel, {'flowfield','template','times','K'});
                    djob.comp{2}.id = struct('space',{{}});% identity matrix
                    if ~labels
                       djob.comp{2}.id.space = {deblank(actGM)}; 
                    else
                       djob.comp{2}.id.space = {deblank(actlbl)};
                    end
            end
            datlas = geg_iydeformation(djob);
            datlasm = spm_vol(datlas.warpedName{1});
            datlasi = spm_read_vols(datlasm);% deformed atlas
            delete(datlasm.fname)
            [~,n,~] = spm_fileparts(datlasm.fname);
            datlasm.fname = fullfile(odir,[n '_' cvt(1:end-4) ext]);
            if save_extras, spm_write_vol(datlasm,datlasi); end
    end
    
    % To avoid GM overlapping (WM/CSF counts in GM), mask the atlas to 
    % GM-counts only. Not necesarry if using PVC segments.
    if ~labels
        datlas1i = datlasi .* (GMi>gmthresh); % create GM segmented atlas
        
        % To introduce WM and CSF onto GTM, the tissue maps are thesholded
        % (or eroded) to avoid overlapping with GM. CSF is always included
        % in the GTM and, if requested, will be zeroed only in the signal
        % vector. 
        newmm = WMm; newcsfm = WMm; 
        if thresh
            WMe = double(WMi>whthr); 
            newmm.fname = fullfile(pth,['thr' namw ext]);
            CSFe = CSFi>whthr;
            newcsfm.fname = fullfile(pth,['thr' namc ext]);
        else
            WMe = geg_erodeMask(WMm,erothresh,erofwhm,1);
            newmm.fname = fullfile(pth,['e' namw ext]);
            CSFe = geg_erodeMask(CSFm,.9,erofwhm,1);
            newcsfm.fname = fullfile(pth,['e' namc ext]);
        end
        if save_extras, spm_write_vol(newmm,WMe); spm_write_vol(newcsfm,CSFe); end
        datlas2i = datlasi .* WMe; % create WM segmented atlas
        datlas3i = datlasi .* CSFe; % create CSF segmented atlas
    else
        datlas1i = datlasi .* GMi; % create GM segmented atlas
        datlas2i = datlasi .* WMi; % create WM segmented atlas
        datlas3i = datlasi .* CSFi; % create CSF segmented atlas
    end
    
    if save_extras, % Save GM/WM/CSF segments non-overlapped atlas
        datlas1m = datlasm; datlas1m.fname = fullfile(pth,['ac1' nam ext]); spm_write_vol(datlas1m,datlas1i);
        datlas2m = datlasm; datlas2m.fname = fullfile(pth,['ac2' nam ext]); spm_write_vol(datlas2m,datlas2i);
        datlas3m = datlasm; datlas3m.fname = fullfile(pth,['ac3' nam ext]); spm_write_vol(datlas3m,datlas3i);
    end
    
    % Make sure there's no overlapping in brain regions
    % Create brn binary image
    brn = plus(datlas1i>0,datlas2i>0); % Joint GM & WM ROIs
    brn = plus(brn>0,datlas3i>0); % add CSF ROIs
    if save_extras,  % Save brn image
        newbrnm = GMm; newbrnm.fname = fullfile(pth,['brn' nam ext]); newbrnm.descrip = 'brain-tissue class';
        spm_write_vol(newbrnm, brn);
    end    
    bov = find(brn(:)>1, 1); % Measure overlapping
    if ~isempty(bov), fprintf(2,'WARNING: overlapping ROIs exist'); end
    if save_extras,  % Create & save background image
        bkgnd = ones(size(GMi)); bkgnd(brn==1) = 0; % background binary image        
        newbgnm = GMm;  newbgnm.fname = fullfile(pth,['bgn' nam ext]); newbgnm.descrip = 'non-tissue class';
        spm_write_vol(newbgnm,bkgnd);
    end
    
    % Joint GM/WM/CSF into a single atlas image
    datlas4i = zeros(size(datlasi));
    datlas4i(datlas1i>0) = datlas1i(datlas1i>0); % add GM ROIs
    % Create different label-IDs for WM
    wmlbls = numel(SubjParcAtlasIds);
    wmRnames = cell(wmlbls,1); wmRids = zeros(wmlbls,1); maxAtls = max(SubjParcAtlasIds);
    for lbl = 1:wmlbls;
        actlbl = SubjParcAtlasIds(lbl);
        datlas2i(round(datlas2i)==(actlbl)) = maxAtls+actlbl;%+100+lbl;
        wmRids(lbl,1) = maxAtls+actlbl;%+100+lbl;
        wmRnames{lbl,1} = ['wm-' SubjParcAtlasNames{lbl}];
    end
    datlas4i(datlas2i>0) = datlas2i(datlas2i>0); maxwm = max(wmRids); % add WM ROIs
    if singleCSF
        datlas4i(datlas3i>0) = maxwm+50; % add CSF as a single ROI
        csfID = maxwm+50;
    else
        datlas4i(datlas3i>0) = datlas3i(datlas3i>0); % add CSF ROIs
    end    
    datlasm.fname = fullfile(odir,[n '_' cvt(1:end-4) '_WM' ext]); datlasm.descrip = 'PETPVE12 - nonOverlaping';
    spm_write_vol(datlasm,datlas4i); 
    
    % Set end time
    fprintf('End time: %s \n',datestr(clock,'local'))
end

% Save PVE corrected values of ROIs for every subject in a txt file
if singleCSF
    allRNames = [SubjParcAtlasNames;'csf';wmRnames]; allRids = [SubjParcAtlasIds; csfID; wmRids];
else
    allRNames = [SubjParcAtlasNames;wmRnames]; allRids = [SubjParcAtlasIds; wmRids];
end
vol_txt = fullfile(odir,[n '_WM.txt']);
fid = fopen(vol_txt, 'w');
for i=1:length(allRNames)
    label = allRNames{i}; labelval = allRids(i);
    fprintf(fid,'%s\t',label);
    fprintf(fid,'%5.4f\n',labelval);
end;
fprintf(fid,'\n');
fclose(fid);

fprintf('\n');