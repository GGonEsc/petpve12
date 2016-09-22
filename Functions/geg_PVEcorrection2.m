function geg_PVEcorrection2(job)
% geg_PVEcorrection2.m
% 
% PVC VOI and Brain MR based: This correction correct for the spillover at
% multiple ROIs of the brain by using the so called geometric transfer
% matrix (GTM). The GTM describes the contribution of each source region to
% the measured regional activities. After GTM has been created searches for
% WM ROIs and use this value to correct PET data via Muller-Gartner
% algorithm. This is often called the modified Muller-Gartner (mMG)
% approach [1]. 
% 
% [1] Rousset OG, Ma Y, Evans AC (1998) Correction for partial volume
% effects in PET: principle and validation. J Nucl Med 39:904-911 
%---------------------------------------------------------------------
% Partial volume correction of PET images (in native space - coregistered)
% corrects PET images for PVE based on Rousset VOI algorithm and 3
% compartment model:
% To calculate true activities in each region, the measured activity in the
% ith region (mi) equals a weighted sum of the true values in all source
% regions, with the weights appearing in the ith row of GTM. 
% mi = GTMi1T1 + GTMi2T2 + ... + GTMiRTR
% Then the estimates of the true activities T were obtained by multiplying
% the vector of measured activities m  [m1, m2 ... mR]by the inverse of the
% GTM.  
% T = GTM-1 * m
% 
% in terms of code the expressions ROI, VOI and parcel are indifferently used
% 
% Author: Gabriel Gonzalez-Escamilla March 2015
%
% If using the included aal atlas, please also reference:
%     Tzourio-Mazoyer N, Landeau B, Papathanassiou D, Crivello F, Etard O,
%     Delcroix N, Mazoyer B, Joliot M. Automated anatomical labeling of
%     activations in SPM using a macroscopic anatomical parcellation of the
%     MNI MRI single-subject brain. Neuroimage. 2002 Jan;15(1):273-89.
%     PUBMED ID: 11771995    
%_______________________________________________________________________
% Copyright (C) 2015
%
% Gabriel Gonzalez-Escamilla
% $Id: geg_PVEcorrection2.m 001 2015-03-24 20:13:03Z $
% 
% 
% rev = '$Rev: 020 $'; % 04-May-2016

if nargin==1
    S   = job.PETdata;
    tim = fieldnames(job.SegImgs);
    tim = tim{1};
    switch tim
        case 'Tsegs'
            N1     = job.SegImgs.Tsegs.tiss1;
            N2     = job.SegImgs.Tsegs.tiss2;
            CSFin = fieldnames(job.PVEopts.CSFsignal.CSFsigIn);
            CSFin = CSFin{1};
            if strncmp(CSFin,'CSFincl',7) % include CSF tissue
                csfmulti = 0;
            elseif strncmp(CSFin,'CSFzeroing',10) % existing CSF ROIs
                csfmulti = 1;
            end
            N3 = job.PVEopts.CSFsignal.CSFcalc.tiss3;
            labels = 0;
        case 'lblsegs'
            N      = job.SegImgs.lblsegs.pvelbl;
            CSFin = fieldnames(job.PVEopts.CSFsignal); % include CSF tissue
            CSFin = CSFin{1};
            if strncmp(CSFin,'CSFincl',7) % include CSF tissue
                csfmulti = 0;
            elseif strncmp(CSFin,'CSFzeroing',10) % existing CSF ROIs
                csfmulti = 1;
            end
            labels = 1;
    end
    PETpsf       = job.PVEopts.fwhm_PSF;
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
    dfn          = fieldnames(job.PVE_AParc_opts);
    dfn          = dfn{1};
    switch dfn
        case 'type1' 
           % 'Parcellation in subject space'
            SubjParcAtlas = job.PVE_AParc_opts.type1.ParcAtlas;
            % As atlas is the same for all subjects it assumes a single txt
            % file containing 2 columns (names and assigned-value)            
            tissmult = false;
        case 'type2'
            % 'Atlas in standard space'
            ParcAtlas               = job.PVE_AParc_opts.type2.ParcAtlas;
            tissmult = true;
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
    S      = spm_select(Inf,'image','Select native coregistered PET images');
    labels = spm_input('Type of segmented image(s)','1','b',{'Tissue_Segments','Labeled_pveMRI'},[0,1]);
    if labels
        N  = spm_select(size(S,1),{'^c0.*\.nii','^p0.*\.nii'},'Select native PVE labeled maps');
    else
        N1 = spm_select(size(S,1),{'^c1.*\.nii','^p1.*\.nii'},'Select native GM maps');
        N2 = spm_select(size(S,1),{'^c2.*\.nii','^p2.*\.nii'},'Select native WM maps');
        N3 = spm_select(size(S,1),{'^c3.*\.nii','^p3.*\.nii'},'Select native CSF maps');
    end
    PETpsf     = spm_input('PET PSF (FWHM in mm)','+1','e',[6,6,6],[1,3]);
    gmthresh   = spm_input('GM threshold','1','e',0,1);
    csfmulti = spm_input('CSF-signal?','+1','b',{'atlas_defined','tissue_map'},[1,0]);
    thresh  = spm_input('Option to include WM/CSF in the GTM?','+1','b',{'threhsold','erode'},[1,0]);
    if ~thresh
        erothresh = geg_petpve12_get_defaults('PVEopts.EroThresh');
        erofwhm   = geg_petpve12_get_defaults('PVEopts.Erofwhm');
    else
        whthr  = geg_petpve12_get_defaults('PVEopts.GTMtissthr');
    end
    dfn = spm_input('Pacellation atlas','+1','b',{'AtlasNative','AtlasStandard'},{'type1','type2'});
    dfn = dfn{1};
    switch dfn
        case 'type1'
            % 'Parcellation in subject space'
            ParcAtlasdes  = spm_select([1,1],'^*\.txt','Select Atlas descriptor file');
            SubjParcAtlas = spm_select(size(S,1),'^*\.nii','Select ROI atlas for every subject');
            tissmult = false;
        case 'type2'
            % 'Atlas in standard space'
            ParcAtlas        = {spm_select([1,1],'^*\.nii','Select Template Atlas in standard space ')};
            ParcAtlasdes     = spm_select([1,1],'^*\.txt','Select Atlas descriptor file');
            tissmult = true;
            djob.out{1}.pull = struct('fnames',{{}},'savedir',{{}},'interp',{{}},'mask',{{}},'fwhm',{{}});
            djob.out{1}.pull.mask   = 0;
            djob.out{1}.pull.interp = @(val)geg_petpve12_get_defaults('PVEopts.RoiAtlasInterp', val{:});
            djob.out{1}.pull.fwhm   = [0,0,0];            
            deftype = spm_input('Input Fields type','+1','b',{'Deformations','FlowFields'},{'defs','dartel'});
            deftype = deftype{1};
            switch deftype
                case 'defs'
                    InDefs  = spm_select(size(S,1),'^iy_.*\.nii','Select deformation-fields for every subject');
                case 'dartel'
                    Inflowfields = spm_select(size(S,1),'^u_.*\.nii','Select flow-fields for every subject');
                    DRTLtemplate = {''}; 
                    DRTLsteps    = 6;
            end 
    end
end
save_segs = geg_petpve12_get_defaults('PVEopts.GTMsavesegs'); % option to save the extracted tissue compartments from PVC-labeled MRimage
save_extras = geg_petpve12_get_defaults('PVEopts.GTMsaveextras'); % option to save the comprobatory images at different stages
save_GTMs = geg_petpve12_get_defaults('PVEopts.GTMsaveGTM'); % option to save the GTM and the ROI uncorrected/sizes values

% Preparing atlas description (same for all subjects)
% Must acknowledge how many ROIs are and the names.
AtlDesc=importdata(char(ParcAtlasdes)); 
allRNames=AtlDesc.textdata; PAtlasIDs = AtlDesc.data; totPAtlasIDs = numel(PAtlasIDs); 
% Add single CSF label (if requested by the user)
if ~csfmulti
	totPAtlasIDs = totPAtlasIDs+1;
    PAtlasIDs = [PAtlasIDs; max(PAtlasIDs)+5];
    allRNames = [allRNames; 'csf'];
end
% Add background ID
totPAtlasIDs = totPAtlasIDs+1;
PAtlasIDs = [PAtlasIDs; max(PAtlasIDs)+10];
allRNames = [allRNames; 'background'];

% get current release number
A = ver;
for i=1:length(A)
    if strcmp(A(i).Name,'PET Partial Volume Effects correction Toolbox')
        r = str2double(A(i).Version);
    end
end

% Running for every subject
for subj=1:size(S,1)
    
    if exist('r','var')
        fprintf('PETPVE12 r%d: %s\n',r,char(S(subj,:)));
    end

    % Set starting time
    fprintf('Starting time: %s \n',datestr(clock,'local'))
    
    % Loading data
    actPT = char(S(subj,:));% corregistered PET name
    PTm = spm_vol(actPT);% PET id matrix
    [pth,nam,ext] = spm_fileparts(actPT);
    odir = pth;
    if strncmpi(ext,'.nii',4)
        cvt = ['pvc2' nam ext];
    else
        cvt = ['pvc2' nam '.nii'];
    end        
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
    
    % PSF selection
    if PETpsf(1)==PETpsf(2) && PETpsf(1)==PETpsf(3)
        fprintf('Working with isotropic PSF. \n')
        smth=PETpsf(1);
    else
        fprintf('Working with non-isotropic PSF. \n')
        smth=PETpsf;
    end 
    if ~labels, [~,namw,~,~] = spm_fileparts(WMm.fname); [~,namc,~,~] = spm_fileparts(CSFm.fname); end
        
    switch dfn
        case 'type1'
            % 'Parcellation atlas in subject space'
            datlasn = SubjParcAtlas{subj};
            datlasm = spm_vol(datlasn);
            datlasi = spm_read_vols(datlasm);
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
            datlasm.fname = fullfile(odir,[n '_' cvt(5:end-4) ext]);
            spm_write_vol(datlasm,datlasi);
    end
    
    % When deforming the atlas from standard space we need to avoid tissue
    % overlapping (WM/CSF counts in GM), here the atlas is masked to
    % avoid any tissue overlapping. 
    % This is not necesarry if using individual atlas, given that the user
    % might have manually drawn the ROIs or create them by some other
    % subject specific delimitation.
    if tissmult
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
        if save_extras, % Save brain image
            newbrnm = GMm; newbrnm.fname = fullfile(pth,['brn' nam ext]); newbrnm.descrip = 'brain-tissue class';
            spm_write_vol(newbrnm, brn);
        end
        bov = find(brn(:)>1, 1); % Check overlapping
        if ~isempty(bov), fprintf(2,'WARNING: overlapping ROIs exist'); end
        if save_extras, % Create & save background binary image
            bkgnd = ones(size(GMi));
            bkgnd(brn==1) = 0; % background image
            newbgnm = GMm;  newbgnm.fname = fullfile(pth,['bgn' nam ext]); newbgnm.descrip = 'non-tissue class';
            spm_write_vol(newbgnm,bkgnd);
        end
        
        % Joint GM/WM/CSF/BGND into a single atlas image
        datlas4i = zeros(size(datlasi));
        datlas4i(datlas1i>0) = datlas1i(datlas1i>0); % add GM ROIs
        datlas4i(datlas2i>0) = datlas2i(datlas2i>0); % add WM ROIs
        if csfmulti
            datlas4i(datlas3i>0) = datlas3i(datlas3i>0); % add CSF ROIs
        else
            datlas4i(datlas3i>0) = max(PAtlasIDs)+5; % add CSF ROIs            
        end
        spm_write_vol(datlasm,datlas4i); % save the WM/GM/CSF
        datlas4i(datlas4i==0) = max(PAtlasIDs)+10; % add background
        if save_extras, % Save atlas with background
            datlasm.fname = fullfile(odir,[n '_' cvt(5:end-4) '_NO' ext]); datlasm.descrip = 'PETPVE12 - nonOverlaping';
            spm_write_vol(datlasm,datlas4i);
        end
    else        
        datlas4i = datlasi;
        datlas4i(datlas4i==0) = max(PAtlasIDs)+10; % add background
        if save_extras, % Save atlas with background
            datlasm.fname = fullfile(odir,[n '_' cvt(5:end-4) '_NO' ext]); datlasm.descrip = 'PETPVE12 - nonOverlaping';
            spm_write_vol(datlasm,datlas4i);
        end
    end
    
    % make sure dimensions are equal (fix it if not)
    [PTi, PT2] = geg_reslice(GMm,PTm,0);
    % eliminate possible NaNs in PET image
    PTi(isnan(PTi)) = 0;
    % save resliced PET
    [~,n2,~] = spm_fileparts(PTm.fname); PT2.fname = fullfile(odir,['r' n2 ext]); spm_write_vol(PT2,PTi);    
        
    % constructing GTM matrix     
    fprintf('Constructing GTM matrix. \n')    
    [parcOPET, GTM, Rsize] = runGTM(PTi, totPAtlasIDs, PAtlasIDs, datlas4i,smth);
    
    % Get rid of NaNs (even when it should not have any)
    parcOPET(isnan(parcOPET)) = 0;
    GTM(isnan(GTM)) = 0; 
    if save_GTMs,  % save the GTM of each subject
        save(fullfile(odir,['GTM', nam, '.mat']),'GTM','parcOPET');
    end
    
    % To invert the square diagonal matrix, single value decomposition can
    % be used [U,S,V] = svd(X), where U is an m by m real or complex
    % unitary matrix (AKA left eigenvectors), S is an m by n rectangular
    % diagonal matrix with non-negative real numbers on the diagonal (a
    % diagonal matrix whose diagonal elements are the singular values of
    % M), and V* (the conjugate transpose of V, or simply the transpose of
    % V if V is real) is an n by n real or complex unitary matrix (columns
    % are the eigenvectors of the M transposed * M matrix - AKA right
    % eigenvectors).
    % This mean that you simple invert each element of the matrix: 
    %         M^(-1) = (USV)^(-1) = (V)^(-1)S^(-1)U^(-1) = VS^(-1)
    % 
    % [U,s,V] = svd(GTM);
    % invGTM = V * inv(s) * (U'); % == invGTM=inv(V)*inv(s)* inv(U)';
    % 
    % This is the Pen-Rose pseudoinverse of a matrix for linear squares
    % problem solution when the matrix is geometric so:
    %   invGTM = pinv(GTM);% pinv uses the sdv to compute the inverse
    %
    % PVE corrected ROI PET values (observed PET activity multiplied by the
    % inverse of the GTM): 
    %   tPET = invGTM * (parcOPET);
    % 
    % The least-squares problem can be solved using backslash in Matlab. In
    % the particular situation of a geometric matrix the lsqlin function
    % does exactly the same as pinv. lsqnonneg avoids negative elements:
    tPET = lsqnonneg(GTM,parcOPET);
        
    if save_GTMs, % Save PVE raw values of ROIs for every subject in a txt file
        vol_txt = fullfile(odir,[nam, '_labels.txt']);
        fid = fopen(vol_txt, 'w');
        for i=1:length(parcOPET)
            labelval = allRNames{i};
            fprintf(fid,'%s\t',labelval);
        end;
        fprintf(fid,'\n');
        for i=1:length(parcOPET)
            labelval = parcOPET(i);
            fprintf(fid,'%5.4f\t',labelval);
        end;
        fclose(fid);
    end
    
    % Save PVE corrected values of ROIs for every subject in a txt file
    vol_txt = fullfile(odir,['pvc', nam, '_labels.txt']);
    fid = fopen(vol_txt, 'w');    
    for i=1:length(tPET)
        labelval = allRNames{i};
        fprintf(fid,'%s\t',labelval);
    end;
    fprintf(fid,'\n');
    for i=1:length(tPET)
        labelval = tPET(i);
        fprintf(fid,'%5.4f\t',labelval);
    end;
    fclose(fid);
    
    if save_GTMs, % Save ROI sizes for every subject in a txt file
        vol_txt = fullfile(odir,['Rsizes', nam, '_labels.txt']);
        fid = fopen(vol_txt, 'w');
        for i=1:length(Rsize)
            labelval = allRNames{i};
            fprintf(fid,'%s\t',labelval);
        end;
        fprintf(fid,'\n');
        for i=1:length(Rsize)
            labelval = Rsize(i);
            fprintf(fid,'%5.4f\t',labelval);
        end;
        fclose(fid);
    end
    
    % Save Atlas-Parcellated image containing GTM corrected values
    datlas5i = zeros(size(datlasi));
    for roi=1:numel(PAtlasIDs)
        actparcID = PAtlasIDs(roi);
        datlas5i(round(datlas4i)==(actparcID)) = tPET(roi);
    end    
    oparc = datlasm; oparc.descrip = 'PETPVE12 - pvcGTM';
    oparc.fname = fullfile(odir,['pvc', nam, '_labels.nii']);
    spm_write_vol(oparc,datlas5i);
    
    % Set end time
    fprintf('End time: %s \n',datestr(clock,'local'))
    
    % clear to save memory
    clear parcOPET tPET invGTM GTM brn    
    
    fprintf('\n');
end

%------------------------------------------------------------------------
%========================================================================
% SUB-FUNCTIONS
%========================================================================
%------------------------------------------------------------------------
function [ObsPET, GTMatrix, ROIsizes] = runGTM(PTi, totPAtlasIDs, PAtlasIDs, datlas4i,smth)
% This function creates the GTM

%-Start progress plot
spm_progress_bar('Init',numel(PAtlasIDs),'Constructing GTM matrix','ROIs completed');

% Pre-Define variables
ObsPET = zeros(totPAtlasIDs,1); ROIsizes = zeros(1,totPAtlasIDs);
GTMatrix = zeros(totPAtlasIDs,totPAtlasIDs);

for parc1=1:totPAtlasIDs
    % Extract ROIi from Atlas to determine the vector with observed
    % uptake values and convolve this region
    actparcID = PAtlasIDs(parc1);
    
    % Identify voxels of actual ROI
    aparcid = round(datlas4i)==(actparcID); % is logical
    aparcid = double(aparcid); % now is numerical
    
    % Count the number of voxels in actual ROIi
    ROIsizes(1,parc1) = nnz(aparcid(~isnan(aparcid)));
    
    % Create and open convolved region
    caparci = zeros(size(aparcid));%This is needed to create the convolved image
    spm_smooth(aparcid,caparci,smth,0);% convolve with PSF
    
    % Create vector with observed ROI activity
    cmask=PTi(aparcid==1 & aparcid~=0);
    ObsPET(parc1,1) = mean(cmask);
    
    for parc2=1:totPAtlasIDs
        % Extract ROIj from Atlas to determine the overlapping with ROIi
        actparc2ID = PAtlasIDs(parc2);
        
        % identify voxels of actual roi
        aparc2i = round(datlas4i)==(actparc2ID); % logical indexing
        aparc2i = double(aparc2i);
        
        % Conform the GTM
        % Compute fractional contribution of parc1 (ROIi) on parc2 (ROIj)
        % this is the overlapping of convolved parc1 on original parc2
        aparc2i(isnan(aparc2i)) = 0; caparci(isnan(caparci)) = 0; aparc2i(isnan(aparc2i)) = 0;
        GTMatrix(parc2,parc1) = (sum(sum(sum(aparc2i.*caparci)))) / nnz(aparc2i);
        
        clear aparc2m xY1 xY2 xY3 aparc2i
    end
    spm_progress_bar('Set',parc1);
    clear caparci aparci
end
spm_progress_bar('Clear')
