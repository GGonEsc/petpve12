function geg_PVEcorrection(job)
% geg_PVEcorrection.m
% 
% PVC Brain MR based : This correction is based on the assumption that
% white matter uptake is homogeneous. All brain pixels are classified as
% white matter (WM) or grey matter (GM) and sorted into respective
% segments. Based on these segments and the assumed PET resolution the
% spill-out from WM to GM can be estimated and subtracted. Similarly, the
% spill-out from GM to the surroundings can be estimated and compensated
% for. The result is a grey matter image with corrected activity values in
% all pixels. This method was introduced by Muller-Gartner et al. [1]. 
% Given a brain PET and an anatomical MRI of a patient, the implementation
% allows the user performing the segmentation and apply the
% Muller-Gartner PVC in a fully automated way.  
% 
% [1]. Muller-Gartner HW, Links JM, Prince JL, Bryan RN, McVeigh E, Leal
% JP, Davatzikos C, Frost JJ. Measurement of radiotracer concentration in
% brain gray matter using positron emission tomography: MRI-based
% correction for partial volume effects. J Cereb Blood Flow Metab.
% 1992;12(4):571-83. 
%---------------------------------------------------------------------
% Partial volume correction of PET images (in native space - coregistered)
% corrects PET images for PVE based on Muller-Gartner algorithm and 3
% compartment model:
% PET-PVC = (PET-WMsignal*s(PSF)WM-CSFsignal*s(PSF)CSF)*GMmask(thr)/s(PSF)GM
% s(PSF): smoothed by point-spread function
% Usage:
% Choose threshold (thr) for GM mask.
% 
% Based on the idea of Michel Grothe 08 November 2013 $ to perform everything 
% using the image calculator (slice-based), instead of using the whole image.
% Author: Gabriel Gonzalez-Escamilla March 2015
%
%   % To apply the PVE correction method the MR image is decomposed (i.e.,
    % segmented) into its constituent tissue types, and a mathematical
    % model of the scanner resolution effects is applied to each in order
    % to account for them so that the radioactivity concentration in the
    % non-uniform tissue may be estimated. 
    % 
    % The tissue types are distinct (masks) and do not overlap. It considers 
    % the area being imaged as being composed of a small number of
    % different tissue types; for example, gray matter (GM), white matter 
    % (WM), and cerebral spinal fluid (CSF) - with all but one (GM) of
    % the types being uniform with respect to radioactivity concentration.
    % 
    % Use regions of interest to mark structures of presumed uniform
    % radioactivity concentration in the tissues of interest, as well as to
    % mark adjacent background areas. For each such region j, create a
    % synthetic image setting pixel values to unity in the given region and
    % zero elsewhere. 
    % 
    % The convolution operator and the spatially invariant PSF of PET is
    % aproximated by a gaussian function. i.e. smoothing the compartments
    % by the FWHM of the PSF.
    % Blur the image by convolution with a Gaussian kernel with a FWHM that
    % emulates the spatial resolution of the particular PET scanner.
    % (see spm_smooth).
    % 
    % 
%_______________________________________________________________________
% Copyright (C) 2015
%
% Gabriel Gonzalez-Escamilla
% $Id: geg_PVEcorrection.m 001 2015-03-02 10:16:23Z $
% 
% 
% rev = '$Rev: 021 $'; % 30-October-2015

if nargin == 1
    S = job.PETdata;
    tim = fieldnames(job.SegImgs);
    tim = tim{1};
    switch tim
        case 'Tsegs'
            N1 = job.SegImgs.Tsegs.tiss1;
            N2 = job.SegImgs.Tsegs.tiss2;
            CSFin = fieldnames(job.PVEopts.CSFsignal);
            CSFin = CSFin{1};
            if strncmp(CSFin,'CSFcalc',7)
                N3 = job.PVEopts.CSFsignal.CSFcalc.tiss3;
                csfzeroing = 0;
            elseif strncmp(CSFin,'CSFzeroing',10)
                csfzeroing = 1;
            end
            labels = 0; 
        case 'lblsegs'
            N = job.SegImgs.lblsegs.pvelbl;
            CSFin = fieldnames(job.PVEopts.CSFsignal);
            CSFin = CSFin{1};
            if strncmp(CSFin,'CSFcalc',7)
                csfzeroing = 0;
            elseif strncmp(CSFin,'CSFzeroing',10)
                csfzeroing = 1;
            end
            labels = 1;
    end
	gmthresh = job.PVEopts.gmthresh;
    PETpsf = job.PVEopts.fwhm_PSF;
    save_conv = job.PVEopts.TissConv;    
    dfn = fieldnames(job.PVE_Const_opts);
    dfn = dfn{1};
    switch dfn
        case 'type1' 
            % 'User pre-estimated values'
            AWM=job.PVE_Const_opts.type1.WMv;
            ACSF=job.PVE_Const_opts.type1.CSFv;
        case 'type2'
            % 'Deep WM Binary Mask in standard space'
            WMmasks = job.PVE_Const_opts.type2.WMmasks;
            djob.out{1}.pull = struct('fnames',{{}},'savedir',{{}},'interp',{{}},'mask',{{}},'fwhm',{{}});
            djob.out{1}.pull.mask   = 0;
            djob.out{1}.pull.interp = 0;
            djob.out{1}.pull.fwhm   = [0 0 0];
            deftype = fieldnames(job.PVE_Const_opts.type2.InDefs);
            deftype = deftype{1};
            switch deftype
                case 'defs'
                    InDefs = job.PVE_Const_opts.type2.InDefs.defs;
                case 'dartel'
                    Inflowfields = job.PVE_Const_opts.type2.InDefs.dartel.flowfield;
                    DRTLtemplate = {''}; 
                    DRTLsteps    = job.PVE_Const_opts.type2.InDefs.dartel.k;
            end      
        case 'type3'
            % 'Thresholding WM/CSF tissue segments'
            wmcsfthresh = job.PVE_Const_opts.type3.wmcsfthresh;
        case 'type4'
            % 'Eroding the WM/CSF tissue segments'
            erothresh = job.PVE_Const_opts.type4.erothresh;
        case 'type5'
            % 'Using values from GTM (mMG)'
            GTMfiles = job.PVE_Const_opts.type5.GTMtxt;
        otherwise
            fprintf('Not a recognized option!\n')
            return;
    end
end

if nargin < 1
    if strcmp(spm('ver'),'SPM2')
        fprintf(2,'NOT currently available for this SPM version. \n');
    else
        S = spm_select(Inf,'image','Select native coregistered PET images');
        labels = spm_input('Type of segmented image(s)','1','b',{'Tissue_Segments','Tissue-labeled'},[0,1]);
        if labels
            N  = spm_select(size(S,1),{'^c0.*\.nii','^p0.*\.nii'},'Select native tissue-labeled image');
        else
            N1 = spm_select(size(S,1),{'^c1.*\.nii','^p1.*\.nii'},'Select native GM maps');
            N2 = spm_select(size(S,1),{'^c2.*\.nii','^p2.*\.nii'},'Select native WM maps');
            N3 = spm_select(size(S,1),{'^c3.*\.nii','^p3.*\.nii'},'Select native CSF maps');
        end
    end
    gmthresh = spm_input('GM threshold','1','e',0,1);
    csfzeroing = spm_input('CSF-zeroing?','+1','b',{'yes','no'},[1,0]);
    PETpsf = spm_input('PET PSF (FWHM in mm)','+1','e',[6,6,6],[1,3]);
    save_conv = spm_input('Save convolved tissues?','+1','b',{'yes','no'},[1,0]);    
    constActMthds = {'User pre-estimated WM&CSF values ','Deep WM Binary Mask in standard space ',...
                'Deep WM Binary Masks in subject''s space ','Thresholding WM/CSF tissue segments ',...
                'Eroding the WM/CSF tissue segments'};
    constActSel = spm_input('WM/CSF constant activity','+1','m',constActMthds);
    constActTypes = {'type1','type2','type3','type4','type5'};
    dfn = constActTypes{constActSel};    
    switch dfn
        case 'type1';
            AWM = spm_input('WM activity vector','+1','e','',[1,Inf]);
            ACSF = spm_input('WM activity vector','+1','e','',[1,Inf]);
        case 'type2';
            WMmasks = {spm_select([1,1],'^*\.nii','Select deep WM map in standard space ')};
            djob.out{1}.pull = struct('fnames',{{}},'savedir',{{}},'interp',{{}},'mask',{{}},'fwhm',{{}});
            djob.out{1}.pull.mask   = 0;
            djob.out{1}.pull.interp = 0;
            djob.out{1}.pull.fwhm   = [0,0,0];
            deftype = spm_input('Input Fields type','+1','b',{'Deformations','FlowFields'},{'defs','dartel'});
            deftype = deftype{1};
            InDefs1  = spm_select(size(S,1),'^*\.nii','Select deformation- or flow-fields for every subject');            
            switch deftype
                case 'defs'
                    InDefs = InDefs1;
                case 'dartel'
                    Inflowfields = InDefs1;
                    DRTLtemplate = {''}; 
                    DRTLsteps    = 6;
            end
        case 'type3';
            wmcsfthresh = spm_input('WM/CSF threshold','+1','e',0,1);
            
        case 'type4';
            erothresh = spm_input('WM erode threshold','+1','e',0,1);  
        case 'type5';
           GTMfiles = spm_select([1,inf],'^*\.txt','Select GTM output txt file(s)');
    end
end
erofwhm = geg_petpve12_get_defaults('PVEopts.Erofwhm');
save_segs = 1;

% Running for every subject
for subj=1:size(S,1)
    
    % get current release number
    A = ver;
    for i=1:length(A)
        if strcmp(A(i).Name,'PET Partial Volume Effects correction Toolbox')
            r = str2double(A(i).Version);
        end
    end    
    if exist('r','var')
        fprintf('PETPVE12 r%d: %s\n',r,char(S(subj,:)));
    end
    % Set starting time
    fprintf('Starting time: %s \n',datestr(clock,'local'))
    
    % Loading data
    actPT = char(S(subj,:));% corregistered PET name
    PTm = spm_vol(actPT);% PET id matrix
    PTi = spm_read_vols(PTm);% PET image
    PTi(isnan(PTi)) = 0; %making sure there are no NaN's in PET
    [pth,nam,ext] = spm_fileparts(actPT);
    odir = pth;
    if strncmpi(ext,'.nii',4)
        cvt = ['pvc' nam ext];
    else
        cvt = ['pvc' nam '.nii'];
    end    
    if ~labels
        actGM = char(N1(subj,:));% unsmoothed GM name
        actWM = char(N2(subj,:));% unsmoothed WM name
        GMm = spm_vol(actGM);% GM header
        GMi = spm_read_vols(GMm);% GM image
        WMm = spm_vol(actWM);% WM id header
        WMi = spm_read_vols(WMm);% WM image
        if ~csfzeroing
            actCSF = char(N3(subj,:));% unsmoothed CSF name
            CSFm = spm_vol(actCSF);% CSF id header
            CSFi = spm_read_vols(CSFm);% CSF image
        end
    else
        actlbl = char(N(subj,:));% unsmoothed labeled image name
        [p,n,e] = spm_fileparts(actlbl);
        lblm = spm_vol(actlbl);
        lbli = spm_read_vols(lblm);% unsmoothed labeled image
        GMi = single(round(lbli)==2);
        WMi = single(round(lbli)==3);
        CSFi = single(round(lbli)==1);
        GMm = lblm;
        GMm.fname = fullfile(p,['c1' n(3:end),e]);
        WMm = lblm;
        WMm.fname = fullfile(p,['c2' n(3:end),e]);
        if save_segs
            spm_write_vol(GMm,GMi);
            spm_write_vol(WMm,WMi);
        end
        if ~csfzeroing
            CSFm = lblm;
            CSFm.fname = fullfile(p,['c3' n(3:end),e]);
            if save_segs, spm_write_vol(CSFm,CSFi); end
        end
    end
        
    % convolving data
    if length(PETpsf)==1
        fprintf('Convolving tissue segments with isotropic PSF. \n')
        smth=PETpsf(1);
    else
        if PETpsf(1)==PETpsf(2) && PETpsf(1)==PETpsf(3)
            fprintf('Convolving tissue segments with isotropic PSF of %s mm. \n',num2str(PETpsf(1)))
            smth=PETpsf;
        else
            fprintf('Convolving tissue segments with non-isotropic PSF of %s mm. \n',num2str(PETpsf))
            smth=PETpsf;
        end
    end
    [pth,nam,ext,~] = spm_fileparts(GMm.fname);
    QGM = fullfile(pth,['s' nam,ext]); 
    cGMm = GMm; cGMm.fname = QGM; cGMm.descrip = 'smoothed';
    cGMi = zeros(size(GMi));
    spm_smooth(GMm,cGMi,smth,0); % create PSF convolved GM image
    %
    [pth,nam,ext,~] = spm_fileparts(WMm.fname);
    QWM = fullfile(pth,['s' nam,ext]); 
    cWMm = WMm; cWMm.fname = QWM; cWMm.descrip = 'smoothed';
    cWMi = zeros(size(WMi));
    spm_smooth(WMm,cWMi,smth,0);% create PSF convolved WM image
    %
    if ~csfzeroing
        [pth,nam,ext,~] = spm_fileparts(CSFm.fname);
        QCSF = fullfile(pth,['s' nam,ext]);
        cCSFm = CSFm; cCSFm.fname = QCSF; cCSFm.descrip = 'smoothed';
        cCSFi = zeros(size(CSFi));
        spm_smooth(CSFm,cCSFi,smth,0);% create PSF convolved CSF image
    end
    %
    % it is mandatory to save the convolved images to account for dimension differences
    spm_write_vol(cGMm,cGMi);
    spm_write_vol(cWMm,cWMi);
    if ~csfzeroing, spm_write_vol(cCSFm,cCSFi); end
    
    % Computation of "Activity" in WM/CSF
    switch dfn
        case 'type1'
            % 'User estimated values'
            x = AWM(subj);% PET constant activity in WM
            y = ACSF(subj);% PET constant activity in CSF
        case 'type2'
            % 'Deep WM Binary Masks in standard space'
            fprintf('Deforming mask to %s,. \n',actPT)
            djob.out{1}.pull.fnames = WMmasks;
            djob.out{1}.pull.saveusr = odir;
            switch deftype
                case 'defs'
                    djob.comp{1}.def = InDefs(subj); % Name of the deformation field (iy_*.nii)
                case 'dartel'
                    djob.comp{1}.dartel = struct('flowfield',{{}},'template',{{}},'times','K',{{}});
                    djob.comp{1}.dartel.flowfield = Inflowfields(subj);% Name of the flow_field (u_*.nii)
                    djob.comp{1}.dartel.template  = DRTLtemplate;
                    djob.comp{1}.dartel.times     = [0 1];
                    djob.comp{1}.dartel.K         = DRTLsteps;
                    djob.comp{1}.dartel = orderfields(djob.comp{1}.dartel, {'flowfield','template','times','K'});
                    djob.comp{2}.id = struct('space',{{}});% identity matrix
                    if ~labels
                       djob.comp{2}.id.space = {deblank(actGM)}; 
                    else
                       djob.comp{2}.id.space = {deblank(actlbl)};
                    end
            end
            dmask = geg_iydeformation(djob);
            dmaskm = spm_vol(dmask.warpedName{1});
            dmaski = spm_read_vols(dmaskm);% deformed WM mask
            delete(dmaskm.fname)
            [~,n,~] = spm_fileparts(dmaskm.fname);
            dmaskm.fname = fullfile(odir,[n '_' cvt(5:end-4) ext]);
            spm_write_vol(dmaskm,dmaski);
            %mpos = dmaski~=0; AWM = PTi(mpos);
            %x = mean(AWM);% PET constant activity in WM
            AWM = geg_calc_constants(PTm,dmaskm,0);
            x = AWM;% PET constant activity in WM
            if csfzeroing
                y = 0;% PET constant activity in CSF
            else
                %CSFe = geg_erodeMask(CSFm,.99,2,1);
                %ACSF = PTi.*CSFe;
                %y = mean(ACSF~=0);% PET constant activity in CSF
                geg_erodeMask(CSFm,.99,2,1);
                [epth,enam,eext] = spm_fileparts(CSFm.fname);
                eCSFm = CSFm; eCSFm.fname = fullfile(epth,['e' enam, eext]);% '_eroded'
                ACSF = geg_calc_constants(PTm,CSFm,0);
                y = ACSF;% PET constant activity in CSF
            end
        case 'type3'
            % 'Thresholding WM/CSF tissue segments' 
            AWM = geg_calc_constants(PTm,WMm,wmcsfthresh);
            x = AWM;% PET constant activity in WM
            if csfzeroing
                y = 0;% PET constant activity in CSF
            else
                ACSF = geg_calc_constants(PTm,CSFm,wmcsfthresh);
                y = ACSF;% PET constant activity in CSF
            end            
        case 'type4'            
            % In case you want to erode the WM or WM roi
            fprintf('Eroding WM tissue for %s,. \n',actPT)
            %WMe = geg_erodeMask(WMm,erothresh,erofwhm,1);            
            %AWM = PTi(WMe>0);
            %x = mean(AWM);% PET constant activity in WM
            geg_erodeMask(WMm,erothresh,erofwhm,1);
            [epth,enam,eext] = spm_fileparts(WMm.fname);
            eWMm = WMm; eWMm.fname = fullfile(epth,['e' enam, eext]);% '_eroded'
            AWM = geg_calc_constants(PTm,eWMm,0);
            x = ACSF;% PET constant activity in CSF
            if csfzeroing
                y = 0;% PET constant activity in CSF
            else
                %CSFe = geg_erodeMask(CSFm,.9,erofwhm,1);
                %ACSF = PTi(CSFe~=0);
                %y = mean(ACSF);% PET constant activity in CSF
                geg_erodeMask(CSFm,.9,erofwhm,1);
                [epth,enam,eext] = spm_fileparts(CSFm.fname);
                eCSFm = CSFm; eCSFm.fname = fullfile(epth,['e' enam, eext]);% '_eroded'
                ACSF = geg_calc_constants(PTm,CSFm,0);
                y = ACSF;% PET constant activity in CSF
            end
        case 'type5'
            gtmVals = importdata(GTMfiles{subj});
            x = gtmVals.data(end-1);
            y = gtmVals.data(end);
    end
    
    fprintf('Partial volume correction as Muller-Gartner on %s \n',actPT)
    % consider any possible problem with image dimensions (i.e. coregistered, but not resliced already)
    cGMm1 = struct2cell(cGMm); cWMm1 = struct2cell(cWMm); 
    PTm1 = struct2cell(PTm); GMm1 = struct2cell(GMm);
    if ~csfzeroing 
        cCSFm1 = struct2cell(cCSFm);
        Vi = [cGMm1, cWMm1, cCSFm1, PTm1, GMm1];% new order of inputs to always match MRI space
    else
        Vi = [cGMm1, cWMm1, PTm1, GMm1];% new order of inputs to always match MRI space
    end
    Headings = {'fname'; 'dim'; 'dt'; 'pinfo'; 'mat'; 'n'; 'descrip'; 'private'};
    Vi = cell2struct(Vi, Headings, 1);
    Vo = PTm.fname; 
    interp = 1;%1=trilinear; 0=nearest;
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
        n = numel(Vi);
        pvec = zeros(Vo.dim(1:3));
        if labels && ~csfzeroing
            % use total gm mask, as it will not have overlapping with the WM or CSF 
            f = '(i4 - (x*i2) - (y*i3)).*(i5)./i1';%'(PTi - (x*cWMi) - (y*cCSFi)).*(GMi)./cGMi';
        elseif labels && csfzeroing
            % use total gm mask, as it will not have overlapping with the WM or CSF 
            f = '(i3 - (x*i2)).*(i4)./i1';%'(PTi - (x*cWMi)).*(GMi)./cGMi';
        elseif ~labels && ~csfzeroing
            % apply GM-threshold to avoid tissue overlapping
            f = '(i4 - (x*i2) - (y*i3)).*(i5>gmthresh)./i1';%'(PTi - (x*cWMi) - (y*cCSFi)).*(GMi>gmthresh)./cGMi';
        elseif ~labels && csfzeroing
            % apply GM-threshold to avoid tissue overlapping
            f = '(i3 - (x*i2)).*(i4>gmthresh)./i1';%'(PTi - (x*cWMi) - (y*cCSFi)).*(GMi>gmthresh)./cGMi';
        end        
        
        %-Start progress plot
        spm_progress_bar('Init',Vo.dim(3),f,'planes completed');
        
        %-Loop over planes computing result Y
        for p = 1:Vo.dim(3)
            B = spm_matrix([0 0 -p 0 0 0 1 1 1]);
            for i = 1:n
                M = inv(B * inv(Vo.mat) * Vi(i).mat);
                d = spm_slice_vol(Vi(i), M, Vo.dim(1:2), [interp,NaN]);
                eval(['i',num2str(i),'=d;']);
            end            
            try
                eval(['Yp = ' f ';']);
            catch ME1
                % Get last segment of the error message identifier.
                idSegLast = regexp(ME1.identifier, '(?<=:)\w+$', 'match');
                fprintf(2,['\n can''t evaluate %s because %s \n',f,idSegLast]);
            end
            if prod(Vo.dim(1:2)) ~= numel(Yp)
                error(['"',f,'" produced incompatible image.']); end
            pvec(:,:,p) = reshape(Yp,Vo.dim(1:2));
            spm_progress_bar('Set',p);
        end
        spm_progress_bar('Clear')
    else
        % here images do not need to be resliced
        if labels && ~csfzeroing
            % use total gm mask, as it will not have overlapping with the WM or CSF 
            pvec = (PTi - (x*cWMi) - (y*cCSFi)).*(GMi)./cGMi;
        elseif labels && csfzeroing
            pvec = (PTi - (x*cWMi)).*(GMi)./cGMi;
        elseif ~labels && ~csfzeroing
            % apply GM-threshold to avoid tissue overlapping
            pvec = (PTi - (x*cWMi) - (y*cCSFi)).*(GMi>gmthresh)./cGMi;
        elseif ~labels && csfzeroing
            pvec = (PTi - (x*cWMi)).*(GMi>gmthresh)./cGMi;
        end        
    end
    
    % writing data 
    fprintf('saving: %s \n',cvt)
    if exist('Vo','var')
        oPT = Vo;
        oPT.fname   = fullfile(odir, cvt);
        oPT.descrip = 'PETPVE12 - pvcMG';
    else
        nfname=fullfile(odir, cvt);
        oPT = PTm;
        oPT.fname   = nfname;
        oPT.descrip = 'PETPVE12 - pvcMG';
    end
    pvec(isnan(pvec)) = 0; % removing NaNs
    spm_write_vol(oPT,pvec);
    
    % save the convolved tissues
    if ~save_conv
        delete(QGM, QWM)
        if ~csfzeroing, delete(QCSF); end
    end
    clear PTi GMi cGMi cWMi cCSFi Vi Vo sts str refstr oPT pvec
    
    % Set end time
    fprintf('End time: %s \n',datestr(clock,'local'))
    
    fprintf('\n')
end

