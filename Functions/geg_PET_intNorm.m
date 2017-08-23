function geg_PET_intNorm(job)
% This function performs the scaling of tracer uptake to a reference
% region, this is ussually called as intensity normalization of PET images.
% An ideal reference region should not be affected by brain pathology and
% should be easy to image/analyse. 
% However, the choice of the appropriate reference region is especially is
% a fundamental methodological issue and is particullar problematic in
% subjects with neurodegenerative disorders who show early metabolic and
% perfusion deficits. 
%_______________________________________________________________________
% Copyright (C) 2015
%
% Gabriel Gonzalez-Escamilla
% $Id: geg_PET_intNorm.m 001 2015-04-09 12:04:36Z $
% 

% rev = '$Rev: 005 $'; % 23-August-2017

if nargin == 1
    S = job.data;
    dfn = fieldnames(job.PETnorm_opts);
    dfn = dfn{1};
    switch dfn
        case 'type1' 
           % 'Parcellation in subject space'
            SubjParcAtlas = job.PETnorm_opts.type1.RefMask;
            % As atlas is the same for all subjects it assumes a single txt
            % file containing 2 columns (names and assigned-value)
            isvol = true;
        case 'type2'
            % 'Atlas in standard space'
            djob.out{1}.pull = struct('fnames',{{}},'savedir',{{}},'interp',{{}},'mask',{{}},'fwhm',{{}});
            djob.out{1}.pull.mask   = 0;
            djob.out{1}.pull.interp = 0;
            djob.out{1}.pull.fwhm   = [0 0 0];
            deftype = fieldnames(job.PETnorm_opts.type2.InDefs);
            deftype = deftype{1};
            isvol = true;
        case 'type3'    
            % User defined values
            SubjsRefAct = job.PETnorm_opts.type3.RefVals;
            isvol = false;
        otherwise
            fprintf('Nothing else implemented!\n')
            return;
    end
    save_subj_defROI = 0;
else
    error('Not yet implemented')
end

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
    
    % Loading data
    actPT = char(S(subj,:));% corregistered PET name
    PTm = spm_vol(actPT);% PET id matrix
    PTi = spm_read_vols(PTm);% PET image
    [pth,nam,ext] = spm_fileparts(actPT);
    odir = pth;
    if strncmpi(ext,'.nii',4)
        cvt = ['n' nam ext];
    else
        cvt = ['n' nam '.nii'];
    end
    
    switch dfn
        case 'type1'
            % 'ROI in subject space'
            rmask = SubjParcAtlas{subj};
            rmaskm = spm_vol(rmask);
            % rmaski = spm_read_vols(rmaskm);
        case 'type2'
            % 'ROI atlas in standard space'
            fprintf('Deforming template atlas to %s,. \n',actPT)
            djob.out{1}.pull.fnames = job.PETnorm_opts.type2.RefMask;
            djob.out{1}.pull.savedir.saveusr = {odir};
            switch deftype
                case 'defs'
                    djob.comp{1}.def = job.PETnorm_opts.type2.InDefs.defs(subj); % Name of the deformation field (iy_*.nii)
                case 'dartel'
                    djob.comp{1}.dartel = job.PETnorm_opts.type2.InDefs.dartel;
                    djob.comp{1}.dartel = rmfield(djob.comp{1}.dartel,'flowfield');
                    djob.comp{1}.dartel.flowfield = job.PETnorm_opts.type2.InDefs.dartel.flowfield(subj);% Name of the flow_field (u_*.nii)
                    djob.comp{1}.dartel.times = [1 0];% always perform 'Backward';
                    djob.comp{1}.dartel = orderfields(djob.comp{1}.dartel, {'flowfield','template','times','K'});
                    djob.comp{2}.id = struct('space',{{}});% identity matrix
                    djob.comp{2}.id.space = {deblank(actPT)}; 
            end
            drmask = geg_iydeformation(djob);
            rmaskm = spm_vol(drmask.warpedName{1});
            % rmaski = spm_read_vols(rmaskm);% deformed mask
            delete(rmaskm.fname)
            if save_subj_defROI
                [~,n,~] = spm_fileparts(rmaskm.fname);
                rmaskm.fname = fullfile(odir,[n '_' cvt(5:end-4) ext]);
                spm_write_vol(rmaskm,rmaski);
            end
        case 'type3'  
            maskact = SubjsRefAct(subj);
    end
    
    if isvol
        % make sure dimensions are equal (fix it if not)
        [rmaski,~] = geg_reslice(PTm,rmaskm,0);
    end
    
    % compute reference activity
    if isvol
        maskact = PTi.*rmaski;
        maskact = mean(maskact(:));
    end
    normPET = PTi ./ maskact;
    
    % writing data 
    % Save intensity corrected PET image
    fprintf('saving: %s \n', cvt)
    % create a new header to avoid strange outputs
    oPT.Nt      = nifti;
    oPT.Nt.dat  = file_array(fullfile(odir,cvt),...
        PTm.dim(1:3),...
        [spm_type('int32') spm_platform('bigend')],...
        0,1/255,0);
    oPT.Nt.mat  = PTm.mat;
    oPT.Nt.mat0 = PTm.mat;
    oPT.Nt.descrip = 'Intensity Normalized';
    create(oPT.Nt);
    ind  = PTm.n;
    oPT.Nt.dat(:,:,:,ind(1),ind(2))=double(normPET);
    
    fprintf('\n')
end
