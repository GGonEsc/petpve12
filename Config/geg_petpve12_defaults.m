function geg_petpve12_defaults
% Sets the defaults for PET-PVE toolbox
% FORMAT geg_petpve12_defaults
%_______________________________________________________________________
%
% This file is intended to be customised for the site.
%
% Care must be taken when modifying this file
%_______________________________________________________________________
% $Id: geg_petpve12_defaults.m 004 2015-03-03 08:58:11Z $

%rev = '$Rev: 020 $'; 04-May-2016

global petpve

% Estimation options
%=======================================================================
petpve.opts.tpm       = {fullfile(spm('dir'),'tpm','TPM.nii')}; % TPM.nii
petpve.opts.ngaus     = [2 2 2 3 4 2];  % Gaussians per class
petpve.opts.affreg    = 'mni';    % Affine regularisation for segmentation
petpve.opts.warpreg   = [0 0.001 0.5 0.05 0.2];  % Warping regularisation
petpve.opts.biasreg   = 0.0001; % Bias regularisation
petpve.opts.biasfwhm  = 60;   % Bias FWHM
petpve.opts.samp      = 3;    % Sampling distance

% Writing options
%=======================================================================
% segmentations:
%   native    0/1   (none/yes)
petpve.output.bias.native  = 1; % bias corrected (m0) 
petpve.output.label.native = 1; % labeled image (c0*)
petpve.output.GM.native    = 1;  % GM (c1*)
petpve.output.WM.native    = 1;  % WM (c2*)
petpve.output.CSF.native   = 1; % CSF (c3*)
petpve.output.BM.native    = 1; % binary brain mask (bm*)
petpve.output.SS.native    = [0 1];% skull-strip image [T1 Bias]
% [1 0]; % Original T1 skull-stripped (ss0*)
% [0 1]; % Bias corrected skull-stripped (ss1*)

petpve.output.warps        = [0 1];% order is [forward inverse]

% Extended writing options
%=======================================================================
petpve.extopts.darteltpm   = {fullfile(spm('dir'),'toolbox','petpve12','Atlases','IXI550','Template_1_IXI550_MNI152.nii')};% Template
petpve.extopts.print     = 1; % Display and print results

% bias correction options
%=======================================================================
petpve.bias.nits_bias    = 8;
petpve.bias.biasfwhm     = 60;
petpve.bias.biasreg      = 1e-6;
petpve.bias.lmreg        = 1e-6;

% Skull-stripping using GM/WM tissues options
%=======================================================================
petpve.SSopts.morph      = 1; % use morphological operations for skull-stripping
petpve.SSopts.finalmask  = 0; % Generate a more dilated mask to fill remaining holes
petpve.SSopts.thresh     = 0.5; % Threshold for brain mask from precomputed GM/WM segments

% PVE correction options (MG/mMG & GTM)
%=======================================================================
petpve.PVEopts.PSF            = [6 6 6]; % PET point spread function
petpve.PVEopts.GMthr          = 0.5; % Gray matter threshold to MG output
petpve.PVEopts.WMCSFthr       = 0.9; % Threshold to compute PET signal in WM/CSF
petpve.PVEopts.GTMtissthr     = 0.5; % Threshold to produce non-overlapping tissues for GTM inclusion
petpve.PVEopts.CSFzeroing     = 0; % consider CSF signal equal to zero, in GTM only affects the signal vector not the GTM
petpve.PVEopts.TissConv       = 0; % Save the convolved tissue segments
petpve.PVEopts.EroThresh      = 0.99; % Threshold for WM tissue erosion
petpve.PVEopts.Erofwhm        = 2;% Gaussian weighting function for tissue erosion
petpve.PVEopts.RoiAtlas       = {fullfile(spm('dir'),'toolbox','petpve12','Atlases','Desikan-Killiany_MNI_SPM12.nii')};
petpve.PVEopts.RoiAtlasDesc   = {fullfile(spm('dir'),'toolbox','petpve12','Atlases','Desikan-Killiany_MNI_SPM12.txt')};
petpve.PVEopts.RoiAtlasInterp = 0; % Nearest neighbour
petpve.PVEopts.GTMsavesegs    = 1; % Saves a bunch of user useless images
petpve.PVEopts.GTMsaveextras  = 0; % Saves a bunch of user useless images
petpve.PVEopts.GTMsaveGTM     = 1; % Saves the GTM, uncorrected observed ROI activity and ROI sizes text files

% PET intensity-normalization options
%=======================================================================
petpve.PETnormOpts.MNIroi = {fullfile(spm('dir'),'toolbox','petpve12','Atlases','Desikan-Killiany_MNI_cerebellum.nii')};

% PET WM-ROIs creation options
%=======================================================================
petpve.WMrois.RoiAtlas       = {fullfile(spm('dir'),'toolbox','petpve12','Atlases','AAL_MNI_SPM12.nii')};
petpve.WMrois.RoiAtlasDesc   = {fullfile(spm('dir'),'toolbox','petpve12','Atlases','AAL_MNI_SPM12.txt')};
petpve.WMrois.savesegs    = 1; % Saves a bunch of user useless images
petpve.WMrois.saveextras  = 0; % Saves a bunch of user useless images

% expert options for segmentation and Skull-Strippping
%=======================================================================
petpve.extopts.cleanup     = 1;    % Cleanup: 0 - no; 1 - light; 2 -thorough
petpve.extopts.finalmask   = 1;    % Final masking: 0 - no; 1 - yes
petpve.extopts.gcut        = 1;    % Skull-stripping with graph-cut: 0 - no; 1 - yes
petpve.extopts.kmeans      = 1;    % segmentation initialization: 0 - new segment; 1 - Kmeans
petpve.extopts.mrf         = 0.3;  % MRF weighting 0.15
petpve.extopts.sanlm       = 2;    % use SANLM filter: 0 - no SANLM; 1 - SANLM with single-threading; 2 - SANLM with multi-threading
petpve.extopts.bias_fwhm   = 60;   % FWHM of Kmeans internal bias correction
petpve.extopts.histeq_deep = 1;    % weighting of local histogram equalization: 0 - no; 1 - full weighting (not recommended)
petpve.extopts.histeq_mask = {fullfile(spm('dir'),'toolbox','petpve12','Atlases','IXI550','histeq_mask.nii')};

% realign options (longitudinal alignment)
%=======================================================================
petpve.realign.halfway    = 1; % use halfway registration: 0 - no; 1 - yes
petpve.realign.weight     = 1; % weight registration with inverse std: 0 - no; 1 - yes
petpve.realign.ignore_mat = 0; % ignore exisiting positional information: 0 - no; 1 - yes

% Coregistration defaults
%==========================================================================
petpve.coreg.estimate.cost_fun = 'nmi';
petpve.coreg.estimate.sep      = [4 2];
petpve.coreg.estimate.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
petpve.coreg.estimate.fwhm     = [7 7];
petpve.coreg.write.interp      = 4;
petpve.coreg.write.wrap        = [0 0 0];
petpve.coreg.write.mask        = 0;
petpve.coreg.write.prefix      = 'r';