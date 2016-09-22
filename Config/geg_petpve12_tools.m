function tools = geg_petpve12_tools
% wrapper for calling PET-PVC extra utilities
%
%_______________________________________________________________________
% Gabriel González-Escamilla
% $Id: geg_petpvc12_tools.m 001 2015-03-08 16:20:04Z $

% rev = '$Rev: 005 $'; 13-April-2015

%_______________________________________________________________________

data      = cfg_files;
data.tag  = 'data';
data.name = 'Volumes';
data.filter  = 'image';
data.ufilter = '.*';
data.num     = [1 Inf];

%========================================================================
%========================================================================
%------------------------------------------------------------------------
% Check sample homogeneity using covariance
%------------------------------------------------------------------------
%========================================================================
%========================================================================
% if any square has the same value as other, then the image is repeated (i.e. is the same subject) 
data.help = {[...
'Select all images. Images have to be in the same orientation with same voxel size and dimensions ',...
'(e.g. normalized images)']};

c         = cfg_entry;
c.tag     = 'c';
c.name    = 'Vector';
c.help    = {'Vector of nuisance values'};
c.strtype = 'e';
c.num     = [Inf 1];

nuisance       = cfg_branch;
nuisance.tag   = 'nuisance';
nuisance.name  = 'Nuisance';
nuisance.val   = {c};
nuisance.help  = {'Add a nuisance parameter to be removed from data'};

slice      = cfg_entry;
slice.tag  = 'slice';
slice.name = 'Show slice (in mm)?';
slice.strtype = 'e';
slice.num  = [1 1];
slice.val  = {0};
slice.help = {'Choose slice in mm.'};

gap      = cfg_entry;
gap.tag  = 'gap';
gap.name = 'Gap to skip slices';
gap.strtype = 'e';
gap.num  = [1 1];
gap.val  = {5};
gap.help = {['',...
'To speed up calculations you can define that only every x slice the covariance is estimated.']};

scale      = cfg_menu;
scale.tag  = 'scale';
scale.name = 'Proportional scaling?';
scale.labels = {'no','yes'};
scale.values = {0 1};
scale.val    = {0};
scale.help   = {[...
'This option should be only used if image intensity is not scaled (e.g. T1 images) ',...
'or if images have to be scaled during statistical analysis (e.g. modulated images).']};

transform         = cfg_repeat;
transform.tag     = 'transform';
transform.name    = 'Nuisance';
transform.help    = {'This option allows for the specification of nuisance effects to be removed from the data. ',...
'A potential nuisance parameter can be age. In this case the variance explained by age will be removed prior to ',...
'the calculation of the covariance.'};
transform.values  = {nuisance};
transform.num     = [0 Inf];

check_cov      = cfg_exbranch;
check_cov.tag  = 'check_cov';
check_cov.name = 'Check sample homogeneity using covariance';
check_cov.val  = {data,scale,slice,gap,transform};
check_cov.prog = @geg_check_cov;
check_cov.help = {[...
'If you have a reasonable sample size artefacts are easily overseen. In order to identify images with poor image quality ',...
'or even artefacts you can use this function. Images have to be in the same orientation with same voxel size and dimension ',...
'(e.g. normalized images). The idea of this tool is to check the covariance of all files across the sample.'],...
'',[...
'The covariance is calculated between all images and the mean for each image is plotted using a boxplot and the indicated ',...
'filenames. The smaller the mean covariance the more deviant is this image from the sample mean. ',...
'In the plot outliers from ',...
'the sample are usually isolated from the majority of images which are clustered around the sample mean. The mean ',...
'covariance is plotted at the y-axis and the x-axis reflects the image order. Images are plotted from left to right which is helpful if ',...
'you have selected the images in the order of different sub-groups. Furthermore this is also useful for fMRI images which can be ',...
'also used with this tool. The proportional scaling option should be only used if image intensity is not scaled (e.g. T1 images) ',...
'or if images have to be scaled during statistical analysis (e.g. modulated images).']};

%========================================================================
%========================================================================
%------------------------------------------------------------------------
% Estimate coregistration parameters for MRI-PET data
%------------------------------------------------------------------------
%========================================================================
%========================================================================
%--------------------------------------------------------------------------
% ref Reference Images
refs         = cfg_files;
refs.tag     = 'refs';
refs.name    = 'Reference Images';
refs.help    = {'This is the image that is assumed to remain stationary (sometimes known as the target or template image), normally a structural image.'};
refs.filter  = 'image';
refs.ufilter = '.*';
refs.num     = [1 Inf];
%--------------------------------------------------------------------------
% source Source Images
sources         = cfg_files;
sources.tag     = 'sources';
sources.name    = 'Source Images';
sources.help    = {'This is the image that is jiggled about (moved) to best match the reference. One image per reference is expected. Normally, the functional image'};
sources.filter  = 'image';
sources.ufilter = '.*';
sources.num     = [1 Inf];
%--------------------------------------------------------------------------
% other Other Images
others         = cfg_files;
others.tag     = 'others';
others.name    = 'Other Images';
others.val     = {{''}};
others.help    = {'These are any images that need to remain in alignment with the source image. Only one image per reference is expected'};
others.filter  = 'image';
others.ufilter = '.*';
others.num     = [0 Inf];
%--------------------------------------------------------------------------
% cost_fun Objective Function
cost_fun         = cfg_menu;
cost_fun.tag     = 'cost_fun';
cost_fun.name    = 'Objective Function';
cost_fun.help    = {'Registration involves finding parameters that either maximise or minimise some objective function. For inter-modal registration, use Mutual Information/* \cite{collignon95,wells96}*/, Normalised Mutual Information/* \cite{studholme99}*/, or Entropy Correlation Coefficient/* \cite{maes97}*/.For within modality, you could also use Normalised Cross Correlation.'};
cost_fun.labels  = {'Mutual Information'
                    'Normalised Mutual Information'
                    'Entropy Correlation Coefficient'
                    'Normalised Cross Correlation'}';
cost_fun.values  = {'mi'
                    'nmi'
                    'ecc'
                    'ncc'}';
cost_fun.def     = @(val)geg_petpve12_get_defaults('coreg.estimate.cost_fun', val{:});
%--------------------------------------------------------------------------
% sep Separation
sep         = cfg_entry;
sep.tag     = 'sep';
sep.name    = 'Separation';
sep.help    = {'The average distance between sampled points (in mm).  Can be a vector to allow a coarse registration followed by increasingly fine ones.'};
sep.strtype = 'r';
sep.num     = [1 Inf];
sep.def     = @(val)geg_petpve12_get_defaults('coreg.estimate.sep', val{:});
%--------------------------------------------------------------------------
% tol Tolerances
tol         = cfg_entry;
tol.tag     = 'tol';
tol.name    = 'Tolerances';
tol.help    = {'The accuracy for each parameter.  Iterations stop when differences between successive estimates are less than the required tolerance.'};
tol.strtype = 'r';
tol.num     = [1 12];
tol.def     = @(val)geg_petpve12_get_defaults('coreg.estimate.tol', val{:});
%--------------------------------------------------------------------------
% fwhm Histogram Smoothing
fwhm         = cfg_entry;
fwhm.tag     = 'fwhm';
fwhm.name    = 'Histogram Smoothing';
fwhm.help    = {'Gaussian smoothing to apply to the 256x256 joint histogram. Other information theoretic coregistration methods use fewer bins, but Gaussian smoothing seems to be more elegant.'};
fwhm.strtype = 'r';
fwhm.num     = [1 2];
fwhm.def     = @(val)geg_petpve12_get_defaults('coreg.estimate.fwhm', val{:});
%--------------------------------------------------------------------------
% eoptions Estimation Options
eoptions         = cfg_branch;
eoptions.tag     = 'eoptions';
eoptions.name    = 'Estimation Options';
eoptions.val     = {cost_fun sep tol fwhm};
eoptions.help    = {'Various registration options, which are passed to the Powell optimisation algorithm/* \cite{press92}*/.'};
%--------------------------------------------------------------------------
% estimate Coreg: Estimate
coreg_estimate         = cfg_exbranch;
coreg_estimate.tag     = 'coreg_estimate';
coreg_estimate.name    = 'Coregister: Estimate';
coreg_estimate.val     = {refs sources others eoptions};
coreg_estimate.help    = {'The registration method used here is based on work by Collignon et al/* \cite{collignon95}*/. The original interpolation method described in this paper has been changed in order to give a smoother cost function.  The images are also smoothed slightly, as is the histogram.  This is all in order to make the cost function as smooth as possible, to give faster convergence and less chance of local minima.'
                    ''
                    'At the end of coregistration, the voxel-to-voxel affine transformation matrix is displayed, along with the histograms for the images in the original orientations, and the final orientations.  The registered images are displayed at the bottom.'
                    ''
                    'Registration parameters are stored in the headers of the "source" and the "other" images.'}';
coreg_estimate.prog = @geg_run_coreg;
coreg_estimate.vout = @vout_estimate;

%========================================================================
%========================================================================
%------------------------------------------------------------------------
% Estimate coregistration parameters for MRI-PET data then reslice
%------------------------------------------------------------------------
%========================================================================
%========================================================================
% source Source Image
sources.help    = {'This is the image that is jiggled about to best match the reference. One image per reference is expected. ',...
    'These images are resliced to the same dimensions, voxel sizes, orientation etc as the space defining image.'};
%--------------------------------------------------------------------------
% interp Interpolation
interp         = cfg_menu;
interp.tag     = 'interp';
interp.name    = 'Interpolation';
interp.help    = {'The method by which the images are sampled when being written in a different space. Nearest Neighbour is fastest, but not normally recommended. It can be useful for re-orienting images while preserving the original intensities (e.g. an image consisting of labels). Trilinear Interpolation is OK for PET, or realigned and re-sliced fMRI. If subject movement (from an fMRI time series) is included in the transformations then it may be better to use a higher degree approach. Note that higher degree B-spline interpolation/* \cite{thevenaz00a,unser93a,unser93b}*/ is slower because it uses more neighbours.'};
interp.labels  = {'Nearest neighbour'
                  'Trilinear'
                  '2nd Degree B-Spline'
                  '3rd Degree B-Spline'
                  '4th Degree B-Spline'
                  '5th Degree B-Spline'
                  '6th Degree B-Spline'
                  '7th Degree B-Spline'}';
interp.values  = {0 1 2 3 4 5 6 7};
interp.def     = @(val)geg_petpve12_get_defaults('coreg.write.interp', val{:});
%--------------------------------------------------------------------------
% wrap Wrapping
wrap         = cfg_menu;
wrap.tag     = 'wrap';
wrap.name    = 'Wrapping';
wrap.help    = {'These are typically:'
                '    No wrapping - for PET or images that have already been spatially transformed.'
                '    Wrap in Y   - for (un-resliced) MRI where phase encoding is in the Y direction (voxel space).'}';
wrap.labels  = {'No wrap'
                'Wrap X'
                'Wrap Y'
                'Wrap X & Y'
                'Wrap Z'
                'Wrap X & Z'
                'Wrap Y & Z'
                'Wrap X, Y & Z'}';
wrap.values  = {[0 0 0] [1 0 0] [0 1 0] [1 1 0] [0 0 1] [1 0 1] [0 1 1]...
               [1 1 1]};
wrap.def     = @(val)geg_petpve12_get_defaults('coreg.write.wrap', val{:});
%--------------------------------------------------------------------------
% mask Masking
mask         = cfg_menu;
mask.tag     = 'mask';
mask.name    = 'Masking';
mask.help    = {'Because of subject motion, different images are likely to have different patterns of zeros from where it was not possible to sample data. With masking enabled, the program searches through the whole time series looking for voxels which need to be sampled from outside the original images. Where this occurs, that voxel is set to zero for the whole set of images (unless the image format can represent NaN, in which case NaNs are used where possible).'};
mask.labels  = {'Mask images'
                'Dont mask images'}';
mask.values  = {1 0};
mask.def     = @(val)geg_petpve12_get_defaults('coreg.write.mask', val{:});
%--------------------------------------------------------------------------
% prefix Filename Prefix
prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Filename Prefix';
prefix.help    = {'Specify the string to be prepended to the filenames of the resliced image file(s). Default prefix is ''r''.'};
prefix.strtype = 's';
prefix.num     = [1 Inf];
prefix.def     = @(val)geg_petpve12_get_defaults('coreg.write.prefix', val{:});
%--------------------------------------------------------------------------
% roptions Reslice Options
roptions         = cfg_branch;
roptions.tag     = 'roptions';
roptions.name    = 'Reslice Options';
roptions.val     = {interp wrap mask prefix};
roptions.help    = {'Various reslicing options.'};
%--------------------------------------------------------------------------
% estwrite Coreg: Estimate & Reslice
coreg_estwrite      = cfg_exbranch;
coreg_estwrite.tag  = 'coreg_estwrite';
coreg_estwrite.name = 'Coregister: Estimate & Reslice';
coreg_estwrite.val  = {refs sources others eoptions roptions };
coreg_estwrite.help = {'The registration method used here is based on work by Collignon et al/* \cite{collignon95}*/. The original interpolation method described in this paper has been changed in order to give a smoother cost function.  The images are also smoothed slightly, as is the histogram.  This is all in order to make the cost function as smooth as possible, to give faster convergence and less chance of local minima.'
                 ''
                 'At the end of coregistration, the voxel-to-voxel affine transformation matrix is displayed, along with the histograms for the images in the original orientations, and the final orientations.  The registered images are displayed at the bottom.'
                 ''
                 'Registration parameters are stored in the headers of the "source" and the "other" images. These images are also resliced to match the source image voxel-for-voxel. The resliced images are named the same as the originals except that they are prefixed by ''r''.'}';
coreg_estwrite.prog = @geg_run_coreg;
coreg_estwrite.vout = @vout_estwrite;

%========================================================================
%========================================================================
%------------------------------------------------------------------------
% Check that the headers of an image are the same (centered at same 0,0,0)
%------------------------------------------------------------------------
%========================================================================
%========================================================================
refs         = cfg_files;
refs.tag     = 'refs';
refs.name    = 'Reference images';
refs.help    = {'Select Reference images. '};
refs.filter  = {'image'};
refs.ufilter = '.*';
refs.num     = [1 Inf];

srcs         = cfg_files;
srcs.tag     = 'srcs';
srcs.name    = 'Source images';
srcs.help    = {'Select images to check match with reference. '};
srcs.filter  = {'image'};
srcs.ufilter = '.*';
srcs.num     = [1 Inf];

ChkHdrs      = cfg_exbranch;
ChkHdrs.tag = 'ChkHdrs';
ChkHdrs.name = 'PETPVE12: check PET-MRI corregister';
ChkHdrs.val = {refs,srcs};
ChkHdrs.prog   = @geg_macthNIfTI_headers;
ChkHdrs.help   = {[...
    ' writes all the rc*.nii files necessary to perform DARTEL alignment', ...
    '. ']};

%========================================================================
%========================================================================
%------------------------------------------------------------------------
% PET image intensity normalization
%------------------------------------------------------------------------
%========================================================================
%========================================================================
data.help = {[...
'Select PET data (e.g. PVEc-PET images) for processing. ',...
'This assumes that there is one scan for each subject. ']};

RefMask         = cfg_files;
RefMask.tag     = 'RefMask';
RefMask.name    = 'Reference region masks';
RefMask.help    = {'Select a binary mask of the reference region for each subject. They must all have the same image dimensions, orientation, voxel size etc, as their corresponding subject.'};
RefMask.filter  = {'image'};
RefMask.ufilter = '.*';
RefMask.num     = [1 Inf];

type1      = cfg_branch;
type1.tag  = 'type1';
type1.name = 'Individual masks of reference region';
type1.val  = {RefMask};
type1.help = {['It uses a mask in subject space (i.e. hand drawn/native space segmented) to automatically compute the constant "activity" in wihtin this area, ',...
    'and use this value for intensity normalization. ',...
    'Assumes the CSF concentration is zero. ']};
% ---------------------------------------------------------------------
RefMask.help    = {'Select one binary reference mask in standard space. It will be deformed using the input flow_fields/deformation_matrices to match each subject.'};
RefMask.num     = [1 1];
RefMask.def     = @(val)geg_petpve12_get_defaults('PETnormOpts.MNIroi', val{:});

defs        = cfg_files;
defs.tag    = 'defs';
defs.name   = 'Inverse Deformation Field';
defs.filter = '.*iy_.*\.nii$';
defs.num    = [1 inf];
defs.help   = {[...
'Deformations can be thought of as vector fields, and represented ',...
'by three-volume images.  In SPM, deformation fields are saved in ',...
'NIfTI format, with dimensions xdim x ydim x zdim x 1 x 3. ',...
'Each voxel contains the x, y and z mm coordinates of where the deformation points. ']};

ffield         = files('Flow field','flowfield','nifti',[1 inf]);
ffield.ufilter = '^u_.*';
ffield.help    = {...
    ['The flow field stores the deformation information. '...
     'The same field can be used for both forward or backward deformations '...
     '(or even, in principle, half-way or exaggerated deformations). '...
     'One image per subject is expected. ']};

K      = mnu('Time Steps','K',...
        {'1','2','4','8','16','32','64','128','256','512'},...
        {0,1,2,3,4,5,6,7,8,9});
K.val  = {6};
K.help = {...
    ['The number of time points used for solving the '...
     'partial differential equations.  A single time point would be '...
     'equivalent to a small deformation model. '...
     'Smaller values allow faster computations, '...
     'but are less accurate in terms '...
     'of inverse consistency and may result in the one-to-one mapping '...
     'breaking down.']};
 
template        = cfg_files;
template.tag    = 'template';
template.name   = 'Dartel Template';
template.filter = 'nifti';
template.num    = [0 1];
template.val    = {{''}};
template.help   = {...
['Select the final Template file generated by Dartel. This will be affine '...
 'registered with a TPM file, such that the resulting spatially normalised '...
 'images are closely aligned to MNI space. Leave empty if you do not wish to '...
 'incorporate a transform to MNI space '...
 '(ie just click ''done'' on the file selector, without selecting any images).']};

drtl      = branch('Dartel flow','dartel',{ffield,K,template});
drtl.help = {'Imported Dartel flow field.'};

InDefs        = cfg_choice;
InDefs.tag    = 'InDefs';
InDefs.name   = 'Deformation input';
InDefs.values = {drtl,defs};
InDefs.help   = {'This option allows for the specification of the field images to deform from MNI to subject''s space. '};
InDefs.val    = {defs};

type2      = cfg_branch;
type2.tag  = 'type2';
type2.name = 'Standard space mask of reference region';
type2.val  = {RefMask, InDefs};
type2.help = {['It uses a user specified iy_* matrices to transform the specified mask (i.e. cerebellum cortex) in standard space to subject''s space. ',...
    'And then compute the constant "activity" in that mask, that will be used to perform the intensity normalization of corresponding PET image. ',...
    'By default it assumes the CSF concentration is zero. ']};
% ---------------------------------------------------------------------
Refvals         = cfg_entry;
Refvals.tag     = 'RefVals';
Refvals.name    = 'Reference values vector';
Refvals.help    = {'Vector of "activity" values in reference region.'};
Refvals.strtype = 'r';
Refvals.num     = [Inf 1];

type3      = cfg_branch;
type3.tag  = 'type3';
type3.name = 'User estimated values';
type3.val  = {Refvals};
type3.help = {'It uses pre-computed constant "activity" values from a reference region, entered by the user as single vector. '};
% ---------------------------------------------------------------------

PETnorm_opts        = cfg_choice;
PETnorm_opts.tag    = 'PETnorm_opts';
PETnorm_opts.name   = 'Reference activity';
PETnorm_opts.values = {type1, type2, type3};
PETnorm_opts.help   = {'This option allows for the specification of the method user for the WM/CSF "activity", assumed to be a constant value. Pre-estimated values can be entered.'};
PETnorm_opts.val    = {type2};

PETnorm      = cfg_exbranch;
PETnorm.tag  = 'PETnorm';
PETnorm.name = 'PETPVE12: PET intensity normalization';
PETnorm.val  = {data,PETnorm_opts};
PETnorm.prog = @geg_PET_intNorm;
PETnorm.help = {[...
    'Computes the mean value of PET activit in the selected mask, ', ...
    'and then divides the activity at each voxel by this value. ']};

%========================================================================
%========================================================================
%------------------------------------------------------------------------
% The toolbox tools definitions (petpve.tools.)
%------------------------------------------------------------------------
%========================================================================
%========================================================================
tools        = cfg_choice;
tools.name   = 'Tools';
tools.tag    = 'tools';
tools.values = {check_cov,ChkHdrs,coreg_estimate,coreg_estwrite,PETnorm};

return

%========================================================================
%========================================================================
% Sub-Functions
%________________________________________________________________________
%________________________________________________________________________
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%________________________________________________________________________
%________________________________________________________________________

%==========================================================================
function dep = vout_estimate(~)
dep(1)            = cfg_dep;
dep(1).sname      = 'Coregistered Images';
dep(1).src_output = substruct('.','cfiles');
dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
dep(2)            = cfg_dep;
dep(2).sname      = 'Coregistration Matrix';
dep(2).src_output = substruct('.','M');
dep(2).tgt_spec   = cfg_findspec({{'strtype','r'}});
%________________________________________________________________________
%------------------------------------------------------------------------
function dep = vout_estwrite(job)
depe = vout_estimate(job);
depc = vout_reslice(job);
dep = [depe depc];
%________________________________________________________________________
%------------------------------------------------------------------------
function dep = vout_reslice(~)
dep(1)            = cfg_dep;
dep(1).sname      = 'Resliced Images';
dep(1).src_output = substruct('.','rfiles');
dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
%________________________________________________________________________
%========================================================================
function files_item = files(name, tag, fltr, num)
files_item        = cfg_files;
files_item.name   = name;
files_item.tag    = tag;
files_item.filter = fltr;
files_item.num    = num;
function menu_item = mnu(name, tag, labels, values)
menu_item        = cfg_menu;
menu_item.name   = name;
menu_item.tag    = tag;
menu_item.labels = labels;
menu_item.values = values;
function branch_item = branch(name, tag, val)
branch_item      = cfg_branch;
branch_item.name = name;
branch_item.tag  = tag;
branch_item.val  = val;