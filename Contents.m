% PET Partial Volume Effects correction Toolbox
% Version  010  (PETPVC12)  23-March-2015
% _________________________________________________________________________
%  http://-...-/
% _________________________________________________________________________
% Partial Volume Effects correction toolbox for PET images (SPM12)
% PET - Partial Volume Effects correction Toolbox (PETPVE12)
% Copyright (C) 2015 Gabriel González-Escamilla gabriellbk@gmail.com
% 
% =========================================================================
% Description
% =========================================================================
% This toolbox is a collection of extensions to the segmentation algorithm 
% of SPM12 (Wellcome Department of Cognitive Neurology) to provide Partial
% volume Effects correction algorithms. It was developed by Gabriel
% Gonzalez-Escamilla, with the colaboration of Michel Grothe, Catharina 
% Lange and Ralph Buchert and is available to the scientific community
% under the terms of the GNU General Public License.
%
% General files
%   PET-PVEc_12.man              - notes on petpve12 toolbox
%   Contents.m                   - this file
%   spm_petpve12.m               - toolbox wrapper to call functions
%   tbx_cfg_petpve12.m           - configuration file for jobs
%
% petpve12 configurations
% 
%   cg_petpve12_debug.m              - print debug information for SPM8 and petpve12
%   cg_petpve12_defaults.m           - sets the defaults for petpve12
%   cg_petpve12_get_defaults.m       - defaults for petpve12
%   cg_petpve12_tools.m              - wrapper for calling toolbox utilities
%   cg_petpve12_update.m             - check for new updates
%   cg_petpve12_spatial.m            - wrapper to perform segmentation and/or PVEc mehods
%
% Utility and functions
% 
% The subfolder 'Functions' contains all the function that make work this
% toolbox.
% 
% Mex- and c-functions
% All come from the VBM8 toolbox
%
% Templates/Images
% VBM8 IXI550 images + AAL atlas
