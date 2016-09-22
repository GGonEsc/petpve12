function out = geg_run_coreg(job)
% geg_run_coreg.m
%
% Helps with coregister process.
% SPM only takes one image per subject (one reference, 1 source).
% This wraper allws the user to make all computations over multiple
% subjects at once.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%_______________________________________________________________________
% Copyright (C) 2015
% 
% Gabriel Gonzalez-Escamilla
% based on spm_run_coreg.m
% $Id: geg_PVEcorrection.m 001 2015-09-30 12:58:30Z $
% 
% 
% rev = '$Rev: 002 $'; % 01-October-2015


if nargin == 1
    T = job.refs;
    S = job.sources;
    O = job.others;    
end

for subj=1:size(T,1)
    actref = T(subj);
    actsource = S(subj);
    if isempty(O{:}),  actother = {}; else  actother = O(subj); end
    if ~exist('actother','var') || isempty(actother), actother = {}; end
    PO = [actsource(:); actother(:)];
    PO = spm_select('expand',PO);
    
    %-Coregister
    %--------------------------------------------------------------------------
    if isfield(job,'eoptions')
        x  = spm_coreg(char(actref), char(actsource), job.eoptions);
        
        M  = spm_matrix(x);
        MM = zeros(4,4,numel(PO));
        for j=1:numel(PO)
            MM(:,:,j) = spm_get_space(PO{j});
        end
        for j=1:numel(PO)
            spm_get_space(PO{j}, M\MM(:,:,j));
        end
    end
    
    %-Reslice
    %--------------------------------------------------------------------------
    if isfield(job,'roptions')
        P            = char(actref{:},actsource{:},actother{:});
        flags.mask   = job.roptions.mask;
        flags.mean   = 0;
        flags.interp = job.roptions.interp;
        flags.which  = 1;
        flags.wrap   = job.roptions.wrap;
        flags.prefix = job.roptions.prefix;
        
        spm_reslice(P, flags);
    end
    
    %-Dependencies
    %--------------------------------------------------------------------------
    if isfield(job,'eoptions')
        out.cfiles   = PO;
        out.M        = M;
    end
    if isfield(job,'roptions')
        out.rfiles   = spm_file(PO, 'prefix',job.roptions.prefix);
    end
    
end