function geg_PETframes_avg(job)
% This function performs the average of PET images frames.
% Dynamic PET images can be processed and will result in a corrected
% dynamic series. However, for the matching (soregister) and PVEc steps, a
% static PET image showing anatomical information is required. 
% In the case of a static scan, as in the example above, this step is not
% required.
%_______________________________________________________________________
% Copyright (C) 2015
%
% Gabriel Gonzalez-Escamilla
% $Id: geg_PET_mean.m 001 2015-04-21 14:12:36Z $
% 

% rev = '$Rev: 001 $'; 

if nargin == 1
    S = job.data;    
else
    S = spm_select(Inf,'image','Select dynamic PET images');    
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
    
    % Loading 4D data
    actPT = char(S(subj,:));% corregistered PET name
    PTm = spm_vol(actPT);% PET id matrix
    PTi = spm_read_vols(PTm);% PET image
    [pth,nam,ext] = spm_fileparts(actPT);
    if ndims(PTi)==4
        odir = pth;
        if strncmpi(ext,'.nii',4)
            cvt = ['a' nam ext];
        else
            cvt = ['a' nam '.nii'];
        end
        
        % average frames
        % I might allow the user to specify frames to average
        avgfrm = mean(PTi,4);
        
        % writing data
        % Save averaged PET image
        fprintf('saving: %s \n', cvt)
        oPT = PTm;
        oPT.fname = fullfile(odir, cvt);
        spm_write_vol(oPT,avgfrm);
    else
        fprintf('Image %s is already static \n',actPT)
    end
    
    fprintf('\n')
end