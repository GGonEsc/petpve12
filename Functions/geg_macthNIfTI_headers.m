function geg_macthNIfTI_headers(job)
% After any corregistration process it makes sure that two images are
% centered at the same 0,0,0
% 
%---------------------------------------------------------------------
% Author: Gabriel González-Escamilla 25-March-2015 17:38:03Z $
%

%rev = '$Rev: 003 $'; 08-April-2015

if nargin == 1   
    S1 = job.refs;
    S2 = job.srcs;
else    
    % Get files
    S1 = spm_select(Inf,'image','Select Reference image');    
    S2 = spm_select(Inf,'image','Select Source image');    
end

h = warndlg('This function does not perform coregistration or reslicing, only verifies that the two headers are coregistered, and equals image origin','!! Warning !!');
uiwait(h)

if size(S1,1)~=size(S2,1)
    fprintf('Not an equal number of files: %s references & %s images to check match \n',num2str(length(S1)),num2str(length(S2)))
   return 
end

% Running for every subject
for subj=1:size(S1,1)
    % % /0/ load original images
    V1=spm_vol(char(S1(subj,:)));
    V2=spm_vol(char(S2(subj,:)));
    disp(['Reference: ' V1.fname]);
    disp(['Target: ' V2.fname]);
    
    % Check that the PETi matches with the MRI (at least on the headers)
    if any(V2.dim~=V1.dim) && sum(any(V2.mat~=V1.mat))~=0
        fprintf(2,['The input source %s appears to be not in corregitration with reference data, ',...
            'please use the coregister function and then come back again.\n'],V2.fname);
        return;
    elseif all(V2.dim==V1.dim) && sum(any(V2.mat~=V1.mat))~=0
        fprintf(['The input pet %s appears to be not completely matched with reference data, ',...
            'I''ll fix that!!!.\n'],V2.fname)
        
        % Load files
        %[Yref,~] = spm_read_vols(V1);%XYZref
        [YSourc,~] = spm_read_vols(V2);%XYZtarg
        
        % writing data
        filename=V2.fname;
        [pathstr,name,ext,versn] = spm_fileparts(filename);
        disp(['saving: ' 'm' name ext versn])
        ext='.nii';
        nfilename=[pathstr '\m' name ext];
        V3=V2;
        V3.fname=nfilename;
        V3.mat=V1.mat;% sustituir la matriz del header
        spm_write_vol(V3,YSourc);
        fprintf('\n')
    else
        fprintf('Everything is in order. Not doing anything\n')
    end
end