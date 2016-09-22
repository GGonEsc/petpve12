function geg_calc_atlasROImeans_bothNativ()
% This computes everything on a slice-basis intead of image-basis, which is
% supposed to be more correct. It assumes that the parcelled-atlas and the
% functional volumes are in the same space or at least coregistered.
% It outputs an xls file with the values of each ROI (columns) and each
% subject (rows). 
% 
% Saves computed values to a xls file.
% 
% Based on Ashburners script to measure a ROI volume
% Modified by M. Grothe to include only positive values
% Adapted & Integrated by G. Gonzalez-Escamilla
%_______________________________________________________________________
% $Id: geg_calc_atlasROImeans.m 001 2015-10-05 22:44:53Z $
% 
% 
% rev = '$Rev: 002 $'; % 07-October-2015

% data
S  = spm_select(inf,'image','Select PET images');
atlas  = spm_select(size(S,1),'image','Select template atlas for every subject');
AtlasDesc = spm_select([1,1],'^*\.txt','Select Atlas descriptor file');
opth = spm_select([1,1],'Dir','Select output directory');

AtlDesc=importdata(char(AtlasDesc));
ParcAtlasNames=AtlDesc.textdata;
ParcAtlasIds=AtlDesc.data;
signal = zeros(size(S,1),numel(ParcAtlasIds));
average = zeros(size(S,1),numel(ParcAtlasIds));
snames = cell(size(S,1),1);

nam = 'template_avgVals';
tots_xls = fullfile(opth,nam);
xlswrite(tots_xls,['SubjID', ParcAtlasNames'],'ROIavg','A1')

% Running for every subject
for subj=1:size(S,1)
    % Loading data
    actPT = char(S(subj,:));% PET name
    [~,snam,~] = spm_fileparts(actPT);
    PTm = spm_vol(actPT);% PET header    
    snames{subj,1} = snam;
    
    % 'Template atlas in subject space'
    actAtl = char(atlas(subj,:));% PET name
    atlasm = spm_vol(actAtl);
    
    % Dimensions should be equal, then I'll not check it
    datlasi = spm_read_vols(atlasm);
    PTi = spm_read_vols(PTm);
    
    fprintf('Estimating Regions: ');
    for k=1:numel(ParcAtlasIds)
        actp = ParcAtlasIds(k);
        fprintf('%d,',actp) % progression indicator
        img = PTi;
        
        % Identify actual parcel
        actparc = round(datlasi)==(actp);% logical
        actparc = double(actparc);% numerical
        
        % Computing mean values
        img  = img.*actparc;
        img(~isfinite(img)) = 0;
        tot = sum(img(img>0)); %sums only positive voxels
        img(img>0)=1; %sets positive voxels to 1
        vox = sum(img(img>0)); %counts only positive voxels by summing 1s
        signal(subj,k) = tot;
        % divide sum of voxel-values by total number of contributing voxels
        % (i.e. voxels with positive value in (PET) image)
        average(subj,k) = signal(subj,k)/vox;
    end
    fprintf('\n')
end

fprintf('\n')
% print to a xls file
xlswrite(tots_xls,snames,'ROIavg','A2')
xlswrite(tots_xls,average,'ROIavg','B2')
disp('Done')