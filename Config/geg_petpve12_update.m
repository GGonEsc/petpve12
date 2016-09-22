function geg_petpve12_update(update)
% check for new updates
%
% FORMAT geg_petpve12_update(update)
% update - allow installation of update
% 
% This function will connect itself to the server, compare the
% version number of the updates with the one of the VBM8 installation 
% currently in the MATLAB path and will display the outcome.
%_______________________________________________________________________
% Gabriel González-Escamilla
% $Id: geg_petpve12_update.m 001 2015-04-13 10:54:23Z $

% rev = '$Rev: 001 $';

if nargin == 0
  update = 0;
end

% get current release number
r = 0;
A = ver;
for i=1:length(A)
  if strcmp(A(i).Name,'PET Partial Volume Effects correction Toolbox')
    r = str2double(A(i).Version);
  end
end

if update
    error('Not implemented yet.');
end

% get new release number (must be included in the zip file)
url = 'http://---/petpve12/';
if usejava('jvm')
  [s,sts] = urlread(url);
  if ~sts
    fprintf('Cannot access %s. Please check your proxy and/or firewall to allow access.\n. You can download your update at %s\n',url,url); 
    return
  end
else
  fprintf('Please enable Java (JVM) to use update function.\n. You can download your update at %s\n',url); 
  return
end

n = regexp(s,'petpve12_r(\d.*?)\.zip','tokens'); % this should be the structure of the file
if isempty(n)
  fprintf('There are no new releases available yet.\n');
  return;
else
  % get largest release number
  rnew = [];
  for i=1:length(n)
    rnew = [rnew str2double(n{i})];
  end
  rnew = max(rnew);
end

if rnew > r
  fprintf('A new version of PETPVE12 is available on: %s\n',url);
  fprintf('Your version: %d - New version: %d\n',r,rnew);
  if ~update
    fprintf('In order to update use Toolbox|PETPVE12|Check for updates\n');%,r,rnew
  end

  if update
    d = fullfile(spm('Dir'),'toolbox'); 
    overwrite = spm_input('Update',1,'m','Do not update|Download zip-file only|Overwrite old PETPVE12 installation',[-1 0 1],3);
    switch overwrite
    case 1
      try
        % list mex-files and delete these files to prevent that old
        % compiled files are used
        mexfiles = dir(fullfile(d,'petpve12','*.mex*'));% but now they have been moved into directories!!!
        for i=1:length(mexfiles)
          name = fullfile(d,'vbm8',mexfiles(i).name);
          spm_unlink(name);
        end
        fprintf('Download PETPVE12\n');
        s = unzip([url sprintf('petpve12_r%d.zip',rnew)], d);% should match with the name of the zip file
        fprintf('%d files have been updated.\nSPM should be restarted.\n',numel(s));
        restart = spm_input('Restart SPM',1,'m','no|yes',[0 1],2);
        if restart
          rehash
          toolbox_path_cache
          eval('spm pet');
        end
      catch
        fprintf('Update failed: check file permissions. Download zip-file only.\n');
        web([url sprintf('petpve12_r%d.zip',rnew)],'-browser');
        fprintf('Unzip file to %s\n',d);
      end
    case 0
      web([url sprintf('petpve12_r%d.zip',rnew)],'-browser');
      fprintf('Unzip file to %s\n',d);
    end
  end
elseif update
  fprintf('You already have the newest version %d.\n',r);
end