function petpve = tbx_cfg_petpve12
% Configuration file for jobs
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Department of Imaging Neuroscience

% based on John Ashburners version of
% tbx_cfg_preproc8.m
%
% Gabriel Gonzalez-Escamilla
% $Id: tbx_cfg_petpve12.m 429 2015-03-08 20:47:21Z $

% rev = '$Rev: 003 $';

% Check if is the complete toolbox on the path
mpath=path;
posit=strfind(mpath,';');
if isempty(posit), posit=strfind(mpath,':'); end
inPath=cell(1);
for i=1:length(posit)
    if i==1
        inPath{i}=mpath(1:posit(i)-1);
    elseif i==length(posit)
        inPath{i}=mpath(posit(i)+1:end);
    else
        inPath{i}=mpath(posit(i)+1:posit(i+1)-1);
    end    
end
tbxDir=fullfile(spm('Dir'),'toolbox','petpve12');
tbxscnt=dir(tbxDir);
for i=3:length(tbxscnt)
    tbxsfldr=fullfile(tbxDir,tbxscnt(i).name);
    if isdir(tbxsfldr) && ~isempty(find(strcmp(inPath,tbxsfldr), 1))
    
    elseif isdir(tbxsfldr) && isempty(find(strcmp(inPath,tbxsfldr), 1))
        addpath(tbxsfldr)
    end
end

addpath(fileparts(which(mfilename)));

%_______________________________________________________________________
%------------------------------------------------------------------------
tools = geg_petpve12_tools;
spatial = geg_petpve12_spatial;
%------------------------------------------------------------------------
petpve  = cfg_choice;
petpve.name = 'petpve';
petpve.tag  = 'petpve';
petpve.values = {tools,spatial};
%------------------------------------------------------------------------
%------------------------------------------------------------------------
