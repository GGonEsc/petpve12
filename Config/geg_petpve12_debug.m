function geg_petpve12_debug
%geg_pvc12_debug	print debug information for SPM12 and PVC12
%
% FORMAT geg_petpve12_debug
%
%__________________________________________________________________________
% Gabriel Gonzalez-Escamilla
% $Id: geg_petpve12_debug.m 001 2015-03-08 15:43:40Z $

rev = '$Rev: 001 $';

% print last error
fprintf('\nLast error message:\n');
fprintf('-------------------------------------------------------------------------------------\n');
fprintf('-------------------------------------------------------------------------------------\n');
try
	er = lasterror;
	fprintf('%s\n',er.message);
	if isfield(er,'stack')
		for i=1:length(er.stack)
			fprintf('%s at line %g\n',char(er.stack(i).file),er.stack(i).line);
		end
	end
catch
	fprintf('%s\n',lasterr);
end

fprintf('-------------------------------------------------------------------------------------\n');
fprintf('-------------------------------------------------------------------------------------\n');

fprintf('\nVersion information:\n');
fprintf('-------------------------------------------------------------------------------------\n');

ver

