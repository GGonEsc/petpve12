function geg_centerOfMass
% use center-of-mass (COM) to roughly correct for differences in the
% position between image and template
% The center of mass is always close to the AC
%
% This script comes from VBM8 toolbox

S = spm_select(Inf,'image','select images to set the new origin');
VF = spm_vol(S);
n = size(S,1);

% pre-estimated COM of MNI template
com_reference = [0 -20 -30];
% com_reference = [0 -20 -15];% closer to the AC

for i=1:n
  fprintf('Setting orig to center-of-mass for %s\n',VF(i).fname);
  Affine = eye(4);
  vol = spm_read_vols(VF(i));
  avg = mean(vol(:));
  avg = mean(vol(find(vol>avg)));
  
	% don't use background values in COM calculation
	[x,y,z] = ind2sub(size(vol),find(vol>avg));
	com = VF(i).mat(1:3,:)*[mean(x) mean(y) mean(z) 1]';
	com = com';

	M = spm_get_space(VF(i).fname);
	Affine(1:3,4) = (com - com_reference)';
  spm_get_space(VF(i).fname,Affine\M);
end
fprintf('Finished\n')