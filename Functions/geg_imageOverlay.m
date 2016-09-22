function geg_imageOverlay

% background image
imgs = spm_select(1, 'image', 'Select Background Image');

geg_image('Init',imgs)
geg_image('AddColorOverlay')

fprintf('\n')
