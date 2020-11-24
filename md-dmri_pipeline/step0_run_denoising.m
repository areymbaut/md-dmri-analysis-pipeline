function step0_run_denoising(input_parameters)

default_window_size = 5;

%% Load information from the input_parameters structure
data_path = input_parameters.data_directory;
data_file = input_parameters.data_file;
mask_file = input_parameters.mask_file;

[data, h] = mdm_nii_read(fullfile(data_path, data_file));
data = double(data);

%% Set window (triplet of odd integers, works for isotropic-resolution data)
sz = size(data);
if length(sz) == 4 % 4D diffusion volume
    if all(sz(1:3) > default_window_size)
        window = default_window_size*[1 1 1];
    else
        a = min([default_window_size, sz(1)]);
        b = min([default_window_size, sz(2)]);
        c = min([default_window_size, sz(3)]);
        if mod(a,2) == 0, a = a-1; end
        if mod(b,2) == 0, b = b-1; end
        if mod(c,2) == 0, c = c-1; end
        window = [a b c];
    end
elseif length(sz) == 3 % 3D diffusion volume
    if all(sz(1:2) > default_window_size)
        window = default_window_size*[1 1];
    else
        a = min([default_window_size, sz(1)]);
        b = min([default_window_size, sz(2)]);
        if mod(a,2) == 0, a = a-1; end
        if mod(b,2) == 0, b = b-1; end
        window = [a b];
    end
end

%% Run denoising
if strcmp(mask_file,'')
    [data_denoised, noise_std_map_before, noise_std_map_after, nb_principal_components_map] = MP_PCA_denoise(data, window);
else
    mask = mdm_nii_read(fullfile(data_path, mask_file));
    mask = logical(mask);
    [data_denoised, noise_std_map_before, noise_std_map_after, nb_principal_components_map] = MP_PCA_denoise(data, window, mask);
end

%% Save images
% Denoised images
mdm_nii_write(data_denoised, fullfile(data_path, 'data_dn.nii.gz'), h);

% Map of estimated noise standard deviation variance before denoising
mdm_nii_write(noise_std_map_before, fullfile(data_path, 'noise_std_before.nii.gz'), h);

% Map of estimated noise standard deviation variance after denoising
mdm_nii_write(noise_std_map_after, fullfile(data_path, 'noise_std_after.nii.gz'), h);

% Number of detected signal principal components
mdm_nii_write(nb_principal_components_map, fullfile(data_path, 'nb_principal_components.nii.gz'), h);
