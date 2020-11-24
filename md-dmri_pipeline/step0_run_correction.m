function step0_run_correction(input_parameters)

%% Load information from the input_parameters structure
data_path = input_parameters.data_directory;
data_file = input_parameters.data_file;
xps_file = input_parameters.xps_file;

%% Prepare options
opt = mdm_opt();

if input_parameters.do_rotate_bvec
    opt.mdm.mec.do_rotate_bvec = 1;
else
    opt.mdm.mec.do_rotate_bvec = 0;
end
 
s.nii_fn = fullfile(data_path, data_file);
s.xps = mdm_xps_load(fullfile(data_path, xps_file));
output_path = data_path;
 
% Write the elastix parameter file
p = elastix_p_affine(200);
p_fn = elastix_p_write(p, fullfile(data_path, 'p.txt'));
  
% Run an extrapolation-based registration
mdm_s_mec(s, p_fn, output_path, opt);
% s_registered.nii_fn;
% s_registered.xps;