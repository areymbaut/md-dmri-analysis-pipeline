function step0_run_covariance(input_parameters)

%% Load information from the input_parameters structure
method = 'dtd_covariance';
data_path = input_parameters.data_directory;
covariance_path = input_parameters.covariance_directory;
data_file = input_parameters.data_file;
mask_file = input_parameters.mask_file;
xps_file = input_parameters.xps_file;

%% Prepare signal
s.nii_fn = fullfile(data_path, data_file);
if strcmp(mask_file,'')
    s.mask_fn = fullfile(data_path, 'data_mask.nii.gz');
    input_parameters.mask_file = 'data_mask.nii.gz';
else
    s.mask_fn = fullfile(data_path, mask_file);
end
s.xps = mdm_xps_load(fullfile(data_path, xps_file));

%% Prepare options
opt = mdm_opt();
eval(['opt = ' method '_opt(opt);'])

if strcmp(mask_file,'')
    opt.do_mask = 1;
else
    opt.do_mask = 0;
end
opt.do_data2fit = 1;
opt.do_datafit2chisq = 0;
opt.do_fit2param = 1;
opt.do_param2nii = 1;
opt.do_m2pdf = 0;   

opt.verbose = 1;
opt.do_overwrite = 1;

%% Bootstrap
msf_mkdir(covariance_path);
paths.nii_path = fullfile(covariance_path, 'nii_maps');
paths.mfs_fn = fullfile(covariance_path, 'mfs.mat');
paths.chisq_fn = fullfile(covariance_path, 'chisq.mat');
paths.ind_fn = fullfile(covariance_path, 'ind.mat');
paths.dps_fn = fullfile(covariance_path, 'dps.mat');
paths = mdm_paths(paths);
dtd_covariance_pipe(s, paths, opt); % parfor enabled in mio_opt.m
