function step1_run_bootstrap(input_parameters)

%% Load information from the input_parameters structure
method = input_parameters.inversion_method;
data_path = input_parameters.data_directory;
bootstrap_path = input_parameters.bootstrap_directory;
nb_inversions = input_parameters.nb_MC_inversions;
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
s.xps = mdm_xps_load(fullfile(data_path, xps_file)); %#ok<STRNU>

%% Prepare options
opt = mdm_opt();
eval(['opt = ' method '_opt(opt);'])

if strcmp(mask_file,'')
    opt.do_mask = 1;
else
    opt.do_mask = 0;
end
opt.do_data2fit = 1;
opt.do_bootstrap = 1;
opt.do_shuffle = 0;
opt.do_datafit2chisq = 0;
opt.do_fit2param = 0;
opt.do_param2nii = 0;
opt.do_m2pdf = 0;

opt.verbose = 1;
opt.do_overwrite = 1;

% opt.dtd.dmin = 1/max(s.xps.b);
% opt.dtd.dmax = 1/min(s.xps.b);

%% Bootstrap
msf_mkdir(bootstrap_path);

for nbs = 1:nb_inversions
    bs_o = fullfile(bootstrap_path,num2str(nbs));
    if isfolder(bs_o)
        fprintf(['Bootstrap ' num2str(nbs) ' already run\n'])
    else
        fprintf([num2str(nbs) '\n'])
        msf_mkdir(bs_o);
        
        paths.nii_path = fullfile(bs_o, 'nii_maps');
        paths.mfs_fn = fullfile(bs_o, 'mfs.mat');
        paths.chisq_fn = fullfile(bs_o, 'chisq.mat');
        paths.ind_fn = fullfile(bs_o, 'ind.mat');
        paths.dps_fn = fullfile(bs_o, 'dps.mat');
        paths = mdm_paths(paths);
        
        eval([method '_pipe(s, paths, opt);']) % parfor enabled in mio_opt.m
    end
end