function global_inversion_pipeline(input_parameters)

%% PREPARATION

% Set path to MRtrix if needed
PATH = getenv('PATH');
if ~contains(PATH, 'mrtrix')
    setenv('PATH', [PATH ':' input_parameters.mrtrix_directory])
end

% Transfer data
initial_directory = input_parameters.initial_directory;
data_directory = input_parameters.data_directory;
if ~exist(data_directory, 'dir')
    msf_mkdir(data_directory)
    copyfile(fullfile(initial_directory, input_parameters.data_file), data_directory)
    copyfile(fullfile(initial_directory, input_parameters.xps_file), data_directory)
    if ~strcmp(input_parameters.mask_file,'')
        copyfile(fullfile(initial_directory, input_parameters.mask_file), data_directory)
    end
end

% Start parallel pool if not already open
if isempty(gcp('nocreate'))
    if isfinite(input_parameters.nb_parallel_workers)
        parpool('local', input_parameters.nb_parallel_workers);
    else
        parpool('local');
    end
end

%% FILTER OUT XPS FIELDS
xps_fields_to_filter_out = {'protocol_dimensions'};
if logical(length(xps_fields_to_filter_out))
   xps = mdm_xps_load(fullfile(data_directory, input_parameters.xps_file));
   for m = 1:length(xps_fields_to_filter_out)
       field = xps_fields_to_filter_out{m};
       if isfield(xps,field)
           xps = rmfield(xps, field);
       end
   end
   mdm_xps_save(xps, fullfile(data_directory, input_parameters.xps_file))
end

%% RENAME XPS FIELDS (if saturation recovery experiment)
xps_fields_initial = {'ts'};
xps_fields_final = {'tr'};
if logical(length(xps_fields_initial))
   xps = mdm_xps_load(fullfile(data_directory, input_parameters.xps_file));
   for m = 1:length(xps_fields_initial)
       old_field = xps_fields_initial{m};
       new_field = xps_fields_final{m};
       if isfield(xps, old_field)
           [xps.(new_field)] = xps.(old_field);
           xps = rmfield(xps, old_field);
       end
   end
   mdm_xps_save(xps, fullfile(data_directory, input_parameters.xps_file))
end

%% DATA CORRECTION
if input_parameters.do_data_correction 
    if ~isfile(fullfile(data_directory,'data_mc.nii.gz'))
        fprintf('\nSTEP 0: Running correction via extrapolation-based references...\n')
        t_start = tic;
        step0_run_correction(input_parameters)
        toc(t_start)
    else
        fprintf('\nSTEP 0 ALREADY PERFORMED: data corrected using extrapolation-based references...\n')
    end
    input_parameters.data_file = 'data_mc.nii.gz'; % Adapt for name change upon correction ("mc" = "motion-corrected" + eddy-corrected)
    input_parameters.xps_file = 'data_mc_xps.mat';
end

%% ANALYSIS (DTD_COVARIANCE)
clearvars -except input_parameters parallel_pool
if input_parameters.do_covariance && strcmp(input_parameters.inversion_method, 'dtd') && ~exist(input_parameters.covariance_directory, 'dir')
    fprintf('\nSTEP 0: Running covariance tensor approximation...\n')
    t_start = tic;
    step0_run_covariance(input_parameters)
    toc(t_start)
elseif input_parameters.do_covariance && strcmp(input_parameters.inversion_method, 'dtd')
    fprintf('\nSTEP 0 ALREADY PERFORMED: data processed using the covariance tensor approximation...\n')
end
    
%% ANALYSIS (MONTE CARLO REALIZATIONS)
clearvars -except input_parameters parallel_pool
bootstrap_path = input_parameters.bootstrap_directory;
method = input_parameters.inversion_method;
if ~exist(fullfile(bootstrap_path, num2str(input_parameters.nb_MC_inversions)), 'dir')
    if input_parameters.nb_MC_inversions == 1
        fprintf(['\nSTEP 1: Running ' method ' processing with ' num2str(input_parameters.nb_MC_inversions) ' Monte Carlo realization...\n'])
    else
        fprintf(['\nSTEP 1: Running ' method ' processing with ' num2str(input_parameters.nb_MC_inversions) ' Monte Carlo realizations...\n'])
    end
    t_start = tic;
    step1_run_bootstrap(input_parameters)
    toc(t_start)
else
    if input_parameters.nb_MC_inversions == 1
        fprintf(['\nSTEP 1 ALREADY PERFORMED: data processed using ' method ' with ' num2str(input_parameters.nb_MC_inversions) ' Monte Carlo realization\n'])
    else
        fprintf(['\nSTEP 1 ALREADY PERFORMED: data processed using ' method ' with ' num2str(input_parameters.nb_MC_inversions) ' Monte Carlo realizations\n'])
    end
end

%% FROM MONTE CARLO REALIZATIONS TO PARAMETER MAPS AND ODFs
clearvars -except input_parameters parallel_pool
pdf_path = input_parameters.parameter_maps_directory;
if ~isfile(fullfile(pdf_path, 'technicolor_maps.pdf')) || input_parameters.repeat_colormaps
    fprintf('\nSTEP 2: Computing parameter maps and ODFs...\n')
    t_start = tic;
    if strcmp(input_parameters.inversion_method, 'dtr2r1d')
        step2_bootstrap_to_maps_ODFs_6D(input_parameters)
    else
        step2_bootstrap_to_maps_ODFs(input_parameters)
    end
    toc(t_start)
else
    fprintf('\nSTEP 2 ALREADY PERFORMED: parameter maps and ODFs computed\n')
end

%% SH-ODF PEAKS
clearvars -except input_parameters parallel_pool
odf_path = input_parameters.odf_directory;
if ~exist(fullfile(odf_path,'nb_ODF_peaks.nii.gz'), 'file')
    fprintf('\nSTEP 3: Computing ODF peaks...\n')
    t_start = tic;
    if strcmp(input_parameters.inversion_method, 'dtr2r1d')
        step3_get_odf_peaks_and_metrics_6D(input_parameters)
    else
        step3_get_odf_peaks_and_metrics(input_parameters)
    end
    toc(t_start)
else
    fprintf('\nSTEP 3 ALREADY PERFORMED: ODF peaks computed\n')
end

%% CLUSTERING
if input_parameters.do_clustering
    clearvars -except input_parameters parallel_pool
    clustering_path = input_parameters.clustering_directory;
    if ~exist(clustering_path, 'dir')
        fprintf('\nSTEP 4: Performing Monte-Carlo density-peak clustering...\n')
        t_start = tic;
        if strcmp(input_parameters.inversion_method, 'dtr2r1d')
            step4_density_peak_clustering_6D(input_parameters)
        else
            step4_density_peak_clustering(input_parameters)
        end
        toc(t_start)
    else
        fprintf('\nSTEP 4 ALREADY PERFORMED: Monte-Carlo density-peak clustering performed\n')
    end
end


%% THE END
clearvars
fprintf('\n')
delete(gcp('nocreate'));
fprintf('ANALYSIS FINISHED\n')

end