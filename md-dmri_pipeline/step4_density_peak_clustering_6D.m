function step4_density_peak_clustering_6D(input_parameters)
% Clustering of local fiber orientations

%% Load information from the input_parameters structure
% GENERAL INFORMATION
NBS = input_parameters.nb_MC_inversions;
method = input_parameters.inversion_method;
data_directory = input_parameters.data_directory;
bootstrap_directory = input_parameters.bootstrap_directory;
odf_directory = input_parameters.odf_directory;
clustering_directory = input_parameters.clustering_directory;
clustering_figure_path = fullfile(clustering_directory, 'fiber_clustering_figures');
mask_file = input_parameters.mask_file;
dir_flips = input_parameters.dir_flips;

nb_parameters = 7; % Should be equal to the inversion dimension (physical parameters) + 1 (configuration weight)

% CLUSTERING PARAMETERS
threshold_d_delta = input_parameters.threshold_d_delta_thin;
struct_clustering = input_parameters.structure_clustering;
do_plot = struct_clustering.do_plot;
max_nb_clusters = struct_clustering.max_nb_clusters;

% NUMBER OF ODF PEAKS
[nb_peaks_map, nifti_header] = mdm_nii_read(fullfile(odf_directory,'nb_peaks_maximal.nii.gz'));

%% STORE MONTE CARLO REALIZATIONS
data_file = load(fullfile(bootstrap_directory,'1','mfs.mat'));
[Nx, Ny, Nz, dimension] = size(data_file.mfs.m);
ind_dpar = 2:nb_parameters:dimension;
ind_dperp = 3:nb_parameters:dimension;
ind_theta = 4:nb_parameters:dimension;
ind_phi = 5:nb_parameters:dimension;
ind_r2 = 6:nb_parameters:dimension;
ind_r1 = 7:nb_parameters:dimension;
ind_w = (nb_parameters+1):nb_parameters:dimension;
% nifti_header = data_file.mfs.nii_h;

nb_nonzero_nodes_map = zeros([Nx,Ny,Nz,NBS]);
S0_list = zeros(Nx, Ny, Nz, NBS);

% Merge all bootstrap solutions
all_bs_indices_thin = [];
all_dpar_thin = [];
all_dperp_thin = [];
all_theta_thin = [];
all_phi_thin = [];
all_w_thin = [];
all_r1_thin = [];
all_r2_thin = [];

for n = 1:NBS
    data_file = load(fullfile(bootstrap_directory,num2str(n),'mfs.mat'));
    nb_nonzero_nodes_map(:,:,:,n) = data_file.mfs.m(:,:,:,1); % Voxel-wise map of the number of non-zero nodes
    
    dpar = data_file.mfs.m(:,:,:,ind_dpar);
    dperp = data_file.mfs.m(:,:,:,ind_dperp);
    theta = data_file.mfs.m(:,:,:,ind_theta);
    phi = data_file.mfs.m(:,:,:,ind_phi);
    w = data_file.mfs.m(:,:,:,ind_w);
    S0_list(:,:,:,n) = sum(w,4);
    w = w./sum(w,4); % All solutions are normalized within the voxel
    r1 = data_file.mfs.m(:,:,:,ind_r1);
    r2 = data_file.mfs.m(:,:,:,ind_r2);
    
    % Filter out not thin solutions
    d_delta = (dpar - dperp)./(dpar + 2.*dperp);
    d_delta(~isfinite(d_delta)) = 0;
    
    ind_relevant = logical(squeeze((w > 0).*(d_delta > threshold_d_delta).*(r1 > 0).*(r2 > 0)));
    dpar(~ind_relevant) = NaN; % NaNs are necessary to keep the number of recorded elements constant for all voxels, for the later concatenations (see below)
    dperp(~ind_relevant) = NaN;
    theta(~ind_relevant) = NaN;
    phi(~ind_relevant) = NaN;
    w(~ind_relevant) = NaN; 
    r1(~ind_relevant) = NaN;
    r2(~ind_relevant) = NaN;
%     w = w/sum(w(ind_relevant)); % Solutions are normalized within the thin bin for each bootstrap
    
    all_bs_indices_thin = cat(4, all_bs_indices_thin, n*ones(size(dpar))); % 'These' concatenations (see above)
    all_dpar_thin = cat(4, all_dpar_thin, dpar);
    all_dperp_thin = cat(4, all_dperp_thin, dperp);
    all_theta_thin = cat(4, all_theta_thin, theta);
    all_phi_thin = cat(4, all_phi_thin, phi);
    all_w_thin = cat(4, all_w_thin, w);   
    all_r1_thin = cat(4, all_r1_thin, r1);
    all_r2_thin = cat(4, all_r2_thin, r2);
end

%% CLUSTERING INPUTS AND S0 MAP
msf_mkdir(clustering_directory)
if do_plot
    msf_mkdir(clustering_figure_path)
end
save(fullfile(clustering_directory, 'input_clustering.mat'), '-struct', 'struct_clustering');
mdm_nii_write(median(S0_list,4), fullfile(clustering_directory, 'S0_median.nii.gz'), nifti_header);
mdm_nii_write(iqr(S0_list,4), fullfile(clustering_directory, 'S0_iqr.nii.gz'), nifti_header);

%% FIBER MASK
if strcmp(mask_file, '')
    mask = data_file.mfs.mask;
else
    mask = mdm_nii_read(fullfile(data_directory, mask_file));
end
peak_mask = double(nb_peaks_map > 0);
node_mask = double(sum(nb_nonzero_nodes_map,4) > 0);
combined_mask = (mask > 0).*(peak_mask > 0).*(node_mask > 0);

mdm_nii_write(mask, [clustering_directory '/mask_data.nii.gz'], nifti_header);
mdm_nii_write(node_mask, [clustering_directory '/mask_node.nii.gz'], nifti_header);
mdm_nii_write(peak_mask, [clustering_directory '/mask_peaks.nii.gz'], nifti_header);
mdm_nii_write(combined_mask, [clustering_directory '/mask_combined.nii.gz'], nifti_header);

%% LOOP OVER VOXELS

% Decrease dimensionality to lower memory cost
[~, ~, ~, size_all] = size(all_bs_indices_thin);
all_bs_indices_thin = reshape(all_bs_indices_thin, Nx*Ny*Nz, size_all);
all_dpar_thin = reshape(all_dpar_thin, Nx*Ny*Nz, size_all);
all_dperp_thin = reshape(all_dperp_thin, Nx*Ny*Nz, size_all);
all_theta_thin = reshape(all_theta_thin, Nx*Ny*Nz, size_all);
all_phi_thin = reshape(all_phi_thin, Nx*Ny*Nz, size_all);
all_r1_thin = reshape(all_r1_thin, Nx*Ny*Nz, size_all);
all_r2_thin = reshape(all_r2_thin, Nx*Ny*Nz, size_all);
nb_peaks_map = reshape(nb_peaks_map, Nx*Ny*Nz, 1);
combined_mask = reshape(combined_mask, Nx*Ny*Nz, 1)';

% Retain only non-masked-out voxels.
all_bs_indices_thin = all_bs_indices_thin(combined_mask>0, :);
all_dpar_thin = all_dpar_thin(combined_mask>0, :);
all_dperp_thin = all_dperp_thin(combined_mask>0, :);
all_theta_thin = all_theta_thin(combined_mask>0, :);
all_phi_thin = all_phi_thin(combined_mask>0, :);
all_r1_thin = all_r1_thin(combined_mask>0, :);
all_r2_thin = all_r2_thin(combined_mask>0, :);
nb_peaks_map = nb_peaks_map(combined_mask>0);
[N, ~] = size(all_bs_indices_thin);
         
map_cluster_metrics = cell(N, 1);
cluster_median_orientations = single(zeros(N, 3*max_nb_clusters));
cluster_median_orientations_weighted_by_w = single(zeros(N, 3*max_nb_clusters));
mask_repeated_clustering = zeros(N, 1);
mask_odd_clustering = zeros(N, 1);

parfor v = 1:N
    %% FIBER CLUSTERING
    bs_indices_fibers = squeeze(all_bs_indices_thin(v,:))';
    dpar_fibers = squeeze(all_dpar_thin(v,:))';
    dperp_fibers = squeeze(all_dperp_thin(v,:))';
    theta_fibers = squeeze(all_theta_thin(v,:))';
    phi_fibers = squeeze(all_phi_thin(v,:))';
    w_fibers = squeeze(all_w_thin(v,:))';
    
    ind_relevant = ~isnan(theta_fibers);
    bs_indices_fibers = bs_indices_fibers(ind_relevant);
    dpar_fibers = dpar_fibers(ind_relevant);
    dperp_fibers = dperp_fibers(ind_relevant);
    theta_fibers = theta_fibers(ind_relevant);
    phi_fibers = phi_fibers(ind_relevant);
    w_fibers = w_fibers(ind_relevant); % Retains the normalization over all solutions
    
    x_fibers = (-1)^dir_flips(1).*sin(theta_fibers).*cos(phi_fibers);
    y_fibers = (-1)^dir_flips(2).*sin(theta_fibers).*sin(phi_fibers);
    z_fibers = (-1)^dir_flips(3).*cos(theta_fibers);
    
    r1_fibers = squeeze(all_r1_thin(v,:))';
    r1_fibers = r1_fibers(ind_relevant);
    t1_fibers = msf_notfinite2zero(1./r1_fibers);
    
    r2_fibers = squeeze(all_r2_thin(v,:))';
    r2_fibers = r2_fibers(ind_relevant);
    t2_fibers = msf_notfinite2zero(1./r2_fibers);
    
    if length(w_fibers) < NBS % Not even one thin-bin solution per bootstrap realization in average
        mask_odd_clustering(v) = 1;
        continue
    end
    
    w_fibers_clustering = w_fibers/sum(w_fibers); % Normalization over the ensemble of bootstrap thin-bin solutions, for clustering only
    optimal_K_initial = nb_peaks_map(v); % Use the number of ODF peaks to set the optimal number of clusters
    [optimal_K, ind_fibers_all_clusters, ind_fibers_out_of_clustering, rho] = density_peaks_clustering_improved(optimal_K_initial, theta_fibers, phi_fibers, w_fibers_clustering, struct_clustering);
    
    % Book-keeping of special cases
    if optimal_K ~= optimal_K_initial
        mask_repeated_clustering(v) = 1;
    end
    
    if optimal_K == -1
        mask_odd_clustering(v) = 1;
        continue
    end
    
    % Flip the points so that they end up in geometrical clusters (applies to the rest of the clustering procedure)
    for m = 1:optimal_K
        ind_in_cluster = ind_fibers_all_clusters{m};
        
        % Use the point of maximal density in the cluster as reference
        [~, i_max] = max(rho(ind_in_cluster));
        ind = ind_in_cluster(i_max);
        if z_fibers(ind) < 0
            x_fibers(ind) = -x_fibers(ind);
            y_fibers(ind) = -y_fibers(ind);
            z_fibers(ind) = -z_fibers(ind);
        end
        x_ref = x_fibers(ind);
        y_ref = y_fibers(ind);
        z_ref = z_fibers(ind);
        
        for i = 1:length(ind_in_cluster)
            ind = ind_in_cluster(i);
            if dot([x_fibers(ind), y_fibers(ind), z_fibers(ind)], [x_ref, y_ref, z_ref]) < 0
                x_fibers(ind) = -x_fibers(ind);
                y_fibers(ind) = -y_fibers(ind);
                z_fibers(ind) = -z_fibers(ind);
            end
        end
    end
    
    % Cluster visualization
    if do_plot
        average_x_clusters = zeros(1,optimal_K);
        average_y_clusters = zeros(1,optimal_K);
        average_z_clusters = zeros(1,optimal_K);
        for m = 1:optimal_K
            ind_in_cluster = ind_fibers_all_clusters{m};
            w_fibers_cluster = w_fibers(ind_in_cluster);
            w_fibers_cluster = w_fibers_cluster/sum(w_fibers_cluster); % Normalization, just to compute the cluster medians
            
            x_fibers_cluster = x_fibers(ind_in_cluster);
            y_fibers_cluster = y_fibers(ind_in_cluster);
            z_fibers_cluster = z_fibers(ind_in_cluster);
            
            average_x_clusters(m) = sum(w_fibers_cluster.*x_fibers_cluster)/sum(w_fibers_cluster); % No need to normalize these
            average_y_clusters(m) = sum(w_fibers_cluster.*y_fibers_cluster)/sum(w_fibers_cluster);
            average_z_clusters(m) = sum(w_fibers_cluster.*z_fibers_cluster)/sum(w_fibers_cluster);
        end
        colors_cluster = cluster_orientational_colors(average_x_clusters, average_y_clusters, average_z_clusters);
        
        plot_fiber_sphere(x_fibers, y_fibers, z_fibers, w_fibers_clustering, ind_fibers_all_clusters, ind_fibers_out_of_clustering, optimal_K, colors_cluster, clustering_figure_path, vx, vy, vz)
    end
    
    %% CLUSTER-WISE CLASSIFICATION OF THE BOOTSTRAP SOLUTIONS
    x_clusters = cell([1, optimal_K]);
    y_clusters = cell([1, optimal_K]);
    z_clusters = cell([1, optimal_K]);
    dpar_clusters = cell([1, optimal_K]);
    dperp_clusters = cell([1, optimal_K]);
    theta_clusters = cell([1, optimal_K]);
    phi_clusters = cell([1, optimal_K]);
    w_clusters = cell([1, optimal_K]);
    r1_clusters = cell([1, optimal_K]);
    t1_clusters = cell([1, optimal_K]);
    r2_clusters = cell([1, optimal_K]);
    t2_clusters = cell([1, optimal_K]);
    
    geomedian_x_clusters = zeros([1, optimal_K]);
    geomedian_y_clusters = zeros([1, optimal_K]);
    geomedian_z_clusters = zeros([1, optimal_K]);
    
    median_orientations = single(zeros([1, 3*max_nb_clusters]));
    median_orientations_weighted_by_w = single(zeros([1, 3*max_nb_clusters]));
    
    for m = 1:optimal_K
        ind_in_cluster = ind_fibers_all_clusters{m};
        bs_ind_in_cluster = bs_indices_fibers(ind_in_cluster);
        for n = 1:NBS
            ind_in_cluster_and_bootstrap = ind_in_cluster(bs_ind_in_cluster == n);
            if ~isempty(ind_in_cluster_and_bootstrap)
                w_relevant = w_fibers(ind_in_cluster_and_bootstrap);
                x_relevant = x_fibers(ind_in_cluster_and_bootstrap);
                y_relevant = y_fibers(ind_in_cluster_and_bootstrap);
                z_relevant = z_fibers(ind_in_cluster_and_bootstrap);
                dpar_relevant = dpar_fibers(ind_in_cluster_and_bootstrap);
                dperp_relevant = dperp_fibers(ind_in_cluster_and_bootstrap);
                
                % To ensure that the averages of x, y and z be on the unit sphere
                average_x = sum(w_relevant.*x_relevant)/sum(w_relevant);
                average_y = sum(w_relevant.*y_relevant)/sum(w_relevant);
                average_z = sum(w_relevant.*z_relevant)/sum(w_relevant);
                list_average = [average_x, average_y, average_z]/norm([average_x, average_y, average_z]);
                [average_theta, average_phi] = cartesian2spherical_unit_sphere(list_average(1), list_average(2), list_average(3));
                
                w_clusters{m} = [w_clusters{m} sum(w_relevant)]; % Still weighted by the overall thin bin fraction
                x_clusters{m} = [x_clusters{m} list_average(1)];
                y_clusters{m} = [y_clusters{m} list_average(2)];
                z_clusters{m} = [z_clusters{m} list_average(3)];
                theta_clusters{m} = [theta_clusters{m} average_theta];
                phi_clusters{m} = [phi_clusters{m} average_phi];
                dpar_clusters{m} = [dpar_clusters{m} sum(w_relevant.*dpar_relevant)/sum(w_relevant)];
                dperp_clusters{m} = [dperp_clusters{m} sum(w_relevant.*dperp_relevant)/sum(w_relevant)];
                
                r1_relevant = r1_fibers(ind_in_cluster_and_bootstrap);
                t1_relevant = t1_fibers(ind_in_cluster_and_bootstrap);
                r1_clusters{m} = [r1_clusters{m} sum(w_relevant.*r1_relevant)/sum(w_relevant)];
                t1_clusters{m} = [t1_clusters{m} sum(w_relevant.*t1_relevant)/sum(w_relevant)];
                
                r2_relevant = r2_fibers(ind_in_cluster_and_bootstrap);
                t2_relevant = t2_fibers(ind_in_cluster_and_bootstrap);
                r2_clusters{m} = [r2_clusters{m} sum(w_relevant.*r2_relevant)/sum(w_relevant)];
                t2_clusters{m} = [t2_clusters{m} sum(w_relevant.*t2_relevant)/sum(w_relevant)];
            end
        end
        
        %% Find the weighted geometric median of the cluster points to extract a median orientation for each cluster
        [~, average_phi] = cartesian2spherical_unit_sphere(mean(x_clusters{m}), mean(y_clusters{m}), mean(z_clusters{m}));
        
        bounds_theta = [min(theta_clusters{m}) max(theta_clusters{m})];
        bounds_phi = [min(phi_clusters{m}) max(phi_clusters{m})];
        if ~(average_phi > bounds_phi(1)) && ~(average_phi < bounds_phi(2)) % Flip if the cluster is around the (x > 0) xz half-plane
            bounds_phi = [bounds_phi(2) bounds_phi(1)+2*pi];
        end
        
        delta_angle = 1*pi/180;
        possible_geomedian_theta = linspace(bounds_theta(1), bounds_theta(2), ceil(diff(bounds_theta)/delta_angle));
        possible_geomedian_phi = linspace(bounds_phi(1), bounds_phi(2), ceil(diff(bounds_phi)/delta_angle));
        
        if isempty(possible_geomedian_theta) || isempty(possible_geomedian_phi)
            mask_odd_clustering(v) = 1;
            break
        end
        
        total_distance = zeros([length(possible_geomedian_theta) length(possible_geomedian_phi)]);
        for i = 1:length(possible_geomedian_theta)
            for j = 1:length(possible_geomedian_phi)
                total_distance(i,j) = sum(w_clusters{m}.*acos(abs(cosine_angular_difference(possible_geomedian_theta(i),possible_geomedian_phi(j), theta_clusters{m}, phi_clusters{m}))));
            end
        end
        
        [row, col] = find(ismember(total_distance, min(total_distance(:))));
        geomedian_theta = possible_geomedian_theta(row(1));
        geomedian_phi = possible_geomedian_phi(col(1));
        
        geomedian_x_clusters(m) = sin(geomedian_theta)*cos(geomedian_phi);
        geomedian_y_clusters(m) = sin(geomedian_theta)*sin(geomedian_phi);
        geomedian_z_clusters(m) = cos(geomedian_theta);
        
        median_orientations(3*(m-1)+1:3*(m-1)+3) = single([geomedian_x_clusters(m), geomedian_y_clusters(m), geomedian_z_clusters(m)]);
        median_orientations_weighted_by_w(3*(m-1)+1:3*(m-1)+3) = single([geomedian_x_clusters(m), geomedian_y_clusters(m), geomedian_z_clusters(m)]*median(w_clusters{m}));
    end
    cluster_median_orientations(v,:) = median_orientations;
    cluster_median_orientations_weighted_by_w(v,:) = median_orientations_weighted_by_w;
    
    % Necessary because of the isempty(possible_geomedian_theta) || isempty(possible_geomedian_phi) condition
    if mask_odd_clustering(v) == 1
        continue
    end
    
    %% Store cluster metrics
    cluster_metrics = struct;
    cluster_metrics.n = optimal_K;
    cluster_metrics.x = x_clusters;
    cluster_metrics.y = y_clusters;
    cluster_metrics.z = z_clusters;
    cluster_metrics.dpar = dpar_clusters;
    cluster_metrics.dperp = dperp_clusters;
    cluster_metrics.theta = theta_clusters;
    cluster_metrics.phi = phi_clusters;
    cluster_metrics.w = w_clusters;
    cluster_metrics.median_orientations = median_orientations;
    cluster_metrics.median_orientations_weighted_by_w = median_orientations_weighted_by_w;
    cluster_metrics.r1 = r1_clusters;
    cluster_metrics.t1 = t1_clusters;
    cluster_metrics.r2 = r2_clusters;
    cluster_metrics.t2 = t2_clusters;
    map_cluster_metrics{v} = cluster_metrics;
    
    % Cluster visualization
    if do_plot
        colors_cluster = cluster_orientational_colors(geomedian_x_clusters, geomedian_y_clusters, geomedian_z_clusters);
        plot_fiber_sphere_no_outlier(x_clusters, y_clusters, z_clusters, w_clusters, optimal_K, colors_cluster, clustering_figure_path, vx, vy, vz)
    end
end

dummy_map = cell([Nx*Ny*Nz, 1]);
dummy_ind = 1:Nx*Ny*Nz;
dummy_ind = dummy_ind(combined_mask > 0);
for v = 1:N
    dummy_map{dummy_ind(v)} = map_cluster_metrics{v};
end
map_cluster_metrics = reshape(dummy_map, Nx, Ny, Nz);

dummy = zeros(Nx*Ny*Nz, 1);
dummy(:, combined_mask>0) = mask_repeated_clustering;
mask_repeated_clustering = reshape(dummy, Nx, Ny, Nz);

dummy = zeros(Nx*Ny*Nz, 1);
dummy(:, combined_mask>0) = mask_odd_clustering;
mask_odd_clustering = reshape(dummy, Nx, Ny, Nz);

dummy = zeros(Nx*Ny*Nz, 3*max_nb_clusters);
dummy(:, combined_mask>0) = cluster_median_orientations;
cluster_median_orientations = single(reshape(dummy, Nx, Ny, Nz, 3*max_nb_clusters));

dummy = zeros(Nx*Ny*Nz, 3*max_nb_clusters);
dummy(:, combined_mask>0) = cluster_median_orientations_weighted_by_w;
cluster_median_orientations_weighted_by_w = single(reshape(dummy, Nx, Ny, Nz, 3*max_nb_clusters));

%% CREATING MASKS and MAPS
mdm_nii_write(mask_repeated_clustering, fullfile(clustering_directory, 'mask_repeated_clustering.nii.gz'), nifti_header);
mdm_nii_write(mask_odd_clustering, fullfile(clustering_directory, 'mask_odd_clustering.nii.gz'), nifti_header);
mdm_nii_write(cluster_median_orientations, fullfile(clustering_directory, 'cluster_median_orientations.nii.gz'), nifti_header);
mdm_nii_write(cluster_median_orientations_weighted_by_w, fullfile(clustering_directory, 'cluster_median_orientations_weighted_by_w.nii.gz'), nifti_header);
save(fullfile(clustering_directory, 'cluster_metrics.mat'), 'map_cluster_metrics');

end