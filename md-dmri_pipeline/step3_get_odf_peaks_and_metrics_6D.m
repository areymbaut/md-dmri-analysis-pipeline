function step3_get_odf_peaks_and_metrics_6D(input_parameters)

%% Load information from the input_parameters structure
method = input_parameters.inversion_method;
max_n_peaks = input_parameters.max_nb_odf_peaks;
data_directory = input_parameters.data_directory;
bootstrap_directory = input_parameters.bootstrap_directory;
odf_directory = input_parameters.odf_directory;
mask_file = input_parameters.mask_file;
nb_MeshNodes = input_parameters.nb_mesh_nodes;
threshold_w = input_parameters.threshold_ODF_w;

mfs = mdm_mfs_load(fullfile(bootstrap_directory, '1', 'mfs.mat'));
h = mfs.nii_h;
if strcmp(mask_file, '')
    mask = mfs.mask;
    mask_fn = fullfile(data_directory, 'data_mask.nii.gz');
    mdm_nii_write(double(mask), mask_fn, h);
    input_parameters.mask_file = 'data_mask.nii.gz';
else
    mask_fn = fullfile(data_directory, input_parameters.mask_file);
    mask = mdm_nii_read(mask_fn);
end
 
opt = mdm_opt(); %#ok<NASGU>
eval(['opt = ' method '_opt(opt);'])

load(fullfile(odf_directory, ['odf_bsmedian_' num2str(nb_MeshNodes)]), 'odf_bsmedian'); 

%% Use the discrete ODFs as spherical functions (SF) to fit spherical harmonics (SH) coefficients
if ~strcmp(input_parameters.mrtrix_directory, '')
    % SH fitting
    struct_sf_to_sh.x = odf_bsmedian.x;
    struct_sf_to_sh.y = odf_bsmedian.y;
    struct_sf_to_sh.z = odf_bsmedian.z;
    struct_sf_to_sh.w = odf_bsmedian.w_bin{1};
    struct_sf_to_sh.diso = odf_bsmedian.diso_bin{1};
    struct_sf_to_sh.sqddelta = odf_bsmedian.sqddelta_bin{1};
    struct_sf_to_sh.r1 = odf_bsmedian.r1_bin{1};
    struct_sf_to_sh.t1 = odf_bsmedian.t1_bin{1};
    struct_sf_to_sh.r2 = odf_bsmedian.r2_bin{1};
    struct_sf_to_sh.t2 = odf_bsmedian.t2_bin{1};
    sf_to_sh(struct_sf_to_sh, 8, [0 1 0], mask_fn, odf_directory, method)
    
    % Get SH peaks (requires MRtrix)
    SH_coeffs_fn = fullfile(odf_directory, 'SH_coeffs_ODF_normalized.nii.gz');
    SH_peaks_fn = fullfile(odf_directory, 'SH_peaks.nii.gz');
    command = ['sh2peaks -force -num ' num2str(max_n_peaks) ' -threshold ' num2str(threshold_w) ' -mask ' mask_fn ' ' SH_coeffs_fn ' ' SH_peaks_fn];
    [stat, res] = system(command);
    
    % Compute the map of the number of SH peaks
    abs_sh_peaks = abs(msf_notfinite2zero(mdm_nii_read(SH_peaks_fn)));
    [Nx, Ny, Nz, ~] = size(odf_bsmedian.w_bin{1});
    nb_sh_peaks_map = zeros([Nx,Ny,Nz]);
    for i = 1:max_n_peaks
        ind = (3*(i-1)+1):(3*(i-1)+3);
        is_peak_map = double(logical(sum(abs_sh_peaks(:,:,:,ind),4)));
        nb_sh_peaks_map = nb_sh_peaks_map + is_peak_map;
    end
    mdm_nii_write(floor(nb_sh_peaks_map), fullfile(odf_directory, 'nb_SH_peaks.nii.gz'), h);
else
    [Nx, Ny, Nz, ~] = size(odf_bsmedian.w_bin{1});
    nb_sh_peaks_map = zeros([Nx,Ny,Nz]);
end

%% Work on the discrete ODFs
odf_norm = 0.3*max(odf_bsmedian.w_bin{1}(:));

% Decrease dimensionality to lower memory cost
nb_sh_peaks_map = reshape(nb_sh_peaks_map, Nx*Ny*Nz, 1);
odf_bsmedian.w_bin = reshape(odf_bsmedian.w_bin{1}, Nx*Ny*Nz, nb_MeshNodes);
odf_bsmedian.diso_bin = reshape(odf_bsmedian.diso_bin{1}, Nx*Ny*Nz, nb_MeshNodes);
odf_bsmedian.sqddelta_bin = reshape(odf_bsmedian.sqddelta_bin{1}, Nx*Ny*Nz, nb_MeshNodes);
odf_bsmedian.r1_bin = reshape(odf_bsmedian.r1_bin{1}, Nx*Ny*Nz, nb_MeshNodes);
odf_bsmedian.t1_bin = reshape(odf_bsmedian.t1_bin{1}, Nx*Ny*Nz, nb_MeshNodes);
odf_bsmedian.r2_bin = reshape(odf_bsmedian.r2_bin{1}, Nx*Ny*Nz, nb_MeshNodes);
odf_bsmedian.t2_bin = reshape(odf_bsmedian.t2_bin{1}, Nx*Ny*Nz, nb_MeshNodes);
mask = reshape(mask, Nx*Ny*Nz, 1)';

% Retain only non-masked-out voxels
nb_sh_peaks_map = nb_sh_peaks_map(mask>0);
odf_bsmedian.w_bin = odf_bsmedian.w_bin(mask>0, :);
odf_bsmedian.diso_bin = odf_bsmedian.diso_bin(mask>0, :);
odf_bsmedian.sqddelta_bin = odf_bsmedian.sqddelta_bin(mask>0, :);
odf_bsmedian.r1_bin = odf_bsmedian.r1_bin(mask>0, :);
odf_bsmedian.t1_bin = odf_bsmedian.t1_bin(mask>0, :);
odf_bsmedian.r2_bin = odf_bsmedian.r2_bin(mask>0, :);
odf_bsmedian.t2_bin = odf_bsmedian.r2_bin(mask>0, :);
N = length(nb_sh_peaks_map);

all_peaks = single(zeros(N, 3*max_n_peaks));
all_sddelta_weighted_peaks = single(zeros(N, 3*max_n_peaks));
all_w_weighted_peaks = single(zeros(N, 3*max_n_peaks));

cell_odf_peaks = cell(N, 1);
nb_odf_peaks_map = zeros([N, 1]);
max_nb_peaks_map = zeros([N, 1]);

parfor v = 1:N
    n_sh_peaks = nb_sh_peaks_map(v);
    
    odf_s = struct;
    conn = struct;
    structure = struct;
    
    odf_s.n = odf_bsmedian.n;
    odf_s.x = odf_bsmedian.x;
    odf_s.y = odf_bsmedian.y;
    odf_s.z = odf_bsmedian.z;
    odf_s.c = abs([odf_s.x odf_s.y odf_s.z]);
    odf_s.tri = odf_bsmedian.tri;
    odf_s.w = squeeze(odf_bsmedian.w_bin(v,:))'/odf_norm;
    odf_s.w_normalized = odf_s.w/sum(odf_s.w);
    odf_s.diso = squeeze(odf_bsmedian.diso_bin(v,:))';
    odf_s.sddelta = squeeze(odf_bsmedian.sqddelta_bin(v,:))';
    odf_s.r1 = squeeze(odf_bsmedian.r1_bin(v,:))';
    odf_s.t1 = squeeze(odf_bsmedian.t1_bin(v,:))';
    odf_s.r2 = squeeze(odf_bsmedian.r2_bin(v,:))';
    odf_s.t2 = squeeze(odf_bsmedian.t2_bin(v,:))';
    odf_s.verts = repmat(odf_s.w,[1 3]).*[odf_s.x odf_s.y odf_s.z];
    odf_s.norm_verts = vecnorm(odf_s.verts')';
    
    TR = triangulation(odf_s.tri, odf_s.verts);
    conn.basis = vertexAttachments(TR);
    indx = false(size(conn.basis,1),1);
    for i = 1:size(conn.basis,1)
        conn.tri = conn.basis{i};
        conn.verts = unique([odf_s.tri(conn.tri,1); odf_s.tri(conn.tri,2); odf_s.tri(conn.tri,3)]);
        if all(odf_s.norm_verts(i) >= odf_s.norm_verts(conn.verts))
            indx(i) = 1;
        end
    end
    
    % Filtering out low-probability peaks
    indw = odf_s.w_normalized/max(odf_s.w_normalized) >= threshold_w;
    odf_verts_peaks = odf_s.verts(indx & indw, :);
    diso_peaks = odf_s.diso(indx & indw, :);
    sddelta_peaks = odf_s.sddelta(indx & indw, :);
    r1_peaks = odf_s.r1(indx & indw, :);
    t1_peaks = odf_s.t1(indx & indw, :);
    r2_peaks = odf_s.r2(indx & indw, :);
    t2_peaks = odf_s.t2(indx & indw, :);
    w_peaks = odf_s.w_normalized(indx & indw, :);
    
    odf_s = [];
    
    for i = 1:nnz(indx & indw)
        odf_verts_peaks(i,:) = odf_verts_peaks(i,:)/norm(odf_verts_peaks(i,:));
    end
    
    % Filtering out redundant antipodal points and taking up to n_peaks peaks with highest weights
    minimal_angular_separation = odf_bsmedian.minimal_angular_separation;
    relevant_ind = [];
    checked_ind = [];
    for i = 1:nnz(indx & indw)-1
        if ~ismember(i, checked_ind)
            for j = i+1:nnz(indx & indw)
                if abs(dot(odf_verts_peaks(i,:),odf_verts_peaks(j,:))) > cos(2*minimal_angular_separation)
                    if w_peaks(i) >= w_peaks(j)
                        relevant_ind = [relevant_ind i];
                    else
                        relevant_ind = [relevant_ind j];
                    end
                    checked_ind = [checked_ind j];
                end
            end
        end
    end
    
    odf_verts_peaks = odf_verts_peaks(relevant_ind, :);
    diso_peaks = diso_peaks(relevant_ind);
    sddelta_peaks = sddelta_peaks(relevant_ind);
    r1_peaks = r1_peaks(relevant_ind);
    t1_peaks = t1_peaks(relevant_ind);
    r2_peaks = r2_peaks(relevant_ind);
    t2_peaks = t2_peaks(relevant_ind);
    w_peaks = w_peaks(relevant_ind);
    
    norm_odf_verts_peaks = sqrt(odf_verts_peaks(:,1).^2+odf_verts_peaks(:,2).^2 + odf_verts_peaks(:,3).^2);
    [~, indx] = sort(norm_odf_verts_peaks,'descend');
    
    norm_peaks = double.empty;
    n_peaks = max_n_peaks;
    if numel(indx) < max_n_peaks
        odf_verts_peaks(1:numel(indx),:) = odf_verts_peaks(indx,:);
        norm_peaks(1:numel(indx),:) = norm_odf_verts_peaks(indx);
        diso_peaks = diso_peaks(indx);
        sddelta_peaks = sddelta_peaks(indx);
        r1_peaks = r1_peaks(indx);
        t1_peaks = t1_peaks(indx);
        r2_peaks = r2_peaks(indx);
        t2_peaks = t2_peaks(indx);
        w_peaks = w_peaks(indx);
        n_peaks = numel(indx);
    else
        odf_verts_peaks = odf_verts_peaks(indx(1:n_peaks),:);
        norm_peaks = norm_odf_verts_peaks(indx(1:n_peaks));
        diso_peaks = diso_peaks(indx(1:n_peaks));
        sddelta_peaks = sddelta_peaks(indx(1:n_peaks));
        r1_peaks = r1_peaks(indx(1:n_peaks));
        t1_peaks = t1_peaks(indx(1:n_peaks));
        r2_peaks = r2_peaks(indx(1:n_peaks));
        t2_peaks = t2_peaks(indx(1:n_peaks));
        w_peaks = w_peaks(indx(1:n_peaks));
    end
    
    nb_odf_peaks_map(v) = n_peaks;
    max_nb_peaks_map(v) = max(n_sh_peaks, n_peaks);
    
    % Necessary for parallel computing
    orientation_peaks = odf_verts_peaks;
    v1 = zeros([1, 3*max_n_peaks]);
    v2 = zeros([1, 3*max_n_peaks]);
    v3 = zeros([1, 3*max_n_peaks]);
    for i = 1:n_peaks
        orientation_peaks(i,:) = orientation_peaks(i,:)/norm_peaks(i); % Peak normalization
        ind = (3*(i-1)+1):(3*(i-1)+3);
        v1(ind) = orientation_peaks(i,:); % Peak
        v2(ind) = sddelta_peaks(i)*orientation_peaks(i,:); % Peak weigthted by its sddelta
        v3(ind) = w_peaks(i)*orientation_peaks(i,:); % Peak weigthted by its w
    end
    all_peaks(v, :) = single(v1); % Peak
    all_sddelta_weighted_peaks(v, :) = single(v2); % Peak weigthted by its sddelta
    all_w_weighted_peaks(v, :) = single(v3); % Peak weigthted by its w
    
    structure.n_peaks = n_peaks;
    structure.odf_norm = odf_norm;
    structure.vert_peaks = odf_verts_peaks;
    structure.orientation_peaks = orientation_peaks;
    structure.diso_peaks = diso_peaks;
    structure.sddelta_peaks = sddelta_peaks;
    structure.r1_peaks = r1_peaks;
    structure.t1_peaks = t1_peaks;
    structure.r2_peaks = r2_peaks;
    structure.t2_peaks = t2_peaks;
    structure.w_peaks = w_peaks;
    cell_odf_peaks{v} = structure;
end


%% Revert dimensionality reduction
dummy_map = zeros(Nx*Ny*Nz, 1);
dummy_map(mask>0) = nb_odf_peaks_map;
nb_odf_peaks_map = reshape(dummy_map, Nx, Ny, Nz);

dummy_map = zeros(Nx*Ny*Nz, 1);
dummy_map(mask>0) = max_nb_peaks_map;
max_nb_peaks_map = reshape(dummy_map, Nx, Ny, Nz);

dummy_map = cell([Nx*Ny*Nz, 1]);
dummy_ind = 1:Nx*Ny*Nz;
dummy_ind = dummy_ind(mask > 0);
for v = 1:N
    dummy_map{dummy_ind(v)} = cell_odf_peaks{v};
end
cell_odf_peaks = reshape(dummy_map, Nx, Ny, Nz);

dummy_peak = zeros(Nx*Ny*Nz, 3*max_n_peaks);
dummy_peak(mask>0, :) = all_peaks;
all_peaks = single(reshape(dummy_peak, Nx, Ny, Nz, 3*max_n_peaks));

dummy_peak = zeros(Nx*Ny*Nz, 3*max_n_peaks);
dummy_peak(mask>0, :) = all_sddelta_weighted_peaks;
all_sddelta_weighted_peaks = single(reshape(dummy_peak, Nx, Ny, Nz, 3*max_n_peaks));

dummy_peak = zeros(Nx*Ny*Nz, 3*max_n_peaks);
dummy_peak(mask>0, :) = all_w_weighted_peaks;
all_w_weighted_peaks = single(reshape(dummy_peak, Nx, Ny, Nz, 3*max_n_peaks));

%% Final saving 
mdm_nii_write(nb_odf_peaks_map, fullfile(odf_directory, 'nb_ODF_peaks.nii.gz'), h);
mdm_nii_write(max_nb_peaks_map, fullfile(odf_directory, 'nb_peaks_maximal.nii.gz'), h);

mdm_nii_write(all_peaks, fullfile(odf_directory, 'ODF_peaks.nii.gz'), h);
mdm_nii_write(all_sddelta_weighted_peaks, fullfile(odf_directory, 'ODF_peaks_weighted_by_sddelta.nii.gz'), h);
mdm_nii_write(all_w_weighted_peaks, fullfile(odf_directory, 'ODF_peaks_weighted_by_w.nii.gz'), h);
save(fullfile(odf_directory,'ODF_peaks_metrics.mat'), 'cell_odf_peaks');

end