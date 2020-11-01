function step2_bootstrap_to_maps_ODFs(input_parameters)

%% Load information from the input_parameters structure
method = input_parameters.inversion_method;
nb_MC_inversions = input_parameters.nb_MC_inversions;
data_directory = input_parameters.data_directory;
bootstrap_directory = input_parameters.bootstrap_directory;
parameter_maps_directory = input_parameters.parameter_maps_directory;
odf_directory = input_parameters.odf_directory;
framework_directory = input_parameters.framework_directory;
threshold_d_iso_big = input_parameters.threshold_d_iso_big;
threshold_d_delta_thin = input_parameters.threshold_d_delta_thin;
threshold_weight_thin = input_parameters.threshold_weight_thin;
mask_file = input_parameters.mask_file;
dir_flips = input_parameters.dir_flips;
nb_MeshNodes = input_parameters.nb_mesh_nodes;

FLAG_relaxation = false;
if strcmp(method,'dtd')
    nb_dimension = 4;
    Nparams = 16;
elseif strcmp(method,'dtr2d')
    nb_dimension = 5;
    Nparams = 23;
    FLAG_relaxation = true;
    char_relaxation = 'r2';
    char_relaxation2 = 't2';
    char_plot_relaxation = 'R_2';
    char_plot_relaxation2 = 'T_2';
elseif strcmp(method,'dtr1d')
    nb_dimension = 5;
    Nparams = 23;
    FLAG_relaxation = true;
    char_relaxation = 'r1';
    char_relaxation2 = 't1';
    char_plot_relaxation = 'R_1';
    char_plot_relaxation2 = 'T_1';
end

%% Define bins [big, thick, thin]
opt = mdm_opt(); %#ok<NASGU>
eval(['opt = ' method '_opt(opt);'])
eval(['opt_bis = opt.' method ';'])
opt_bis.bin_disomin = [threshold_d_iso_big 0 0]*1e-9; 
opt_bis.bin_disomax = [5 threshold_d_iso_big threshold_d_iso_big]*1e-9;
opt_bis.bin_sqddeltamin = [0 0 threshold_d_delta_thin^2]; % Thin bin has d_delta^2 > 0.5^2 = 0.25
opt_bis.bin_sqddeltamax = [1 threshold_d_delta_thin^2 1]; 
Nbins = numel(opt_bis.bin_disomax);
n_thin_bin = find(opt_bis.bin_sqddeltamin > 0, 1, 'first');

% diso_mid = 0.65;
% dratio = 4;
% sddelta = ((dratio-1)/(2+dratio))^2;
% opt_bis.bin_disomin = [diso_mid 0 0]*1e-9; 
% opt_bis.bin_disomax = [5 diso_mid 5]*1e-9;
% opt_bis.bin_sqddeltamin = [0 0 sddelta];
% opt_bis.bin_sqddeltamax = [sddelta sddelta 1]; 
% Nbins = numel(opt_bis.bin_disomax);
% n_thin_bin = 3;

%% Colormap bounds
diso_clim = input_parameters.clim.diso_clim;
msqddelta_clim = input_parameters.clim.msqddelta_clim;
vdiso_clim = input_parameters.clim.vdiso_clim;
vsqddelta_clim = input_parameters.clim.vsqddelta_clim;
cvdisosqddelta_clim = input_parameters.clim.cvdisosqddelta_clim;

if strcmp(method, 'dtr2d')
    mr_clim = input_parameters.clim.mr_clim;
    vr_clim = input_parameters.clim.vr_clim;
    cvdisor_clim = input_parameters.clim.cvdisor_clim;
    cvsqddeltar_clim = input_parameters.clim.cvsqddeltar_clim;    
    mt_clim = input_parameters.clim.mt_clim;
    mt_clim_blackandwhite = input_parameters.clim.mt_clim_blackandwhite;
    vt_clim = input_parameters.clim.vt_clim;
    cvdisot_clim = input_parameters.clim.cvdisot_clim;
    cvsqddeltat_clim = input_parameters.clim.cvsqddeltat_clim; 
elseif strcmp(method, 'dtr1d')
    mr_clim = input_parameters.clim.mr_clim;
    vr_clim = input_parameters.clim.vr_clim;
    cvdisor_clim = input_parameters.clim.cvdisor_clim;
    cvsqddeltar_clim = input_parameters.clim.cvsqddeltar_clim;    
    mt_clim = input_parameters.clim.mt_clim;
    mt_clim_blackandwhite = input_parameters.clim.mt_clim_blackandwhite;
    vt_clim = input_parameters.clim.vt_clim;
    cvdisot_clim = input_parameters.clim.cvdisot_clim;
    cvsqddeltat_clim = input_parameters.clim.cvsqddeltat_clim;  
end

%% Get basic info
mfs = mdm_mfs_load(fullfile(bootstrap_directory, '1', 'mfs.mat'));
m = mfs.m;
nii_h = mfs.nii_h;
if strcmp(mask_file, '')
    mask = mfs.mask;
else
    mask = mdm_nii_read(fullfile(data_directory, mask_file));
end
sz = size(m);
Nx = sz(1);
Ny = sz(2);
Nz = sz(3);
ind = false(sz(4),1);
ind(2:nb_dimension+1:end) = 1;
% nn = (sz(4)-1)/(nb_dimension+1);

%% Prepare median ODFs
angles_path = fullfile(framework_directory,'tools','uvec','repulsion_angles_tri');
angles = load(fullfile(angles_path, strcat(num2str(nb_MeshNodes), '.mat')));

% Find the minimal angular separation in the grid to set kappa
[nb_triangles, ~] = size(angles.tri);
minimal_angular_distances = zeros([1 3*nb_triangles]);
for i = 1:nb_triangles
    ind_triangle = (3*(i-1)+1):(3*(i-1)+3);
    ind_points = squeeze(angles.tri(i, :));
    ind_1 = [ind_points(1) ind_points(1) ind_points(2)];
    ind_2 = [ind_points(2) ind_points(3) ind_points(3)];
    minimal_angular_distances_on_triangles = zeros([1 3]);
    for j = 1:3
        minimal_angular_distances_on_triangles(j) = acos(abs(cos(angles.theta(ind_1(j))).*cos(angles.theta(ind_2(j))) + sin(angles.theta(ind_1(j))).*sin(angles.theta(ind_2(j))).*cos(angles.phi(ind_1(j))-angles.phi(ind_2(j)))));
    end
    minimal_angular_distances(ind_triangle) = minimal_angular_distances_on_triangles;
end
minimal_angular_separation = median(minimal_angular_distances);
factor = 1.5;
standard_deviation = factor*minimal_angular_separation;
kappa = 1/(2*standard_deviation^2);

odf_smooth.n = nb_MeshNodes;
odf_smooth.x = sin(angles.theta).*cos(angles.phi);
odf_smooth.y = sin(angles.theta).*sin(angles.phi);
odf_smooth.z = cos(angles.theta);
odf_smooth.tri = angles.tri;
odf_smooth.kappa = kappa; 

odf_bsmedian.n = nb_MeshNodes;
odf_bsmedian.x = sin(angles.theta).*cos(angles.phi);
odf_bsmedian.y = sin(angles.theta).*sin(angles.phi);
odf_bsmedian.z = cos(angles.theta);
odf_bsmedian.tri = angles.tri;
odf_bsmedian.minimal_angular_separation = minimal_angular_separation;
odf_bsmedian.kappa = odf_smooth.kappa;
odf_bsmedian.standard_deviation = standard_deviation;

%% Pre-read all of bootstrap solutions
m_global = zeros([nb_MC_inversions, Nx, Ny, Nz, sz(4)]);
parfor nBS = 1:nb_MC_inversions
    mfs = mdm_mfs_load(fullfile(fullfile(bootstrap_directory,num2str(nBS)), 'mfs.mat'));
    m_global(nBS,:,:,:,:) = mfs.m;
end

% Decrease dimensionality to lower memory cost
[~, ~, ~, ~, size_m] = size(m_global);
m_global = reshape(m_global, nb_MC_inversions, Nx*Ny*Nz, size_m);
mask = reshape(mask, Nx*Ny*Nz, 1)';

% Retain only non-masked-out voxels.
m_global = m_global(:, mask>0, :);
[~, N, ~] = size(m_global);

%% Compute medians across bootstrap solutions
s0 = zeros(N, 1);
mdiso = zeros(N, 1);
msqddelta = zeros(N, 1);
vdiso = zeros(N, 1);
vsqddelta = zeros(N, 1);
cvdisosqddelta = zeros(N, 1);

f_bin = zeros(Nbins, N);
mdxx_bin = zeros(Nbins, N);
mdyy_bin = zeros(Nbins, N);
mdzz_bin = zeros(Nbins, N);
mdiso_bin = zeros(Nbins, N);
msqddelta_bin = zeros(Nbins, N);

odf_w = zeros(N, nb_MeshNodes);
odf_diso = zeros(N, nb_MeshNodes);
odf_sqddelta = zeros(N, nb_MeshNodes);

if FLAG_relaxation
    mr = zeros(N, 1);
    mt = zeros(N, 1);
    vr = zeros(N, 1);
    vt = zeros(N, 1);
    cvdisor = zeros(N, 1);
    cvdisot = zeros(N, 1);
    cvsqddeltar = zeros(N, 1);
    cvsqddeltat = zeros(N, 1);

    mr_bin = zeros(Nbins, N);
    mt_bin = zeros(Nbins, N);
    
    odf_r = zeros(N, nb_MeshNodes);
    odf_t = zeros(N, nb_MeshNodes);
end

if FLAG_relaxation
    parfor v = 1:N      
        temp_s0 = zeros(1, nb_MC_inversions);
        temp_mdiso = zeros(1, nb_MC_inversions);
        temp_msqddelta = zeros(1, nb_MC_inversions);
        temp_mr = zeros(1, nb_MC_inversions);
        temp_mt = zeros(1, nb_MC_inversions);
        temp_vdiso = zeros(1, nb_MC_inversions);
        temp_vsqddelta = zeros(1, nb_MC_inversions);
        temp_vr = zeros(1, nb_MC_inversions);
        temp_vt = zeros(1, nb_MC_inversions);
        temp_cvdisosqddelta = zeros(1, nb_MC_inversions);
        temp_cvdisor = zeros(1, nb_MC_inversions);
        temp_cvdisot = zeros(1, nb_MC_inversions);
        temp_cvsqddeltar = zeros(1, nb_MC_inversions);
        temp_cvsqddeltat = zeros(1, nb_MC_inversions);
        
        temp_f_bin = zeros(Nbins, nb_MC_inversions);
        temp_mdxx_bin = zeros(Nbins, nb_MC_inversions);
        temp_mdyy_bin = zeros(Nbins, nb_MC_inversions);
        temp_mdzz_bin = zeros(Nbins, nb_MC_inversions);
        temp_mdiso_bin = zeros(Nbins, nb_MC_inversions);
        temp_msqddelta_bin = zeros(Nbins, nb_MC_inversions);
        temp_mr_bin = zeros(Nbins, nb_MC_inversions);
        temp_mt_bin = zeros(Nbins, nb_MC_inversions);
        
        temp_odf_w = zeros(nb_MeshNodes, nb_MC_inversions);
        temp_odf_diso = zeros(nb_MeshNodes, nb_MC_inversions);
        temp_odf_sqddelta = zeros(nb_MeshNodes, nb_MC_inversions);
        temp_odf_r = zeros(nb_MeshNodes, nb_MC_inversions);
        temp_odf_t = zeros(nb_MeshNodes, nb_MC_inversions);
        
        for nBS = 1:nb_MC_inversions
            m = squeeze(m_global(nBS,v,:));
            dpar = m(circshift(ind,0,1));
            dperp = m(circshift(ind,1,1));
            theta = m(circshift(ind,2,1));
            phi = m(circshift(ind,3,1));
            w = m(circshift(ind,nb_dimension,1));
            r = m(circshift(ind,4,1));
            t = msf_notfinite2zero(1./r);
            if strcmp(method,'dtr2d') && isfield(opt_bis,'r2extrap') == 1
                w = w.*exp(-r*opt_bis.r2extrap);
            end
            
            dpar = dpar(w>0);
            dperp = dperp(w>0);
            theta = theta(w>0);
            phi = phi(w>0);
            w = w(w>0);
            r = r(w>0);
            t = t(w>0);
            
            diso = (dpar + 2*dperp)/3;
            daniso = (dpar - dperp)/3;
            ddelta = msf_notfinite2zero(daniso./diso);
            sqdaniso = daniso.^2;
            sqddelta = ddelta.^2;
            dratio = msf_notfinite2zero(dpar./dperp);
            [dxx,dyy,dzz,dxy,dxz,dyz] = dtd_pars2elements(dpar,dperp,theta,phi);
            
            dtds = struct('w',w,'dpar',dpar,'dperp',dperp,'theta',theta,'phi',phi,'diso',diso,'daniso',daniso,'ddelta',ddelta,...
                'sqdaniso',sqdaniso,'sqddelta',sqddelta,'dratio',dratio,'dxx',dxx,'dyy',dyy,'dzz',dzz,'dxy',dxy,'dxz',dxz,'dyz',dyz,'r',r,'t',t);
            
            f = w./sum(w);
            
            temp_s0(nBS) = sum(w);
            temp_mdiso(nBS) = sum(f.*diso);
            temp_msqddelta(nBS) = sum(f.*sqddelta);
            temp_mr(nBS) = sum(f.*r);
            temp_mt(nBS) = sum(f.*t);
            temp_vdiso(nBS) = sum(f.*(diso-temp_mdiso(nBS)).^2);
            temp_vsqddelta(nBS) = sum(f.*(sqddelta-temp_msqddelta(nBS)).^2);
            temp_vr(nBS) = sum(f.*(r-temp_mr(nBS)).^2);
            temp_vt(nBS) = sum(f.*(t-temp_mt(nBS)).^2);
            temp_cvdisosqddelta(nBS) = sum(f.*(diso-temp_mdiso(nBS)).*(sqddelta-temp_msqddelta(nBS)));
            temp_cvdisor(nBS) = sum(f.*(diso-temp_mdiso(nBS)).*(r-temp_mr(nBS)));
            temp_cvdisot(nBS) = sum(f.*(diso-temp_mdiso(nBS)).*(t-temp_mt(nBS)));
            temp_cvsqddeltar(nBS) = sum(f.*(sqddelta-temp_msqddelta(nBS)).*(r-temp_mr(nBS)));
            temp_cvsqddeltat(nBS) = sum(f.*(sqddelta-temp_msqddelta(nBS)).*(t-temp_mt(nBS)));
            
            for nbin = 1:Nbins
                ind_bin = false([length(f) 4]);
                ind_bin(:,1) = diso >= opt_bis.bin_disomin(nbin);
                ind_bin(:,2) = diso <= opt_bis.bin_disomax(nbin);
                ind_bin(:,3) = sqddelta >= opt_bis.bin_sqddeltamin(nbin);
                ind_bin(:,4) = sqddelta <= opt_bis.bin_sqddeltamax(nbin);
                ind_bin = all(ind_bin,2);
                
                if nnz(ind_bin) == 0
                    temp_f_bin(nbin,nBS) = 0;
                    temp_mdxx_bin(nbin,nBS) = NaN;
                    temp_mdyy_bin(nbin,nBS) = NaN;
                    temp_mdzz_bin(nbin,nBS) = NaN;
                    temp_mdiso_bin(nbin,nBS) = NaN;
                    temp_msqddelta_bin(nbin,nBS) = NaN;
                    temp_mr_bin(nbin,nBS) = NaN;
                    temp_mt_bin(nbin,nBS) = NaN;
                else
                    temp_f_bin(nbin,nBS) = sum(f.*ind_bin);
                    temp_mdxx_bin(nbin,nBS) = sum(f.*dxx.*ind_bin)/sum(f.*ind_bin);
                    temp_mdyy_bin(nbin,nBS) = sum(f.*dyy.*ind_bin)/sum(f.*ind_bin);
                    temp_mdzz_bin(nbin,nBS) = sum(f.*dzz.*ind_bin)/sum(f.*ind_bin);
                    temp_mdiso_bin(nbin,nBS) = sum(f.*diso.*ind_bin)/sum(f.*ind_bin);
                    temp_msqddelta_bin(nbin,nBS) = sum(f.*sqddelta.*ind_bin)/sum(f.*ind_bin);
                    temp_mr_bin(nbin,nBS) = sum(f.*r.*ind_bin)/sum(f.*ind_bin);
                    temp_mt_bin(nbin,nBS) = sum(f.*t.*ind_bin)/sum(f.*ind_bin);
                end
                
                dtds_temp = dtds;
                dtds_temp.w = dtds.w.*ind_bin;
                
                if nbin == n_thin_bin
                    if nnz(ind_bin) == 0 % No bin solution found for that specific bootstrap realization
                        temp_odf_w(:,nBS) = zeros([1, nb_MeshNodes]);
                        temp_odf_diso(:,nBS) = NaN([1, nb_MeshNodes]);
                        temp_odf_sqddelta(:,nBS) = NaN([1, nb_MeshNodes]);
                        temp_odf_r(:,nBS) = NaN([1, nb_MeshNodes]);
                        temp_odf_t(:,nBS) = NaN([1, nb_MeshNodes]);
                    else
                        
                        % Discrete ODF
                        odf_discrete = struct;
                        odf_discrete.x = (-1)^dir_flips(1).*sin(dtds_temp.theta).*cos(dtds_temp.phi);
                        odf_discrete.y = (-1)^dir_flips(2).*sin(dtds_temp.theta).*sin(dtds_temp.phi);
                        odf_discrete.z = (-1)^dir_flips(3).*cos(dtds_temp.theta);
                        odf_discrete.w = dtds_temp.w;
                        odf_discrete.diso = dtds_temp.diso;
                        odf_discrete.sqddelta = dtds_temp.sqddelta;
                        odf_discrete.r = dtds_temp.r;
                        odf_discrete.t = dtds_temp.t;
                        
                        % Projection of the discrete ODF onto the smooth grid
                        odf_vox = dist_podf_discrete2smooth(odf_discrete, odf_smooth);
                        temp_odf_w(:,nBS) = odf_vox.w(:);
                        temp_odf_diso(:,nBS) = odf_vox.diso(:);
                        temp_odf_sqddelta(:,nBS) = odf_vox.sqddelta(:);
                        temp_odf_r(:,nBS) = odf_vox.r(:);
                        temp_odf_t(:,nBS) = odf_vox.t(:);
                    end
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        s0(v) = msf_notfinite2zero(nanmedian(temp_s0));
        mdiso(v) = msf_notfinite2zero(nanmedian(temp_mdiso));
        msqddelta(v) = msf_notfinite2zero(nanmedian(temp_msqddelta));
        mr(v) = msf_notfinite2zero(nanmedian(temp_mr));
        mt(v) = msf_notfinite2zero(nanmedian(temp_mt));
        vdiso(v) = msf_notfinite2zero(nanmedian(temp_vdiso));
        vsqddelta(v) = msf_notfinite2zero(nanmedian(temp_vsqddelta));
        vr(v) = msf_notfinite2zero(nanmedian(temp_vr));
        vt(v) = msf_notfinite2zero(nanmedian(temp_vt));
        cvdisosqddelta(v) = msf_notfinite2zero(nanmedian(temp_cvdisosqddelta));
        cvdisor(v) = msf_notfinite2zero(nanmedian(temp_cvdisor));
        cvdisot(v) = msf_notfinite2zero(nanmedian(temp_cvdisot));
        cvsqddeltar(v) = msf_notfinite2zero(nanmedian(temp_cvsqddeltar));
        cvsqddeltat(v) = msf_notfinite2zero(nanmedian(temp_cvsqddeltat));
        
        list_f_bin = zeros(Nbins, 1);
        list_mdiso_bin = zeros(Nbins, 1);
        list_msqddelta_bin = zeros(Nbins, 1);
        list_mdxx_bin = zeros(Nbins, 1);
        list_mdyy_bin = zeros(Nbins, 1);
        list_mdzz_bin = zeros(Nbins, 1);
        list_mr_bin = zeros(Nbins, 1);
        list_mt_bin = zeros(Nbins, 1);
        for nbin = 1:Nbins
            list_f_bin(nbin) = msf_notfinite2zero(nanmedian(squeeze(temp_f_bin(nbin,:))));
            list_mdiso_bin(nbin) = msf_notfinite2zero(nanmedian(squeeze(temp_mdiso_bin(nbin,:))));
            list_msqddelta_bin(nbin) = msf_notfinite2zero(nanmedian(squeeze(temp_msqddelta_bin(nbin,:))));
            list_mdxx_bin(nbin) = msf_notfinite2zero(nanmedian(squeeze(temp_mdxx_bin(nbin,:))));
            list_mdyy_bin(nbin) = msf_notfinite2zero(nanmedian(squeeze(temp_mdyy_bin(nbin,:))));
            list_mdzz_bin(nbin) = msf_notfinite2zero(nanmedian(squeeze(temp_mdzz_bin(nbin,:))));
            list_mr_bin(nbin) = msf_notfinite2zero(nanmedian(squeeze(temp_mr_bin(nbin,:))));
            list_mt_bin(nbin) = msf_notfinite2zero(nanmedian(squeeze(temp_mt_bin(nbin,:))));
        end
        f_bin(:,v) = list_f_bin;
        mdiso_bin(:,v) = list_mdiso_bin;
        msqddelta_bin(:,v) = list_msqddelta_bin;
        mdxx_bin(:,v) = list_mdxx_bin;
        mdyy_bin(:,v) = list_mdyy_bin;
        mdzz_bin(:,v) = list_mdzz_bin;
        mr_bin(:,v) = list_mr_bin;
        mt_bin(:,v) = list_mt_bin;
        
        odf_w(v,:) = msf_notfinite2zero(nanmedian(temp_odf_w,2));
        odf_diso(v,:) = msf_notfinite2zero(nanmedian(temp_odf_diso,2));
        odf_sqddelta(v,:) = msf_notfinite2zero(nanmedian(temp_odf_sqddelta,2));
        odf_r(v,:) = msf_notfinite2zero(nanmedian(temp_odf_r,2));
        odf_t(v,:) = msf_notfinite2zero(nanmedian(temp_odf_t,2));
    end
else
    parfor v =  1:N     
        temp_s0 = zeros(1, nb_MC_inversions);
        temp_mdiso = zeros(1, nb_MC_inversions);
        temp_msqddelta = zeros(1, nb_MC_inversions);
        temp_vdiso = zeros(1, nb_MC_inversions);
        temp_vsqddelta = zeros(1, nb_MC_inversions);
        temp_cvdisosqddelta = zeros(1, nb_MC_inversions);
        
        temp_f_bin = zeros(Nbins, nb_MC_inversions);
        temp_mdxx_bin = zeros(Nbins, nb_MC_inversions);
        temp_mdyy_bin = zeros(Nbins, nb_MC_inversions);
        temp_mdzz_bin = zeros(Nbins, nb_MC_inversions);
        temp_mdiso_bin = zeros(Nbins, nb_MC_inversions);
        temp_msqddelta_bin = zeros(Nbins, nb_MC_inversions);
        
        temp_odf_w = zeros(nb_MeshNodes, nb_MC_inversions);
        temp_odf_diso = zeros(nb_MeshNodes, nb_MC_inversions);
        temp_odf_sqddelta = zeros(nb_MeshNodes, nb_MC_inversions);
        
        for nBS = 1:nb_MC_inversions
            m = squeeze(m_global(nBS,v,:));
            dpar = m(circshift(ind,0,1));
            dperp = m(circshift(ind,1,1));
            theta = m(circshift(ind,2,1));
            phi = m(circshift(ind,3,1));
            w = m(circshift(ind,nb_dimension,1));
            
            dpar = dpar(w>0);
            dperp = dperp(w>0);
            theta = theta(w>0);
            phi = phi(w>0);
            w = w(w>0);
            
            diso = (dpar + 2*dperp)/3;
            daniso = (dpar - dperp)/3;
            ddelta = msf_notfinite2zero(daniso./diso);
            sqdaniso = daniso.^2;
            sqddelta = ddelta.^2;
            dratio = msf_notfinite2zero(dpar./dperp);
            [dxx,dyy,dzz,dxy,dxz,dyz] = dtd_pars2elements(dpar,dperp,theta,phi);
            
            dtds = struct('w',w,'dpar',dpar,'dperp',dperp,'theta',theta,'phi',phi,'diso',diso,'daniso',daniso,'ddelta',ddelta,...
                'sqdaniso',sqdaniso,'sqddelta',sqddelta,'dratio',dratio,'dxx',dxx,'dyy',dyy,'dzz',dzz,'dxy',dxy,'dxz',dxz,'dyz',dyz);
            
            f = w./sum(w);
            
            temp_s0(nBS) = sum(w);
            temp_mdiso(nBS) = sum(f.*diso);
            temp_msqddelta(nBS) = sum(f.*sqddelta);
            temp_vdiso(nBS) = sum(f.*(diso-temp_mdiso(nBS)).^2);
            temp_vsqddelta(nBS) = sum(f.*(sqddelta-temp_msqddelta(nBS)).^2);
            temp_cvdisosqddelta(nBS) = sum(f.*(diso-temp_mdiso(nBS)).*(sqddelta-temp_msqddelta(nBS)));
            
            for nbin = 1:Nbins
                ind_bin = false([length(f) 4]);
                ind_bin(:,1) = diso >= opt_bis.bin_disomin(nbin);
                ind_bin(:,2) = diso <= opt_bis.bin_disomax(nbin);
                ind_bin(:,3) = sqddelta >= opt_bis.bin_sqddeltamin(nbin);
                ind_bin(:,4) = sqddelta <= opt_bis.bin_sqddeltamax(nbin);
                ind_bin = all(ind_bin,2);
                
                if nnz(ind_bin) == 0
                    temp_f_bin(nbin,nBS) = 0;
                    temp_mdxx_bin(nbin,nBS) = NaN;
                    temp_mdyy_bin(nbin,nBS) = NaN;
                    temp_mdzz_bin(nbin,nBS) = NaN;
                    temp_mdiso_bin(nbin,nBS) = NaN;
                    temp_msqddelta_bin(nbin,nBS) = NaN;
                else
                    temp_f_bin(nbin,nBS) = sum(f.*ind_bin);
                    temp_mdxx_bin(nbin,nBS) = sum(f.*dxx.*ind_bin)/sum(f.*ind_bin);
                    temp_mdyy_bin(nbin,nBS) = sum(f.*dyy.*ind_bin)/sum(f.*ind_bin);
                    temp_mdzz_bin(nbin,nBS) = sum(f.*dzz.*ind_bin)/sum(f.*ind_bin);
                    temp_mdiso_bin(nbin,nBS) = sum(f.*diso.*ind_bin)/sum(f.*ind_bin);
                    temp_msqddelta_bin(nbin,nBS) = sum(f.*sqddelta.*ind_bin)/sum(f.*ind_bin);
                end
                
                dtds_temp = dtds;
                dtds_temp.w = dtds.w.*ind_bin;
                
                if nbin == n_thin_bin
                    if nnz(ind_bin) == 0 % No bin solution found for that specific bootstrap realization
                        temp_odf_w(:,nBS) = zeros([1, nb_MeshNodes]);
                        temp_odf_diso(:,nBS) = NaN([1, nb_MeshNodes]);
                        temp_odf_sqddelta(:,nBS) = NaN([1, nb_MeshNodes]);
                    else
                        
                        % Discrete ODF
                        odf_discrete = struct;
                        odf_discrete.x = (-1)^dir_flips(1).*sin(dtds_temp.theta).*cos(dtds_temp.phi);
                        odf_discrete.y = (-1)^dir_flips(2).*sin(dtds_temp.theta).*sin(dtds_temp.phi);
                        odf_discrete.z = (-1)^dir_flips(3).*cos(dtds_temp.theta);
                        odf_discrete.w = dtds_temp.w;
                        odf_discrete.diso = dtds_temp.diso;
                        odf_discrete.sqddelta = dtds_temp.sqddelta;
                        
                        % Projection of the discrete ODF onto the smooth grid
                        odf_vox = dist_podf_discrete2smooth(odf_discrete, odf_smooth);
                        temp_odf_w(:,nBS) = odf_vox.w(:);
                        temp_odf_diso(:,nBS) = odf_vox.diso(:);
                        temp_odf_sqddelta(:,nBS) = odf_vox.sqddelta(:);
                    end
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        s0(v) = msf_notfinite2zero(nanmedian(temp_s0));
        mdiso(v) = msf_notfinite2zero(nanmedian(temp_mdiso));
        msqddelta(v) = msf_notfinite2zero(nanmedian(temp_msqddelta));
        vdiso(v) = msf_notfinite2zero(nanmedian(temp_vdiso));
        vsqddelta(v) = msf_notfinite2zero(nanmedian(temp_vsqddelta));
        cvdisosqddelta(v) = msf_notfinite2zero(nanmedian(temp_cvdisosqddelta));
        
        list_f_bin = zeros(Nbins, 1);
        list_mdiso_bin = zeros(Nbins, 1);
        list_msqddelta_bin = zeros(Nbins, 1);
        list_mdxx_bin = zeros(Nbins, 1);
        list_mdyy_bin = zeros(Nbins, 1);
        list_mdzz_bin = zeros(Nbins, 1);
        for nbin = 1:Nbins
            list_f_bin(nbin) = msf_notfinite2zero(nanmedian(squeeze(temp_f_bin(nbin,:))));
            list_mdiso_bin(nbin) = msf_notfinite2zero(nanmedian(squeeze(temp_mdiso_bin(nbin,:))));
            list_msqddelta_bin(nbin) = msf_notfinite2zero(nanmedian(squeeze(temp_msqddelta_bin(nbin,:))));
            list_mdxx_bin(nbin) = msf_notfinite2zero(nanmedian(squeeze(temp_mdxx_bin(nbin,:))));
            list_mdyy_bin(nbin) = msf_notfinite2zero(nanmedian(squeeze(temp_mdyy_bin(nbin,:))));
            list_mdzz_bin(nbin) = msf_notfinite2zero(nanmedian(squeeze(temp_mdzz_bin(nbin,:))));
        end
        f_bin(:,v) = list_f_bin;
        mdiso_bin(:,v) = list_mdiso_bin;
        msqddelta_bin(:,v) = list_msqddelta_bin;
        mdxx_bin(:,v) = list_mdxx_bin;
        mdyy_bin(:,v) = list_mdyy_bin;
        mdzz_bin(:,v) = list_mdzz_bin;
        
        odf_w(v,:) = msf_notfinite2zero(nanmedian(temp_odf_w,2));
        odf_diso(v,:) = msf_notfinite2zero(nanmedian(temp_odf_diso,2));
        odf_sqddelta(v,:) = msf_notfinite2zero(nanmedian(temp_odf_sqddelta,2));
    end
end

%% Revert dimensionality reduction
dummy_map = zeros(Nx*Ny*Nz, 1);
dummy_map(mask>0) = s0;
s0 = reshape(dummy_map, Nx, Ny, Nz);

dummy_map = zeros(Nx*Ny*Nz, 1);
dummy_map(mask>0) = mdiso;
mdiso = reshape(dummy_map, Nx, Ny, Nz);

dummy_map = zeros(Nx*Ny*Nz, 1);
dummy_map(mask>0) = msqddelta;
msqddelta = reshape(dummy_map, Nx, Ny, Nz);

dummy_map = zeros(Nx*Ny*Nz, 1);
dummy_map(mask>0) = vdiso;
vdiso = reshape(dummy_map, Nx, Ny, Nz);

dummy_map = zeros(Nx*Ny*Nz, 1);
dummy_map(mask>0) = vsqddelta;
vsqddelta = reshape(dummy_map, Nx, Ny, Nz);

dummy_map = zeros(Nx*Ny*Nz, 1);
dummy_map(mask>0) = cvdisosqddelta;
cvdisosqddelta = reshape(dummy_map, Nx, Ny, Nz);

if FLAG_relaxation
    dummy_map = zeros(Nx*Ny*Nz, 1);
    dummy_map(mask>0) = mr;
    mr = reshape(dummy_map, Nx, Ny, Nz);
    
    dummy_map = zeros(Nx*Ny*Nz, 1);
    dummy_map(mask>0) = mt;
    mt = reshape(dummy_map, Nx, Ny, Nz);
    
    dummy_map = zeros(Nx*Ny*Nz, 1);
    dummy_map(mask>0) = vr;
    vr = reshape(dummy_map, Nx, Ny, Nz);
    
    dummy_map = zeros(Nx*Ny*Nz, 1);
    dummy_map(mask>0) = vt;
    vt = reshape(dummy_map, Nx, Ny, Nz);
    
    dummy_map = zeros(Nx*Ny*Nz, 1);
    dummy_map(mask>0) = cvdisor;
    cvdisor = reshape(dummy_map, Nx, Ny, Nz);
    
    dummy_map = zeros(Nx*Ny*Nz, 1);
    dummy_map(mask>0) = cvdisot;
    cvdisot = reshape(dummy_map, Nx, Ny, Nz);
    
    dummy_map = zeros(Nx*Ny*Nz, 1);
    dummy_map(mask>0) = cvsqddeltar;
    cvsqddeltar = reshape(dummy_map, Nx, Ny, Nz);
    
    dummy_map = zeros(Nx*Ny*Nz, 1);
    dummy_map(mask>0) = cvsqddeltat;
    cvsqddeltat = reshape(dummy_map, Nx, Ny, Nz);
end

dummy_bin = zeros(Nbins, Nx*Ny*Nz);
dummy_bin(:, mask>0) = f_bin;
f_bin = reshape(dummy_bin, Nbins, Nx, Ny, Nz);

dummy_bin = zeros(Nbins, Nx*Ny*Nz);
dummy_bin(:, mask>0) = mdiso_bin;
mdiso_bin = reshape(dummy_bin, Nbins, Nx, Ny, Nz);

dummy_bin = zeros(Nbins, Nx*Ny*Nz);
dummy_bin(:, mask>0) = msqddelta_bin;
msqddelta_bin = reshape(dummy_bin, Nbins, Nx, Ny, Nz);

dummy_bin = zeros(Nbins, Nx*Ny*Nz);
dummy_bin(:, mask>0) = mdxx_bin;
mdxx_bin = reshape(dummy_bin, Nbins, Nx, Ny, Nz);

dummy_bin = zeros(Nbins, Nx*Ny*Nz);
dummy_bin(:, mask>0) = mdyy_bin;
mdyy_bin = reshape(dummy_bin, Nbins, Nx, Ny, Nz);

dummy_bin = zeros(Nbins, Nx*Ny*Nz);
dummy_bin(:, mask>0) = mdzz_bin;
mdzz_bin = reshape(dummy_bin, Nbins, Nx, Ny, Nz);

if FLAG_relaxation
    dummy_bin = zeros(Nbins, Nx*Ny*Nz);
    dummy_bin(:, mask>0) = mr_bin;
    mr_bin = reshape(dummy_bin, Nbins, Nx, Ny, Nz);
    
    dummy_bin = zeros(Nbins, Nx*Ny*Nz);
    dummy_bin(:, mask>0) = mt_bin;
    mt_bin = reshape(dummy_bin, Nbins, Nx, Ny, Nz);
end


dummy_odf = zeros(Nx*Ny*Nz, nb_MeshNodes);
dummy_odf(mask>0, :) = odf_w;
odf_w = reshape(dummy_odf, Nx, Ny, Nz, nb_MeshNodes);

dummy_odf = zeros(Nx*Ny*Nz, nb_MeshNodes);
dummy_odf(mask>0, :) = odf_diso;
odf_diso = reshape(dummy_odf, Nx, Ny, Nz, nb_MeshNodes);

dummy_odf = zeros(Nx*Ny*Nz, nb_MeshNodes);
dummy_odf(mask>0, :) = odf_sqddelta;
odf_sqddelta = reshape(dummy_odf, Nx, Ny, Nz, nb_MeshNodes);

if FLAG_relaxation
    dummy_odf = zeros(Nx*Ny*Nz, nb_MeshNodes);
    dummy_odf(mask>0, :) = odf_r;
    odf_r = reshape(dummy_odf, Nx, Ny, Nz, nb_MeshNodes);
    
    dummy_odf = zeros(Nx*Ny*Nz, nb_MeshNodes);
    dummy_odf(mask>0, :) = odf_t;
    odf_t = reshape(dummy_odf, Nx, Ny, Nz, nb_MeshNodes);
end

f_comp1 = squeeze(f_bin(1,:,:,:));
f_comp2 = squeeze(f_bin(2,:,:,:));
f_comp3 = squeeze(f_bin(3,:,:,:));
mdiso_comp1 = squeeze(mdiso_bin(1,:,:,:));
mdiso_comp2 = squeeze(mdiso_bin(2,:,:,:));
mdiso_comp3 = squeeze(mdiso_bin(3,:,:,:));
msqddelta_comp1 = squeeze(msqddelta_bin(1,:,:,:));
msqddelta_comp2 = squeeze(msqddelta_bin(2,:,:,:));
msqddelta_comp3 = squeeze(msqddelta_bin(3,:,:,:));
mdxx_comp1 = squeeze(mdxx_bin(1,:,:,:));
mdxx_comp2 = squeeze(mdxx_bin(2,:,:,:));
mdxx_comp3 = squeeze(mdxx_bin(3,:,:,:));
mdyy_comp1 = squeeze(mdyy_bin(1,:,:,:));
mdyy_comp2 = squeeze(mdyy_bin(2,:,:,:));
mdyy_comp3 = squeeze(mdyy_bin(3,:,:,:));
mdzz_comp1 = squeeze(mdzz_bin(1,:,:,:));
mdzz_comp2 = squeeze(mdzz_bin(2,:,:,:));
mdzz_comp3 = squeeze(mdzz_bin(3,:,:,:));

if FLAG_relaxation
   mr_comp1 = squeeze(mr_bin(1,:,:,:));
   mr_comp2 = squeeze(mr_bin(2,:,:,:));
   mr_comp3 = squeeze(mr_bin(3,:,:,:));
   mt_comp1 = squeeze(mt_bin(1,:,:,:));
   mt_comp2 = squeeze(mt_bin(2,:,:,:));
   mt_comp3 = squeeze(mt_bin(3,:,:,:));
end

%% Save ODFs
odf_bsmedian.w_bin{1} = odf_w;
odf_bsmedian.diso_bin{1} = odf_diso;
odf_bsmedian.sqddelta_bin{1} = odf_sqddelta;
if strcmp(method,'dtr2d')
    odf_bsmedian.r2_bin{1} = odf_r;
    odf_bsmedian.t2_bin{1} = odf_t;
elseif strcmp(method,'dtr1d')
    odf_bsmedian.r1_bin{1} = odf_r;
    odf_bsmedian.t1_bin{1} = odf_t;
end

odf_bsmedian_fn = fullfile(odf_directory,strcat('odf_bsmedian_',num2str(nb_MeshNodes)));
msf_mkdir(fileparts(odf_bsmedian_fn));
save(odf_bsmedian_fn, 'odf_bsmedian','-v7.3');

%% Write nifti maps
smax = max(s0(:));

if strcmp(mask_file, '')
    mfs = mdm_mfs_load(fullfile(bootstrap_directory, '1', 'mfs.mat'));
    mask = mfs.mask;
else
    mask = mdm_nii_read(fullfile(data_directory, mask_file));
end
mask = double(mask);

mdm_nii_write(double(mask),fullfile(parameter_maps_directory,strcat(method, '_mask.nii.gz')),nii_h,0);
mdm_nii_write(mask.*s0,fullfile(parameter_maps_directory,strcat(method, '_s0.nii.gz')),nii_h,0);
mdm_nii_write(mask.*mdiso,fullfile(parameter_maps_directory,strcat(method, '_mdiso.nii.gz')),nii_h,0);
mdm_nii_write(mask.*msqddelta,fullfile(parameter_maps_directory,strcat(method, '_msqddelta.nii.gz')),nii_h,0);
mdm_nii_write(mask.*vdiso,fullfile(parameter_maps_directory,strcat(method, '_vdiso.nii.gz')),nii_h,0);
mdm_nii_write(mask.*vsqddelta,fullfile(parameter_maps_directory,strcat(method, '_vsqddelta.nii.gz')),nii_h,0);
mdm_nii_write(mask.*cvdisosqddelta,fullfile(parameter_maps_directory,strcat(method, '_cvdisosqddelta.nii.gz')),nii_h,0);
if FLAG_relaxation
    mdm_nii_write(mask.*mr,fullfile(parameter_maps_directory,strcat(method, '_m', char_relaxation, '.nii.gz')),nii_h,0);
    mdm_nii_write(mask.*mt,fullfile(parameter_maps_directory,strcat(method, '_m', char_relaxation2, '.nii.gz')),nii_h,0);
    mdm_nii_write(mask.*vr,fullfile(parameter_maps_directory,strcat(method, '_v', char_relaxation, '.nii.gz')),nii_h,0);
    mdm_nii_write(mask.*vt,fullfile(parameter_maps_directory,strcat(method, '_v', char_relaxation2, '.nii.gz')),nii_h,0);
    mdm_nii_write(mask.*cvsqddeltar,fullfile(parameter_maps_directory,strcat(method, '_cvsqddelta', char_relaxation, '.nii.gz')),nii_h,0);
    mdm_nii_write(mask.*cvsqddeltat,fullfile(parameter_maps_directory,strcat(method, '_cvsqddelta', char_relaxation2, '.nii.gz')),nii_h,0);
    mdm_nii_write(mask.*cvdisor,fullfile(parameter_maps_directory,strcat(method, '_cvdiso', char_relaxation, '.nii.gz')),nii_h,0);
    mdm_nii_write(mask.*cvdisot,fullfile(parameter_maps_directory,strcat(method, '_cvdiso', char_relaxation2, '.nii.gz')),nii_h,0);
end

mdm_nii_write(mask.*f_comp3, fullfile(parameter_maps_directory,strcat(method, '_fraction_thin.nii.gz')),nii_h,0); %#ok<*USENS>
mdm_nii_write(mask.*f_comp2, fullfile(parameter_maps_directory,strcat(method, '_fraction_thick.nii.gz')),nii_h,0);
mdm_nii_write(mask.*f_comp1, fullfile(parameter_maps_directory,strcat(method, '_fraction_big.nii.gz')),nii_h,0);
mdm_nii_write(mask.*(f_comp1 + f_comp2), fullfile(parameter_maps_directory,strcat(method, '_fraction_not_thin.nii.gz')),nii_h,0);
mdm_nii_write(double(mask.*f_comp3 >= threshold_weight_thin), fullfile(parameter_maps_directory,strcat(method, '_thin_mask.nii.gz')),nii_h,0);

clear dxx dyy dzz daniso danison sdaniso sqddelta


%% Plot parameter maps

pixaspect = nii_h.pixdim(3)/nii_h.pixdim(2);

f = figure('Visible', 'off');

if length(size(mask)) == 2 % One slice
    nk = 1;
    Nslices = 1;
else
    is_something = mask(:,:,:) > 0;
    is_something_z = logical(squeeze(sum(squeeze(sum(is_something,1)),1)));
    min_z = find(is_something_z, 1, 'first');
    max_z = find(is_something_z, 1, 'last');
    nk = min_z:max_z;
    Nslices = numel(nk);
end

papersize = 3*[Nparams (Nslices+1)*pixaspect];
fontsize_text = 15;

position.dbottom = 1/(Nslices+1);
dleft = 1/Nparams;
position.height = 1.01*position.dbottom;
position.width = 1.01*dleft;
position.left = 0;

%% s0
im3d = s0;
clim = 0.4*smax*[0 1];
axh_v = mplot_slicescolumn(cat(3,mask(:,:,nk).*im3d(:,:,nk),zeros([sz(1),sz(2)])),position,clim);
text(axh_v(end), 0.5, 0.5, '$S_0$', 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')

%% mdiso
clear im3d
im3d = mdiso;
clim = diso_clim;
position.left = position.left + dleft;
axh_v = mplot_slicescolumn(cat(3,mask(:,:,nk).*im3d(:,:,nk),zeros([sz(1),sz(2)])),position,clim);
text(axh_v(end), 0.5, 0.5, '$\mathrm{E}[D_\mathrm{iso}]$', 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')

%% msqddelta
clear im3d
im3d = msqddelta;
clim = msqddelta_clim;
position.left = position.left + dleft;
axh_v = mplot_slicescolumn(cat(3,mask(:,:,nk).*im3d(:,:,nk),zeros([sz(1),sz(2)])),position,clim);
text(axh_v(end), 0.5, 0.5, '$\mathrm{E}[D_\Delta^2]$', 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')

if FLAG_relaxation
    %% mr
    clear im3d
    im3d = mr;
    clim = mr_clim;
    position.left = position.left + dleft;
    axh_v = mplot_slicescolumn(cat(3,mask(:,:,nk).*im3d(:,:,nk),zeros([sz(1),sz(2)])),position,clim);
    text(axh_v(end), 0.5, 0.5, ['$\mathrm{E}[' char_plot_relaxation ']$'], 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')
end

%% vdiso
clear im3d
im3d = vdiso;
clim = vdiso_clim;
position.left = position.left + dleft;
axh_v = mplot_slicescolumn(cat(3,mask(:,:,nk).*im3d(:,:,nk),zeros([sz(1),sz(2)])),position,clim);
text(axh_v(end), 0.5, 0.5, '$\mathrm{V}[D_\mathrm{iso}]$', 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')

%% vsqddelta
clear im3d
im3d = vsqddelta; 
clim = vsqddelta_clim;
position.left = position.left + dleft;
axh_v = mplot_slicescolumn(cat(3,mask(:,:,nk).*im3d(:,:,nk),zeros([sz(1),sz(2)])),position,clim);
text(axh_v(end), 0.5, 0.5, '$\mathrm{V}[D_\Delta^2]$', 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')

if FLAG_relaxation
    % vr
    clear im3d
    im3d = vr;
    clim = vr_clim;
    position.left = position.left + dleft;
    axh_v = mplot_slicescolumn(cat(3,mask(:,:,nk).*im3d(:,:,nk),zeros([sz(1),sz(2)])),position,clim);
    text(axh_v(end), 0.5, 0.5, ['$\mathrm{V}[' char_plot_relaxation ']$'], 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')
end

%% cvdisosqddelta
clear im3d
im3d = cvdisosqddelta; 
clim = cvdisosqddelta_clim;
position.left = position.left + dleft;
axh_v = mplot_slicescolumn(cat(3,mask(:,:,nk).*im3d(:,:,nk),zeros([sz(1),sz(2)])),position,clim);
for n = 1:numel(axh_v), colormap(axh_v(n),mplot_cmaphotcold(64)), end
text(axh_v(end), 0.5, 0.5, '$\mathrm{C}[D_\mathrm{iso},D_\Delta^2]$', 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')
    
im3d = mask.*im3d;
im3d(im3d < clim(1)) = clim(1);
im3d(im3d > clim(2)) = clim(2);
im3d_ind = ceil((im3d - clim(1))/(clim(2) - clim(1))*size(mplot_cmaphotcold(64),1));
I = zeros([3 msf_size(im3d,3)]);
for k = nk
    RGB = ind2rgb(im3d_ind(:,:,k),mplot_cmaphotcold(64));
    I(1,:,:,k) = RGB(:,:,1);
    I(2,:,:,k) = RGB(:,:,2);
    I(3,:,:,k) = RGB(:,:,3);
end
I(isnan(I)) = 0;
I(isinf(I)) = 0;
I( I(:) > 1 ) = 1;
I( I(:) < 0 ) = 0;
mdm_nii_write(255*I, fullfile(parameter_maps_directory, strcat(method, '_cvdisosqddelta_rgb.nii.gz')), nii_h, 1);

if FLAG_relaxation
    %% cvdisor
    clear im3d
    im3d = cvdisor; 
    clim = cvdisor_clim;
    position.left = position.left + dleft;
    axh_v = mplot_slicescolumn(cat(3,mask(:,:,nk).*im3d(:,:,nk),zeros([sz(1),sz(2)])),position,clim);
    for n = 1:numel(axh_v), colormap(axh_v(n),mplot_cmaphotcold(64)), end
    text(axh_v(end), 0.5, 0.5, ['$\mathrm{C}[D_\mathrm{iso},' char_plot_relaxation ']$'], 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')

    im3d = mask.*im3d;
    im3d(im3d < clim(1)) = clim(1);
    im3d(im3d > clim(2)) = clim(2);
    im3d_ind = ceil((im3d - clim(1))/(clim(2) - clim(1))*size(mplot_cmaphotcold(64),1));
    I = zeros([3 msf_size(im3d,3)]);
    for k = nk
        RGB = ind2rgb(im3d_ind(:,:,k),mplot_cmaphotcold(64));
        I(1,:,:,k) = RGB(:,:,1);
        I(2,:,:,k) = RGB(:,:,2);
        I(3,:,:,k) = RGB(:,:,3);
    end
    I(isnan(I)) = 0;
    I(isinf(I)) = 0;
    I( I(:) > 1 ) = 1;
    I( I(:) < 0 ) = 0;
    mdm_nii_write(255*I, fullfile(parameter_maps_directory, strcat(method, '_cvdisor_rgb.nii.gz')), nii_h, 1);
    
    %% cvsqddeltar
    clear im3d
    im3d = cvsqddeltar; 
    clim = cvsqddeltar_clim;
    position.left = position.left + dleft;
    axh_v = mplot_slicescolumn(cat(3,mask(:,:,nk).*im3d(:,:,nk),zeros([sz(1),sz(2)])),position,clim);
    for n = 1:numel(axh_v), colormap(axh_v(n),mplot_cmaphotcold(64)), end
    text(axh_v(end), 0.5, 0.5, ['$\mathrm{C}[D_\Delta^2,' char_plot_relaxation ']$'], 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')
    
    im3d = mask.*im3d;
    im3d(im3d < clim(1)) = clim(1);
    im3d(im3d > clim(2)) = clim(2);
    im3d_ind = ceil((im3d - clim(1))/(clim(2) - clim(1))*size(mplot_cmaphotcold(64),1));
    I = zeros([3 msf_size(im3d,3)]);
    for k = nk
        RGB = ind2rgb(im3d_ind(:,:,k),mplot_cmaphotcold(64));
        I(1,:,:,k) = RGB(:,:,1);
        I(2,:,:,k) = RGB(:,:,2);
        I(3,:,:,k) = RGB(:,:,3);
    end
    I(isnan(I)) = 0;
    I(isinf(I)) = 0;
    I( I(:) > 1 ) = 1;
    I( I(:) < 0 ) = 0;
    mdm_nii_write(255*I, fullfile(parameter_maps_directory, strcat(method, '_cvsqddeltar_rgb.nii.gz')), nii_h, 1);
end

%% diso_bins
clear im3d
cind = (cat(3,mdiso_comp3(:,:,nk),zeros([sz(1),sz(2)]))-min(diso_clim))/(max(diso_clim)-min(diso_clim));
im3d = dist_cind2rgb_jet(cind);
im3d.bright = cat(3,mask(:,:,nk).*f_comp3(:,:,nk),zeros([sz(1),sz(2)]));
clim = [0 1];
position.left = position.left + dleft;
axh_v = mplot_slicescolumn(im3d,position,clim);
text(axh_v(end), 0.5, 0.5, '$\mathrm{E}[D_\mathrm{iso}]_{\mathrm{thin}}$', 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')

cind = (mdiso_comp3-min(diso_clim))/(max(diso_clim)-min(diso_clim));
c = dist_cind2rgb_jet(cind);
c.bright = mask.*f_comp3;
I = zeros([3 msf_size(c.r,3)]);
I(1,:,:,:) = c.bright.*c.r;
I(2,:,:,:) = c.bright.*c.g;
I(3,:,:,:) = c.bright.*c.b;
I(isnan(I)) = 0;
I(isinf(I)) = 0;
I(I(:) > 1) = 1;
I(I(:) < 0) = 0;
mdm_nii_write(uint8(255*I), fullfile(parameter_maps_directory, strcat(method, '_mdiso_thin_rgb.nii.gz')), nii_h, 1);

clear im3d
cind = (cat(3,mdiso_comp2(:,:,nk),zeros([sz(1),sz(2)]))-min(diso_clim))/(max(diso_clim)-min(diso_clim));
im3d = dist_cind2rgb_jet(cind);
im3d.bright = cat(3,mask(:,:,nk).*f_comp2(:,:,nk),zeros([sz(1),sz(2)]));
clim = [0 1];
position.left = position.left + dleft;
axh_v = mplot_slicescolumn(im3d,position,clim);
text(axh_v(end), 0.5, 0.5, '$\mathrm{E}[D_\mathrm{iso}]_{\mathrm{thick}}$', 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')

cind = (mdiso_comp2-min(diso_clim))/(max(diso_clim)-min(diso_clim));
c = dist_cind2rgb_jet(cind);
c.bright = mask.*f_comp2;
I = zeros([3 msf_size(c.r,3)]);
I(1,:,:,:) = c.bright.*c.r;
I(2,:,:,:) = c.bright.*c.g;
I(3,:,:,:) = c.bright.*c.b;
I(isnan(I)) = 0;
I(isinf(I)) = 0;
I(I(:) > 1) = 1;
I(I(:) < 0) = 0;
mdm_nii_write(uint8(255*I), fullfile(parameter_maps_directory, strcat(method, '_mdiso_thick_rgb.nii.gz')), nii_h, 1);

clear im3d
cind = (cat(3,mdiso_comp1(:,:,nk),zeros([sz(1),sz(2)]))-min(diso_clim))/(max(diso_clim)-min(diso_clim));
im3d = dist_cind2rgb_jet(cind);
im3d.bright = cat(3,mask(:,:,nk).*f_comp1(:,:,nk),zeros([sz(1),sz(2)]));
clim = [0 1];
position.left = position.left + dleft;
axh_v = mplot_slicescolumn(im3d,position,clim);
text(axh_v(end), 0.5, 0.5, '$\mathrm{E}[D_\mathrm{iso}]_{\mathrm{big}}$', 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')

cind = (mdiso_comp1-min(diso_clim))/(max(diso_clim)-min(diso_clim));
c = dist_cind2rgb_jet(cind);
c.bright = mask.*f_comp1;
I = zeros([3 msf_size(c.r,3)]);
I(1,:,:,:) = c.bright.*c.r;
I(2,:,:,:) = c.bright.*c.g;
I(3,:,:,:) = c.bright.*c.b;
I(isnan(I)) = 0;
I(isinf(I)) = 0;
I(I(:) > 1) = 1;
I(I(:) < 0) = 0;
mdm_nii_write(uint8(255*I), fullfile(parameter_maps_directory, strcat(method, '_mdiso_big_rgb.nii.gz')), nii_h, 1);

%% msqddelta_bins
clear im3d
cind = (cat(3,msqddelta_comp3(:,:,nk),zeros([sz(1),sz(2)]))-min(msqddelta_clim))/(max(msqddelta_clim)-min(msqddelta_clim));
im3d = dist_cind2rgb_jet(cind);
im3d.bright = cat(3,mask(:,:,nk).*f_comp3(:,:,nk),zeros([sz(1),sz(2)]));
clim = [0 1];
position.left = position.left + dleft;
axh_v = mplot_slicescolumn(im3d,position,clim);
text(axh_v(end), 0.5, 0.5, '$\mathrm{E}[D_\Delta^2]_{\mathrm{thin}}$', 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')

cind = (msqddelta_comp3-min(msqddelta_clim))/(max(msqddelta_clim)-min(msqddelta_clim));
c = dist_cind2rgb_jet(cind);
c.bright = mask.*f_comp3;
I = zeros([3 msf_size(c.r,3)]);
I(1,:,:,:) = c.bright.*c.r;
I(2,:,:,:) = c.bright.*c.g;
I(3,:,:,:) = c.bright.*c.b;
I(isnan(I)) = 0;
I(isinf(I)) = 0;
I(I(:) > 1) = 1;
I(I(:) < 0) = 0;
mdm_nii_write(uint8(255*I), fullfile(parameter_maps_directory, strcat(method, '_msqddelta_thin_rgb.nii.gz')), nii_h, 1);

clear im3d
cind = (cat(3,msqddelta_comp2(:,:,nk),zeros([sz(1),sz(2)]))-min(msqddelta_clim))/(max(msqddelta_clim)-min(msqddelta_clim));
im3d = dist_cind2rgb_jet(cind);
im3d.bright = cat(3,mask(:,:,nk).*f_comp2(:,:,nk),zeros([sz(1),sz(2)]));
clim = [0 1];
position.left = position.left + dleft;
axh_v = mplot_slicescolumn(im3d,position,clim);
text(axh_v(end), 0.5, 0.5, '$\mathrm{E}[D_\Delta^2]_{\mathrm{thick}}$', 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')

cind = (msqddelta_comp2-min(msqddelta_clim))/(max(msqddelta_clim)-min(msqddelta_clim));
c = dist_cind2rgb_jet(cind);
c.bright = mask.*f_comp2;
I = zeros([3 msf_size(c.r,3)]);
I(1,:,:,:) = c.bright.*c.r;
I(2,:,:,:) = c.bright.*c.g;
I(3,:,:,:) = c.bright.*c.b;
I(isnan(I)) = 0;
I(isinf(I)) = 0;
I(I(:) > 1) = 1;
I(I(:) < 0) = 0;
mdm_nii_write(uint8(255*I), fullfile(parameter_maps_directory, strcat(method, '_msqddelta_thick_rgb.nii.gz')), nii_h, 1);

clear im3d
cind = (cat(3,msqddelta_comp1(:,:,nk),zeros([sz(1),sz(2)]))-min(msqddelta_clim))/(max(msqddelta_clim)-min(msqddelta_clim));
im3d = dist_cind2rgb_jet(cind);
im3d.bright = cat(3,mask(:,:,nk).*f_comp1(:,:,nk),zeros([sz(1),sz(2)]));
clim = [0 1];
position.left = position.left + dleft;
axh_v = mplot_slicescolumn(im3d,position,clim);
text(axh_v(end), 0.5, 0.5, '$\mathrm{E}[D_\Delta^2]_{\mathrm{big}}$', 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')

cind = (msqddelta_comp1-min(msqddelta_clim))/(max(msqddelta_clim)-min(msqddelta_clim));
c = dist_cind2rgb_jet(cind);
c.bright = mask.*f_comp1;
I = zeros([3 msf_size(c.r,3)]);
I(1,:,:,:) = c.bright.*c.r;
I(2,:,:,:) = c.bright.*c.g;
I(3,:,:,:) = c.bright.*c.b;
I(isnan(I)) = 0;
I(isinf(I)) = 0;
I(I(:) > 1) = 1;
I(I(:) < 0) = 0;
mdm_nii_write(uint8(255*I), fullfile(parameter_maps_directory, strcat(method, '_msqddelta_big_rgb.nii.gz')), nii_h, 1);

%% mr_bins
if FLAG_relaxation
    clear im3d
    cind = (cat(3,mr_comp3(:,:,nk),zeros([sz(1),sz(2)]))-min(mr_clim))/(max(mr_clim)-min(mr_clim));
    im3d = dist_cind2rgb_jet(cind);
    im3d.bright = cat(3,mask(:,:,nk).*f_comp3(:,:,nk),zeros([sz(1),sz(2)]));
    clim = [0 1];
    position.left = position.left + dleft;
    axh_v = mplot_slicescolumn(im3d,position,clim);
    text(axh_v(end), 0.5, 0.5, ['$\mathrm{E}[' char_plot_relaxation ']_{\mathrm{thin}}$'], 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')
    
    cind = (mr_comp3-min(mr_clim))/(max(mr_clim)-min(mr_clim));
    c = dist_cind2rgb_jet(cind);
    c.bright = mask.*f_comp3;
    I = zeros([3 msf_size(c.r,3)]);
    I(1,:,:,:) = c.bright.*c.r;
    I(2,:,:,:) = c.bright.*c.g;
    I(3,:,:,:) = c.bright.*c.b;
    I(isnan(I)) = 0;
    I(isinf(I)) = 0;
    I(I(:) > 1) = 1;
    I(I(:) < 0) = 0;
    mdm_nii_write(uint8(255*I), fullfile(parameter_maps_directory, strcat(method, '_m', char_relaxation, '_thin_rgb.nii.gz')), nii_h, 1);
    
    clear im3d
    cind = (cat(3,mr_comp2(:,:,nk),zeros([sz(1),sz(2)]))-min(mr_clim))/(max(mr_clim)-min(mr_clim));
    im3d = dist_cind2rgb_jet(cind);
    im3d.bright = cat(3,mask(:,:,nk).*f_comp2(:,:,nk),zeros([sz(1),sz(2)]));
    clim = [0 1];
    position.left = position.left + dleft;
    axh_v = mplot_slicescolumn(im3d,position,clim);
    text(axh_v(end), 0.5, 0.5, ['$\mathrm{E}[' char_plot_relaxation ']_{\mathrm{thick}}$'], 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')
    
    cind = (mr_comp2-min(mr_clim))/(max(mr_clim)-min(mr_clim));
    c = dist_cind2rgb_jet(cind);
    c.bright = mask.*f_comp2;
    I = zeros([3 msf_size(c.r,3)]);
    I(1,:,:,:) = c.bright.*c.r;
    I(2,:,:,:) = c.bright.*c.g;
    I(3,:,:,:) = c.bright.*c.b;
    I(isnan(I)) = 0;
    I(isinf(I)) = 0;
    I(I(:) > 1) = 1;
    I(I(:) < 0) = 0;
    mdm_nii_write(uint8(255*I), fullfile(parameter_maps_directory, strcat(method, '_m', char_relaxation, '_thick_rgb.nii.gz')), nii_h, 1);
    
    clear im3d
    cind = (cat(3,mr_comp1(:,:,nk),zeros([sz(1),sz(2)]))-min(mr_clim))/(max(mr_clim)-min(mr_clim));
    im3d = dist_cind2rgb_jet(cind);
    im3d.bright = cat(3,mask(:,:,nk).*f_comp1(:,:,nk),zeros([sz(1),sz(2)]));
    clim = [0 1];
    position.left = position.left + dleft;
    axh_v = mplot_slicescolumn(im3d,position,clim);
    text(axh_v(end), 0.5, 0.5, ['$\mathrm{E}[' char_plot_relaxation ']_{\mathrm{big}}$'], 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')
    
    cind = (mr_comp1-min(mr_clim))/(max(mr_clim)-min(mr_clim));
    c = dist_cind2rgb_jet(cind);
    c.bright = mask.*f_comp1;
    I = zeros([3 msf_size(c.r,3)]);
    I(1,:,:,:) = c.bright.*c.r;
    I(2,:,:,:) = c.bright.*c.g;
    I(3,:,:,:) = c.bright.*c.b;
    I(isnan(I)) = 0;
    I(isinf(I)) = 0;
    I(I(:) > 1) = 1;
    I(I(:) < 0) = 0;
    mdm_nii_write(uint8(255*I), fullfile(parameter_maps_directory, strcat(method, '_m', char_relaxation, '_big_rgb.nii.gz')), nii_h, 1);
end

%% fractions_bins
clear im3d
im3d.r = cat(3,mdxx_comp3(:,:,nk),zeros([sz(1),sz(2)]));
im3d.g = cat(3,mdyy_comp3(:,:,nk),zeros([sz(1),sz(2)]));
im3d.b = cat(3,mdzz_comp3(:,:,nk),zeros([sz(1),sz(2)]));
im3d.bright = cat(3,mask(:,:,nk).*f_comp3(:,:,nk),zeros([sz(1),sz(2)]));
clim = [0 1];
position.left = position.left + dleft;
axh_v = mplot_slicescolumn(im3d,position,clim);
text(axh_v(end), 0.5, 0.5, '$f_{\mathrm{color,thin}}$', 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')

norm_md = sqrt(mdxx_comp3.^2 + mdyy_comp3.^2 + mdzz_comp3.^2);
c.r = squeeze(abs(mdxx_comp3)./norm_md)./mask;
c.g = squeeze(abs(mdyy_comp3)./norm_md)./mask;
c.b = squeeze(abs(mdzz_comp3)./norm_md)./mask;
c.bright = mask.*f_comp3;
I = zeros([3 msf_size(c.r,3)]);
I(1,:,:,:) = c.bright.*c.r;
I(2,:,:,:) = c.bright.*c.g;
I(3,:,:,:) = c.bright.*c.b;
I(isnan(I)) = 0;
I(isinf(I)) = 0;
I(I(:) > 1) = 1;
I(I(:) < 0) = 0;
mdm_nii_write(uint8(255*I), fullfile(parameter_maps_directory, strcat(method, '_fraction_thin_rgb.nii.gz')), nii_h, 1);

clear im3d
im3d.r = cat(3,mdxx_comp2(:,:,nk),zeros([sz(1),sz(2)]));
im3d.g = cat(3,mdyy_comp2(:,:,nk),zeros([sz(1),sz(2)]));
im3d.b = cat(3,mdzz_comp2(:,:,nk),zeros([sz(1),sz(2)]));
im3d.bright = cat(3,mask(:,:,nk).*f_comp2(:,:,nk),zeros([sz(1),sz(2)]));
clim = [0 1];
position.left = position.left + dleft;
axh_v = mplot_slicescolumn(im3d,position,clim);
text(axh_v(end), 0.5, 0.5, '$f_{\mathrm{color,thick}}$', 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')

norm_md = sqrt(mdxx_comp2.^2 + mdyy_comp2.^2 + mdzz_comp2.^2);
c.r = squeeze(abs(mdxx_comp2)./norm_md)./mask;
c.g = squeeze(abs(mdyy_comp2)./norm_md)./mask;
c.b = squeeze(abs(mdzz_comp2)./norm_md)./mask;
c.bright = mask.*f_comp2;
I = zeros([3 msf_size(c.r,3)]);
I(1,:,:,:) = c.bright.*c.r;
I(2,:,:,:) = c.bright.*c.g;
I(3,:,:,:) = c.bright.*c.b;
I(isnan(I)) = 0;
I(isinf(I)) = 0;
I(I(:) > 1) = 1;
I(I(:) < 0) = 0;
mdm_nii_write(uint8(255*I), fullfile(parameter_maps_directory, strcat(method, '_fraction_thick_rgb.nii.gz')), nii_h, 1);

clear im3d
im3d.r = cat(3,mdxx_comp1(:,:,nk),zeros([sz(1),sz(2)]));
im3d.g = cat(3,mdyy_comp1(:,:,nk),zeros([sz(1),sz(2)]));
im3d.b = cat(3,mdzz_comp1(:,:,nk),zeros([sz(1),sz(2)]));
im3d.bright = cat(3,mask(:,:,nk).*f_comp1(:,:,nk),zeros([sz(1),sz(2)]));
clim = [0 1];
position.left = position.left + dleft;
axh_v = mplot_slicescolumn(im3d,position,clim);
text(axh_v(end), 0.5, 0.5, '$f_{\mathrm{color,big}}$', 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')

norm_md = sqrt(mdxx_comp1.^2 + mdyy_comp1.^2 + mdzz_comp1.^2);
c.r = squeeze(abs(mdxx_comp1)./norm_md)./mask;
c.g = squeeze(abs(mdyy_comp1)./norm_md)./mask;
c.b = squeeze(abs(mdzz_comp1)./norm_md)./mask;
c.bright = mask.*f_comp1;
I = zeros([3 msf_size(c.r,3)]);
I(1,:,:,:) = c.bright.*c.r;
I(2,:,:,:) = c.bright.*c.g;
I(3,:,:,:) = c.bright.*c.b;
I(isnan(I)) = 0;
I(isinf(I)) = 0;
I(I(:) > 1) = 1;
I(I(:) < 0) = 0;
mdm_nii_write(uint8(255*I), fullfile(parameter_maps_directory, strcat(method, '_fraction_big_rgb.nii.gz')), nii_h, 1);

%% Segmentation
clear im3d
im3d.r = cat(3,f_comp3(:,:,nk),zeros([sz(1),sz(2)]));
im3d.g = cat(3,f_comp2(:,:,nk),zeros([sz(1),sz(2)]));
im3d.b = cat(3,f_comp1(:,:,nk),zeros([sz(1),sz(2)]));
im3d.bright = cat(3,mask(:,:,nk).*(f_comp1(:,:,nk) + f_comp2(:,:,nk) + f_comp3(:,:,nk)),zeros([sz(1),sz(2)]));
clim = [0 1];
position.left = position.left + dleft;
axh_v = mplot_slicescolumn(im3d,position,clim);
text(axh_v(end), 0.5, 0.5, 'segment.', 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')

norm_md = ones(size(f_comp1));
c.r = squeeze(abs(f_comp3)./norm_md)./mask;
c.g = squeeze(abs(f_comp2)./norm_md)./mask;
c.b = squeeze(abs(f_comp1)./norm_md)./mask;
c.bright = mask.*(f_comp1 + f_comp2 + f_comp3);
I = zeros([3 msf_size(c.r,3)]);
I(1,:,:,:) = c.bright.*c.r;
I(2,:,:,:) = c.bright.*c.g;
I(3,:,:,:) = c.bright.*c.b;
I(isnan(I)) = 0;
I(isinf(I)) = 0;
I(I(:) > 1) = 1;
I(I(:) < 0) = 0;
mdm_nii_write(uint8(255*I), fullfile(parameter_maps_directory, strcat(method, '_segmentation_rgb.nii.gz')), nii_h, 1);

set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize);

msf_mkdir(parameter_maps_directory);
eval(['print ' parameter_maps_directory '/technicolor_maps -loose -dpdf'])
clf(f)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
if FLAG_relaxation
    f = figure('Visible', 'off');
    
    position.dbottom = 1/(Nslices+1);
    dleft = 1/Nparams;
    position.height = 1.01*position.dbottom;
    position.width = 1.01*dleft;
    position.left = 0;
    
    % s0
    im3d = s0;
    clim = 0.4*smax*[0 1];
    axh_v = mplot_slicescolumn(cat(3,mask(:,:,nk).*im3d(:,:,nk),zeros([sz(1),sz(2)])),position,clim);
    text(axh_v(end), 0.5, 0.5, '$S_0$', 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')
    
    % mdiso
    clear im3d
    im3d = mdiso;
    clim = diso_clim;
    position.left = position.left + dleft;
    axh_v = mplot_slicescolumn(cat(3,mask(:,:,nk).*im3d(:,:,nk),zeros([sz(1),sz(2)])),position,clim);
    text(axh_v(end), 0.5, 0.5, '$\mathrm{E}[D_\mathrm{iso}]$', 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')
    
    % msqddelta
    clear im3d
    im3d = msqddelta;
    clim = msqddelta_clim;
    position.left = position.left + dleft;
    axh_v = mplot_slicescolumn(cat(3,mask(:,:,nk).*im3d(:,:,nk), zeros([sz(1),sz(2)])),position,clim);
    text(axh_v(end), 0.5, 0.5, '$\mathrm{E}[D_\Delta^2]$', 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')
   
    % mt
    clear im3d
    im3d = mt;
    clim = mt_clim_blackandwhite;
    position.left = position.left + dleft;
    axh_v = mplot_slicescolumn(cat(3,mask(:,:,nk).*im3d(:,:,nk), zeros([sz(1),sz(2)])),position,clim);
    text(axh_v(end), 0.5, 0.5, ['$\mathrm{E}[' char_plot_relaxation2 ']$'], 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')
    
    % vdiso
    clear im3d
    im3d = vdiso;
    clim = vdiso_clim;
    position.left = position.left + dleft;
    axh_v = mplot_slicescolumn(cat(3,mask(:,:,nk).*im3d(:,:,nk),zeros([sz(1),sz(2)])),position,clim);
    text(axh_v(end), 0.5, 0.5, '$\mathrm{V}[D_\mathrm{iso}]$', 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')
    
    % vsqddelta
    clear im3d
    im3d = vsqddelta;
    clim = vsqddelta_clim;
    position.left = position.left + dleft;
    axh_v = mplot_slicescolumn(cat(3,mask(:,:,nk).*im3d(:,:,nk),zeros([sz(1),sz(2)])),position,clim);
    text(axh_v(end), 0.5, 0.5, '$\mathrm{V}[D_\Delta^2]$', 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')
    
    % vt
    clear im3d
    im3d = vt;
    clim = vt_clim;
    position.left = position.left + dleft;
    axh_v = mplot_slicescolumn(cat(3,mask(:,:,nk).*im3d(:,:,nk),zeros([sz(1),sz(2)])),position,clim);
    text(axh_v(end), 0.5, 0.5, ['$\mathrm{V}[' char_plot_relaxation2 ']$'], 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')
    
    % cvdisosqddelta
    clear im3d
    im3d = cvdisosqddelta; 
    clim = cvdisosqddelta_clim;
    position.left = position.left + dleft;
    axh_v = mplot_slicescolumn(cat(3,mask(:,:,nk).*im3d(:,:,nk),zeros([sz(1),sz(2)])),position,clim);
    for n = 1:numel(axh_v), colormap(axh_v(n),mplot_cmaphotcold(64)), end
    text(axh_v(end), 0.5, 0.5, '$\mathrm{C}[D_\mathrm{iso},D_\Delta^2]$', 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')
   
    % cvdisot
    clear im3d
    im3d = cvdisot; 
    clim = cvdisot_clim;
    position.left = position.left + dleft;
    axh_v = mplot_slicescolumn(cat(3,mask(:,:,nk).*im3d(:,:,nk),zeros([sz(1),sz(2)])),position,clim);
    for n = 1:numel(axh_v), colormap(axh_v(n),mplot_cmaphotcold(64)), end
    text(axh_v(end), 0.5, 0.5, ['$\mathrm{C}[D_\mathrm{iso},' char_plot_relaxation2 ']$'], 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')
    
    im3d = mask.*im3d;
    im3d(im3d < clim(1)) = clim(1);
    im3d(im3d > clim(2)) = clim(2);
    im3d_ind = ceil((im3d - clim(1))/(clim(2) - clim(1))*size(mplot_cmaphotcold(64),1));
    I = zeros([3 msf_size(im3d,3)]);
    for k = nk
        RGB = ind2rgb(im3d_ind(:,:,k),mplot_cmaphotcold(64));
        I(1,:,:,k) = RGB(:,:,1);
        I(2,:,:,k) = RGB(:,:,2);
        I(3,:,:,k) = RGB(:,:,3);
    end
    I(isnan(I)) = 0;
    I(isinf(I)) = 0;
    I( I(:) > 1 ) = 1;
    I( I(:) < 0 ) = 0;
    mdm_nii_write(255*I, fullfile(parameter_maps_directory, strcat(method, '_cvdisot_rgb.nii.gz')), nii_h, 1);
    
    % cvsqddeltat
    clear im3d
    im3d = cvsqddeltat; 
    clim = cvsqddeltat_clim;
    position.left = position.left + dleft;
    axh_v = mplot_slicescolumn(cat(3,mask(:,:,nk).*im3d(:,:,nk),zeros([sz(1),sz(2)])),position,clim);
    for n = 1:numel(axh_v), colormap(axh_v(n),mplot_cmaphotcold(64)), end
    text(axh_v(end), 0.5, 0.5, ['$\mathrm{C}[D_\Delta^2,' char_plot_relaxation2 ']$'], 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')
    
    im3d = mask.*im3d;
    im3d(im3d < clim(1)) = clim(1);
    im3d(im3d > clim(2)) = clim(2);
    im3d_ind = ceil((im3d - clim(1))/(clim(2) - clim(1))*size(mplot_cmaphotcold(64),1));
    I = zeros([3 msf_size(im3d,3)]);
    for k = nk
        RGB = ind2rgb(im3d_ind(:,:,k),mplot_cmaphotcold(64));
        I(1,:,:,k) = RGB(:,:,1);
        I(2,:,:,k) = RGB(:,:,2);
        I(3,:,:,k) = RGB(:,:,3);
    end
    I(isnan(I)) = 0;
    I(isinf(I)) = 0;
    I( I(:) > 1 ) = 1;
    I( I(:) < 0 ) = 0;
    mdm_nii_write(255*I, fullfile(parameter_maps_directory, strcat(method, '_cvsqddeltat_rgb.nii.gz')), nii_h, 1);
    
    clear im3d
    cind = (cat(3,mdiso_comp3(:,:,nk),zeros([sz(1),sz(2)]))-min(diso_clim))/(max(diso_clim)-min(diso_clim));
    im3d = dist_cind2rgb_jet(cind);
    im3d.bright = cat(3,mask(:,:,nk).*f_comp3(:,:,nk),zeros([sz(1),sz(2)]));
    clim = [0 1];
    position.left = position.left + dleft;
    axh_v = mplot_slicescolumn(im3d,position,clim);
    text(axh_v(end), 0.5, 0.5, '$\mathrm{E}[D_\mathrm{iso}]_{\mathrm{thin}}$', 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')
        
    clear im3d
    cind = (cat(3,mdiso_comp2(:,:,nk),zeros([sz(1),sz(2)]))-min(diso_clim))/(max(diso_clim)-min(diso_clim));
    im3d = dist_cind2rgb_jet(cind);
    im3d.bright = cat(3,mask(:,:,nk).*f_comp2(:,:,nk),zeros([sz(1),sz(2)]));
    clim = [0 1];
    position.left = position.left + dleft;
    axh_v = mplot_slicescolumn(im3d,position,clim);
    text(axh_v(end), 0.5, 0.5, '$\mathrm{E}[D_\mathrm{iso}]_{\mathrm{thick}}$', 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')
    
    clear im3d
    cind = (cat(3,mdiso_comp1(:,:,nk),zeros([sz(1),sz(2)]))-min(diso_clim))/(max(diso_clim)-min(diso_clim));
    im3d = dist_cind2rgb_jet(cind);
    im3d.bright = cat(3,mask(:,:,nk).*f_comp1(:,:,nk),zeros([sz(1),sz(2)]));
    clim = [0 1];
    position.left = position.left + dleft;
    axh_v = mplot_slicescolumn(im3d,position,clim);
    text(axh_v(end), 0.5, 0.5, '$\mathrm{E}[D_\mathrm{iso}]_{\mathrm{big}}$', 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')
    
    clear im3d
    cind = (cat(3,msqddelta_comp3(:,:,nk),zeros([sz(1),sz(2)]))-min(msqddelta_clim))/(max(msqddelta_clim)-min(msqddelta_clim));
    im3d = dist_cind2rgb_jet(cind);
    im3d.bright = cat(3,mask(:,:,nk).*f_comp3(:,:,nk),zeros([sz(1),sz(2)]));
    clim = [0 1];
    position.left = position.left + dleft;
    axh_v = mplot_slicescolumn(im3d,position,clim);
    text(axh_v(end), 0.5, 0.5, '$\mathrm{E}[D_\Delta^2]_{\mathrm{thin}}$', 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')
    
    clear im3d
    cind = (cat(3,msqddelta_comp2(:,:,nk),zeros([sz(1),sz(2)]))-min(msqddelta_clim))/(max(msqddelta_clim)-min(msqddelta_clim));
    im3d = dist_cind2rgb_jet(cind);
    im3d.bright = cat(3,mask(:,:,nk).*f_comp2(:,:,nk),zeros([sz(1),sz(2)]));
    clim = [0 1];
    position.left = position.left + dleft;
    axh_v = mplot_slicescolumn(im3d,position,clim);
    text(axh_v(end), 0.5, 0.5, '$\mathrm{E}[D_\Delta^2]_{\mathrm{thick}}$', 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')
    
    clear im3d
    cind = (cat(3,msqddelta_comp1(:,:,nk),zeros([sz(1),sz(2)]))-min(msqddelta_clim))/(max(msqddelta_clim)-min(msqddelta_clim));
    im3d = dist_cind2rgb_jet(cind);
    im3d.bright = cat(3,mask(:,:,nk).*f_comp1(:,:,nk),zeros([sz(1),sz(2)]));
    clim = [0 1];
    position.left = position.left + dleft;
    axh_v = mplot_slicescolumn(im3d,position,clim);
    text(axh_v(end), 0.5, 0.5, '$\mathrm{E}[D_\Delta^2]_{\mathrm{big}}$', 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')
    
    clear im3d
    cind = (cat(3,mt_comp3(:,:,nk),zeros([sz(1),sz(2)]))-min(mt_clim))/(max(mt_clim)-min(mt_clim));
    im3d = dist_cind2rgb_jet(cind);
    im3d.bright = cat(3,mask(:,:,nk).*f_comp3(:,:,nk),zeros([sz(1),sz(2)]));
    clim = [0 1];
    position.left = position.left + dleft;
    axh_v = mplot_slicescolumn(im3d,position,clim);
    text(axh_v(end), 0.5, 0.5, ['$\mathrm{E}[' char_plot_relaxation2 ']_{\mathrm{thin}}$'], 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')
    
    cind = (mt_comp3-min(mt_clim))/(max(mt_clim)-min(mt_clim));
    c = dist_cind2rgb_jet(cind);
    c.bright = mask.*f_comp3;
    I = zeros([3 msf_size(c.r,3)]);
    I(1,:,:,:) = c.bright.*c.r;
    I(2,:,:,:) = c.bright.*c.g;
    I(3,:,:,:) = c.bright.*c.b;
    I(isnan(I)) = 0;
    I(isinf(I)) = 0;
    I(I(:) > 1) = 1;
    I(I(:) < 0) = 0;
    mdm_nii_write(uint8(255*I), fullfile(parameter_maps_directory, strcat(method, '_m', char_relaxation2, '_thin_rgb.nii.gz')), nii_h, 1);
    
    clear im3d
    cind = (cat(3,mt_comp2(:,:,nk),zeros([sz(1),sz(2)]))-min(mt_clim))/(max(mt_clim)-min(mt_clim));
    im3d = dist_cind2rgb_jet(cind);
    im3d.bright = cat(3,mask(:,:,nk).*f_comp2(:,:,nk),zeros([sz(1),sz(2)]));
    clim = [0 1];
    position.left = position.left + dleft;
    axh_v = mplot_slicescolumn(im3d,position,clim);
    text(axh_v(end), 0.5, 0.5, ['$\mathrm{E}[' char_plot_relaxation2 ']_{\mathrm{thick}}$'], 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')
    
    cind = (mt_comp2-min(mt_clim))/(max(mt_clim)-min(mt_clim));
    c = dist_cind2rgb_jet(cind);
    c.bright = mask.*f_comp2;
    I = zeros([3 msf_size(c.r,3)]);
    I(1,:,:,:) = c.bright.*c.r;
    I(2,:,:,:) = c.bright.*c.g;
    I(3,:,:,:) = c.bright.*c.b;
    I(isnan(I)) = 0;
    I(isinf(I)) = 0;
    I(I(:) > 1) = 1;
    I(I(:) < 0) = 0;
    mdm_nii_write(uint8(255*I), fullfile(parameter_maps_directory, strcat(method, '_m', char_relaxation2, '_thick_rgb.nii.gz')), nii_h, 1);
    
    clear im3d
    cind = (cat(3,mt_comp1(:,:,nk),zeros([sz(1),sz(2)]))-min(mt_clim))/(max(mt_clim)-min(mt_clim));
    im3d = dist_cind2rgb_jet(cind);
    im3d.bright = cat(3,mask(:,:,nk).*f_comp1(:,:,nk),zeros([sz(1),sz(2)]));
    clim = [0 1];
    position.left = position.left + dleft;
    axh_v = mplot_slicescolumn(im3d,position,clim);
    text(axh_v(end), 0.5, 0.5, ['$\mathrm{E}[' char_plot_relaxation2 ']_{\mathrm{big}}$'], 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')
    
    cind = (mt_comp1-min(mt_clim))/(max(mt_clim)-min(mt_clim));
    c = dist_cind2rgb_jet(cind);
    c.bright = mask.*f_comp1;
    I = zeros([3 msf_size(c.r,3)]);
    I(1,:,:,:) = c.bright.*c.r;
    I(2,:,:,:) = c.bright.*c.g;
    I(3,:,:,:) = c.bright.*c.b;
    I(isnan(I)) = 0;
    I(isinf(I)) = 0;
    I(I(:) > 1) = 1;
    I(I(:) < 0) = 0;
    mdm_nii_write(uint8(255*I), fullfile(parameter_maps_directory, strcat(method, '_m', char_relaxation2, '_big_rgb.nii.gz')), nii_h, 1);
    
    clear im3d
    im3d.r = cat(3,mdxx_comp3(:,:,nk),zeros([sz(1),sz(2)]));
    im3d.g = cat(3,mdyy_comp3(:,:,nk),zeros([sz(1),sz(2)]));
    im3d.b = cat(3,mdzz_comp3(:,:,nk),zeros([sz(1),sz(2)]));
    im3d.bright = cat(3,mask(:,:,nk).*f_comp3(:,:,nk),zeros([sz(1),sz(2)]));
    clim = [0 1];
    position.left = position.left + dleft;
    axh_v = mplot_slicescolumn(im3d,position,clim);
    text(axh_v(end), 0.5, 0.5, '$f_{\mathrm{color,thin}}$', 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')
    
    clear im3d
    im3d.r = cat(3,mdxx_comp2(:,:,nk),zeros([sz(1),sz(2)]));
    im3d.g = cat(3,mdyy_comp2(:,:,nk),zeros([sz(1),sz(2)]));
    im3d.b = cat(3,mdzz_comp2(:,:,nk),zeros([sz(1),sz(2)]));
    im3d.bright = cat(3,mask(:,:,nk).*f_comp2(:,:,nk),zeros([sz(1),sz(2)]));
    clim = [0 1];
    position.left = position.left + dleft;
    axh_v = mplot_slicescolumn(im3d,position,clim);
    text(axh_v(end), 0.5, 0.5, '$f_{\mathrm{color,thick}}$', 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')
    
    clear im3d
    im3d.r = cat(3,mdxx_comp1(:,:,nk),zeros([sz(1),sz(2)]));
    im3d.g = cat(3,mdyy_comp1(:,:,nk),zeros([sz(1),sz(2)]));
    im3d.b = cat(3,mdzz_comp1(:,:,nk),zeros([sz(1),sz(2)]));
    im3d.bright = cat(3,mask(:,:,nk).*f_comp1(:,:,nk),zeros([sz(1),sz(2)]));
    clim = [0 1];
    position.left = position.left + dleft;
    axh_v = mplot_slicescolumn(im3d,position,clim);
    text(axh_v(end), 0.5, 0.5, '$f_{\mathrm{color,big}}$', 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')    
    
    clear im3d
    im3d.r = cat(3,f_comp3(:,:,nk),zeros([sz(1),sz(2)]));
    im3d.g = cat(3,f_comp2(:,:,nk),zeros([sz(1),sz(2)]));
    im3d.b = cat(3,f_comp1(:,:,nk),zeros([sz(1),sz(2)]));
    im3d.bright = cat(3,mask(:,:,nk).*(f_comp1(:,:,nk) + f_comp2(:,:,nk) + f_comp3(:,:,nk)),zeros([sz(1),sz(2)]));
    clim = [0 1];
    position.left = position.left + dleft;
    axh_v = mplot_slicescolumn(im3d,position,clim);
    text(axh_v(end), 0.5, 0.5, 'segment.', 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')
    
    set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize);
    
    msf_mkdir(parameter_maps_directory);
    eval(['print ' parameter_maps_directory '/technicolor_maps2 -loose -dpdf'])
    clf(f)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mdxx = sum(dxx.*w,4)./s0;
% mdyy = sum(dyy.*w,4)./s0;
% mdzz = sum(dzz.*w,4)./s0;
% mdxy = sum(dxy.*w,4)./s0;
% mdxz = sum(dxz.*w,4)./s0;
% mdyz = sum(dyz.*w,4)./s0;
% norm_md = sqrt(mdxx.^2 + mdyy.^2 + mdzz.^2);

% FA and color FA
% [~, E_shear, ~] = tm_1x21_iso();
% reg = 0.0001;
% sz_md = size(mdxx);
% FA = zeros(sz_md);
% for vx = 1:sz_md(1)
%     for vy = 1:sz_md(2)
%         for vz = 1:sz_md(3)
%             if mask(vx,vy,vz)
%                 dt_1x6 = [mdxx(vx,vy,vz),mdyy(vx,vy,vz),mdzz(vx,vy,vz),sqrt(2)*mdxy(vx,vy,vz),sqrt(2)*mdxz(vx,vy,vz),sqrt(2)*mdyz(vx,vy,vz)]*1e9;
%                 dt2_1x21 = tm_1x6_to_1x21(dt_1x6);
%                 V_shear2 = tm_inner(dt2_1x21, E_shear);
%                 MD = tm_md(dt_1x6);
%                 FA_value = tm_fa(max(MD, reg), V_shear2);
%                 FA(vx,vy,vz) = mio_min_max_cut(FA_value, -1, 2);
%                 %FA(vx,vy,vz) = tm_fa([mdxx(vx,vy,vz),mdyy(vx,vy,vz),mdzz(vx,vy,vz),sqrt(2)*mdxy(vx,vy,vz),sqrt(2)*mdxz(vx,vy,vz),sqrt(2)*mdyz(vx,vy,vz)]);
%             end
%         end
%     end
% end
% 
% c.bright = FA;
% c.r = squeeze(abs(mdxx)./norm_md)./mask;
% c.g = squeeze(abs(mdyy)./norm_md)./mask;
% c.b = squeeze(abs(mdzz)./norm_md)./mask;
% I = zeros([3 msf_size(c.r,3)]);
% I(1,:,:,:) = c.bright.*c.r;
% I(2,:,:,:) = c.bright.*c.g;
% I(3,:,:,:) = c.bright.*c.b;
% I(isnan(I)) = 0;
% I(isinf(I)) = 0;
% I(I(:) > 1) = 1;
% I(I(:) < 0) = 0;
% mdm_nii_write(uint8(255*I), fullfile(parameter_maps_directory, strcat(method, '_FA_rgb.nii.gz')), nii_h, 1);

%Normalized measures
% vdison = vdiso./mdiso.^2; vdison(isfinite(vdison)~=1) = 0;
% vsdanisonn = vsdanison./msdanison.^2; vsdanisonn(isfinite(vsdanisonn)~=1) = 0;
% cvdisosdanisonn = cvdisosdanison./(mdiso.*msdanison); cvdisosdanisonn(isfinite(cvdisosdanisonn)~=1) = 0;
% cvsdanisonr2n = cvsdanisonr2./(msdanison.*mr2); cvsdanisonr2n(isfinite(cvsdanisonr2n)~=1) = 0;
% cvdisor2n = cvdisor2./(mdiso.*mr2); cvdisor2n(isfinite(cvdisor2)~=1) = 0;

% color FA
% clear im3d
% im3d.r = cat(3,mdxx(:,:,nk),zeros([sz(1),sz(2)]));
% im3d.g = cat(3,mdyy(:,:,nk),zeros([sz(1),sz(2)]));
% im3d.b = cat(3,mdzz(:,:,nk),zeros([sz(1),sz(2)]));
% im3d.bright = cat(3,mask(:,:,nk).*FA(:,:,nk),zeros([sz(1),sz(2)]));
% clim = FA_clim;
% position.left = position.left + dleft;
% axh_v = mplot_slicescolumn(im3d,position,clim);
% text(axh_v(end), 0.5, 0.5, 'FA$_{\mathrm{color}}$', 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')
% clear im3d

%     % FA
%     im3d = FA;
%     clim = FA_clim;
%     position.left = position.left + dleft;
%     axh_v = mplot_slicescolumn(cat(3,mask(:,:,nk).*im3d(:,:,nk),zeros([sz(1),sz(2)])),position,clim);
%     text(axh_v(end), 0.5, 0.5, 'FA', 'Units', 'Normalized', 'Fontsize', fontsize_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.9290 0.6940 0.1250], 'Interpreter', 'latex')

% cind = (cvdisosdanison-min(clim))/(max(clim)-min(clim));
% cmap = mplot_cmaphotcold(64);
% N = size(cmap,1);
% cind(cind > 1) = 1;
% cind(cind < 0) = 0;
% c.bright = mask.*cind;
% cind = ceil(N*cind);
% cind(cind == 0) = 1;
% c.r = zeros(size(sz));
% c.g = zeros(size(sz));
% c.b = zeros(size(sz));
% for i = 1:sz(1)
%     for j = 1:sz(2)
%         for k = 1:sz(3)
%             color = cmap(cind(i,j,k),:);
%             c.r(i,j,k) = color(1);
%             c.g(i,j,k) = color(2);
%             c.b(i,j,k) = color(3);
%         end
%     end
% end
% I = zeros([3 msf_size(c.r,3)]);
% I(1,:,:,:) = c.bright.*c.r;
% I(2,:,:,:) = c.bright.*c.g;
% I(3,:,:,:) = c.bright.*c.b;
% I(isnan(I)) = 0;
% I(isinf(I)) = 0;
% I(I(:) > 1) = 1;
% I(I(:) < 0) = 0;
% mdm_nii_write(uint8(255*I), fullfile(parameter_maps_directory, strcat(method, '_cvdisosdanison_rgb.nii.gz')), nii_h, 1);

% clear im3d
% cind = (mdiso_comp4(:,:,nk)-min(diso_clim))/(max(diso_clim)-min(diso_clim));
% im3d = dist_cind2rgb_jet(cind);
% im3d.bright = mask(:,:,nk).*f_comp4(:,:,nk);
% clim = [0 1];
% position.left = position.left + dleft;
% axh_v = mplot_slicescolumn(im3d,position,clim);
