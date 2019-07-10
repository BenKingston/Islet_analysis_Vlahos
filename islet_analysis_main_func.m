



function [vessel_radius_img, cell_dist_vessel_img, cell_ins_int_img, ind_cell_data, tissue_data] = islet_analysis_main_func(pre_nuclei,post_nuclei,post_vessels,insulin_ch,save_dir,sample_name)

shortfile = sample_name;
display(['Analyzing ' shortfile])
tic
%%%Define the tissue boundary using the nuclei channel
    dapi = pre_nuclei;
    
    fill_tissue = fillholes(dapi>40);
    open_tissue = opening(fill_tissue,4,'elliptic');
    erosion_tissue = erosion(open_tissue,3,'elliptic');
    gauss_tissue = gaussf(erosion_tissue,1);
    tissue_area = uint8(gauss_tissue>0);
    
    num_slices = size(tissue_area,3);

%%%%Load in Ilastik segmented image files and insulin channel
nuclei_post_processed = post_nuclei;
vessels_post_processed = post_vessels;

%%%%Load insulin channel and trim to tissue boundary
ins = insulin_ch;
tissue_area_bin = tissue_area>0;
ins_crop = ins.*uint16(tissue_area_bin);

%%%%Analysis of individual cells
nuclabel = label(nuclei_post_processed>0,1,25,0);
nuc_dilated_label = dip_growregions(nuclabel,[],[],1,2,'low_first');
nuc_dilated_label = uint32(nuc_dilated_label);
All_tissue_stats = regionprops3(tissue_area_bin,ins_crop,'Centroid','MeanIntensity','Volume','SurfaceArea');

all_cells_int_stats = regionprops3(nuc_dilated_label,ins_crop,'MeanIntensity','VoxelIdxList');
all_cells_int_stats_struct = table2struct(all_cells_int_stats);

vess_thresh = vessels_post_processed>0;
dt_vessels = bwdist(vess_thresh);
all_cells_dist_stats = regionprops3(nuc_dilated_label,dt_vessels,'MeanIntensity','VoxelIdxList');
all_cells_dist_stats_struct = table2struct(all_cells_dist_stats);

max_cell_in_tissue = max(nuc_dilated_label(:));
min_cell_in_tissue = min(nuc_dilated_label(:))+1;

new_all_nuc_int =  uint32(nuc_dilated_label);     
    
    for b = min_cell_in_tissue:max_cell_in_tissue
            new_all_nuc_int(all_cells_int_stats_struct(b).VoxelIdxList) = all_cells_int_stats_struct(b).MeanIntensity;
    end

new_all_nuc_dist =  uint32(nuc_dilated_label);     

    for b = min_cell_in_tissue:max_cell_in_tissue
            new_all_nuc_dist(all_cells_dist_stats_struct(b).VoxelIdxList) = all_cells_dist_stats_struct(b).MeanIntensity;
    end

new_all_nuc_int_crop = new_all_nuc_int.*uint32(tissue_area_bin);    
new_all_nuc_dist_crop = new_all_nuc_dist.*uint32(tissue_area_bin);      
nuclei_dilated_post_processed_crop = nuc_dilated_label.*uint32(tissue_area_bin);      
all_cells_int_stats.Properties.VariableNames = {'Int_VoxelIDs' 'MeanInsIntensity'};
all_cells_dist_stats.Properties.VariableNames = {'Dist_VoxelIDs' 'MeanDisttovessel_px'};
all_cell_stats = [all_cells_int_stats.MeanInsIntensity all_cells_dist_stats.MeanDisttovessel_px];

%%%%Vessel Diameter analysis
vess_neg = vess_thresh~=1;
vess_neg_dia = bwdist(vess_neg);
vess_dia = vess_neg_dia + vess_thresh;

vessel_vol_total = sum(vessels_post_processed(:));
sum_vess_dist = sum(vess_dia(:));

avg_vess_dia = sum_vess_dist/vessel_vol_total;

%%%%Vessel volume analysis
tissuearea_vol = sum(tissue_area(:));
vessel_vol = sum(vessels_post_processed(:));
vessel_percent = (vessel_vol/tissuearea_vol)*100;

table = array2table([tissuearea_vol vessel_vol vessel_percent avg_vess_dia]);
table.Properties.VariableNames = {'Totaltissuevolpx' 'Totalvesselvolpx' 'PercentVolisVessel' 'Avg_vess_diameter_px'};

cd(save_dir)
table_name = strcat(shortfile,'_total vol_tissue_vess.csv');
writetable(table,table_name);
 
%%%%%Final analyzed images and tables
vessel_radius_img = vess_dia;
cell_dist_vessel_img = new_all_nuc_dist_crop;
cell_ins_int_img = new_all_nuc_int_crop;
ind_cell_data = all_cell_stats;
tissue_data = table;

%%%%%Write analyzed image files and tables to whole image data and
%%%%%individual cell data
cd(save_dir)
table_name = strcat(shortfile,'_total vol_tissue_vess.csv');
writetable(table,table_name);

all_cell_name = strcat(shortfile,'all_cells_stats.csv');
writetable(array2table(all_cell_stats),all_cell_name);

met_mean_np_int_name = strcat(shortfile,'_mean Ins Intensity','.tif');
met_mean_dist_name = strcat(shortfile,'_mean distance to vessel','.tif');
met_nuclei_label_name = strcat(shortfile,'_labeled cells','.tif');
vessel_distimg_name = strcat(shortfile,'_radius_vessel','.tif');

num_slices = size(ins_crop,3);

imwrite(uint16(new_all_nuc_int_crop(:,:,1)),met_mean_np_int_name);
         
  for p = 2:num_slices
          imwrite(uint16(new_all_nuc_int_crop(:,:,p)),met_mean_np_int_name, 'WriteMode','append');
  end

imwrite(uint16(new_all_nuc_dist_crop(:,:,1)),met_mean_dist_name);
         
  for p = 2:num_slices
          imwrite(uint16(new_all_nuc_dist_crop(:,:,p)),met_mean_dist_name, 'WriteMode','append');
  end

imwrite(uint16(nuclei_dilated_post_processed_crop(:,:,1)),met_nuclei_label_name);
         
  for p = 2:num_slices
          imwrite(uint16(nuclei_dilated_post_processed_crop(:,:,p)),met_nuclei_label_name, 'WriteMode','append');
  end
  
  imwrite(uint16(vess_dia(:,:,1)),vessel_distimg_name);
         
  for p = 2:num_slices
          imwrite(uint16(vess_dia(:,:,p)),vessel_distimg_name, 'WriteMode','append');
  end
 
toc
end

