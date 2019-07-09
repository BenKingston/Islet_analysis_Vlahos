


function [post_nuclei, post_nuclei_dilate, post_vessels] = post_process_islet_func(pre_nuclei,seg_nuclei,seg_vessels,save_dir,sample_name)

shortfile = sample_name;
display(['Post-processing ' shortfile])

    dapi = pre_nuclei;

    fill_tissue = fillholes(dapi>40);
    open_tissue = opening(fill_tissue,4,'elliptic');
    erosion_tissue = erosion(open_tissue,3,'elliptic');
    gauss_tissue = gaussf(erosion_tissue,1);
    tissue_area = uint8(gauss_tissue>0);
    
    num_slices = size(tissue_area,3);

%%%%Load in Ilastik segmented image files
nuclei_seg = seg_nuclei;
vessels_seg = seg_vessels;

%%%%%Post-processing of Ilastik segmented nuclei channel
nuclei_seg_bin = nuclei_seg==3;
tissue_area_bin = tissue_area>0;
nuclei_seg_crop = (nuclei_seg_bin.*tissue_area_bin)>0;


threshnuc = nuclei_seg_crop>0;
    threshnuc3 = opening(threshnuc,2);
    dt_threshnuc = dt(threshnuc3);
     
    seeds = maxima(dt_threshnuc,2,0);
    seeds2 = dilation(seeds,2)>0;
    image_out = waterseed(seeds2,max(dt_threshnuc)-dt_threshnuc,1,0,0);
    threshnuc4 = (uint8(threshnuc3))>0;
    threshnuc4(image_out) = false; 
 
nuclabel = label(threshnuc4>0,1,7,10000);
nuc_dilated_label = dip_growregions(nuclabel,[],[],1,2,'low_first');
 final_nuc = nuclabel;
final_dilated_nuc = nuc_dilated_label;

nuclei_processed = uint32(final_nuc);
nuclei_dilated_processed = uint32(final_dilated_nuc);
%%%%%Post-processing of Ilastik segmented blood vessel channel
vessel_seg_bin = vessels_seg==3;
vessel_seg_crop = vessel_seg_bin.*tissue_area_bin;

vessels_processed = vessel_seg_crop>0;

%%%%%Final segmentation images
post_nuclei = nuclei_processed;
post_vessels = vessels_processed;
post_nuclei_dilate = nuclei_dilated_processed;
%%%%%Write post-processed files

nuclei_processed_name = strcat(shortfile,'_post_processed_nuclei.tif');
nuclei_dilated_processed_name = strcat(shortfile,'_post_processed_dialted_nuclei.tif');
vessels_processed_name = strcat(shortfile,'_post_processed_vessels.tif');

num_slices = size(vessels_seg,3);

cd(save_dir)
imwrite(uint16(nuclei_processed(:,:,1)),nuclei_processed_name);
         
  for p = 2:num_slices
          imwrite(uint16(nuclei_processed(:,:,p)),nuclei_processed_name, 'WriteMode','append');
  end

imwrite(uint16(nuclei_dilated_processed(:,:,1)),nuclei_dilated_processed_name);
         
  for p = 2:num_slices
          imwrite(uint16(nuclei_dilated_processed(:,:,p)),nuclei_dilated_processed_name, 'WriteMode','append');
  end  

 imwrite(uint8(vessels_processed(:,:,1)),vessels_processed_name);
         
    for p = 2:num_slices
            imwrite(uint8(vessels_processed(:,:,p)),vessels_processed_name, 'WriteMode','append');
    end   

end


