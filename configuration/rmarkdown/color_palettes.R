timepoint_palette <-  c('0200' = '#CBD2CA', '0600' = '#B8DECC', '1000' = '#7d80fa', '1400' = '#f46d43', '1800' = '#B284BE', '2200' = '#797979')
zt_palette <-  c('20' = '#CBD2CA', '0' = '#B8DECC', '4' = '#7d80fa', '8' = '#f46d43', '12' = '#B284BE', '16' = '#797979')
sample_type_palette <-  c('ctc_single' = 'firebrick', 'ctc_cluster' = 'blue', 'ctc_cluster_wbc'  = 'black')
sample_type_palette_legend <-  c('Single CTCs' = 'firebrick', 'CTC-Clusters' = 'blue', 'CTC-WBC Clusters'  = 'black')

timepoint_sample_type_palette <- c(
  'active_ctc_single' = '#5e5fcc',
  'active_ctc_cluster' = '#009ffa',
  'active_ctc_cluster_wbc' = '#00e5f8',
  'resting_ctc_single' = '#f60239',
  'resting_ctc_cluster' = '#fd6763',
  'resting_ctc_cluster_wbc' = '#ffc4c4'

)

timepoint_palette <- c(timepoint_palette,
                       'active' = timepoint_sample_type_palette[['active_ctc_single']],
                       'resting' = timepoint_sample_type_palette[['resting_ctc_single']])


timepoint_sample_type_legend_palette <- timepoint_sample_type_palette
names(timepoint_sample_type_legend_palette) <- c(
  'Active phase Single CTCs',
  'Active phase CTC-cluster',
  'Active phase WBC-CTC Clusters',
  'Rest phase Single CTCs',
  'Rest phase CTC-cluster',
  'Rest phase WBC-CTC Clusters'
)

timepoint_sample_type_legend_palette_2 <- timepoint_sample_type_palette
names(timepoint_sample_type_legend_palette_2) <- c(
  'Active phase Single CTCs',
  'Active phase CTC Clusters',
  'Active phase CTC-WBC Clusters',
  'Rest phase Single CTCs',
  'Rest phase CTC Clusters',
  'Rest phase CTC-WBC Clusters'
)



zt_sample_type_legend_palette <- timepoint_sample_type_palette
names(zt_sample_type_legend_palette) <- c(
  'ZT16 Single CTCs',
  'ZT16 CTC-Clusters',
  'ZT16 CTC-WBC Clusters',
  'ZT4 Single CTCs',
  'ZT4 CTC-Clusters',
  'ZT4 CTC-WBC Clusters'
)
