# SISAL_d13C_extract.R
#
# This script extracts the SISAL database and systematically filters
# the d18O to sites with at least 4,000 years of coverage out of the last 12,000 years
#
# The logic to apply the filter:
#   1. Query the MySQL database for all samples with ages <= 12000 BP(1950)
#      i.   The 'ages' are taken original_chronology.interp_age
#      ii.  The entity cannot be 'superseded'. Only the 'current' or 
#           'current, partially modified' entities are to be included
#      iii. gaps and hiatuses are to be included as these will be needed to
#           determine the 4,000 years coverage filter
#
#   2. Loop through each site from the queried table from (1.)
#      i.   if the site only has one entity:
#           a. if the entity has no hiatuses or gaps
#               - calculate the age range. if they are over
#                 4,000 years, keep the site. 
#               - Save the age range in a summary table as 'coverage'
#           b. if the entity has hiatuses or gaps
#               - align the gaps/hiatuses with the samples (based on
#                 either the depth_sample or sample_id)
#               - group the samples into groups separated by the 
#                 gaps or hiatuses and calculate the age range for
#                 each group. Sum up the age range and if they are
#                 over 4,000 years, keep the site.
#               - Save the sum of the age range in a summary table as 
#                 'coverage'
#     ii.   if the site has more than one entity:
#           a. Group samples based on entity/group using the hiatuses
#              and gaps
#               - Loop through each entity/group
#                   * find min and max ages for each entity/group
#               - Calculate the age range based on the min and max ages
#                 of each entity, taking into account the intersection of
#                 max and min ages of different entities. 
#                   * Loop through the entity/group with the youngest min 
#                     age and look for the entity/group with the next 
#                     youngest min age and whether the min age is greater 
#                     than the max age of the first one. If so, then compare 
#                     the max age of the two entity/group and take the older 
#                     max age from the two. Loop through all entity/group 
#                     until there are no more entities which intersects with 
#                     the first entity/group. Calculate the age range and 
#                     store for this entities/groups combination. Exclude 
#                     all the entity/group that intersects with the first 
#                     one from the next search. Do the same for the next 
#                     entity/group with the youngest min age that is still 
#                     left.
#               - Sum up all the age range for all the entity/group. If this 
#                 is over 4,000, keep the site.
#               - Save the sum of the age range in a summary table as 
#                 'coverage'
#



get_site_coverage <- function(threshold){
  
  # 1. Query the SISAL database ------------------------------------------####
  query <- paste("SELECT * FROM 
(SELECT site.site_id, site.site_name, site.latitude, site.longitude, site.elevation, 
entity.entity_name, entity.entity_status, entity.corresponding_current, 
sample.*, gap, hiatus, d18O_measurement, d18O_precision, 
interp_age 
FROM site JOIN entity USING(site_id) 
JOIN sample USING(entity_id) 
LEFT JOIN hiatus USING(sample_id) 
LEFT JOIN gap USING(sample_id) 
LEFT JOIN original_chronology USING(sample_id) 
LEFT JOIN d18O USING(sample_id) 
WHERE entity_id IN (SELECT distinct(entity_id) FROM sample JOIN d18O USING(sample_id))) t 
WHERE (interp_age <= 12000 OR gap IS NOT NULL OR hiatus IS NOT NULL) 
AND (entity_status != 'superseded');", sep = " ")
  
  dt <- dbGetQuery(mydb, query)
  
  # 2. Perform filter on table queried from (1.)--------------------------####
  # Initiate an empty table to store the data from sites selected
  
  dt_site <- data.frame(site_id = integer(), coverage = double())
  dt_site <- c()
  
  # Loop through each site.
  for (i in unique(dt[['site_id']])){
    dt_subset <- dt[dt[['site_id']] == i,]
    num_ent <- length(unique(dt_subset[['entity_id']]))
    # If they are all gaps/hiatuses, these are first removed from the analysis
    all_gap_hiatus <- all((dt_subset[['gap']] == 'G') | (dt_subset[['hiatus']] == 'H'))
    if (is.na(all_gap_hiatus)){
      all_gap_hiatus <- F
      }
    if (all_gap_hiatus){
      next
      }
    if (num_ent <= 0){
      print(paste('This is impossible, site_id = ', i, ' does not have an entity', sep = ''))
      } else if (num_ent == 1){
        # 2.i. site only has one entity
        if (any((dt_subset[['gap']] == 'G') | (dt_subset[['hiatus']] == 'H'), na.rm = T)){
          # 2.i.b.
          # sort dt_subset based on sample_id
          dt_subset <- dt_subset[order(dt_subset$sample_id),]
          # sort dt_subset based on depths
          dt_subset <- dt_subset[order(dt_subset$depth_sample),]
          row.names(dt_subset) <- NULL # reset the index
          # group the samples 
          grp_ctr <- 1
          for (k in 1:dim(dt_subset)[1]){
            if (is.na(dt_subset[k,'gap']) & is.na(dt_subset[k,'hiatus'])){
              dt_subset$grp[k] = toString(grp_ctr)
              } else {
                dt_subset$grp[k] = NA
                grp_ctr = grp_ctr + 1
              }
            }
          dt_sum <- subset(dt_subset, !is.na(grp)) %>% group_by(grp) %>% summarise(coverage = diff(range(interp_age, na.rm = T)), min_age = min(interp_age, na.rm = T), max_age = max(interp_age, na.rm = T))
          coverage <- sum(dt_sum[['coverage']][is.finite(dt_sum[['coverage']])], na.rm = T)
          dt_site <- rbind(dt_site, c(i, coverage))
          
          } else {
            # 2.i.a. There are no gaps or hiatuses
            coverage <- diff(range(dt_subset$interp_age, na.rm = T))
            dt_site <- rbind(dt_site, c(i, coverage))
            }
        
        } else if (num_ent > 1){
          # 2.ii.
          # 2.ii.a.
          grp_ctr <- 0
          ent_ctr <- 0
          for (j in unique(dt_subset[['entity_id']])){
            grp_ctr <- grp_ctr + 1
            ent_ctr <- ent_ctr + 1
            dt_subset_ent <- dt_subset[dt_subset[['entity_id']] == j,]
            # sort dt_subset_ent based on sample_id and depths
            dt_subset_ent <- dt_subset_ent[order(dt_subset_ent$sample_id, dt_subset_ent$depth_sample),]
            row.names(dt_subset_ent) <- NULL # reset the index
            for (k in 1:dim(dt_subset_ent)[1]){
              if (is.na(dt_subset_ent[k,'gap']) & is.na(dt_subset_ent[k,'hiatus'])){
                dt_subset_ent$grp[k] = toString(grp_ctr)
                } else {
                  grp_ctr = grp_ctr + 1
                  dt_subset_ent$grp[k] = NA
                }
              }
            if (ent_ctr == 1){
              dt_ent_grp_str <- dt_subset_ent
              } else {
                dt_ent_grp_str <- rbind(dt_ent_grp_str, dt_subset_ent)
              }
            }
          
          # group the samples into groups
          dt_sum <- dt_ent_grp_str %>% group_by(grp) %>% summarise(min_age = min(interp_age, na.rm = T), max_age = max(interp_age, na.rm = T))
          
          if (is.na(dt_sum[nrow(dt_sum), 'grp']) == T){ 
            dt_sum <- dt_sum[-nrow(dt_sum),]
            } # remove any na/inf row that may result grom hiatuses/gaps that do not occur within any Hol entities
          
          # If still multiple groups, group overlapping entities together, non-overlapping seperately
          if (nrow(dt_sum) > 1){
            dt_sum <- dt_sum[order(dt_sum$min_age),]
            row.names(dt_sum) <- NULL # reset the index
            str_minmax <- c()
            min_age <- dt_sum[1, 'min_age']
            max_age <- dt_sum[1, 'max_age']
            for (k in 1:dim(dt_sum)[1]){
              if (k == dim(dt_sum)[1]){ # if this is the end, store it in str_minmax
                if (dt_sum[k, 'min_age'] <= max_age){
                  if (dt_sum[k, 'max_age'] > max_age){ # if entity is overlapping
                    max_age <- dt_sum[k, 'max_age']
                    str_minmax <- rbind(str_minmax, c(min_age, max_age))
                    } else if (dt_sum[k, 'max_age'] <= max_age){ # if grp is entirely encompassed by k=1 entity
                      str_minmax <- rbind(str_minmax, c(min_age, max_age))
                      }
                  } else if (dt_sum[k, 'min_age'] > max_age){ # if final grp is non-overlapping with previous group
                    str_minmax <- rbind(str_minmax, c(min_age, max_age)) #save last grp/cluster of grps
                    min_age <- dt_sum[k, 'min_age']
                    max_age <- dt_sum[k, 'max_age']
                    str_minmax <- rbind(str_minmax, c(min_age, max_age)) # final entity as new group
                    }
                } else if (dt_sum[k, 'min_age'] <= max_age){ #if grp is overlapping with previous
                  if (dt_sum[k, 'max_age'] > max_age){
                    max_age <- dt_sum[k, 'max_age']
                    }
                  } else {
                    str_minmax <- rbind(str_minmax, c(min_age, max_age)) # if grp is non-overlapping, save previous, set new min and max for this grp
                    min_age <- dt_sum[k, 'min_age']
                    max_age <- dt_sum[k, 'max_age']
                  }
              }
            coverage <- sum(as.numeric(str_minmax[,2]) - as.numeric(str_minmax[,1]))
            dt_site <- rbind(dt_site, c(i, coverage))
            } else {
              coverage <- dt_sum$max_age- dt_sum$min_age
              dt_site <- rbind(dt_site, c(i, coverage))
            }
        }
  }
  dt_site <- data.frame(dt_site)
  colnames(dt_site) <- c('site_id', 'site_coverage')
  site_sum_dt <- dbGetQuery(mydb, 'SELECT site_name, site_id, latitude, longitude, elevation FROM site;')
  dt_site_merge <- merge(dt_site, site_sum_dt)
  dt_site_merge <- dt_site_merge %>% filter(site_coverage >= threshold)
}





