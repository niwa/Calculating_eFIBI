calculate.efibi<-function(Data,models,key,seed=1){
  ###########################################################################
  # Check for duplicated taxa within any sample --------------------------
  ###########################################################################
  Data %>%
    dplyr::count(site.id, date, species_code) %>%
    dplyr::filter(n > 1)->duplicates
  
  if (nrow(duplicates) > 0) {
    stop(
      "Duplicated taxa detected within some samples. Please ensure species codes are unique within each sample.",
      call. = FALSE
    )
  }
  
  ###########################################################################
  # Calculate observed richness within each sub-metric-----------------------
  ###########################################################################
  message("Calculating observed sub-metric richness")
  Data%>%
    group_by(site.id,date,genus)%>%
    mutate(
      has_species = any(rank == "species")
    )%>%
    ungroup()%>%
    #Fish are counted towards richness if identified to species level, or, when identified only to genus level, if they are the only detection of that genus in the sample.
    mutate(counted = case_when(
      rank == "genus" & !has_species ~ 1,
      rank == "species" ~ 1,
      TRUE ~ 0
    ))%>%
    #assign a sub-metric to each fish taxon
    left_join(key, by= "species_code")%>%#
    #taxa not present in the key do not belong to any sub-metric 
    mutate(across(
      c(native, benthic_riffle, benthic_pool,
        pelagic_pool, intolerant, non_native),
      ~tidyr::replace_na(.x, 0)
    ))%>%
    #calculate the number of taxa within each sub-metric for each sample
    group_by(site.id,date,elevation,distance)%>%
    summarise(across(
      all_of(c("native","benthic_riffle","benthic_pool","pelagic_pool","intolerant","non_native")),
      ~ sum(.x[counted == 1], na.rm = TRUE)
    ),.groups = "drop_last")%>%
    ungroup()%>%
    mutate(p.native=case_when(native+non_native==0~0,
                              TRUE~native/(native+non_native)))->Data.richness
  
  ############################################################################
  # Calculate maximum species richness at each site in Data ----------------
  ############################################################################
  #maximum species richness for each submetric are estimated as the 95th quantile of the predictive posterior distribution of the associated model
  message("Estimating maximum richness")
  set.seed(seed)
  Data.richness$m.n=predict(models$Native,newdata=Data.richness,probs=0.95,robust=TRUE,
                            ndraws=1200)[,"Q95"]
  Data.richness$m.br=predict(models$Benthic_riffle,newdata=Data.richness,probs=0.95,robust=TRUE,
                             ndraws=1200)[,"Q95"]
  Data.richness$m.bp=predict(models$Benthic_pool,newdata=Data.richness,probs=0.95,robust=TRUE,
                             ndraws=1200)[,"Q95"]
  Data.richness$m.pp=predict(models$Pelagic_pool,newdata=Data.richness,probs=0.95,robust=TRUE,
                             ndraws=1200)[,"Q95"]
  Data.richness$m.i=predict(models$Intolerant,newdata=Data.richness,probs=0.95,robust=TRUE,
                            ndraws=1200)[,"Q95"]
  ###########################################################################
  # Calculate eFIBI ---------------------------------------------------------
  ###########################################################################
  message("Calculating eFIBI")
  Data.richness%>%
    #calculate sub-metric scores
    mutate(native_score=case_when(native>(m.n*2/3)~5*2,
                                  native>=(m.n*1/3)~3*2,
                                  native<(m.n*1/3)~1*2),
           benthic_riffle_score=case_when(benthic_riffle>(m.br*2/3)~5*2,
                                          benthic_riffle>=(m.br*1/3)~3*2,
                                          benthic_riffle<(m.br*1/3)~1*2),
           benthic_pool_score=case_when(benthic_pool>(m.bp*2/3)~5*2,
                                        benthic_pool>=(m.bp*1/3)~3*2,
                                        benthic_pool<(m.bp*1/3)~1*2),
           pelagic_pool_score=case_when(pelagic_pool>(m.pp*2/3)~5*2,
                                        pelagic_pool>=(m.pp*1/3)~3*2,
                                        pelagic_pool<(m.pp*1/3)~1*2),
           intolerant_score=case_when(intolerant>(m.i*2/3)~5*2,
                                      intolerant>=(m.i*1/3)~3*2,
                                      intolerant<(m.i*1/3)~1*2),
           prop_native_score=case_when(p.native>2/3~5*2,
                                       p.native>=1/3~3*2,
                                       p.native<1/3~1*2),
           #calculate eFIBI as the sum of submetric scores
           eFIBI=native_score+
             benthic_riffle_score+
             benthic_pool_score+
             pelagic_pool_score+
             intolerant_score+
             prop_native_score)%>%
    dplyr::select(-m.n,-m.br,-m.bp,-m.pp,-m.i,-elevation,-distance)->result
  
  return(result)
}