#----------------------- Continuous PFDA
##---------------------- Uadjusted
PFDA_inflam_lm <- data.frame(OlinkID = NA_character_, Value = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:353){
  
  s_lm <- (lm(as.formula(paste0(protein_in_inflammation[i], "~ PFDA_Aug21")), 
              data = BioMe_proteome_PFAS_wide))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = coef(s_lm)
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  PFDA_inflam_lm <- rbind(PFDA_inflam_lm, c(protein_in_inflammation[i], as.numeric(r.est[2, c(1, 2, 3, 4)])))
}



PFDA_inflam_lm <- PFDA_inflam_lm[-1,]
PFDA_inflam_lm$z.value <- as.numeric(PFDA_inflam_lm$z.value)
PFDA_inflam_lm$p.value <- as.numeric(PFDA_inflam_lm$p.value)


q <-qvalue::qvalue(as.numeric(PFDA_inflam_lm$p.value), lambda=0)
PFDA_inflam_lm$q.value <-  q$qvalues


PFDA_inflam_unlm_results<- PFDA_inflam_lm %>% 
  left_join(BioMe_proteome_PFAS_inflammation[,c("OlinkID", "Protein_name")], by="OlinkID") %>% 
  distinct()


write.csv(PFDA_inflam_unlm_results,
          "~/Projects/BioMe/proteome/input/exwas/exwas_PFDA_inflam_unlm_cont.csv",
          row.names = F)



##---------------------- Adjusted
PFDA_inflam_lm <- data.frame(OlinkID = NA_character_, Value = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:353){
  
  s_lm <- (glm(as.formula(paste0("status~", protein_in_inflammation[i], "+ self_reported_race + age_at_enrollment + sex + bmi_at_enrollment+ c_date_enrl")), 
              data = BioMe_proteome_PFAS_wide, family = binomial()))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = coef(s_lm)
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  PFDA_inflam_lm <- rbind(PFDA_inflam_lm, c(protein_in_inflammation[i], as.numeric(r.est[2, c(1, 2, 3, 4)])))
}



PFDA_inflam_lm <- PFDA_inflam_lm[-1,]
PFDA_inflam_lm$z.value <- as.numeric(PFDA_inflam_lm$z.value)
PFDA_inflam_lm$p.value <- as.numeric(PFDA_inflam_lm$p.value)


q <-qvalue::qvalue(as.numeric(PFDA_inflam_lm$p.value), lambda=0)
PFDA_inflam_lm$q.value <-  q$qvalues


PFDA_inflam_adlm_results<- PFDA_inflam_lm %>% 
  left_join(BioMe_proteome_PFAS_inflammation[,c("OlinkID", "Protein_name")], by="OlinkID") %>% 
  distinct()


write.csv(PFDA_inflam_adlm_results,
          "~/Projects/BioMe/proteome/input/exwas/exwas_T2D_inflam_adlm_cont.csv",
          row.names = F)

