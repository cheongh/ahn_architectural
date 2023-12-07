library(malevnc)
library(neuprintr)
library(nat.nblast)
library(nat)

msahn_bodyids = c(13926, 12536)
manc_t2_bodyids = neuprint_list2df(neuprint_fetch_custom(cypher="MATCH (a:Neuron) WHERE a.somaNeuromere = 'T2' AND (a.upstream > 25 OR a.downstream > 25) AND a.status IN ['Traced','Sensory Anchor','Soma Anchor','Anchor','Primary Anchor','Cervical Anchor','PRT Orphan'] RETURN a.bodyId AS bodyid", timeout=2000))

msahn_swcs = neuprint_read_skeletons(msahn_bodyids)
manc_t2_swcs = neuprint_read_neurons(manc_t2_bodyids, OmitFailures = TRUE)

#rescale to um for nblast
for (i in 1:length(msahn_swcs)) msahn_swcs[[i]]$d[c("X", "Y", "Z")] = msahn_swcs[[i]]$d[c("X", "Y", "Z")]*8/1000
for (i in 1:length(manc_t2_swcs)) manc_t2_swcs[[i]]$d[c("X", "Y", "Z")] = manc_t2_swcs[[i]]$d[c("X", "Y", "Z")]*8/1000

manc_t2_swcs_original = manc_t2_swcs
msahn_swcs = prune_twigs(msahn_swcs, twig_length = 10)
manc_t2_swcs = prune_twigs(manc_t2_swcs, twig_length = 10)

#need to rescale to microns BEFORE doing symmetric_manc
msahn_swcs_mirror = symmetric_manc(msahn_swcs, mirror = TRUE)
msahn_swcs = symmetric_manc(msahn_swcs, mirror = FALSE)
manc_t2_swcs = symmetric_manc(manc_t2_swcs, mirror = FALSE)


msahn_dotprops_mirror = dotprops(msahn_swcs_mirror, resample=1, k=5)
msahn_dotprops = dotprops(msahn_swcs, resample=1, k=5)
manc_t2_dotprops = dotprops(manc_t2_swcs, resample=1, k=5)

nblast_df = nlapply(msahn_dotprops, FUN = function(x) {temp = nblast(x, target = manc_t2_dotprops, version = 2, normalised = TRUE);
	sort(temp, decreasing = TRUE)[2:51]})
#not actually df, create the df properly below
nblast_output_df = data.frame(segID = names(nblast_df))
nblast_output_df[paste0("match_", 1:50)] = t(sapply(nblast_df, FUN = function(x) names(x)))
nblast_output_df[paste0("norm_score_", 1:50)] = t(sapply(nblast_df, FUN = function(x) x))

google_output_ss = "https://docs.google.com/spreadsheets/d/1jbapg_QSySZggY69PTnh43K9A5L3hhzIu1anPa10Oac/edit#gid=0"
sheet_write(nblast_output_df, google_output_ss, "MsAHN_nblast_MANC")

nblast_df = nlapply(msahn_dotprops_mirror, FUN = function(x) {temp = nblast(x, target = manc_t2_dotprops, version = 2, normalised = TRUE);
	sort(temp, decreasing = TRUE)[2:51]})
#not actually df, create the df properly below
nblast_output_df = data.frame(segID = names(nblast_df))
nblast_output_df[paste0("match_", 1:50)] = t(sapply(nblast_df, FUN = function(x) names(x)))
nblast_output_df[paste0("norm_score_", 1:50)] = t(sapply(nblast_df, FUN = function(x) x))

google_output_ss = "https://docs.google.com/spreadsheets/d/1jbapg_QSySZggY69PTnh43K9A5L3hhzIu1anPa10Oac/edit#gid=0"
sheet_write(nblast_output_df, google_output_ss, "MsAHN_mirrored_nblast_MANC")





mtahn_bodyids = c(42819, 11003)
manc_t2_bodyids = neuprint_list2df(neuprint_fetch_custom(cypher="MATCH (a:Neuron) WHERE a.somaNeuromere = 'T3' AND (a.upstream > 25 OR a.downstream > 25) AND a.status IN ['Traced','Sensory Anchor','Soma Anchor','Anchor','Primary Anchor','Cervical Anchor','PRT Orphan'] RETURN a.bodyId AS bodyid", timeout=2000))

mtahn_swcs = neuprint_read_skeletons(mtahn_bodyids)
manc_t3_swcs = neuprint_read_neurons(manc_t2_bodyids, OmitFailures = TRUE)

#rescale to um for nblast
for (i in 1:length(mtahn_swcs)) mtahn_swcs[[i]]$d[c("X", "Y", "Z")] = mtahn_swcs[[i]]$d[c("X", "Y", "Z")]*8/1000
for (i in 1:length(manc_t3_swcs)) manc_t3_swcs[[i]]$d[c("X", "Y", "Z")] = manc_t3_swcs[[i]]$d[c("X", "Y", "Z")]*8/1000

manc_t3_swcs_original = manc_t3_swcs
mtahn_swcs = prune_twigs(mtahn_swcs, twig_length = 10)
manc_t3_swcs = prune_twigs(manc_t3_swcs, twig_length = 10)

#need to rescale to microns BEFORE doing symmetric_manc
mtahn_swcs_mirror = symmetric_manc(mtahn_swcs, mirror = TRUE)
mtahn_swcs = symmetric_manc(mtahn_swcs, mirror = FALSE)
manc_t3_swcs = symmetric_manc(manc_t3_swcs, mirror = FALSE)


mtahn_dotprops_mirror = dotprops(mtahn_swcs_mirror, resample=1, k=5)
mtahn_dotprops = dotprops(mtahn_swcs, resample=1, k=5)
manc_t2_dotprops = dotprops(manc_t3_swcs, resample=1, k=5)

nblast_df = nlapply(mtahn_dotprops, FUN = function(x) {temp = nblast(x, target = manc_t2_dotprops, version = 2, normalised = TRUE);
	sort(temp, decreasing = TRUE)[2:51]})
#not actually df, create the df properly below
nblast_output_df = data.frame(segID = names(nblast_df))
nblast_output_df[paste0("match_", 1:50)] = t(sapply(nblast_df, FUN = function(x) names(x)))
nblast_output_df[paste0("norm_score_", 1:50)] = t(sapply(nblast_df, FUN = function(x) x))

google_output_ss = "https://docs.google.com/spreadsheets/d/1jbapg_QSySZggY69PTnh43K9A5L3hhzIu1anPa10Oac/edit#gid=0"
sheet_write(nblast_output_df, google_output_ss, "MtAHN_nblast_MANC")

nblast_df = nlapply(mtahn_dotprops_mirror, FUN = function(x) {temp = nblast(x, target = manc_t2_dotprops, version = 2, normalised = TRUE);
	sort(temp, decreasing = TRUE)[2:51]})
#not actually df, create the df properly below
nblast_output_df = data.frame(segID = names(nblast_df))
nblast_output_df[paste0("match_", 1:50)] = t(sapply(nblast_df, FUN = function(x) names(x)))
nblast_output_df[paste0("norm_score_", 1:50)] = t(sapply(nblast_df, FUN = function(x) x))

google_output_ss = "https://docs.google.com/spreadsheets/d/1jbapg_QSySZggY69PTnh43K9A5L3hhzIu1anPa10Oac/edit#gid=0"
sheet_write(nblast_output_df, google_output_ss, "MtAHN_mirrored_nblast_MANC")