library(fancr)
library(malevnc)
library(googlesheets4)
library(nat)
library(neuprintr)
library(nat.nblast)

google_output_ss = "https://docs.google.com/spreadsheets/d/1jbapg_QSySZggY69PTnh43K9A5L3hhzIu1anPa10Oac/edit#gid=0"

swcs = read.neurons(paths = "fanc_skels/dn", pattern = "\\d+\\.swc", format = "swc")
xformed_swcs = transform_fanc2manc(swcs, inverse = F)

google_dn_type_ss = "https://docs.google.com/spreadsheets/d/1diedvCNkyBPIPDxPTBOldKVBzLAL6aY7hph7C6YlUc8/edit#gid=0"
google_dn_type_ss = read_sheet(google_dn_type_ss, sheet = 1)
google_dn_type_ss = google_dn_type_ss[google_dn_type_ss$class1 == "DN",]
dn_swcs = neuprint_read_skeletons(google_dn_type_ss$bodyid)
#rescale dn_swcs to be right size from voxel units (8nm) to nm
for (i in 1:length(dn_swcs)) {
	dn_swcs[[i]]$d$X = dn_swcs[[i]]$d$X * 8
	dn_swcs[[i]]$d$Y = dn_swcs[[i]]$d$Y * 8
	dn_swcs[[i]]$d$Z = dn_swcs[[i]]$d$Z * 8
}

pruned_xformed_swcs = prune_twigs(xformed_swcs, twig_length = 10000)
pruned_dn_swcs = prune_twigs(dn_swcs, twig_length = 10000)

#rescale to um for nblast
for (i in 1:length(pruned_xformed_swcs)) pruned_xformed_swcs[[i]]$d[c("X", "Y", "Z")] = pruned_xformed_swcs[[i]]$d[c("X", "Y", "Z")]/1000
for (i in 1:length(pruned_dn_swcs)) pruned_dn_swcs[[i]]$d[c("X", "Y", "Z")] = pruned_dn_swcs[[i]]$d[c("X", "Y", "Z")]/1000

fanc_dotprops = dotprops(pruned_xformed_swcs, resample=1, k=5)
dn_dotprops = dotprops(pruned_dn_swcs, resample=1, k=5)

nblast_df = nlapply(fanc_dotprops, FUN = function(x) {temp = nblast(x, target = append(dn_dotprops, neuronlist(x)), version = 2, normalised = TRUE);
	sort(temp, decreasing = TRUE)[2:11]})
#not actually df, create the df properly below
nblast_output_df = data.frame(segID = names(nblast_df))
nblast_output_df[paste0("match_", 1:10)] = t(sapply(nblast_df, FUN = function(x) names(x)))
nblast_output_df[paste0("norm_score_", 1:10)] = t(sapply(nblast_df, FUN = function(x) x))
nblast_output_df[paste0("cell_type_", 1:10)] = apply(nblast_output_df[paste0("match_", 1:10)], MARGIN = 2, FUN = function(x) sapply(x, FUN = function(y)
	ifelse(!is.na(google_dn_type_ss$ourtype[match(y, google_dn_type_ss$bodyid)]),
	google_dn_type_ss$ourtype[match(y, google_dn_type_ss$bodyid)],
	google_dn_type_ss$ourmatch[match(y, google_dn_type_ss$bodyid)])))
	
sheet_write(nblast_output_df, google_output_ss)