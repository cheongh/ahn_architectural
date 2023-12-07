library(fafbseg)
library(gplots)
library(viridis)
library(dplyr)
library(tidyr)
library(dendsort)
library(ggplot2)
library(RColorBrewer)
library(svglite)
library(igraph)

mtahns =  c("MtAHN-L"="720575940614269393", "MtAHN-R"="720575940622346876")
msahns = c("MsAHN-R"="720575940630175276", "MsAHN-L"="720575940626130469")
#Tect INs
type_1 = c(23510,22289,22472,101206,21307,23628,20774,18519,157981,20590,18048,14502,16970,19315,18370,14625,32328,15337,23485,18657,18517,20537,15586,36805,22521,21381,23462,24945,153911,23191,23638,13060,22220,24442,155196,20586,16638,15788,17055,15113,101682,17945)

#color lookup table for vertices
graph_lut = c(MsAHN = "#ffe099", MtAHN = "#ffff64", SN = "#db58d7", ASN = "#ffc9b4", AN = "#9632fa", IN = "#ff6666", DN = "#5ce7e1", MN = "#1910c7", EN = "#80c8ff", AEN = "#b9b4ff", ND = "#aaaaaa", VCN = "#fb9fbc", VPN = "#6f97f2")
superclass_abbrv = c("central" = "IN", "visual_projection" = "VPN",  "descending" = "DN", "optic"="OLN", "motor" = "MN", "visual_centrifugal" = "VCN", "sensory" = "SN", "ascending" = "AN")

#fast method to do sum matrix rows and column by group
combine_matrix_rows_cols_by_group =function(adj, groups=NA, row_groups=NA, col_groups=NA, sparse = FALSE) {
	if (all(is.na(groups)) & (all(is.na(row_groups))|all(is.na(col_groups)))) stop("Either groups or both row_groups and col_groups need to be specified")
	if (all(is.na(row_groups))|all(is.na(col_groups))) i = j = as.character(groups) 
	else {
		i = as.character(row_groups)
		j = as.character(col_groups)
	}
	i[is.na(i)] = 'NA'
	j[is.na(j)] = 'NA'
	i = Matrix::sparse.model.matrix(~i+0,transpose=T)
	rownames(i) = substr(rownames(i),2,1000)
	j = Matrix::sparse.model.matrix(~j+0)
	colnames(j) = substr(colnames(j),2,1000)
	adj = i %*% adj %*% j
	if (sparse) adj else as.matrix(adj)
}

squarify = function(adj) {
	row_col_names = unique(c(rownames(adj), colnames(adj)))
	new_adj = matrix(0, length(row_col_names), length(row_col_names), dimnames=list(row_col_names,row_col_names))
	new_adj[rownames(adj), colnames(adj)] = adj
	return(new_adj)
}

#prepost is pre or postsynaptic connections here, not sites
get_flywire_syn_neuropil_per_partner = function(root_ids, prepost) {
	all_syns = flywire_connectome_data("syn")
	if (prepost == "POST") {
		collect(filter(all_syns, pre_pt_root_id %in% root_ids))
	} else if (prepost == "PRE") {
		collect(filter(all_syns, post_pt_root_id %in% root_ids))
	} else stop("prepost must be PRE or POST")
}

get_flywire_syn_neuropil = function(root_ids, prepost) {
	syn_per_partner = get_flywire_syn_neuropil_per_partner(root_ids, prepost)
	if (prepost=="PRE") {
		syn_per_partner = aggregate(syn_count~post_pt_root_id+neuropil, data=syn_per_partner, FUN = sum)
		syn_per_partner = tidyr::pivot_wider(syn_per_partner, id_cols = post_pt_root_id, names_from = neuropil, values_from = syn_count)
		row_names = as.character(syn_per_partner$post_pt_root_id)
		syn_per_partner = as.matrix(syn_per_partner[,colnames(syn_per_partner)!="post_pt_root_id",])
		rownames(syn_per_partner) = row_names
		return(syn_per_partner)
	} else {
		syn_per_partner = aggregate(syn_count~pre_pt_root_id+neuropil, data=syn_per_partner, FUN = sum)
		syn_per_partner = tidyr::pivot_wider(syn_per_partner, id_cols = pre_pt_root_id, names_from = neuropil, values_from = syn_count)
		row_names = as.character(syn_per_partner$pre_pt_root_id)
		syn_per_partner = as.matrix(syn_per_partner[,colnames(syn_per_partner)!="pre_pt_root_id",])
		rownames(syn_per_partner) = row_names
		return(syn_per_partner)
	}
}

simple_row_hierarchical_clustering <- function(mat, Rowv = NULL) {
  d=dist(mat,method="euclidean")
  fit = hclust(d, method="ward.D2")
  if (is.null(Rowv)) {
    fit = dendsort(fit)
    return(mat[fit$order,])
  } else {
    fit = reorder(as.dendrogram(fit),Rowv)
    return(mat[labels(fit),])
  }
}

#flywire neuropils
np_order = c("ME","LO",
			"MB_CA",
			"LAL","GA",#LX
			"AOTU","AVLP","PVLP","PLP","WED",#VLNP
			"LH",
			"SLP", "SIP",#SNP
			"CRE","SCL","ICL","IB", #INP
			"VES","EPA","GOR","SPS","IPS",#VMNP
			"AMMC","SAD","CAN","FLA",#PENP
			"GNG")
super_np = list(OL=c("ME","LO"),
			LX=c("LAL","GA"),#LX
			VLNP=c("AOTU","AVLP","PVLP","PLP","WED"),#VLNP
			SNP=c("SLP", "SIP"),#SNP
			INP=c("CRE","SCL","ICL","IB")) #INP)
all_np_order = c("OL",
			"MB_CA",
			"LX",
			"VLNP",
			"LH",
			"SNP",
			"INP",
			"VES","EPA","GOR","SPS","IPS",#VMNP
			"AMMC","SAD","CAN","FLA",#PENP
			"GNG")
		