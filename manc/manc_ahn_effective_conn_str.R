library(neuprintr)
library(malevnc)
library(RColorBrewer)
library(colorspace)
library(nat)
library(googlesheets4)
library(gplots)
library(viridis)
library(ggplot2)
box::use(mancfuns/han)
source("shared_functions.R")

msahn = c(12536, 13926) #12536--L, 13926--R
names(msahn) = c("MsAHN-R", "MsAHN-L")
mtahn = c(11003, 42819) #11003--L, 42819--R
names(mtahn) = c("MtAHN-R", "MtAHN-L")

np_thres = 0.1

subclass_color_lut = c(nm = "#019948", wm="#ED2690", hm="#4D8AC9", lg="#D9A528", ad="#768798", xm="#6A8B3B", notum="#EDD4B6",ventral_prothoracic="#7A4FA0")

#effective conn str to just sns and mns
mn_sn_dat = neuprintr::neuprint_list2df(neuprintr::neuprint_fetch_custom(cypher=paste0("MATCH (a:Neuron) WHERE a.class IN ['", paste0(names(han$manc_neuron_type_lookup)[(han$manc_neuron_type_lookup %in% c('MN','SN'))], collapse="','"), "'] RETURN a.bodyId AS bodyid, a.class AS class, a.group AS group, a.subclass AS subclass, a.synonyms AS synonyms, a.rootSide AS root_side, a.somaSide AS soma_side, a.exitNerve AS exit_nerve, a.entry_nerve AS entry_nerve"), timeout=2000))
mn_sn_dat$abbrv = han$class2abbrv(mn_sn_dat$class)
mn_sn_dat$subclass[mn_sn_dat$subclass %in% c(NA, "NA", "")] = "unknown"
ahn_conn_str = han$get_partners_by_multistep_conn(inputids = c(msahn, mtahn), outputids = mn_sn_dat$bodyid[mn_sn_dat$abbrv=="MN"], path_length=2, syn_frac_thres = 0, by_group = TRUE, remove_dn_us = TRUE, remove_sn_us = FALSE, return_type = "detailed")

ahn_sn_conn_str = ahn_conn_str[[2]][,colnames(ahn_conn_str[[2]]) %in% mn_sn_dat$subclass[mn_sn_dat$abbrv=="SN"]]
sn_order = c("neck hair plate","ventral prothorax pCO","wing","wing base d.Rad.A CS","wing CS","wing margin bristle","haltere","haltere CS","haltere CS dF2","haltere small diameter","notum","prothoracic leg","prothoracic leg FeCO claw","prothoracic leg FeCO club","prothoracic leg FeCO hook","prothoracic leg hair plate","prothoracic leg TrCS","mesothoracic leg","mesothoracic leg FeCO claw","mesothoracic leg FeCO club","mesothoracic leg FeCO hook","mesothoracic leg TrCS","metathoracic leg","metathoracic leg FeCO claw","metathoracic leg FeCO club","metathoracic leg FeCO hook","metathoracic leg TrCS","abdomen","ad")
sn_color_cats = list(nm="neck hair plate",ventral_prothoracic="ventral prothorax pCO",wm=c("wing","wing base d.Rad.A CS","wing CS","wing margin bristle"),hm=c("haltere","haltere CS","haltere CS dF2","haltere small diameter"),
	notum="notum",lg=c("prothoracic leg","prothoracic leg FeCO claw","prothoracic leg FeCO club","prothoracic leg FeCO hook","prothoracic leg hair plate","prothoracic leg TrCS","mesothoracic leg","mesothoracic leg FeCO claw","mesothoracic leg FeCO club","mesothoracic leg FeCO hook","mesothoracic leg TrCS","metathoracic leg","metathoracic leg FeCO claw","metathoracic leg FeCO club","metathoracic leg FeCO hook","metathoracic leg TrCS"),ad=c("abdomen","ad"))
ahn_sn_conn_str = ahn_sn_conn_str[,sn_order]
sn_color = sapply(colnames(ahn_sn_conn_str), FUN = function(x) names(sn_color_cats)[which(sapply(sn_color_cats, FUN=function(y) x %in% y))])
sn_color = subclass_color_lut[sn_color]
pdf(file=paste0("ahn_sn_indir_conn_str", Sys.Date(),".pdf"), height=(300)/72, width=(600)/72)
print(heatmap.2(as.matrix(ahn_sn_conn_str), col = viridis(100),scale="none",density.info='none',xlab="sensory category",ylab="AHN", #labRow = NA,
	ColSideColors = sn_color, dendrogram='none',Rowv=FALSE,Colv=FALSE,trace='none',sepcolor='black',cexRow = 0.15 + 1/log10(ncol(ahn_sn_conn_str)), margins = c(7, 7)))
dev.off()
plot_dat=data.frame(msahn=ahn_sn_conn_str["12536",], mtahn=ahn_sn_conn_str["11003",], subclass=factor(names(subclass_color_lut)[match(sn_color, subclass_color_lut)], levels=names(subclass_color_lut)))
pdf(file=paste0("ahn_sn_indir_conn_str_scatter_", Sys.Date(),".pdf"), height=(300)/72, width=(400)/72)
ggplot(plot_dat, aes(x=mtahn, y=msahn, color=subclass)) +
	geom_point() +
	theme_classic() +
	scale_color_manual(values=subclass_color_lut) +
	xlab("MtAHN")+ylab("MsAHN")
dev.off()

ahn_mn_conn_str = ahn_conn_str[[2]][,colnames(ahn_conn_str[[2]]) %in% mn_sn_dat$group[mn_sn_dat$abbrv=="MN"]]
temp_mn_subclass = mn_sn_dat$subclass[match(colnames(ahn_mn_conn_str), mn_sn_dat$group)]
temp_mn_subclass[temp_mn_subclass %in% c('fl','ml','hl')] = 'lg'
ahn_mn_conn_str = ahn_mn_conn_str[,order(match(temp_mn_subclass, names(subclass_color_lut)))]
mn_color = mn_sn_dat$subclass[match(colnames(ahn_mn_conn_str), mn_sn_dat$group)]
mn_color[mn_color %in% c('fl','ml','hl')] = 'lg'
mn_color = subclass_color_lut[mn_color]
pdf(file=paste0("ahn_mn_indir_conn_str", Sys.Date(),".pdf"), height=(300)/72, width=(600)/72)
print(heatmap.2(as.matrix(ahn_mn_conn_str), col = viridis(100),scale="none",density.info='none',xlab="MN subclass",ylab="AHN", #labRow = NA,
	ColSideColors = mn_color, dendrogram='none',Rowv=FALSE,Colv=FALSE,trace='none',sepcolor='black',cexRow = 0.15 + 1/log10(ncol(ahn_mn_conn_str)), margins = c(7, 7)))
	legend(x=0.4, y=1, legend=names(subclass_color_lut), fill=subclass_color_lut, col=NA, bg=NA, horiz = F)
dev.off()
plot_dat=data.frame(msahn=ahn_mn_conn_str["12536",], mtahn=ahn_mn_conn_str["11003",], subclass=factor(names(subclass_color_lut)[match(mn_color, subclass_color_lut)], levels=names(subclass_color_lut)))
pdf(file=paste0("ahn_mn_indir_conn_str_scatter_", Sys.Date(),".pdf"), height=(300)/72, width=(400)/72)
ggplot(plot_dat, aes(x=mtahn, y=msahn, color=subclass)) +
	geom_point() +
	theme_classic() +
	scale_color_manual(values=subclass_color_lut) +
	xlab("MtAHN")+ylab("MsAHN")
dev.off()

#to all neurons
clio_manc_data = neuprintr::neuprint_list2df(neuprintr::neuprint_fetch_custom(cypher=paste0("MATCH (a:Neuron) WHERE a.class IN ['", paste0(names(han$manc_neuron_type_lookup)[!(han$manc_neuron_type_lookup %in% c('ND','Glia'))], collapse="','"), "'] RETURN a.bodyId AS bodyid, a.type AS type, a.class AS class, a.group AS group, a.subclass AS subclass, a.rootSide AS root_side, a.somaSide AS soma_side, a.exitNerve AS exit_nerve, a.entry_nerve AS entry_nerve"), timeout=2000))
clio_manc_data$abbrv = han$class2abbrv(clio_manc_data$class)
clio_manc_data$group = ifelse(is.na(clio_manc_data$group), clio_manc_data$bodyid, clio_manc_data$group)
clio_manc_data$type[clio_manc_data$group %in% type_1] = "Tect IN"
ahn_all_conn_str = han$get_partners_by_multistep_conn(inputids = c(msahn, mtahn), path_length=2, syn_frac_thres = 0.0001, by_group = TRUE, remove_dn_us = TRUE, remove_sn_us = TRUE, return_type = "detailed")
rownames(ahn_all_conn_str[[2]]) = c("MsAHN","MtAHN")[match(c(12536,11003), rownames(ahn_all_conn_str[[2]]))]
ahn_all_conn_str[[2]] = ahn_all_conn_str[[2]][, !(colnames(ahn_all_conn_str[[2]]) %in% c(12536,11003))]
top_50_ms = ahn_all_conn_str[[2]]
top_50_ms = top_50_ms[,order(top_50_ms["MsAHN",],decreasing=T)[1:50]]
top_msahn_neuron_class = clio_manc_data$abbrv[match(colnames(top_50_ms), clio_manc_data$group)]
#top_msahn_neuron_class = ifelse(is.na(top_msahn_neuron_class), clio_manc_data$abbrv[match(colnames(top_50_ms), clio_manc_data$bodyid)], top_msahn_neuron_class)
top_50_ms = top_50_ms[,order(top_msahn_neuron_class)]
top_msahn_neuron_class = top_msahn_neuron_class[order(top_msahn_neuron_class)]
temp_cell_count = sapply(colnames(top_50_ms), FUN = function(x) sum(clio_manc_data$group %in% x))
row_groups = sapply(colnames(top_50_ms), FUN=function(x) {
	type = clio_manc_data$type[match(x, clio_manc_data$group)]
	all_groups_of_type = unique(clio_manc_data$group[clio_manc_data$type %in% type])
	if(length(all_groups_of_type)>1) paste(type, x) else type
	})
row_groups = paste0(row_groups,"(",temp_cell_count,")") #add cell count
colnames(top_50_ms) = row_groups
pdf(file=paste0("msahn_manc_ahn_eff_conn_str_top_50_", Sys.Date(),".pdf"), height = 300/72, width = 1200/72)
print(heatmap.2(as.matrix(top_50_ms), col = viridis(100),scale="none",density.info='none',xlab="cell type",ylab="AHN",
	dendrogram='none',Rowv=FALSE,Colv=FALSE,trace='none',sepcolor='black', margins = c(5, 5), ColSideColors = graph_lut[top_msahn_neuron_class]))
legend("topright", legend=names(graph_lut)[names(graph_lut) %in% top_msahn_neuron_class], fill=graph_lut[names(graph_lut) %in% top_msahn_neuron_class], horiz=TRUE)
dev.off()

top_50_mt = ahn_all_conn_str[[2]]
top_50_mt = top_50_mt[,order(top_50_mt["MtAHN",],decreasing=T)[1:50]]
top_mtahn_neuron_class = clio_manc_data$abbrv[match(colnames(top_50_mt), clio_manc_data$group)]
#top_mtahn_neuron_class = ifelse(is.na(top_mtahn_neuron_class), clio_manc_data$abbrv[match(colnames(top_50_mt), clio_manc_data$bodyid)], top_mtahn_neuron_class)
top_50_mt = top_50_mt[,order(top_mtahn_neuron_class)]
top_mtahn_neuron_class = top_mtahn_neuron_class[order(top_mtahn_neuron_class)]
temp_cell_count = sapply(colnames(top_50_mt), FUN = function(x) sum(clio_manc_data$group %in% x))
row_groups = sapply(colnames(top_50_mt), FUN=function(x) {
	type = clio_manc_data$type[match(x, clio_manc_data$group)]
	all_groups_of_type = unique(clio_manc_data$group[clio_manc_data$type %in% type])
	if(length(all_groups_of_type)>1) paste(type, x) else type
	})
row_groups = paste0(row_groups,"(",temp_cell_count,")") #add cell count
colnames(top_50_mt) = row_groups
pdf(file=paste0("mtahn_manc_ahn_eff_conn_str_top_50_", Sys.Date(),".pdf"), height = 300/72, width = 1200/72)
print(heatmap.2(as.matrix(top_50_mt), col = viridis(100),scale="none",density.info='none',xlab="cell type",ylab="AHN",
	dendrogram='none',Rowv=FALSE,Colv=FALSE,trace='none',sepcolor='black', margins = c(5, 5), ColSideColors = graph_lut[top_mtahn_neuron_class]))
legend("topright", legend=names(graph_lut)[names(graph_lut) %in% top_mtahn_neuron_class], fill=graph_lut[names(graph_lut) %in% top_mtahn_neuron_class], horiz=TRUE)
dev.off()


#plot ALL indirect conn values rank-ordered to show long tail
temp_indirect_conn_ms = sort(ahn_all_conn_str[[2]]["MsAHN",ahn_all_conn_str[[2]]["MsAHN",]>0], decreasing=TRUE)
temp_indirect_conn_mt = sort(ahn_all_conn_str[[2]]["MtAHN",ahn_all_conn_str[[2]]["MtAHN",]>0], decreasing=TRUE)
len = max(length(temp_indirect_conn_ms), length(temp_indirect_conn_mt))
temp_dat = data.frame(len=1:len, ms=c(temp_indirect_conn_ms,rep(NA,len-length(temp_indirect_conn_ms))),
	mt=c(temp_indirect_conn_mt,rep(NA,len-length(temp_indirect_conn_mt))))
temp_dat = pivot_longer(temp_dat, cols=c('ms','mt'))
pdf(file=paste0("msahn_mtahn_manc_indirect_conn_str_line_histogram_", Sys.Date(),".pdf"), height = 300/72, width = 300/72)
ggplot(temp_dat, aes(x=len, y=value, color=name)) +
  geom_line() + xlab('Neuron type (rank-ordered)') + ylab('indirect connection strength')
dev.off()


#plot their ROI and morphology
MANC.surf.scale = MANC.surf
MANC.surf.scale$Vertices = MANC.surf.scale$Vertices*1000
neuron_color_lut = brewer.pal(10, "Paired")
for (i in c("MsAHN","MtAHN")) {
	top_ahn_neurons = colnames(ahn_all_conn_str[[2]])[order(ahn_all_conn_str[[2]][i,],decreasing=T)[1:5]]
	top_ahn_neurons = clio_manc_data$bodyid[clio_manc_data$group %in% top_ahn_neurons]
	test_mesh = read_manc_meshes(top_ahn_neurons)
	neuron_colors = clio_manc_data$group[match(top_ahn_neurons, clio_manc_data$bodyid)]
	temp_cell_types = unique(neuron_colors)
	for (j in 1:length(temp_cell_types)) neuron_colors[neuron_colors==temp_cell_types[j]] = colorRampPalette(neuron_color_lut[(j*2-1):(j*2)])(sum(neuron_colors==temp_cell_types[j]))
	
	open3d()
	userMatrix = matrix(c(0.009574261, -0.2694314, -0.96296543, 0, -0.984437406, -0.1714953, 0.03819529, 0, -0.175436258,  0.9476197, -0.26688156, 0, 0, 0, 0, 1), nrow = 4, ncol = 4, byrow = TRUE)
	par3d(windowRect=c(0,0,1920,800),zoom=0.45,userMatrix=userMatrix)
	for (j in 1:length(test_mesh)) wire3d(test_mesh[[j]],col = neuron_colors[j], lit = FALSE)
	wire3d(MANC.surf.scale, alpha = 0.04, col = "grey")
	rgl.snapshot(paste0(i, "_effective_conn_str_morphology_top_5_MANC_", Sys.Date(), ".png"),fmt="png")
	rgl.close()
}

for (i in c("MsAHN","MtAHN")) {
	top_ahn_neurons = colnames(ahn_all_conn_str[[2]])[order(ahn_all_conn_str[[2]][i,],decreasing=T)[1:50]]
	top_ahn_neurons = clio_manc_data$bodyid[clio_manc_data$group %in% top_ahn_neurons]
	names(top_ahn_neurons) = clio_manc_data$group[match(top_ahn_neurons, clio_manc_data$bodyid)]
	
	neuron_roi = as.matrix(neuprint_get_roiInfo(top_ahn_neurons, all_segments = TRUE))
	neuron_roi[is.na(neuron_roi)] = 0
	rownames(neuron_roi) = neuron_roi[,"bodyid"]
	neuron_roi_us = neuron_roi[, colnames(neuron_roi) %in% paste0(unlist(han$manc_roi_groups), ".upstream")]
	col_groups = gsub("\\.upstream", "", colnames(neuron_roi_us))
	col_groups = sapply(col_groups, function(x) if (!(x %in% unlist(han$manc_roi_groups))) x else names(han$manc_roi_groups)[sapply(han$manc_roi_groups, function(y) x %in% y)])#combine super-neuropils
	row_groups = clio_manc_data$group[match(rownames(neuron_roi_us), clio_manc_data$bodyid)]
	names(row_groups) = clio_manc_data$abbrv[match(row_groups, clio_manc_data$group)] #cell class
	row_groups = sapply(row_groups, FUN=function(x) {
		type = clio_manc_data$type[match(x, clio_manc_data$group)]
		all_groups_of_type = unique(clio_manc_data$group[clio_manc_data$type %in% type])
		if(length(all_groups_of_type)>1) paste(type, x) else type
		})
	row_groups = sapply(row_groups, function(x) paste0(x,"(",sum(row_groups %in% x),")")) #add cell count


	neuron_roi_us = combine_matrix_rows_cols_by_group(neuron_roi_us, row_groups=row_groups, col_groups=col_groups)
	neuron_roi_us = neuron_roi_us/rowSums(neuron_roi_us)

	neuron_roi_ds = neuron_roi[, colnames(neuron_roi) %in% paste0(unlist(han$manc_roi_groups), ".downstream")]
	col_groups = gsub("\\.downstream", "", colnames(neuron_roi_ds))
	col_groups = sapply(col_groups, function(x) if (!(x %in% unlist(han$manc_roi_groups))) x else names(han$manc_roi_groups)[sapply(han$manc_roi_groups, function(y) x %in% y)])#combine super-neuropils
	neuron_roi_ds = combine_matrix_rows_cols_by_group(neuron_roi_ds, row_groups=row_groups, col_groups=col_groups)
	neuron_roi_ds = neuron_roi_ds/rowSums(neuron_roi_ds)
	neuron_roi_ds[names(row_groups)[match(rownames(neuron_roi_ds), row_groups)]=="MN", ] = 0 #set ds of MNs to zero because we don't believe most of the syn predictions
	
	temp = cbind(neuron_roi_us, neuron_roi_ds[rownames(neuron_roi_us),])
	temp_row_order = simple_row_hierarchical_clustering(temp)
	rois_to_keep = apply(temp, 2, function(x) any(x>=np_thres))
	rois_to_keep = unique(colnames(temp)[rois_to_keep])

	pdf(file=paste0(i, "_effective_conn_str_roi_post_sites_manc_", Sys.Date(),".pdf"), height = 450/72, width = 1200/72)
	temp = t(neuron_roi_us[rownames(temp_row_order),])
	temp = temp[order(match(rownames(temp), names(han$manc_roi_groups))),]
	print(heatmap.2(temp[rownames(temp) %in% rois_to_keep,], col = viridis(100),scale="none",density.info='none',xlab="ds neuron",ylab="Postsynaptic sites",
		dendrogram='none',Rowv=FALSE,Colv=FALSE,trace='none',sepcolor='black', margins = c(10, 10), breaks=seq(0,1,0.01)))
	dev.off()

	pdf(file=paste0(i, "_effective_conn_str_roi_pre_sites_manc_", Sys.Date(),".pdf"), height = 450/72, width = 1200/72)
	temp = t(neuron_roi_ds[rownames(temp_row_order),])
	temp = temp[order(match(rownames(temp), names(han$manc_roi_groups))),]
	print(heatmap.2(temp[rownames(temp) %in% rois_to_keep,], col = viridis(100),scale="none",density.info='none',xlab="ds neuron",ylab="Presynaptic sites",
		dendrogram='none',Rowv=FALSE,Colv=FALSE,trace='none',sepcolor='black', margins = c(10, 10), breaks=seq(0,1,0.01)))
	dev.off()
}

