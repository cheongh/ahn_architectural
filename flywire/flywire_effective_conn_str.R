library(fafbseg)
library(gplots)
library(viridis)
library(reticulate)
library(dplyr)
library(tidyr)
library(dendsort)
library(elmr)
library(ggplot2)
library(RColorBrewer)
source("shared_functions.R")

choose_segmentation("public-flywire31")

eff_conn_thres = 0.0005 #reasonably gets top ~50 types ds of MsAHN but hardly any of MtAHN #not used
#threshold of how many input synapses a neuron must have to consider effective conn str 'real'
#doing this because some poorly reconstructed neurons appear as top hits
total_syn_thres = 100


all_neurons = flytable_meta()
all_neurons = all_neurons[!(all_neurons$super_class %in% c(NA, "endocrine", "not_a_neuron")),]
schlegel_groups = read.csv("schlegel_2023_sup_3_morphological_groups.csv", colClasses = 'character') #if cell does not have type, use morphological groups as defined in Schlegel 2023 flywire paper
all_neurons$cell_type = ifelse(is.na(all_neurons$cell_type) | all_neurons$cell_type == "",
	schlegel_groups$morphological_group[match(all_neurons$root_id, schlegel_groups$root_630)],
	all_neurons$cell_type)
han_groups = read.csv("MsAHN ds flywire OLC neurons.csv", colClasses = 'character') #han groups for optic lobe centrifugal neurons ds of MsAHNs
all_neurons$cell_type = ifelse(is.na(all_neurons$cell_type) | all_neurons$cell_type == "",
	han_groups$hc_cell_type[match(all_neurons$root_id, han_groups$root_id)],
	all_neurons$cell_type)
#commands used to get this list
# top_msahn_neurons = names(sort(ahn_effective_conn_str["MsAHN",], decreasing=TRUE)[1:1000])
#  temp2 = all_neurons[match(top_msahn_neurons[!(top_msahn_neurons %in% all_neurons$cell_type)],all_neurons$root_id), c('root_id', 'top_nt', 'hemibrain_type', 'super_class', 'cell_class', 'cell_type', 'side', 'ito_lee_hemilineage')]
# temp2 = temp2[temp2$super_class %in% 'visual_centrifugal',]
# temp2$hc_cell_type = temp3$hc_cell_type[match(temp2$root_id, temp3$root_id)]
# temp2$cell_type[temp2$cell_type %in% "OLT1"] = "OLC1"
# write.csv(temp2, "MsAHN ds flywire OLC neurons.csv")
han_groups_2 = read.csv("MtAHN ds flywire OL neurons.csv", colClasses = 'character')
# top_mtahn_neurons = names(sort(ahn_effective_conn_str["MtAHN",], decreasing=TRUE)[1:1000])
# temp2 = all_neurons[match(top_mtahn_neurons[!(top_mtahn_neurons %in% all_neurons$cell_type)],all_neurons$root_id), c('root_id', 'top_nt', 'hemibrain_type', 'super_class', 'cell_class', 'cell_type', 'side', 'ito_lee_hemilineage')]
# temp2 = temp2[temp2$super_class %in% 'optic',]
# write.csv(temp2, "MtAHN ds flywire OL neurons.csv")
all_neurons$cell_type = ifelse(is.na(all_neurons$cell_type) | all_neurons$cell_type == "",
	han_groups_2$hc_cell_type[match(all_neurons$root_id, han_groups_2$root_id)],
	all_neurons$cell_type)
han_groups_3 = read.csv("AHN direct ds flywire.csv", colClasses = 'character')
han_groups_3$hc_cell_type = ifelse(han_groups_3$hc_cell_type=="",han_groups_3$root_id,han_groups_3$hc_cell_type)
all_neurons$cell_type = ifelse(is.na(all_neurons$cell_type) | all_neurons$cell_type == "",
	han_groups_3$hc_cell_type[match(all_neurons$root_id, han_groups_3$root_id)],
	all_neurons$cell_type)
all_neurons$cell_type[all_neurons$root_id %in% mtahns] = "MtAHN"
all_neurons$cell_type[all_neurons$root_id %in% msahns] = "MsAHN"
all_neurons$cell_type = ifelse(is.na(all_neurons$cell_type) | all_neurons$cell_type == "", all_neurons$root_id, all_neurons$cell_type)

all_adj = flywire_adjacency_matrix2(all_neurons$root_id, sparse = TRUE)
#combine neurons by type for neurons we care about (AHNs and DNs)
all_dn_types = unique(all_neurons$cell_type[all_neurons$super_class %in% "descending"])
temp_group = all_neurons$cell_type[match(colnames(all_adj), all_neurons$root_id)]
temp_group = ifelse(!is.na(temp_group), temp_group, colnames(all_adj))
temp_group_class = all_neurons$super_class[match(colnames(all_adj), all_neurons$root_id)]



all_adj = combine_matrix_rows_cols_by_group(all_adj, groups=temp_group, sparse = TRUE)

all_adj_raw = all_adj

#normalize to input fracs and then remove input to sensory and ascending to avoid considering paths through them
all_adj = all_adj %*% Matrix::Diagonal(x = 1/Matrix::colSums(all_adj))
all_adj[,colnames(all_adj) %in% temp_group[temp_group_class %in% "sensory"]] = 0
all_adj[,colnames(all_adj) %in% temp_group[temp_group_class %in% "ascending"]] = 0
#all_adj[,c("MsAHN","MtAHN")] = 0
diag(all_adj) = 0

#ahn_to_in_adj = all_adj[c("MsAHN","MtAHN"),]
#ahn_to_in_adj = all_adj[,colnames(all_adj) %in% all_dn_types]
ahn_effective_conn_str = all_adj[c("MsAHN","MtAHN"),] %*% all_adj
ahn_effective_conn_str = ahn_effective_conn_str[,order(colSums(as.matrix(ahn_effective_conn_str)), decreasing=TRUE)]
per_neuron_input = Matrix::colSums(all_adj_raw)
per_neuron_input = per_neuron_input/sapply(names(per_neuron_input), function(x) sum(all_neurons$cell_type %in% x))
per_neuron_input = per_neuron_input[per_neuron_input>=total_syn_thres]
ahn_effective_conn_str = ahn_effective_conn_str[, colnames(ahn_effective_conn_str) %in% names(per_neuron_input)]

ahn_effective_conn_str_dns = ahn_effective_conn_str[,colnames(ahn_effective_conn_str) %in% all_dn_types]

ahn_effective_conn_str_thres = ahn_effective_conn_str[, apply(ahn_effective_conn_str, 2, max)>=eff_conn_thres]
categories = apply(ahn_effective_conn_str_thres, 2, function(x) if(all(x>=eff_conn_thres)) 2 else if(x[1]>=eff_conn_thres) 1 else 3)
ahn_effective_conn_str_thres = ahn_effective_conn_str_thres[,order(categories)]
categories = apply(ahn_effective_conn_str_thres, 2, function(x) if(all(x>=eff_conn_thres)) 2 else if(x[1]>=eff_conn_thres) 1 else 3)
row_color_lut = c("green", "magenta", "cyan")

pdf(file=paste0("flywire_ahn_effective_conn_str_", Sys.Date(),".pdf"), height = 300/72, width = 1200/72)
print(heatmap.2(as.matrix(ahn_effective_conn_str_thres), col = viridis(100),scale="none",density.info='none',xlab="ds neuron",ylab="AHN",
	dendrogram='none',Rowv=FALSE,Colv=FALSE,trace='none',sepcolor='black', margins = c(5, 5), ColSideColors = row_color_lut[categories]))
dev.off()

	
pdf(file=paste0("flywire_ahn_effective_conn_str_", Sys.Date(),".pdf"), height = 300/72, width = 1200/72)
print(heatmap.2(as.matrix(ahn_effective_conn_str[,1:200]), col = viridis(100),scale="none",density.info='none',xlab="DN type",ylab="AHN",
	dendrogram='none',Rowv=FALSE,Colv=FALSE,trace='none',sepcolor='black', margins = c(5, 5)))
dev.off()

top_msahn_neurons = names(sort(ahn_effective_conn_str["MsAHN",], decreasing=TRUE)[1:50])
top_msahn_neuron_class = all_neurons$super_class[match(top_msahn_neurons, all_neurons$cell_type)]
top_msahn_neuron_class = superclass_abbrv[top_msahn_neuron_class]
top_msahn_neurons = top_msahn_neurons[order(top_msahn_neuron_class)]
top_msahn_neuron_class = top_msahn_neuron_class[order(top_msahn_neuron_class)]
temp_adj = as.matrix(ahn_effective_conn_str[,top_msahn_neurons])
colnames(temp_adj) = sapply(colnames(temp_adj), function(x) paste0(x,"(",sum(all_neurons$cell_type %in% x),")")) #add cell count
pdf(file=paste0("msahn_flywire_ahn_eff_conn_str_top_50_", Sys.Date(),".pdf"), height = 300/72, width = 1200/72)
print(heatmap.2(as.matrix(temp_adj), col = viridis(100),scale="none",density.info='none',xlab="cell type",ylab="AHN",
	dendrogram='none',Rowv=FALSE,Colv=FALSE,trace='none',sepcolor='black', margins = c(5, 5), ColSideColors = graph_lut[top_msahn_neuron_class]))
legend("topright", legend=names(graph_lut)[names(graph_lut) %in% top_msahn_neuron_class], fill=graph_lut[names(graph_lut) %in% top_msahn_neuron_class], horiz=TRUE)
dev.off()

#plot ALL indirect conn values rank-ordered to show long tail
temp_indirect_conn_ms = sort(ahn_effective_conn_str["MsAHN",ahn_effective_conn_str["MsAHN",]>0], decreasing=TRUE)
temp_indirect_conn_mt = sort(ahn_effective_conn_str["MtAHN",ahn_effective_conn_str["MtAHN",]>0], decreasing=TRUE)
len = max(length(temp_indirect_conn_ms), length(temp_indirect_conn_mt))
temp_dat = data.frame(len=1:len, ms=c(temp_indirect_conn_ms,rep(NA,len-length(temp_indirect_conn_ms))),
	mt=c(temp_indirect_conn_mt,rep(NA,len-length(temp_indirect_conn_mt))))
temp_dat = pivot_longer(temp_dat, cols=c('ms','mt'))
pdf(file=paste0("msahn_mtahn_flywire_indirect_conn_str_line_histogram_", Sys.Date(),".pdf"), height = 300/72, width = 300/72)
ggplot(temp_dat, aes(x=len, y=value, color=name)) +
  geom_line() + xlab('Neuron type (rank-ordered)') + ylab('indirect connection strength')
dev.off()

#temp plot graph of connectivity for example neurons from msahn to target
types_of_interest = c("MsAHN", "DNg18_a")
inbetween_types = rownames(all_adj_raw)[all_adj_raw[types_of_interest[1],]>1 & all_adj_raw[,types_of_interest[2]]>1]
temp_adj_2 = all_adj_raw[c(types_of_interest,inbetween_types), c(types_of_interest,inbetween_types)]
temp_adj_2[-1,-2] = 0
temp_graph = graph_from_adjacency_matrix(adjmatrix = temp_adj_2, mode = "directed", weighted = TRUE, diag = FALSE)
tkplot_graph = tkplot(temp_graph, canvas.width = 650, canvas.height = 650,
	edge.arrow.size = 0.8, edge.arrow.width = 1.5, curved = TRUE, edge.label = E(temp_graph)$weight, edge.width = 2,
	vertex.label = V(temp_graph)$name, vertex.color = V(temp_graph)$color,
	margin = c(0,0,0,0), layout = layout_with_graphopt(temp_graph, spring.constant = 5, charge = 0.01, spring.length = 50), vertex.label.family = "sans")
readline(prompt="Adjust graph layout by dragging nodes around in tkplot window.\nPress ENTER once done with adjusting node positions: ")
tkplot_layout = tkplot.getcoords(tkplot_graph)
tkplot.close(tkplot_graph)
#manually tweak layout and get coordinates
#tkplot_layout = tkplot.getcoords(3)
#after saving coords, go back and relabel vertices to proper names
svglite(paste0(paste(types_of_interest,collapse="_"), "_in_between_graph_", Sys.Date(),".svg"), width = 5, height = 5)
print(plot.igraph(temp_graph, edge.arrow.size = 0.4,
	edge.arrow.width = 1.5, curved = TRUE, edge.label = E(temp_graph)$weight, edge.width = log1p(E(temp_graph)$weight)/2,
	margin = c(0,0,0,0), layout = tkplot_layout, vertex.label.family = "sans", vertex.label.cex = 0.65, edge.label.family = "sans", edge.label.cex = 0.6))
dev.off()

types_of_interest = c("MsAHN", "WED130")


#now do the same for mtahn
top_mtahn_neurons = names(sort(ahn_effective_conn_str["MtAHN",], decreasing=TRUE)[1:50])
top_mtahn_neuron_class = all_neurons$super_class[match(top_mtahn_neurons, all_neurons$cell_type)]
top_mtahn_neuron_class = superclass_abbrv[top_mtahn_neuron_class]
top_mtahn_neurons = top_mtahn_neurons[order(top_mtahn_neuron_class)]
top_mtahn_neuron_class = top_mtahn_neuron_class[order(top_mtahn_neuron_class)]
temp_adj = as.matrix(ahn_effective_conn_str[,top_mtahn_neurons])
colnames(temp_adj) = sapply(colnames(temp_adj), function(x) paste0(x,"(",sum(all_neurons$cell_type %in% x),")")) #add cell count
pdf(file=paste0("mtahn_flywire_ahn_eff_conn_str_top_50_", Sys.Date(),".pdf"), height = 300/72, width = 1200/72)
print(heatmap.2(as.matrix(temp_adj), col = viridis(100),scale="none",density.info='none',xlab="cell type",ylab="AHN",
	dendrogram='none',Rowv=FALSE,Colv=FALSE,trace='none',sepcolor='black', margins = c(5, 5), ColSideColors = graph_lut[top_mtahn_neuron_class]))
legend("topright", legend=names(graph_lut)[names(graph_lut) %in% top_mtahn_neuron_class], fill=graph_lut[names(graph_lut) %in% top_mtahn_neuron_class], horiz=TRUE)
dev.off()

#temp plot graph of connectivity for example neurons from mtahn to target
types_of_interest = c("MtAHN", "DNge046")
inbetween_types = rownames(all_adj_raw)[all_adj_raw[types_of_interest[1],]>1 & all_adj_raw[,types_of_interest[2]]>1]
temp_adj_2 = all_adj_raw[c(types_of_interest,inbetween_types), c(types_of_interest,inbetween_types)]
temp_adj_2[-1,-2] = 0
temp_graph = graph_from_adjacency_matrix(adjmatrix = temp_adj_2, mode = "directed", weighted = TRUE, diag = FALSE)
tkplot_graph = tkplot(temp_graph, canvas.width = 650, canvas.height = 650,
	edge.arrow.size = 0.8, edge.arrow.width = 1.5, curved = TRUE, edge.label = E(temp_graph)$weight, edge.width = 2,
	vertex.label = V(temp_graph)$name, vertex.color = V(temp_graph)$color,
	margin = c(0,0,0,0), layout = layout_with_graphopt(temp_graph, spring.constant = 5, charge = 0.01, spring.length = 50), vertex.label.family = "sans")
readline(prompt="Adjust graph layout by dragging nodes around in tkplot window.\nPress ENTER once done with adjusting node positions: ")
tkplot_layout = tkplot.getcoords(tkplot_graph)
tkplot.close(tkplot_graph)
#manually tweak layout and get coordinates
#tkplot_layout = tkplot.getcoords(3)
#after saving coords, go back and relabel vertices to proper names
svglite(paste0(paste(types_of_interest,collapse="_"), "_in_between_graph_", Sys.Date(),".svg"), width = 5, height = 5)
print(plot.igraph(temp_graph, edge.arrow.size = 0.4,
	edge.arrow.width = 1.5, curved = TRUE, edge.label = E(temp_graph)$weight, edge.width = log1p(E(temp_graph)$weight)/2,
	margin = c(0,0,0,0), layout = tkplot_layout, vertex.label.family = "sans", vertex.label.cex = 0.65, edge.label.family = "sans", edge.label.cex = 0.6))
dev.off()

types_of_interest = c("MtAHN", "PVLP046")


#top_msahn_neurons = ahn_effective_conn_str["MsAHN", ahn_effective_conn_str["MsAHN",]>=0.0005]
#top_msahn_neurons = names(top_msahn_neurons)
#top_msahn_neurons_group2id = all_neurons$root_id[all_neurons$cell_type %in% top_msahn_neurons]
#top_msahn_neurons = c(top_msahn_neurons[top_msahn_neurons %in% all_neurons$root_id], top_msahn_neurons_group2id)
#top_msahn_neurons = all_neurons[all_neurons$root_id %in% top_msahn_neurons, c('root_id', 'top_nt', 'hemibrain_type', 'super_class', 'cell_class', 'cell_type', 'side', 'ito_lee_hemilineage')]
#top_msahn_neurons$eff_conn_str = sapply(top_msahn_neurons$root_id, FUN = function(x) if(x %in% colnames(ahn_effective_conn_str)) ahn_effective_conn_str["MsAHN", x] else NA)
#top_msahn_neurons = top_msahn_neurons[order(top_msahn_neurons$eff_conn_str, decreasing=T),]
#write.csv(top_msahn_neurons, "top_indirect_msahn_targets_flywire.csv")
neuron_color_lut = brewer.pal(10, "Paired")

top_msahn_neurons = names(sort(ahn_effective_conn_str["MsAHN",], decreasing=TRUE)[1:5])
top_msahn_neurons_group2id = all_neurons$root_id[all_neurons$cell_type %in% top_msahn_neurons]
top_msahn_neurons = c(top_msahn_neurons[top_msahn_neurons %in% all_neurons$root_id], top_msahn_neurons_group2id)
test_mesh = read_cloudvolume_meshes(top_msahn_neurons)
neuron_colors = all_neurons$cell_type[match(names(test_mesh), all_neurons$root_id)]
temp_cell_types = unique(neuron_colors)
for (i in 1:length(temp_cell_types)) neuron_colors[neuron_colors==temp_cell_types[i]] = colorRampPalette(neuron_color_lut[(i*2-1):(i*2)])(sum(neuron_colors==temp_cell_types[i]))

nview3d("anterior", zoom = 0.4690333, windowRect = c(45L, 45L, 2141L, 1079L))
plot3d(FAFB.surf, color='grey', alpha = 0.1)
plot3d(test_mesh, lit = FALSE, color=neuron_colors)
rgl.snapshot(paste("msahn_effective_con_str_top_5_types_morphology.png", sep = "_"),fmt="png")
rgl.close()

top_mtahn_neurons = names(sort(ahn_effective_conn_str["MtAHN",], decreasing=TRUE)[1:5])
top_mtahn_neurons_group2id = all_neurons$root_id[all_neurons$cell_type %in% top_mtahn_neurons]
top_mtahn_neurons = c(top_mtahn_neurons[top_mtahn_neurons %in% all_neurons$root_id], top_mtahn_neurons_group2id)
test_mesh = read_cloudvolume_meshes(top_mtahn_neurons)
neuron_colors = all_neurons$cell_type[match(names(test_mesh), all_neurons$root_id)]
temp_cell_types = unique(neuron_colors)
for (i in 1:length(temp_cell_types)) neuron_colors[neuron_colors==temp_cell_types[i]] = colorRampPalette(neuron_color_lut[(i*2-1):(i*2)])(sum(neuron_colors==temp_cell_types[i]))

nview3d("anterior", zoom = 0.4690333, windowRect = c(45L, 45L, 2141L, 1079L))
plot3d(FAFB.surf, color='grey', alpha = 0.1)
plot3d(test_mesh, lit = FALSE, color=neuron_colors)
rgl.snapshot(paste("mtahn_effective_con_str_top_5_types_morphology.png", sep = "_"),fmt="png")
rgl.close()

#get neuropil ramification

np_thres = 0.1

top_msahn_neurons = names(sort(ahn_effective_conn_str["MsAHN",], decreasing=TRUE)[1:50])
top_msahn_neurons_group2id = all_neurons$root_id[all_neurons$cell_type %in% top_msahn_neurons]
top_msahn_neurons = c(top_msahn_neurons[top_msahn_neurons %in% all_neurons$root_id], top_msahn_neurons_group2id)
top_mtahn_neurons = names(sort(ahn_effective_conn_str["MtAHN",], decreasing=TRUE)[1:50])
top_mtahn_neurons_group2id = all_neurons$root_id[all_neurons$cell_type %in% top_mtahn_neurons]
top_mtahn_neurons = c(top_mtahn_neurons[top_mtahn_neurons %in% all_neurons$root_id], top_mtahn_neurons_group2id)
neurons_of_interest = list(msahn = top_msahn_neurons, mtahn = top_mtahn_neurons)

for (i in names(neurons_of_interest)) {
	ahn_ds_ds_pre = get_flywire_syn_neuropil(neurons_of_interest[[i]], "PRE")
	ahn_ds_ds_pre[is.na(ahn_ds_ds_pre)] = 0
	col_groups = gsub("_[LR]", "", colnames(ahn_ds_ds_pre)) #combine L and R neuropils
	col_groups = sapply(col_groups, function(x) if (!(x %in% unlist(super_np))) x else names(super_np)[sapply(super_np, function(y) x %in% y)])#combine super-neuropils
	row_groups = all_neurons$cell_type[match(rownames(ahn_ds_ds_pre), all_neurons$root_id)]
	row_groups = sapply(row_groups, function(x) paste0(x,"(",sum(all_neurons$cell_type %in% x),")")) #add cell count
	#ahn_ds_ds_pre = t(apply(t(ahn_ds_ds_pre), 2, function(x) tapply(x, col_groups, sum, na.rm = TRUE)))
	ahn_ds_ds_pre = combine_matrix_rows_cols_by_group(ahn_ds_ds_pre, row_groups=row_groups, col_groups=col_groups)
	ahn_ds_ds_pre = ahn_ds_ds_pre/rowSums(ahn_ds_ds_pre)

	ahn_ds_ds_post = get_flywire_syn_neuropil(neurons_of_interest[[i]], "POST")
	ahn_ds_ds_post[is.na(ahn_ds_ds_post)] = 0
	col_groups = gsub("_[LR]", "", colnames(ahn_ds_ds_post)) #combine L and R neuropils
	col_groups = sapply(col_groups, function(x) if (!(x %in% unlist(super_np))) x else names(super_np)[sapply(super_np, function(y) x %in% y)])#combine super-neuropils
	row_groups = all_neurons$cell_type[match(rownames(ahn_ds_ds_post), all_neurons$root_id)]
	row_groups = sapply(row_groups, function(x) paste0(x,"(",sum(all_neurons$cell_type %in% x),")")) #add cell count
	#ahn_ds_ds_post = t(apply(t(ahn_ds_ds_post), 2, function(x) tapply(x, col_groups, sum, na.rm = TRUE)))
	ahn_ds_ds_post = combine_matrix_rows_cols_by_group(ahn_ds_ds_post, row_groups=row_groups, col_groups=col_groups)
	ahn_ds_ds_post = ahn_ds_ds_post/rowSums(ahn_ds_ds_post)

	temp = cbind(ahn_ds_ds_pre, ahn_ds_ds_post[rownames(ahn_ds_ds_pre),])
	temp_row_order = simple_row_hierarchical_clustering(temp)
	#temp_col_order = simple_row_hierarchical_clustering(t(temp))
	rois_to_keep = apply(temp, 2, function(x) any(x>=np_thres))
	rois_to_keep = unique(colnames(temp)[rois_to_keep])

	pdf(file=paste0(i, "_ds_ds_roi_post_sites_", Sys.Date(),".pdf"), height = 450/72, width = 1200/72)
	temp = t(ahn_ds_ds_pre[rownames(temp_row_order),])
	temp = temp[order(match(rownames(temp), all_np_order)),]
	print(heatmap.2(temp[rownames(temp) %in% rois_to_keep,], col = viridis(100),scale="none",density.info='none',xlab="ds neuron",ylab="Postsynaptic sites",
		dendrogram='none',Rowv=FALSE,Colv=FALSE,trace='none',sepcolor='black', margins = c(10, 10), breaks=seq(0,1,0.01)))
	dev.off()

	pdf(file=paste0(i, "_ds_ds_roi_pre_sites_", Sys.Date(),".pdf"), height = 450/72, width = 1200/72)
	temp = t(ahn_ds_ds_post[rownames(temp_row_order),])
	temp = temp[order(match(rownames(temp), all_np_order)),]
	print(heatmap.2(temp[rownames(temp) %in% rois_to_keep,], col = viridis(100),scale="none",density.info='none',xlab="ds neuron",ylab="Presynaptic sites",
		dendrogram='none',Rowv=FALSE,Colv=FALSE,trace='none',sepcolor='black', margins = c(10, 10), breaks=seq(0,1,0.01)))
	dev.off()
}