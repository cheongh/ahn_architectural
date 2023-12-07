library(fafbseg)
library(gplots)
library(viridis)
library(dplyr)
library(tidyr)
library(dendsort)
library(elmr)
library(ggplot2)
library(RColorBrewer)
library(svglite)
library(igraph)
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

all_adj = flywire_adjacency_matrix2(inputids = c(msahns,mtahns), outputids = all_neurons$root_id, sparse = TRUE)
all_adj = all_adj[, Matrix::colSums(all_adj)>0]
#combine neurons by type for neurons we care about
all_dn_types = unique(all_neurons$cell_type[all_neurons$super_class %in% "descending"])
temp_group = all_neurons$cell_type[match(colnames(all_adj), all_neurons$root_id)]
temp_group = ifelse(!is.na(temp_group), temp_group, colnames(all_adj))
temp_group_class = all_neurons$super_class[match(colnames(all_adj), all_neurons$root_id)]
temp_ahn_groups = ifelse(rownames(all_adj) %in% msahns, "MsAHN", "MtAHN")

all_adj = combine_matrix_rows_cols_by_group(all_adj, row_groups = temp_ahn_groups, col_groups=temp_group, sparse = TRUE)

# temp = all_neurons[all_neurons$root_id %in% colnames(all_adj), c('root_id', 'top_nt', 'hemibrain_type', 'super_class', 'cell_class', 'cell_type', 'side', 'ito_lee_hemilineage')]
# write.csv(temp, "AHN direct ds flywire.csv")

#ahn_to_in_adj = all_adj[c("MsAHN","MtAHN"),]
#ahn_to_in_adj = all_adj[,colnames(all_adj) %in% all_dn_types]
	
# pdf(file=paste0("flywire_ahn_direct_ds_", Sys.Date(),".pdf"), height = 300/72, width = 1200/72)
# print(heatmap.2(as.matrix(all_adj), col = viridis(100),scale="none",density.info='none',xlab="ds neuron",ylab="AHN",
	# dendrogram='none',Rowv=FALSE,Colv=FALSE,trace='none',sepcolor='black', margins = c(5, 5)))
# dev.off()

#make boxplot in same method as for manc and take boxplot outliers
syn_count_by_group =as.data.frame(t(as.matrix(all_adj)))
syn_count_by_group$group = rownames(syn_count_by_group)
syn_count_by_group = pivot_longer(syn_count_by_group, cols=c("MsAHN","MtAHN"), names_to="ahn",values_to="weight")
syn_count_by_group = syn_count_by_group[syn_count_by_group$weight>0,]
syn_count_by_group$abbrv = all_neurons$super_class[match(syn_count_by_group$group, all_neurons$cell_type)]
syn_count_by_group$abbrv = superclass_abbrv[syn_count_by_group$abbrv]
syn_count_by_group$outlier=FALSE
syn_count_by_group$outlier[syn_count_by_group$ahn == "MsAHN"] = (syn_count_by_group$weight[syn_count_by_group$ahn == "MsAHN"] 
	> min(max(syn_count_by_group$weight[syn_count_by_group$ahn == "MsAHN"],na.rm=T), as.numeric(quantile(syn_count_by_group$weight[syn_count_by_group$ahn == "MsAHN"], 0.75)) + (IQR(syn_count_by_group$weight[syn_count_by_group$ahn == "MsAHN"])*1.5)))
syn_count_by_group$outlier[syn_count_by_group$ahn == "MtAHN"] = (syn_count_by_group$weight[syn_count_by_group$ahn == "MtAHN"] 
	> min(max(syn_count_by_group$weight[syn_count_by_group$ahn == "MtAHN"],na.rm=T), as.numeric(quantile(syn_count_by_group$weight[syn_count_by_group$ahn == "MtAHN"], 0.75)) + (IQR(syn_count_by_group$weight[syn_count_by_group$ahn == "MtAHN"])*1.5)))
png(paste0("flywire_MsAHN_MtAHN_POST_syn_count_by_group_boxplot_", Sys.Date(), ".png"), height=300,width=300)
print(ggplot() + 
	geom_boxplot(syn_count_by_group, outlier.shape = NA, mapping = aes(x=ahn, y=weight)) +
	geom_point(data = syn_count_by_group[syn_count_by_group$outlier,], mapping = aes(x=ahn, y=weight, color=abbrv), position = position_jitter(w = 0.2)) +
	scale_color_manual(values=graph_lut[names(graph_lut) %in% syn_count_by_group$abbrv[syn_count_by_group$outlier]]) +
	theme(panel.background = element_blank(), axis.line = element_line(color = "darkgrey")) +
	xlab("AHN type") +
	ylab("Synapse weight per cell type"))
ggsave(file=paste0("flywire_MsAHN_MtAHN_POST_syn_count_by_group_boxplot_", Sys.Date(), ".svg"))
dev.off()

#plot morphology of top 5 neuron types ds of ahns 
#only plot subset of neurons of each type that have conn with ahns, because LB3_2 (top type) consists of many types
neuron_color_lut = brewer.pal(10, "Paired")
all_adj_raw = flywire_adjacency_matrix2(inputids = c(msahns,mtahns), outputids = all_neurons$root_id, sparse = TRUE)
all_adj_raw = all_adj_raw[, Matrix::colSums(all_adj_raw)>0]
rownames(all_adj_raw) = names(c(msahns,mtahns))[match(rownames(all_adj_raw), c(msahns,mtahns))]
for (i in c("MsAHN","MtAHN")) {
	top_ahn_neurons = order(syn_count_by_group$weight[syn_count_by_group$ahn==i], decreasing=TRUE)[1:5]
	top_ahn_neurons = syn_count_by_group$group[syn_count_by_group$ahn==i][top_ahn_neurons]
	top_ahn_neurons_group2id = all_neurons$root_id[all_neurons$cell_type %in% top_ahn_neurons]
	top_ahn_neurons_group2id = top_ahn_neurons_group2id[top_ahn_neurons_group2id %in% colnames(all_adj_raw)]
	top_ahn_neurons = c(top_ahn_neurons[top_ahn_neurons %in% all_neurons$root_id], top_ahn_neurons_group2id)
	test_mesh = read_cloudvolume_meshes(top_ahn_neurons)
	neuron_colors = all_neurons$cell_type[match(names(test_mesh), all_neurons$root_id)]
	temp_cell_types = unique(neuron_colors)
	for (j in 1:length(temp_cell_types)) neuron_colors[neuron_colors==temp_cell_types[j]] = colorRampPalette(neuron_color_lut[(j*2-1):(j*2)])(sum(neuron_colors==temp_cell_types[j]))

	nview3d("anterior", zoom = 0.4690333, windowRect = c(45L, 45L, 2141L, 1079L))
	plot3d(FAFB.surf, color='grey', alpha = 0.1)
	plot3d(test_mesh, lit = FALSE, color=neuron_colors)
	rgl.snapshot(paste(i, "direct_ds_top_5_types_morphology_flywire.png", sep = "_"),fmt="png")
	rgl.close()
}


#get neuropil ramification for top ds partners (boxplot outliers)


np_thres = 0.1

top_msahn_neurons = syn_count_by_group$group[syn_count_by_group$outlier & syn_count_by_group$ahn=="MsAHN"]
top_msahn_neurons = all_neurons$root_id[all_neurons$cell_type %in% top_msahn_neurons]
top_msahn_neurons = top_msahn_neurons[top_msahn_neurons %in% colnames(all_adj_raw)]
names(top_msahn_neurons) = all_neurons$cell_type[match(top_msahn_neurons, all_neurons$root_id)]
top_mtahn_neurons = syn_count_by_group$group[syn_count_by_group$outlier & syn_count_by_group$ahn=="MtAHN"]
top_mtahn_neurons = all_neurons$root_id[all_neurons$cell_type %in% top_mtahn_neurons]
top_mtahn_neurons = top_mtahn_neurons[top_mtahn_neurons %in% colnames(all_adj_raw)]
names(top_mtahn_neurons) = all_neurons$cell_type[match(top_mtahn_neurons, all_neurons$root_id)]
neurons_of_interest = list(msahn = top_msahn_neurons, mtahn = top_mtahn_neurons)

for (i in names(neurons_of_interest)) {
	ahn_ds_ds_pre = get_flywire_syn_neuropil(neurons_of_interest[[i]], "PRE")
	ahn_ds_ds_pre[is.na(ahn_ds_ds_pre)] = 0
	col_groups = gsub("_[LR]", "", colnames(ahn_ds_ds_pre)) #combine L and R neuropils
	col_groups = sapply(col_groups, function(x) if (!(x %in% unlist(super_np))) x else names(super_np)[sapply(super_np, function(y) x %in% y)])#combine super-neuropils
	row_groups = all_neurons$cell_type[match(rownames(ahn_ds_ds_pre), all_neurons$root_id)]
	row_groups = sapply(row_groups, function(x) paste0(x,"(",sum(names(neurons_of_interest[[i]]) %in% x),")")) #add cell count
	#ahn_ds_ds_pre = t(apply(t(ahn_ds_ds_pre), 2, function(x) tapply(x, col_groups, sum, na.rm = TRUE)))
	ahn_ds_ds_pre = combine_matrix_rows_cols_by_group(ahn_ds_ds_pre, row_groups=row_groups, col_groups=col_groups)
	ahn_ds_ds_pre = ahn_ds_ds_pre/rowSums(ahn_ds_ds_pre)

	ahn_ds_ds_post = get_flywire_syn_neuropil(neurons_of_interest[[i]], "POST")
	ahn_ds_ds_post[is.na(ahn_ds_ds_post)] = 0
	col_groups = gsub("_[LR]", "", colnames(ahn_ds_ds_post)) #combine L and R neuropils
	col_groups = sapply(col_groups, function(x) if (!(x %in% unlist(super_np))) x else names(super_np)[sapply(super_np, function(y) x %in% y)])#combine super-neuropils
	row_groups = all_neurons$cell_type[match(rownames(ahn_ds_ds_post), all_neurons$root_id)]
	row_groups = sapply(row_groups, function(x) paste0(x,"(",sum(names(neurons_of_interest[[i]]) %in% x),")")) #add cell count
	#ahn_ds_ds_post = t(apply(t(ahn_ds_ds_post), 2, function(x) tapply(x, col_groups, sum, na.rm = TRUE)))
	ahn_ds_ds_post = combine_matrix_rows_cols_by_group(ahn_ds_ds_post, row_groups=row_groups, col_groups=col_groups)
	ahn_ds_ds_post = ahn_ds_ds_post/rowSums(ahn_ds_ds_post)

	temp = cbind(ahn_ds_ds_pre, ahn_ds_ds_post[rownames(ahn_ds_ds_pre),])
	temp_row_order = simple_row_hierarchical_clustering(temp)
	#temp_col_order = simple_row_hierarchical_clustering(t(temp))
	rois_to_keep = apply(temp, 2, function(x) any(x>=np_thres))
	rois_to_keep = unique(colnames(temp)[rois_to_keep])

	pdf(file=paste0(i, "_ds_roi_post_sites_", Sys.Date(),".pdf"), height = 450/72, width = 1200/72)
	temp = t(ahn_ds_ds_pre[rownames(temp_row_order),])
	temp = temp[order(match(rownames(temp), all_np_order)),]
	print(heatmap.2(temp[rownames(temp) %in% rois_to_keep,], col = viridis(100),scale="none",density.info='none',xlab="ds neuron",ylab="Postsynaptic sites",
		dendrogram='none',Rowv=FALSE,Colv=FALSE,trace='none',sepcolor='black', margins = c(10, 10), breaks=seq(0,1,0.01)))
	dev.off()

	pdf(file=paste0(i, "_ds_roi_pre_sites_", Sys.Date(),".pdf"), height = 450/72, width = 1200/72)
	temp = t(ahn_ds_ds_post[rownames(temp_row_order),])
	temp = temp[order(match(rownames(temp), all_np_order)),]
	print(heatmap.2(temp[rownames(temp) %in% rois_to_keep,], col = viridis(100),scale="none",density.info='none',xlab="ds neuron",ylab="Presynaptic sites",
		dendrogram='none',Rowv=FALSE,Colv=FALSE,trace='none',sepcolor='black', margins = c(10, 10), breaks=seq(0,1,0.01)))
	dev.off()
}















#bar chart of ahns input frac by cell class
all_adj_raw = flywire_adjacency_matrix2(inputids = c(msahns,mtahns), outputids = all_neurons$root_id, sparse = TRUE)
all_adj_raw = all_adj_raw[, Matrix::colSums(all_adj_raw)>0]
all_adj_raw =as.data.frame(t(as.matrix(all_adj_raw)))
colnames(all_adj_raw) = names(c(msahns,mtahns))[match(colnames(all_adj_raw), c(msahns,mtahns))]
all_adj_raw$root_id = rownames(all_adj_raw)
all_adj_raw$super_class = all_neurons$super_class[match(all_adj_raw$root_id, all_neurons$root_id)]
all_adj_raw$abbrv = superclass_abbrv[all_adj_raw$super_class]
all_adj_raw = pivot_longer(all_adj_raw, cols=c("MsAHN-R", "MsAHN-L", "MtAHN-L", "MtAHN-R"), names_to="ahn",values_to="weight")
per_class_summary = aggregate(weight~abbrv+ahn, data=all_adj_raw, FUN=sum)
per_class_summary$abbrv = factor(per_class_summary$abbrv, levels=names(graph_lut)[names(graph_lut) %in% per_class_summary$abbrv])

svglite(paste0("flywire_ahn_ds_synapse_frac_", Sys.Date(),".svg"), width = 6, height = 5)
print(ggplot(per_class_summary, aes(fill=abbrv, y=weight, x=ahn)) + 
	geom_bar(position="fill", stat="identity") +
	xlab("") +
	ylab("% synapses") +
	scale_x_discrete() +
	theme(text = element_text(size=15), axis.line.y.left = element_line(color = "black"), axis.line.x.bottom = element_line(color = "black"), panel.background = element_blank()) +
	scale_fill_manual(values = graph_lut[names(graph_lut) %in% per_class_summary$abbrv]))
dev.off()


#graph of ds neurons binned by connectivity with ahn type and side
#overlay pie chart with output fractions by cell class
all_adj_raw = flywire_adjacency_matrix2(inputids = c(msahns,mtahns), outputids = all_neurons$root_id, sparse = TRUE)
all_adj_raw = all_adj_raw[, Matrix::colSums(all_adj_raw)>0]
all_adj_raw =as.data.frame(t(as.matrix(all_adj_raw)))
colnames(all_adj_raw) = names(c(msahns,mtahns))[match(colnames(all_adj_raw), c(msahns,mtahns))]
all_adj_raw$root_id = rownames(all_adj_raw)
all_adj_raw$super_class = all_neurons$super_class[match(all_adj_raw$root_id, all_neurons$root_id)]
all_adj_raw$abbrv = superclass_abbrv[all_adj_raw$super_class]
all_adj_raw$conn_group = apply(all_adj_raw[,c("MsAHN-L", "MsAHN-R", "MtAHN-L", "MtAHN-R")], 1, FUN = function(x) {
	if (sum(x>0)>2) "all_ahn" else {temp = paste(as.numeric(x>0),collapse="")
		if(temp %in% c("1001","0110")) "all_ahn" else temp} #if L and R ahns of different types are upstream, lump into all ahn group
	})
all_adj_raw$total_weight = rowSums(all_adj_raw[,c("MsAHN-L", "MsAHN-R", "MtAHN-L", "MtAHN-R")]) 
conn_group_summary = aggregate(total_weight~conn_group+abbrv, data=all_adj_raw, FUN=sum)
#sum of cell count per conn_group NOT considering abbrv
conn_group_summary$cell_count = sapply(conn_group_summary$conn_group, FUN=function(x) sum(all_adj_raw$conn_group==x))
#make the graph
adj_graph = combine_matrix_rows_cols_by_group(as.matrix(all_adj_raw[,c("MsAHN-L", "MsAHN-R", "MtAHN-L", "MtAHN-R")]),
	col_groups = c("MsAHN-L", "MsAHN-R", "MtAHN-L", "MtAHN-R"), row_groups = all_adj_raw$conn_group)
adj_graph = t(adj_graph)
adj_graph = squarify(adj_graph)
graph = graph_from_adjacency_matrix(adjmatrix = adj_graph, mode = "directed", weighted = TRUE, diag = FALSE)
classes_of_interest = superclass_abbrv[superclass_abbrv %in% all_adj_raw$abbrv]
pie_values = as.data.frame(pivot_wider(conn_group_summary[,c("conn_group","abbrv","total_weight")], names_from=abbrv, values_from=total_weight, values_fill = 0))
rownames(pie_values) = pie_values$conn_group
temp_vert_attr = as.list(pie_values[match(V(graph)$name,pie_values$conn_group),])
temp_vert_attr = append(temp_vert_attr, list(name=V(graph)$name, cell_count=conn_group_summary$cell_count[match(V(graph)$name, conn_group_summary$conn_group)]))
vertex_attr(graph) = temp_vert_attr
V(graph)$label = ifelse(is.na(V(graph)$cell_count), V(graph)$name, V(graph)$cell_count)
V(graph)$color = ifelse(V(graph)$name %in% names(mtahns), graph_lut["MtAHN"], ifelse(V(graph)$name %in% names(msahns), graph_lut["MsAHN"], NA))
RCy3::createNetworkFromIgraph(graph, title="flywire_ahn_ds")
RCy3::setNodeCustomPieChart(columns=unname(classes_of_interest), colors=unname(graph_lut[classes_of_interest]), style.name = "ahn_paper_pie")




