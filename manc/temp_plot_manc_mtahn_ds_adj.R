library(neuprintr)
box::use(mancfuns/han)
box::use(mancfuns/manc/annotation)
library(gplots)
library(viridis)
library(malevnc)
library(igraph)
library(svglite)

ds_syn_thres = 3
mtahn = c(42819,11003)
#mtahn_and_ds = c(42819,11003,10098,10137,10097,10117,10757,10606,10978,11122,11143,11081,10223,10382,11171,10187,11038,10875,10054,11167,14204,26633,10787,14095,10136,10141,10532,10466,101794,172107)
mtahn_ds = read.csv("syn_count_by_group_outliers.csv")
mtahn_and_ds = neuprint_list2df(neuprint_fetch_custom(cypher=paste0("MATCH (a:Neuron) WHERE a.group IN [", paste(mtahn_ds$group[mtahn_ds$ahn == 'MtAHN'], collapse = ','), "] RETURN a.bodyId AS bodyid"), timeout=2000))
mtahn_and_ds = mtahn_and_ds$bodyid[mtahn_and_ds$bodyid %in% neuprint_connection_table(c(42819,11003), prepost = "POST", threshold = ds_syn_thres)$partner]
all_bodyids = mtahn_and_ds 

clio_manc_data = manc_neuprint_meta(mtahn_and_ds)
clio_manc_data$abbrv = han$class2abbrv(clio_manc_data$class)
clio_manc_data$temp_label = paste(clio_manc_data$abbrv, clio_manc_data$group)
clio_manc_data$temp_label[clio_manc_data$bodyid %in% c(42819,11003)] = "MtAHN"
clio_manc_data$temp_label[clio_manc_data$abbrv == "MN"] = paste(clio_manc_data$type[clio_manc_data$abbrv == "MN"], clio_manc_data$group[clio_manc_data$abbrv == "MN"])
clio_manc_data$temp_label[clio_manc_data$abbrv == "DN"] = ifelse(is.na(annotation$type4id(clio_manc_data$bodyid[clio_manc_data$abbrv == "DN"])),
	paste("DN", clio_manc_data$group[clio_manc_data$abbrv == "DN"]),
	paste(annotation$type4id(clio_manc_data$bodyid[clio_manc_data$abbrv == "DN"]), clio_manc_data$group[clio_manc_data$abbrv == "DN"]))
clio_manc_data$temp_label = sapply(clio_manc_data$temp_label, FUN = function(x) paste0(x, "(", sum(clio_manc_data$temp_label == x), ")"))

names(mtahn_and_ds) = clio_manc_data$temp_label[match(mtahn_and_ds, clio_manc_data$bodyid)]

adj = neuprint_get_adjacency_matrix(mtahn_and_ds, timeout = 2000)
adj = apply(adj, 2, function(x) tapply(x, names(mtahn_and_ds), sum, na.rm = TRUE))
adj = t(apply(t(adj), 2, function(x) tapply(x, names(mtahn_and_ds), sum, na.rm = TRUE)))
#adj_order = order(adj["MtAHN(2)",], decreasing = TRUE)
#adj_order = c(adj_order[length(adj_order)], adj_order[1:(length(adj_order)-1)])
#adj = adj[adj_order, adj_order]
#adj = adj[unique(clio_manc_data$temp_label), unique(clio_manc_data$temp_label)]

d=dist(adj,method="euclidean")
Rowv = rowMeans(adj, na.rm = T)
fit = hclust(d, method="ward.D2")
fit = reorder(as.dendrogram(fit),Rowv)
cluster = adj[labels(fit),labels(fit)]
heatmap.2(cluster, col = viridis(100),scale="none",density.info='none',xlab="Output",ylab="Input",
    dendrogram='none',Rowv=FALSE,Colv=FALSE,trace='none',sepcolor='black', cexRow = 0.22 + 1/log10(ncol(cluster)), cexCol = 0.22 + 1/log10(ncol(cluster)), margins = c(10, 10))
dev.off()

#svglite(filename = paste0("mtahn_adj_matrix_to_partners_by_boxplot_outlier_log10(x+1)_transform_", Sys.Date(),".svg"), height=ceiling((18*nrow(adj)+200)/72), width=ceiling((18*ncol(adj)+200)/72))
png(file=paste0("mtahn_adj_matrix_to_partners_by_boxplot_outlier_log10(x+1)_transform_", Sys.Date(),".png"), height=18*nrow(adj)+200, width=18*ncol(adj)+200)
heatmap.2(log10(adj+1), col = viridis(100),scale="none",density.info='none',xlab="Output",ylab="Input",
    dendrogram='none',Rowv=FALSE,Colv=FALSE,trace='none',sepcolor='black', cexRow = 0.22 + 1/log10(ncol(adj)), cexCol = 0.22 + 1/log10(ncol(adj)), margins = c(10, 10))
dev.off()
write.csv(adj, paste0("mtahn_adj_matrix_to_partners_by_boxplot_outlier_log10(x+1)_transform_", Sys.Date(),".csv"))

svglite(filename = paste0("mtahn_adj_matrix_to_partners_by_boxplot_outlier_log10(x+1)_transform_", Sys.Date(),".svg"), height=ceiling((18*nrow(adj)+200)/72), width=ceiling((18*ncol(adj)+200)/72))
heatmap.2(log10(adj+1), col = viridis(100),scale="none",density.info='none',xlab="Output",ylab="Input",
    dendrogram='none',Rowv=FALSE,Colv=FALSE,trace='none',sepcolor='black', cexRow = 0.22 + 1/log10(ncol(adj)), cexCol = 0.22 + 1/log10(ncol(adj)), margins = c(10, 10))
dev.off()

svglite(filename = paste0("mtahn_adj_matrix_to_partners_by_boxplot_outlier_no_transform_", Sys.Date(),".svg"), height=ceiling((18*nrow(adj)+200)/72), width=ceiling((18*ncol(adj)+200)/72))
heatmap.2(adj, col = viridis(100),scale="none",density.info='none',xlab="Output",ylab="Input",
    dendrogram='none',Rowv=FALSE,Colv=FALSE,trace='none',sepcolor='black', cexRow = 0.22 + 1/log10(ncol(adj)), cexCol = 0.22 + 1/log10(ncol(adj)), margins = c(10, 10))
dev.off()

#han$graph_from_bodyid(mtahn_and_ds, threshold = 100L, by_group = TRUE, save_name = "mtahn_ds_manc_graph")

#export in pajek format
mtahn_graph = graph_from_adjacency_matrix(adj, mode = "directed", weighted = TRUE, diag = FALSE)
mtahn_graph = delete_edges(mtahn_graph, edges = which(E(mtahn_graph)$weight < 100))
graph_lut = c(MsAHN = "#ffe099", MtAHN = "#ffff64", AN = "#9632fa", DN = "#5ce7e1", IN = "#ff6666", MN = "#1910c7", SN = "#db58d7", ASN = "#ffc9b4", EN = "#80c8ff", AEN = "#b9b4ff", ND = "#aaaaaa")
V(mtahn_graph)$color = graph_lut[match(clio_manc_data$abbrv[match(V(mtahn_graph)$name,clio_manc_data$temp_label)] , names(graph_lut))]
V(mtahn_graph)$label = igraph::V(mtahn_graph)$label = sapply(igraph::V(mtahn_graph)$name, FUN = function(x) paste0(stringi::stri_wrap(x, width = 12), collapse = "\n"))
write_graph(mtahn_graph, file = "mtahn_ds_pajek.net", format = "pajek")
write.table(data.frame(index = 1:length(V(mtahn_graph)), name = V(mtahn_graph)$name, ic = "ic", color = V(mtahn_graph)$color), file = "mtahn_ds_pajek_names.txt", sep = " ", row.names = FALSE)

#edge curve determination function. If nodes far, apply gentle curve of 0.2. If near (within 8 positions in ring), apply extreme curve of 1.
ring_graph_edge_curve <-function (graph)
{
    V(graph)$name = 1:length(V(graph))
	n_nodes = length(V(graph))
	graph = as_edgelist(graph)
    difference = as.integer(graph[,2]) - as.integer(graph[,1])
	difference = ifelse(difference > n_nodes/2, 
		-n_nodes + difference, 
		ifelse(difference < -n_nodes/2,
			n_nodes + difference,
			difference))
	curvature = (n_nodes/2 - abs(difference))*sign(difference)/(n_nodes/4)
	# curvature = ifelse(difference > 8, 0.2, ifelse(
		# difference < -8, -0.2, ifelse(
			# difference > 0, 1, ifelse(
				# difference < 0, -1, NA))))
	return(curvature)
}

layout_circle = layout_in_circle(mtahn_graph)
#add more node positions that are mock nodes for line thickness legend
layout_legend = matrix(0, nrow = 10, ncol = 2)
layout_legend[,2] = c(0.8, 0.8, 0.85, 0.85, 0.9, 0.9, 0.95, 0.95, 1,1)
layout_legend[,1] = c(0.8, 1)
layout_circle = rbind(layout_circle, layout_legend)
mtahn_graph_w_legend = add_vertices(mtahn_graph, nv = 10, attr = list(name=c(paste('legend', 1:10))))
mtahn_graph_w_legend = add_edges(mtahn_graph_w_legend, edges = length(V(mtahn_graph))+1:10, attr = list(weight = c(10000, 5000, 1000, 500, 100)))
E(mtahn_graph_w_legend)$curvature = c(ring_graph_edge_curve(mtahn_graph),rep(0,5))

svglite(paste0("mtahn_ds_graph_syn_thres_100_", Sys.Date(),".svg"), width = 10, height = 10)
print(plot.igraph(mtahn_graph, edge.arrow.size = 0.8, vertex.size = 10,
	edge.arrow.width = 1, edge.label = NA, edge.width = E(mtahn_graph)$weight/500, vertex.color = V(mtahn_graph)$color, edge.curved = ring_graph_edge_curve(mtahn_graph),
	margin = c(0,0,0,0), layout = layout_in_circle, vertex.label.family = "sans", vertex.label.cex = 0.65, edge.label.family = "sans", edge.label.cex = 0.6))
dev.off()

svglite(paste0("mtahn_ds_graph_syn_thres_100_w_legend", Sys.Date(),".svg"), width = 10, height = 10)
print(plot.igraph(mtahn_graph_w_legend, edge.arrow.size = 0.8, vertex.size = c(rep(10, length(V(mtahn_graph))), rep(0,10)),
	edge.arrow.width = 1, curved = TRUE, edge.label = NA, edge.width = E(mtahn_graph_w_legend)$weight/500, vertex.color = V(mtahn_graph_w_legend)$color, edge.curved = c(ring_graph_edge_curve(mtahn_graph),rep(0,5)),
	margin = c(0,0,0,0), layout = layout_circle, vertex.label.family = "sans", vertex.label.cex = 0.65, edge.label.family = "sans", edge.label.cex = 0.6))
dev.off()




roi_info = neuprint_get_roiInfo(all_bodyids)

roi_info$sum.pre = neuprint_get_meta(roi_info$bodyid)["pre"]
roi_info[is.na(roi_info)] = 0
roi_groups = list(NTct = c("NTct(UTct-T1)(L)","NTct(UTct-T1)(R)"), WTct = c("WTct(UTct-T2)(L)","WTct(UTct-T2)(R)"), HTct=c("HTct(UTct-T3)(L)","HTct(UTct-T3)(R)"), IntTct="IntTct", LTct="LTct",
	LegNp.T1=c("LegNp(T1)(L)", "LegNp(T1)(R)"), LegNp.T2=c("LegNp(T2)(L)", "LegNp(T2)(R)"), LegNp.T3=c("LegNp(T3)(L)", "LegNp(T3)(R)"),
	Ov=c("Ov(L)", "Ov(R)"), mVAC=c("mVAC(T1)(L)", "mVAC(T1)(R)", "mVAC(T2)(L)", "mVAC(T2)(R)", "mVAC(T3)(L)", "mVAC(T3)(R)"), ANm="ANm")
temp_roi_names = paste0(unlist(roi_groups), ".pre")
roi_info[temp_roi_names[!(temp_roi_names %in% colnames(roi_info))]] = 0
#if is MN, set roi presynaptic site count to zero; we don't believe that these are real
roi_info$abbrv = clio_manc_data$abbrv[match(roi_info$bodyid, clio_manc_data$bodyid)]
roi_info[paste0(names(roi_groups), ".pre")] = sapply(roi_groups, FUN = function(x) rowSums(roi_info[paste0(x, ".pre")])/roi_info$sum.pre)
roi_info[roi_info$abbrv == "MN", paste0(names(roi_groups), ".pre")] = 0
roi_info$group = clio_manc_data$group[match(roi_info$bodyid, clio_manc_data$bodyid)]
roi_info$temp_label = clio_manc_data$temp_label[match(roi_info$bodyid, clio_manc_data$bodyid)]

roi_info$sum.post = neuprint_get_meta(roi_info$bodyid)["post"]
Rowv = rowMeans(roi_info[c('sum.pre','sum.post')], na.rm = T)
temp_roi_names = paste0(unlist(roi_groups), ".post")
roi_info[temp_roi_names[!(temp_roi_names %in% colnames(roi_info))]] = 0
roi_info[paste0(names(roi_groups), ".post")] = sapply(roi_groups, FUN = function(x) rowSums(roi_info[paste0(x, ".post")])/roi_info$sum.post)


roi_info = aggregate(roi_info[c(paste0(names(roi_groups), ".pre"), paste0(names(roi_groups), ".post"))], by = list(roi_info$temp_label), FUN = mean)
colnames(roi_info)[1] =  "temp_label"

d=dist(roi_info[c(paste0(names(roi_groups), ".pre"), paste0(names(roi_groups), ".post"))],method="euclidean")
#Rowv = rowMeans(roi_info[names(roi_groups)], na.rm = T)
#Rowv = rowMeans(roi_info[c('sum.pre','sum.post')], na.rm = T)
fit = hclust(d, method="ward.D2")
fit = reorder(as.dendrogram(fit),Rowv)
cluster = roi_info[labels(fit),]
cluster_mat = as.matrix(cluster[c(paste0(names(roi_groups), ".pre"), paste0(names(roi_groups), ".post"))])
rownames(cluster_mat) = cluster$temp_label
cluster_mat = t(cluster_mat)
#cut dendrogram and make data frame for clusters
#cut_tree=cut(fit,h=3)
#cut_tree_clusters = vector()
#for (i in 1:length(cut_tree$lower)) {
#	cut_tree_clusters = c(cut_tree_clusters, rep(i, length(labels(cut_tree$lower[[i]]))))

pdf(paste0("mtahn_ds_roi_prepost_", Sys.Date(),".pdf"), height=ceiling((18*nrow(cluster_mat)+200)/72), width=ceiling((18*ncol(cluster_mat)+200)/72))
heatmap.2(cluster_mat, col = viridis(100),scale="none",density.info='none',xlab="Neuron group",ylab="ROI input or output", rowsep = 11, 
    dendrogram='none',Rowv=FALSE,Colv=FALSE,trace='none',sepcolor='red', cexRow = 0.22 + 1/log10(ncol(cluster_mat)), cexCol = 0.22 + 1/log10(ncol(cluster_mat)), margins = c(10, 10))
dev.off()