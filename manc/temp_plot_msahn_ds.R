library(neuprintr)
box::use(mancfuns/han)
box::use(mancfuns/manc/annotation)
library(gplots)
library(viridis)
library(malevnc)
library(igraph)
library(svglite)

ds_syn_thres = 3

type_1 = c(23510,22289,22472,101206,21307,23628,20774,18519,157981,20590,18048,14502,16970,19315,18370,14625,32328,15337,23485,18657,18517,20537,15586,36805,22521,21381,23462,24945,153911,23191,23638,13060,22220,24442,155196,20586,16638,15788,17055,15113,101682,17945)
#msahn_ds = c(157981,21880,20590,18048,14502,16970,19315,18370,10592,10892,10014,10088,162240,169914,10361,18309,159325,21713,21600,21430,22078,14017,10669,12330,12080,10733,14898,14820,14638,15208,19239,25868,166245,31108,163424,29662,156982,174427,11917,12262,13488,10749,12971,18928,20331,13228,13586,162114,11776,10054,11167,11998,15092,10402,12076,10546,10699,10740,10755,14521,22149,15555,12978,14174,13303,10601,10632,10724,10539,27477,27132,24562,25102,26453,11182,10764,12906,13206,10200,10433,10812,11682,12216,11681,12175,11915,15345,15270,15701,17039,16059,166277,11089,10833,11663,16559,15783,15098,10178,46457,10543,11294,158815,157206,152548,18188,12092,12130,14258,12797,156765,17814,25859,23022,12486,11628,27493,156144,29130,16158,103576,28296,171347,159171)

msahn = c(12536, 13926)

msahn_ds = read.csv("syn_count_by_group_outliers.csv")
all_bodyids = neuprint_list2df(neuprint_fetch_custom(cypher=paste0("MATCH (a:Neuron) WHERE a.group IN [", paste(msahn_ds$group[msahn_ds$ahn == 'MsAHN'], collapse = ','), "] RETURN a.bodyId AS bodyid"), timeout=2000))
all_bodyids = all_bodyids$bodyid[all_bodyids$bodyid %in% neuprint_connection_table(msahn, prepost = "POST", threshold = ds_syn_thres)$partner]
#all_bodyids = c(msahn, all_bodyids)

#all_bodyids = c(msahn, type_1, msahn_ds)
#all_bodyids = unique(all_bodyids)

clio_manc_data = manc_neuprint_meta(all_bodyids)
clio_manc_data$abbrv = han$class2abbrv(clio_manc_data$class)
#if no group, make bodyid group
clio_manc_data$temp_label = paste(clio_manc_data$abbrv, clio_manc_data$group)
clio_manc_data$temp_label[clio_manc_data$abbrv == "MN"] = paste(clio_manc_data$type[clio_manc_data$abbrv == "MN"], clio_manc_data$group[clio_manc_data$abbrv == "MN"])
clio_manc_data$temp_label[clio_manc_data$abbrv == "DN"] = ifelse(is.na(annotation$type4id(clio_manc_data$bodyid[clio_manc_data$abbrv == "DN"])),
	paste("DN", clio_manc_data$group[clio_manc_data$abbrv == "DN"]),
	paste(annotation$type4id(clio_manc_data$bodyid[clio_manc_data$abbrv == "DN"]), clio_manc_data$group[clio_manc_data$abbrv == "DN"]))
clio_manc_data$temp_label[clio_manc_data$bodyid %in% msahn] = "MsAHN"
clio_manc_data$temp_label[clio_manc_data$bodyid %in% type_1] = paste("Tect IN", clio_manc_data$group[clio_manc_data$bodyid %in% type_1])
clio_manc_data$temp_label = sapply(clio_manc_data$temp_label, FUN = function(x) paste0(x,"(",sum(clio_manc_data$temp_label == x),")"))

adj_all_vnc = neuprint_get_adjacency_matrix(bodyids = clio_manc_data$bodyid, chunksize = 1000L, timeout=2000)
adj_all_vnc_plot_group = clio_manc_data$temp_label[match(rownames(adj_all_vnc), clio_manc_data$bodyid)]

adj_all_vnc = t(apply(t(adj_all_vnc), 2, function(x) tapply(x, adj_all_vnc_plot_group, sum, na.rm = TRUE)))
adj_all_vnc = apply(adj_all_vnc, 2, function(x) tapply(x, adj_all_vnc_plot_group, sum, na.rm = TRUE))

order_by_class = order(clio_manc_data$abbrv[match(colnames(adj_all_vnc), clio_manc_data$temp_label)])
adj_all_vnc = adj_all_vnc[order_by_class, order_by_class]

svglite(paste0("MsAHN_ds_adj_matrix_by_boxplot_outlier_log10(x+1)_transform_",Sys.Date(),".svg"), height=(18*nrow(adj_all_vnc)+200)/72, width=(18*ncol(adj_all_vnc)+200)/72)
heatmap.2(log10(adj_all_vnc+1), col = viridis(100),scale="none",density.info='none',xlab="Output",ylab="Input",
    dendrogram='none',Rowv=FALSE,Colv=FALSE,trace='none',sepcolor='black', cexRow = 0.22 + 1/log10(ncol(adj_all_vnc)), cexCol = 0.22 + 1/log10(ncol(adj_all_vnc)), margins = c(10, 10))
dev.off()

png(paste0("MsAHN_ds_adj_matrix_by_boxplot_outlier_log10(x+1)_transform_",Sys.Date(),".png"), height=(18*nrow(adj_all_vnc)+200), width=(18*ncol(adj_all_vnc)+200))
heatmap.2(log10(adj_all_vnc+1), col = viridis(100),scale="none",density.info='none',xlab="Output",ylab="Input",
    dendrogram='none',Rowv=FALSE,Colv=FALSE,trace='none',sepcolor='black', cexRow = 0.22 + 1/log10(ncol(adj_all_vnc)), cexCol = 0.22 + 1/log10(ncol(adj_all_vnc)), margins = c(10, 10))
dev.off()

svglite(paste0("MsAHN_ds_adj_matrix_by_boxplot_outlier_no_transform_",Sys.Date(),".svg"), height=(18*nrow(adj_all_vnc)+200)/72, width=(18*ncol(adj_all_vnc)+200)/72)
heatmap.2(adj_all_vnc, col = viridis(100),scale="none",density.info='none',xlab="Output",ylab="Input",
    dendrogram='none',Rowv=FALSE,Colv=FALSE,trace='none',sepcolor='black', cexRow = 0.22 + 1/log10(ncol(adj_all_vnc)), cexCol = 0.22 + 1/log10(ncol(adj_all_vnc)), margins = c(10, 10))
dev.off()

#han$graph_from_bodyid(all_bodyids, threshold = 20L, by_group = TRUE, save_name = "msahn_ds_manc_graph")

#plot graph
msahn_graph = graph_from_adjacency_matrix(adj_all_vnc, mode = "directed", weighted = TRUE, diag = FALSE)
msahn_graph = delete_edges(msahn_graph, edges = which(E(msahn_graph)$weight < 100))
graph_lut = c(MsAHN = "#ffe099", MtAHN = "#ffff64", AN = "#9632fa", DN = "#5ce7e1", IN = "#ff6666", MN = "#1910c7", SN = "#db58d7", ASN = "#ffc9b4", EN = "#80c8ff", AEN = "#b9b4ff", ND = "#aaaaaa")
V(msahn_graph)$color = graph_lut[match(clio_manc_data$abbrv[match(V(msahn_graph)$name,clio_manc_data$temp_label)] , names(graph_lut))]
V(msahn_graph)$label = igraph::V(msahn_graph)$label = sapply(igraph::V(msahn_graph)$name, FUN = function(x) paste0(stringi::stri_wrap(x, width = 12), collapse = "\n"))
#write_graph(msahn_graph, file = "msahn_ds_pajek.net", format = "pajek")
#write.table(data.frame(index = 1:length(V(msahn_graph)), name = V(msahn_graph)$name, ic = "ic", color = V(msahn_graph)$color), file = "msahn_ds_pajek_names.txt", sep = " ", row.names = FALSE)

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

layout_circle = layout_in_circle(msahn_graph)
#add more node positions that are mock nodes for line thickness legend
layout_legend = matrix(0, nrow = 10, ncol = 2)
layout_legend[,2] = c(0.8, 0.8, 0.85, 0.85, 0.9, 0.9, 0.95, 0.95, 1,1)
layout_legend[,1] = c(0.8, 1)
layout_circle = rbind(layout_circle, layout_legend)
msahn_graph_w_legend = add_vertices(msahn_graph, nv = 10)
msahn_graph_w_legend = add_edges(msahn_graph_w_legend, edges = length(V(msahn_graph))+1:10, attr = list(weight = c(10000, 5000, 1000, 500, 100)))
svglite(paste0("msahn_ds_graph_syn_thres_100_", Sys.Date(),".svg"), width = 10, height = 10)
print(plot.igraph(msahn_graph, edge.arrow.size = 0.8, vertex.size = 10,
	edge.arrow.width = 1, curved = TRUE, edge.label = NA, edge.width = E(msahn_graph)$weight/500, vertex.color = V(msahn_graph)$color, edge.curved = ring_graph_edge_curve(msahn_graph),
	margin = c(0,0,0,0), layout = layout_in_circle, vertex.label.family = "sans", vertex.label.cex = 0.65, edge.label.family = "sans", edge.label.cex = 0.6))
dev.off()
svglite(paste0("msahn_ds_graph_syn_thres_100_w_legend", Sys.Date(),".svg"), width = 10, height = 10)
print(plot.igraph(msahn_graph_w_legend, edge.arrow.size = 0.8, vertex.size = c(rep(10, length(V(msahn_graph))), rep(0,10)),
	edge.arrow.width = 1, curved = TRUE, edge.label = NA, edge.width = E(msahn_graph_w_legend)$weight/500, vertex.color = V(msahn_graph_w_legend)$color, edge.curved = c(ring_graph_edge_curve(msahn_graph),rep(0,5)),
	margin = c(0,0,0,0), layout = layout_circle, vertex.label.family = "sans", vertex.label.cex = 0.65, edge.label.family = "sans", edge.label.cex = 0.6))
dev.off()

roi_info = neuprint_get_roiInfo(all_bodyids[!(all_bodyids %in% msahn)])

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
#}

pdf(paste0("msahn_ds_roi_prepost_", Sys.Date(),".pdf"), height=ceiling((18*nrow(cluster_mat)+200)/72), width=ceiling((18*ncol(cluster_mat)+200)/72))
heatmap.2(cluster_mat, col = viridis(100),scale="none",density.info='none',xlab="Neuron group",ylab="ROI input or output", rowsep = 11, 
    dendrogram='none',Rowv=FALSE,Colv=FALSE,trace='none',sepcolor='red', cexRow = 0.22 + 1/log10(ncol(cluster_mat)), cexCol = 0.22 + 1/log10(ncol(cluster_mat)), margins = c(10, 10))
dev.off()


