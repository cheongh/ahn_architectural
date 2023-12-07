#analysis of AHNs in MANC
library(malevnc)
library(neuprintr)
library(igraph)
library(ggplot2)
library(RColorBrewer)
library(svglite)
library(googlesheets4)
library(nat)
library(gplots)
library(viridis)

box::use(mancfuns/han)
source("shared_functions.R")

syn_thres_prepost = c(PRE = 10, POST = 3)
np_thres = 0.1

#dn_lut = c(DNg02 = "#3A9E79", DNp54 = "#D15E13", DNp08 = "#7670B1", DNg32 = "#DE2A88", DNp38 = "#6DA62D", "26730" = "#E1AB24", DNg29 = "#A0711C", "17526" = "#666666", "29987" = "#6eb2e6", "DN (summed)" = "#CBCBCB")
dn_lut = c(DNg02 = "#3A9E79", DNxn085 = "#D15E13", DNp08 = "#7670B1", DNg32 = "#DE2A88", DNxn074 = "#6DA62D", "26730" = "#E1AB24", DNg29 = "#A0711C", DNxn155 = "#666666", DNnt007 = "#6eb2e6", "DN (summed)" = "#CBCBCB")

msahn = c(12536, 13926) #12536--L, 13926--R
names(msahn) = c("MsAHN-R", "MsAHN-L")
mtahn = c(11003, 42819) #11003--L, 42819--R
names(mtahn) = c("MtAHN-R", "MtAHN-L")
#manc_dn_all = read.csv("MANC_DN_typing_8-23-21.csv")


#get partners of ahns
for (prepost in c("POST", "PRE")) {
	syn_thres = syn_thres_prepost[prepost]
	ahn_partner = neuprint_connection_table(c(msahn,mtahn), prepost = prepost, threshold = syn_thres)
	clio_manc_data = manc_neuprint_meta(unique(ahn_partner$partner))
	clio_manc_data$type[clio_manc_data$group %in% type_1] = "Tect IN"
	clio_manc_data$group[is.na(clio_manc_data$group)] = clio_manc_data$bodyid[is.na(clio_manc_data$group)]
	clio_manc_data$abbrv = han$class2abbrv(clio_manc_data$class)
	#get syn counts by status
	ahn_partner$status = clio_manc_data$status[match(ahn_partner$partner, clio_manc_data$bodyid)]
	ahn_partner$status[is.na(ahn_partner$status)] = "ND"
	ahn_partner$class = clio_manc_data[match(ahn_partner$partner, clio_manc_data$bodyid), "class"]
	
	ahn_partner_bodyid = data.frame(bodyid = as.character(unique(ahn_partner$partner)))
	ahn_partner_bodyid$type = as.vector(clio_manc_data[match(ahn_partner_bodyid[[1]], clio_manc_data$bodyid), "class"])
	ahn_partner_bodyid$type[is.na(ahn_partner_bodyid$type)] = "ND"
	#change type (this is actually class after changing terminology) in ahn downstream list to abbreviations
	ahn_partner_bodyid$type = han$manc_neuron_type_lookup[match(ahn_partner_bodyid$type, names(han$manc_neuron_type_lookup))]
	ahn_partner_bodyid = ahn_partner_bodyid[!(ahn_partner_bodyid$type %in% c("Glia", "ND", NA)),]

	#get adjacency matrix and set synapses below threshold to zero
	if (prepost == "POST") {
		ahn_partner_adj = neuprint_get_adjacency_matrix(inputids = c(msahn,mtahn), outputids = ahn_partner_bodyid$bodyid)
	} else {
		ahn_partner_adj = neuprint_get_adjacency_matrix(inputids = ahn_partner_bodyid$bodyid, outputids = c(msahn,mtahn))
	}
	ahn_partner_adj[ahn_partner_adj < syn_thres] = 0

	#classify cells on whether they have connectivity to L, R or both sides of each AHN type
	if (prepost == "POST") {
		ahn_partner_bodyid$type_conn = apply(ahn_partner_bodyid, MARGIN = 1, FUN = function(x) {
			ahn_above_thres = (ahn_partner_adj[as.character(c(msahn, mtahn)), x["bodyid"]] > 0)
			if (sum(ahn_above_thres) > 2 | all(ahn_above_thres == c(F,T,T,F)) | all(ahn_above_thres == c(T,F,F,T))) ahn_above_thres = c(T,T,T,T) #if has conn with 3 ahns, or if ahn conn is between left and right of different ahn pairs, lump into conn with all 4 ahns
			paste(x["type"], paste(c("MsR", "MsL", "MtR", "MtL")[ahn_above_thres], collapse = "_"))
			})
	} else {
		ahn_partner_bodyid$type_conn = apply(ahn_partner_bodyid, MARGIN = 1, FUN = function(x) {
			ahn_above_thres = (ahn_partner_adj[x["bodyid"], as.character(c(msahn, mtahn))] > 0)
			if (sum(ahn_above_thres) > 2 | all(ahn_above_thres == c(F,T,T,F)) | all(ahn_above_thres == c(T,F,F,T))) ahn_above_thres = c(T,T,T,T) #if has conn with 3 ahns, or if ahn conn is between left and right of different ahn pairs, lump into conn with all 4 ahns
			paste(x["type"], paste(c("MsR", "MsL", "MtR", "MtL")[ahn_above_thres], collapse = "_"))
			})
	}
	#ahn_partner_bodyid$vertex_label = apply(ahn_partner_bodyid, MARGIN = 1, 
	#	FUN = function(x) paste0(x[2], "(", sum(ahn_partner_bodyid$type_conn == x[3]), ")"))
	ahn_partner_bodyid$vertex_label = apply(ahn_partner_bodyid, MARGIN = 1, 
		FUN = function(x) sum(ahn_partner_bodyid$type_conn == x["type_conn"]))
	ahn_partner_bodyid$vertex_color = graph_lut[ahn_partner_bodyid$type]

	#sum ahn downstream by type and connectivity with AHNs
	ahn_partner_adj_sum = ahn_partner_adj
	if (prepost == "POST") {
		ahn_partner_adj_sum = t(apply(t(ahn_partner_adj_sum), 2, function(x) tapply(x, ahn_partner_bodyid$type_conn, sum, na.rm = TRUE)))
		ahn_partner_adj_sum = squarify(ahn_partner_adj_sum)
	} else {
		ahn_partner_adj_sum = apply(ahn_partner_adj_sum, 2, function(x) tapply(x, ahn_partner_bodyid$type_conn, sum, na.rm = TRUE))
		ahn_partner_adj_sum = squarify(ahn_partner_adj_sum)
	}

	ahn_partner_graph_sum = graph_from_adjacency_matrix(adjmatrix = ahn_partner_adj_sum, mode = "directed", weighted = TRUE, diag = FALSE)
	V(ahn_partner_graph_sum)$label = ahn_partner_bodyid$vertex_label[match(V(ahn_partner_graph_sum)$name, ahn_partner_bodyid$type_conn)]
	V(ahn_partner_graph_sum)$label[match(c(msahn, mtahn), V(ahn_partner_graph_sum)$name)] = names(c(msahn, mtahn))
	V(ahn_partner_graph_sum)$color = ahn_partner_bodyid$vertex_color[match(V(ahn_partner_graph_sum)$name, ahn_partner_bodyid$type_conn)]
	V(ahn_partner_graph_sum)$color[match(c(msahn, mtahn), V(ahn_partner_graph_sum)$name)] = graph_lut[c("MsAHN", "MsAHN", "MtAHN", "MtAHN")]

	tkplot_graph = tkplot(ahn_partner_graph_sum, canvas.width = 650, canvas.height = 650,
		edge.arrow.size = 0.8, edge.arrow.width = 1.5, curved = TRUE, edge.label = E(ahn_partner_graph_sum)$weight, edge.width = 2,
		vertex.label = V(ahn_partner_graph_sum)$name, vertex.color = V(ahn_partner_graph_sum)$color,
		margin = c(0,0,0,0), layout = layout_with_graphopt(ahn_partner_graph_sum, spring.constant = 5, charge = 0.01, spring.length = 50), vertex.label.family = "sans")
	readline(prompt="Adjust graph layout by dragging nodes around in tkplot window.\nPress ENTER once done with adjusting node positions: ")
	tkplot_layout = tkplot.getcoords(tkplot_graph)
	tkplot.close(tkplot_graph)
	#manually tweak layout and get coordinates
	#tkplot_layout = tkplot.getcoords(3)
	#after saving coords, go back and relabel vertices to proper names
	svglite(paste0("ahn_partner_MANC_", prepost, "_syn_thres_", syn_thres, "_", Sys.Date(),".svg"), width = 5, height = 5)
	print(plot.igraph(ahn_partner_graph_sum, edge.arrow.size = 0.4,
		edge.arrow.width = 1.5, curved = TRUE, edge.label = E(ahn_partner_graph_sum)$weight, edge.width = log1p(E(ahn_partner_graph_sum)$weight)/2, vertex.color = V(ahn_partner_graph_sum)$color,
		margin = c(0,0,0,0), layout = tkplot_layout, vertex.label.family = "sans", vertex.label.cex = 0.65, edge.label.family = "sans", edge.label.cex = 0.6))
	dev.off()
	 

	#calculate summary statistics for downstream partners of each AHN
	ahn_partner$abbrv = han$class2abbrv(ahn_partner$class)
	ahn_partner = ahn_partner[!(ahn_partner$status %in% c("Unimportant", "ND")) & !(ahn_partner$abbrv %in% c("Glia","ND")),]
	ahn_partner$ahn = names(c(msahn, mtahn))[match(ahn_partner$bodyid, c(msahn, mtahn))]
	ahn_partner_summary = aggregate(ahn_partner$weight, by = list(ahn_partner$ahn, ahn_partner$abbrv), FUN = sum)
	names(ahn_partner_summary) = c("AHN", "type", "synapse_weight")
	ahn_partner_summary$type = factor(ahn_partner_summary$type, levels = names(graph_lut)[names(graph_lut) %in% ahn_partner_summary$type])

	svglite(paste0("ahn_partner_MANC_synapse_frac_", prepost, "_", syn_thres, "_", Sys.Date(),".svg"), width = 6, height = 5)
	print(ggplot(ahn_partner_summary, aes(fill=type, y=synapse_weight, x=AHN)) + 
		geom_bar(position="fill", stat="identity") +
		xlab("") +
		ylab("% synapses") +
		scale_x_discrete() +
		theme(text = element_text(size=15), axis.line.y.left = element_line(color = "black"), axis.line.x.bottom = element_line(color = "black"), panel.background = element_blank()) +
		scale_fill_manual(values = graph_lut[names(graph_lut) %in% ahn_partner_summary$type]))
	dev.off()
	
	if (prepost == "POST") {
		#make histogram of synapse count by cell type
		#export top 5 cell types downstream of AHNs
		ahn_partner$ahn[ahn_partner$ahn %in% names(msahn)] = "MsAHN"
		ahn_partner$ahn[ahn_partner$ahn %in% names(mtahn)] = "MtAHN"
		ahn_partner$group = clio_manc_data$group[match(ahn_partner$partner, clio_manc_data$bodyid)]
		ahn_partner$group[is.na(ahn_partner$group)] = ahn_partner$partner[is.na(ahn_partner$group)]
		syn_count_by_group = aggregate(weight ~ group + ahn + abbrv, data = ahn_partner, FUN = sum)
		syn_count_by_group$abbrv = factor(syn_count_by_group$abbrv, levels = names(graph_lut)[names(graph_lut) %in% syn_count_by_group$abbrv])
		png(paste0("MANC_MsAHN_POST_syn_count_by_group_", Sys.Date(), ".png"), height=300,width=300)
		print(ggplot(syn_count_by_group[syn_count_by_group$ahn == "MsAHN",], aes(x=weight, fill=abbrv)) + 
			geom_histogram(binwidth = 1, position = "stack") +
			xlab("Synapse count") +
			ylab("Num groups") +
			scale_fill_manual(values=graph_lut[names(graph_lut) %in% syn_count_by_group$abbrv])) # +
			#xlim(0, 200))
		ggsave(file=paste0("MANC_MsAHN_POST_syn_count_by_group_", Sys.Date(), ".svg"))
		dev.off()
		png(paste0("MANC_MtAHN_POST_syn_count_by_group_", Sys.Date(), ".png"), height=300,width=300)
		print(ggplot(syn_count_by_group[syn_count_by_group$ahn == "MtAHN",], aes(x=weight, fill=abbrv)) + 
			geom_histogram(binwidth = 1, position = "stack") +
			xlab("Synapse count") +
			ylab("Num groups") +
			scale_fill_manual(values=graph_lut[names(graph_lut) %in% syn_count_by_group$abbrv]))# +
			#xlim(0, 200))
		ggsave(file=paste0("MANC_MtAHN_POST_syn_count_by_group_", Sys.Date(), ".svg"))
		dev.off()
		
		syn_count_by_group$outlier = FALSE
		syn_count_by_group$outlier[syn_count_by_group$ahn == "MsAHN"] = (syn_count_by_group$weight[syn_count_by_group$ahn == "MsAHN"] 
			> min(max(syn_count_by_group$weight[syn_count_by_group$ahn == "MsAHN"],na.rm=T), as.numeric(quantile(syn_count_by_group$weight[syn_count_by_group$ahn == "MsAHN"], 0.75)) + (IQR(syn_count_by_group$weight[syn_count_by_group$ahn == "MsAHN"])*1.5)))
		syn_count_by_group$outlier[syn_count_by_group$ahn == "MtAHN"] = (syn_count_by_group$weight[syn_count_by_group$ahn == "MtAHN"] 
			> min(max(syn_count_by_group$weight[syn_count_by_group$ahn == "MtAHN"],na.rm=T), as.numeric(quantile(syn_count_by_group$weight[syn_count_by_group$ahn == "MtAHN"], 0.75)) + (IQR(syn_count_by_group$weight[syn_count_by_group$ahn == "MtAHN"])*1.5)))
		png(paste0("MANC_MsAHN_MtAHN_POST_syn_count_by_group_boxplot_", Sys.Date(), ".png"), height=300,width=300)
		print(ggplot() + 
			geom_boxplot(syn_count_by_group, outlier.shape = NA, mapping = aes(x=ahn, y=weight)) +
			geom_point(data = syn_count_by_group[syn_count_by_group$outlier,], mapping = aes(x=ahn, y=weight, color=abbrv), position = position_jitter(w = 0.2)) +
			scale_color_manual(values=graph_lut[names(graph_lut) %in% syn_count_by_group$abbrv[syn_count_by_group$outlier]]) +
			theme(panel.background = element_blank(), axis.line = element_line(color = "darkgrey")) +
			xlab("AHN type") +
			ylab("Synapse weight per group"))
		ggsave(file=paste0("MANC_MsAHN_MtAHN_POST_syn_count_by_group_boxplot_", Sys.Date(), ".svg"))
		dev.off()
		write.csv(syn_count_by_group[syn_count_by_group$outlier,], "syn_count_by_group_outliers.csv")
		
		MANC.surf.scale = MANC.surf
		MANC.surf.scale$Vertices = MANC.surf.scale$Vertices*1000
		neuron_color_lut = brewer.pal(10, "Paired")
		all_adj_raw = neuprint_get_adjacency_matrix(inputids = c(msahn,mtahn), outputids = ahn_partner_bodyid$bodyid)
		rownames(all_adj_raw) = names(c(msahn,mtahn))[match(rownames(all_adj_raw), c(msahn,mtahn))]
		all_adj_raw[all_adj_raw<syn_thres]=0
		all_adj_raw = all_adj_raw[,colSums(all_adj_raw)>=0]
		for (i in c("MsAHN","MtAHN")) {
			top_ahn_neurons = order(syn_count_by_group$weight[syn_count_by_group$ahn==i], decreasing=TRUE)[1:5]
			top_ahn_neurons = syn_count_by_group$group[syn_count_by_group$ahn==i][top_ahn_neurons]
			top_ahn_neurons = clio_manc_data$bodyid[clio_manc_data$group %in% top_ahn_neurons]
			top_ahn_neurons = top_ahn_neurons[top_ahn_neurons %in% colnames(all_adj_raw)]
			test_mesh = read_manc_meshes(top_ahn_neurons)
			neuron_colors = clio_manc_data$group[match(top_ahn_neurons, clio_manc_data$bodyid)]
			temp_cell_types = unique(neuron_colors)
			for (j in 1:length(temp_cell_types)) neuron_colors[neuron_colors==temp_cell_types[j]] = colorRampPalette(neuron_color_lut[(j*2-1):(j*2)])(sum(neuron_colors==temp_cell_types[j]))
			
			open3d()
			userMatrix = matrix(c(0.009574261, -0.2694314, -0.96296543, 0, -0.984437406, -0.1714953, 0.03819529, 0, -0.175436258,  0.9476197, -0.26688156, 0, 0, 0, 0, 1), nrow = 4, ncol = 4, byrow = TRUE)
			par3d(windowRect=c(0,0,1920,800),zoom=0.45,userMatrix=userMatrix)
			for (j in 1:length(test_mesh)) wire3d(test_mesh[[j]],col = neuron_colors[j], lit = FALSE)
			wire3d(MANC.surf.scale, alpha = 0.04, col = "grey")
			rgl.snapshot(paste0(i, "_downstream_morphology_top_5_MANC_", Sys.Date(), ".png"),fmt="png")
			rgl.close()
		}
		#plot input output rois
		for (i in c("MsAHN","MtAHN")) {
			top_ahn_neurons = syn_count_by_group$group[syn_count_by_group$outlier & syn_count_by_group$ahn==i]
			top_ahn_neurons = clio_manc_data$bodyid[clio_manc_data$group %in% top_ahn_neurons]
			top_ahn_neurons = top_ahn_neurons[top_ahn_neurons %in% colnames(all_adj_raw)]
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

			pdf(file=paste0(i, "_ds_roi_post_sites_manc_", Sys.Date(),".pdf"), height = 450/72, width = 1200/72)
			temp = t(neuron_roi_us[rownames(temp_row_order),])
			temp = temp[order(match(rownames(temp), names(han$manc_roi_groups))),]
			print(heatmap.2(temp[rownames(temp) %in% rois_to_keep,], col = viridis(100),scale="none",density.info='none',xlab="ds neuron",ylab="Postsynaptic sites",
				dendrogram='none',Rowv=FALSE,Colv=FALSE,trace='none',sepcolor='black', margins = c(10, 10), breaks=seq(0,1,0.01)))
			dev.off()

			pdf(file=paste0(i, "_ds_roi_pre_sites_manc_", Sys.Date(),".pdf"), height = 450/72, width = 1200/72)
			temp = t(neuron_roi_ds[rownames(temp_row_order),])
			temp = temp[order(match(rownames(temp), names(han$manc_roi_groups))),]
			print(heatmap.2(temp[rownames(temp) %in% rois_to_keep,], col = viridis(100),scale="none",density.info='none',xlab="ds neuron",ylab="Presynaptic sites",
				dendrogram='none',Rowv=FALSE,Colv=FALSE,trace='none',sepcolor='black', margins = c(10, 10), breaks=seq(0,1,0.01)))
			dev.off()
		}
	}
}



#extract only DNs from AHN upstream list, and add identities
manc_dn_all = manc_neuprint_meta("class:descending")
ahn_upstream_dns = ahn_partner[ahn_partner$abbrv == "DN",]
ahn_upstream_dns$dn_identity = manc_dn_all$type[match(ahn_upstream_dns$partner, manc_dn_all$bodyid)]

ahn_upstream_dn_adj = neuprint_get_adjacency_matrix(inputids = unique(ahn_upstream_dns$partner), outputids = c(msahn,mtahn))
rownames(ahn_upstream_dn_adj) = ahn_upstream_dns$dn_identity[match(rownames(ahn_upstream_dn_adj), ahn_upstream_dns$partner)]
rownames(ahn_upstream_dn_adj) = sapply(rownames(ahn_upstream_dn_adj), FUN = function(x) paste0(x,"(",sum(rownames(ahn_upstream_dn_adj) == x), ")"))
colnames(ahn_upstream_dn_adj) = c(names(msahn), names(mtahn))
ahn_upstream_dn_adj[ahn_upstream_dn_adj < syn_thres] = 0

ahn_upstream_dn_adj = apply(ahn_upstream_dn_adj, 2, function(x) tapply(x, rownames(ahn_upstream_dn_adj), sum, na.rm = TRUE))
ahn_upstream_dn_adj = squarify(ahn_upstream_dn_adj)


ahn_upstream_dn_graph = graph_from_adjacency_matrix(adjmatrix = ahn_upstream_dn_adj, mode = "directed", weighted = TRUE, diag = FALSE)

tkplot_graph = tkplot(ahn_upstream_dn_graph, canvas.width = 600, canvas.height = 600,
	edge.arrow.size = 0.8, edge.arrow.width = 1.5, curved = TRUE, edge.label = E(ahn_upstream_dn_graph)$weight, edge.width = 2, vertex.color = "teal",
	margin = c(0,0,0,0), layout = layout_with_graphopt(ahn_upstream_dn_graph, spring.constant = 5, charge = 0.01, spring.length = 50), vertex.label.family = "sans")
readline(prompt="Adjust graph layout by dragging nodes around in tkplot window.\nPress ENTER once done with adjusting node positions: ")
tkplot_layout_upstream_dn = tkplot.getcoords(tkplot_graph)
tkplot.close(tkplot_graph)
#manually tweak layout and get coordinates
#tkplot_layout_upstream_dn = tkplot.getcoords(2)
#after saving coords, go back and relabel vertices to proper names
svglite(paste0("AHN_upstream_DNs_MANC_", Sys.Date(),".svg"), width = 5, height = 5)
print(plot.igraph(ahn_upstream_dn_graph, edge.arrow.size = 0.4,
	edge.arrow.width = 1.5, curved = TRUE, edge.label = E(ahn_upstream_dn_graph)$weight, edge.width = log1p(E(ahn_upstream_dn_graph)$weight)/2, vertex.color = graph_lut["DN"],
	margin = c(0,0,0,0), layout = tkplot_layout_upstream_dn, vertex.label.family = "sans", vertex.label.cex = 0.4, edge.label.family = "sans", edge.label.cex = 0.4))
dev.off()


#calculate summary statistics for upstream partners of each AHN
ahn_upstream_dn_summary = aggregate(ahn_upstream_dns$weight, by = list(ahn_upstream_dns$ahn, ahn_upstream_dns$dn_identity), FUN = sum)
names(ahn_upstream_dn_summary) = c("AHN", "dn_identity", "synapse_weight")
#made this group-to-group connectivity instead of thresholding by individual AHN connectivity
ahn_upstream_dn_summary$AHN_type = ifelse(ahn_upstream_dn_summary$AHN %in% c("MsAHN-L", "MsAHN-R"), "MsAHN", "MtAHN")
ahn_upstream_dn_summary$less_than_5_percent = apply(ahn_upstream_dn_summary, MARGIN = 1, FUN = function(x)
	ifelse(sum(ahn_upstream_dn_summary$synapse_weight[ahn_upstream_dn_summary$dn_identity == x["dn_identity"] & ahn_upstream_dn_summary$AHN_type == x["AHN_type"]]) < 0.05*sum(ahn_upstream_dn_summary$synapse_weight[ahn_upstream_dn_summary$AHN_type == x["AHN_type"]]),
	"DN (summed)", x["dn_identity"]
	))
#new line 3-21-22: if connectivity of a DN with either AHN is above 5%, show DN for both
ahn_upstream_dn_summary$less_than_5_percent = sapply(ahn_upstream_dn_summary$dn_identity, FUN = function(x) {
	temp_eval = ahn_upstream_dn_summary$less_than_5_percent[ahn_upstream_dn_summary$dn_identity == x];
	if (any(temp_eval != "DN (summed)")) return(x) else return("DN (summed)")
	})
ahn_upstream_dn_summary_sum = aggregate(ahn_upstream_dn_summary$synapse_weight, by = list(ahn_upstream_dn_summary$AHN, ahn_upstream_dn_summary$less_than_5_percent), FUN = sum)
colnames(ahn_upstream_dn_summary_sum) = c("AHN", "type", "synapse_weight")
ahn_upstream_dn_summary_sum$AHN_type = ifelse(ahn_upstream_dn_summary_sum$AHN %in% c("MsAHN-L", "MsAHN-R"), "MsAHN", "MtAHN")
dn_type_order = sapply(unique(ahn_upstream_dn_summary_sum$type), FUN = function(x) sum(ahn_upstream_dn_summary_sum$synapse_weight[ahn_upstream_dn_summary_sum$type == x]))
dn_type_order = sort(dn_type_order[names(dn_type_order) != "DN (summed)"], decreasing = TRUE)
ahn_upstream_dn_summary_sum$type = factor(ahn_upstream_dn_summary_sum$type, levels = c(names(dn_type_order), "DN (summed)"))

write.csv(ahn_upstream_dn_summary_sum, paste0("ahn_partner_MANC_upstream_syn_thres_", syn_thres, "_", Sys.Date(),".csv"))

#plot summed upstream DNs
svglite(paste0("AHN_upstream_MANC_DN_synapse_frac_", Sys.Date(),".svg"), width = 6, height = 5)
print(ggplot(ahn_upstream_dn_summary_sum, aes(fill=type, y=synapse_weight, x=AHN)) + 
    geom_bar(position="fill", stat="identity") +
	xlab("") +
	ylab("% synapses") +
	scale_x_discrete() +
	theme(text = element_text(size=15), axis.line.y.left = element_line(color = "black"), axis.line.x.bottom = element_line(color = "black"), panel.background = element_blank()) +
	scale_fill_manual(values = dn_lut[levels(ahn_upstream_dn_summary_sum$type)]))
dev.off()

#plot ROIs of DNs that provide >=5% input
roi_info = ahn_upstream_dn_summary_sum$type[ahn_upstream_dn_summary_sum$type != "DN (summed)"]
roi_info = unique(ahn_upstream_dns$partner[ahn_upstream_dns$dn_identity %in% roi_info])
temp_neuron_info = manc_neuprint_meta(roi_info)
roi_info = neuprint_get_roiInfo(roi_info)
roi_info[is.na(roi_info)] = 0
roi_info$pre = temp_neuron_info$pre[match(roi_info$bodyid, temp_neuron_info$bodyid)]
manc_roi_groups = list(NTct = c("NTct(UTct-T1)(L)","NTct(UTct-T1)(R)"), WTct = c("WTct(UTct-T2)(L)","WTct(UTct-T2)(R)"), HTct=c("HTct(UTct-T3)(L)","HTct(UTct-T3)(R)"), IntTct="IntTct", LTct="LTct", 
	LegNp.T1=c("LegNp(T1)(L)", "LegNp(T1)(R)"), LegNp.T2=c("LegNp(T2)(L)", "LegNp(T2)(R)"), LegNp.T3=c("LegNp(T3)(L)", "LegNp(T3)(R)"),
	Ov=c("Ov(L)", "Ov(R)"), ANm="ANm", mVAC=c("mVAC(T1)(L)", "mVAC(T1)(R)", "mVAC(T2)(L)", "mVAC(T2)(R)", "mVAC(T3)(L)", "mVAC(T3)(R)"))
roi_info[names(manc_roi_groups)] = sapply(manc_roi_groups, FUN = function(x) {
	temp_col_names = paste0(x, ".pre")[paste0(x, ".pre") %in% colnames(roi_info)]
	rowSums(roi_info[temp_col_names], na.rm = TRUE)/roi_info$pre
})
roi_info$type = ahn_upstream_dns$dn_identity[match(roi_info$bodyid, ahn_upstream_dns$partner)]
roi_info_mean = aggregate(cbind(NTct,WTct,HTct,IntTct,LTct,LegNp.T1,LegNp.T2,LegNp.T3,Ov,mVAC,ANm)~type, data = roi_info, FUN = mean)
simple_row_hierarchical_clustering <- function(mat, Rowv = NULL) {
	d=dist(mat,method="euclidean")
	if(is.null(Rowv)) Rowv = Matrix::rowSums(mat)
	fit = hclust(d, method="ward.D2")
	fit = reorder(as.dendrogram(fit),Rowv)
	return(mat[labels(fit),])
	#return(mat[fit$order,])
}
roi_info_mean_ms = roi_info_mean[roi_info_mean$type %in% ahn_upstream_dns$dn_identity[ahn_upstream_dns$ahn %in% c("MsAHN-L","MsAHN-R")],]
rownames(roi_info_mean_ms) = sapply(roi_info_mean_ms$type, FUN = function(x) paste0(x, "(", length(unique(ahn_upstream_dns$partner[ahn_upstream_dns$dn_identity==x])), ")"))
roi_info_mean_ms = as.matrix(roi_info_mean_ms[names(manc_roi_groups)])
roi_info_mean_ms = simple_row_hierarchical_clustering(roi_info_mean_ms)
roi_info_mean_ms = roi_info_mean_ms[nrow(roi_info_mean_ms):1,]

roi_info_mean_mt = roi_info_mean[roi_info_mean$type %in% ahn_upstream_dns$dn_identity[ahn_upstream_dns$ahn %in% c("MtAHN-L","MtAHN-R")],]
rownames(roi_info_mean_mt) = sapply(roi_info_mean_mt$type, FUN = function(x) paste0(x, "(", length(unique(ahn_upstream_dns$partner[ahn_upstream_dns$dn_identity==x])), ")"))
roi_info_mean_mt = as.matrix(roi_info_mean_mt[names(manc_roi_groups)])
roi_info_mean_mt = simple_row_hierarchical_clustering(roi_info_mean_mt)
#roi_info_mean_mt = roi_info_mean_mt[nrow(roi_info_mean_mt):1,]
roi_info_mean_both = cbind(t(roi_info_mean_ms), t(roi_info_mean_mt))
pdf(paste0("ahn_us_dn_roi_info_", Sys.Date(), ".pdf"), height = 400/72, width = 400/72)
print(heatmap.2(roi_info_mean_both, col = viridis(100),scale="none",density.info='none',xlab="DN",ylab="ROI presynaptic sites", colsep = nrow(roi_info_mean_ms), sepcolor = "white", #labRow = NA,
	dendrogram='none',Rowv=FALSE,Colv=FALSE,trace='none',margins = c(9, 9)))
dev.off()



#plot DNs
#msahn upstream dns
msahn_upstream_dn_unique = unique(ahn_upstream_dns$partner[ahn_upstream_dns$ahn == "MsAHN-L" | ahn_upstream_dns$ahn == "MsAHN-R"])
#get rid of "DN (summed)" category for now
dn_mesh_col = dn_lut[ahn_upstream_dns$dn_identity[match(msahn_upstream_dn_unique, ahn_upstream_dns$partner)]]
msahn_upstream_dn_unique = msahn_upstream_dn_unique[!is.na(dn_mesh_col)]
dn_mesh_col = dn_mesh_col[!is.na(dn_mesh_col)]
dn_mesh = read_manc_meshes(msahn_upstream_dn_unique, type = "merged")
#dn_mesh_col = dn_lut[ahn_upstream_dns$dn_identity[match(msahn_upstream_dn_unique, ahn_upstream_dns$partner)]]
#dn_mesh_col[is.na(dn_mesh_col)] = dn_lut["DN (summed)"]
#MANC.tissue.surf.scale = MANC.tissue.surf
#MANC.tissue.surf.scale$Vertices = MANC.tissue.surf.scale$Vertices*1000
MANC.surf.scale = MANC.surf
MANC.surf.scale$Vertices = MANC.surf.scale$Vertices*1000
open3d()
userMatrix = matrix(c(0.009574261, -0.2694314, -0.96296543, 0, -0.984437406, -0.1714953, 0.03819529, 0, -0.175436258,  0.9476197, -0.26688156, 0, 0, 0, 0, 1), nrow = 4, ncol = 4, byrow = TRUE)
par3d(windowRect=c(0,0,1920,800),zoom=0.45,userMatrix=userMatrix)
for (i in 1:length(dn_mesh)) wire3d(dn_mesh[[i]],col = dn_mesh_col[i], lit = FALSE)
#plot3d(MANC.tissue.surf.scale, alpha = 0.1)
wire3d(MANC.surf.scale, alpha = 0.04, col = "grey")
rgl.snapshot(paste0("MANC_MsAHN_upstream_dns_", Sys.Date(), ".png"),fmt="png")
rgl.close()
#mtahn upstream dns
mtahn_upstream_dn_unique = unique(ahn_upstream_dns$partner[ahn_upstream_dns$ahn == "MtAHN-L" | ahn_upstream_dns$ahn == "MtAHN-R"])
dn_mesh_col = dn_lut[ahn_upstream_dns$dn_identity[match(mtahn_upstream_dn_unique, ahn_upstream_dns$partner)]]
mtahn_upstream_dn_unique = mtahn_upstream_dn_unique[!is.na(dn_mesh_col)]
dn_mesh_col = dn_mesh_col[!is.na(dn_mesh_col)]
dn_mesh = read_manc_meshes(mtahn_upstream_dn_unique, type = "merged")
open3d()
par3d(windowRect=c(0,0,1920,800),zoom=0.45,userMatrix=userMatrix)
for (i in 1:length(dn_mesh)) wire3d(dn_mesh[[i]],col = dn_mesh_col[i], lit = FALSE)
#plot3d(MANC.tissue.surf.scale, alpha = 0.1)
wire3d(MANC.surf.scale, alpha = 0.04, col = "grey")
rgl.snapshot(paste0("MANC_MtAHN_upstream_dns_", Sys.Date(), ".png"),fmt="png")
rgl.close()

#make DN graph with summed categories
ahn_upstream_dn_adj_sum = neuprint_get_adjacency_matrix(inputids = unique(ahn_upstream_dns$partner), outputids = c(msahn,mtahn))
ahn_upstream_dn_thres = data.frame(bodyid = rownames(ahn_upstream_dn_adj_sum), type = ahn_upstream_dns$dn_identity[match(rownames(ahn_upstream_dn_adj_sum), ahn_upstream_dns$partner)])
rownames(ahn_upstream_dn_adj_sum) = ahn_upstream_dns$dn_identity[match(rownames(ahn_upstream_dn_adj_sum), ahn_upstream_dns$partner)]
colnames(ahn_upstream_dn_adj_sum) = c(names(msahn), names(mtahn))
ahn_upstream_dn_adj_sum[ahn_upstream_dn_adj_sum < syn_thres] = 0
##determine categories by connectivity with AHNs for entire cell type rather than individual cell
ahn_upstream_dn_thres$conn = sapply(rownames(ahn_upstream_dn_adj_sum), FUN = function(x)
	paste0(c("MsR", "MsL", "MtR", "MtL")[colSums(ahn_upstream_dn_adj_sum[rownames(ahn_upstream_dn_adj_sum) == x, c(names(msahn), names(mtahn)), drop = FALSE]) > 0], collapse = "_"))
ahn_upstream_dn_thres$type = ifelse(ahn_upstream_dn_thres$type %in% ahn_upstream_dn_summary$less_than_5_percent,
	ahn_upstream_dn_thres$type,
	"DN (summed)")
#sapply(rownames(ahn_upstream_dn_adj_sum), FUN = function(x) ifelse(ahn_upstream_dn_summary_sum$[ahn_upstream_dn_summary$dn_identity == x]
#make connectivity categories for determining which DNs to sum
ahn_upstream_dn_thres$type_conn = ifelse(ahn_upstream_dn_thres$type == "DN (summed)",
	paste0("DN_sum_", ahn_upstream_dn_thres$conn),
	rownames(ahn_upstream_dn_adj_sum))
ahn_upstream_dn_thres$vertex_label = sapply(ahn_upstream_dn_thres$type_conn, FUN = function(x) sum(ahn_upstream_dn_thres$type_conn == x))
ahn_upstream_dn_thres$vertex_label = paste0(ahn_upstream_dn_thres$type, "(", ahn_upstream_dn_thres$vertex_label, ")")
ahn_upstream_dn_thres$vertex_color = dn_lut[ahn_upstream_dn_thres$type]
rownames(ahn_upstream_dn_adj_sum) = ahn_upstream_dn_thres$type_conn


ahn_upstream_dn_adj_sum = apply(ahn_upstream_dn_adj_sum, 2, function(x) tapply(x, ahn_upstream_dn_thres$type_conn, sum, na.rm = TRUE))
ahn_upstream_dn_adj_sum = squarify(ahn_upstream_dn_adj_sum)

ahn_upstream_dn_graph_sum = graph_from_adjacency_matrix(adjmatrix = ahn_upstream_dn_adj_sum, mode = "directed", weighted = TRUE, diag = FALSE)
V(ahn_upstream_dn_graph_sum)$label = ahn_upstream_dn_thres$vertex_label[match(V(ahn_upstream_dn_graph_sum)$name, ahn_upstream_dn_thres$type_conn)]
V(ahn_upstream_dn_graph_sum)$label[match(names(c(msahn, mtahn)), V(ahn_upstream_dn_graph_sum)$name)] = c("R", "L", "R", "L")
V(ahn_upstream_dn_graph_sum)$color = ahn_upstream_dn_thres$vertex_color[match(V(ahn_upstream_dn_graph_sum)$name, ahn_upstream_dn_thres$type_conn)]
V(ahn_upstream_dn_graph_sum)$color[match(names(c(msahn, mtahn)), V(ahn_upstream_dn_graph_sum)$name)] = graph_lut[c("MsAHN", "MsAHN", "MtAHN", "MtAHN")]
tkplot_graph = tkplot(ahn_upstream_dn_graph_sum, canvas.width = 600, canvas.height = 600,
	edge.arrow.size = 0.8, edge.arrow.width = 1.5, curved = TRUE, edge.label = E(ahn_upstream_dn_graph_sum)$weight, edge.width = 2, vertex.color = V(ahn_upstream_dn_graph_sum)$color,
	margin = c(0,0,0,0), layout = layout_with_graphopt(ahn_upstream_dn_graph_sum, spring.constant = 5, charge = 0.01, spring.length = 50), vertex.label.family = "sans")
readline(prompt="Adjust graph layout by dragging nodes around in tkplot window.\nPress ENTER once done with adjusting node positions: ")
tkplot_layout_upstream_dn = tkplot.getcoords(tkplot_graph, norm = TRUE)
tkplot.close(tkplot_graph)
svglite(paste0("AHN_upstream_DNs_MANC_summed_", Sys.Date(),".svg"), width = 5, height = 5)
print(plot.igraph(ahn_upstream_dn_graph_sum, edge.arrow.size = 0.4,
	edge.arrow.width = 1.5, curved = TRUE, edge.label = E(ahn_upstream_dn_graph_sum)$weight, edge.width = log1p(E(ahn_upstream_dn_graph_sum)$weight)/2, vertex.color = V(ahn_upstream_dn_graph_sum)$color,
	margin = c(0,0,0,0), layout = tkplot_layout_upstream_dn, vertex.label.family = "sans", vertex.label.cex = 0.65, edge.label.family = "sans", edge.label.cex = 0.6))
dev.off()


#make the above graph but with separate DN sides
ahn_upstream_dn_thres$side = manc_dn_all$rootSide[match(ahn_upstream_dn_thres$bodyid, manc_dn_all$bodyid)]
ahn_upstream_dn_thres$side = gsub("HS", "", ahn_upstream_dn_thres$side)
ahn_upstream_dn_adj_sum = neuprint_get_adjacency_matrix(inputids = unique(ahn_upstream_dns$partner), outputids = c(msahn,mtahn))
ahn_upstream_dn_adj_sum[ahn_upstream_dn_adj_sum < syn_thres] = 0
##determine categories by connectivity with AHNs for entire cell type rather than individual cell
ahn_upstream_dn_thres$type = ifelse(ahn_upstream_dn_thres$type %in% ahn_upstream_dn_summary$less_than_5_percent,
	ahn_upstream_dn_thres$type,
	"DN (summed)")
#make connectivity categories for determining which DNs to sum
ahn_upstream_dn_thres$type_conn = ifelse(ahn_upstream_dn_thres$type == "DN (summed)",
	paste0("DN_sum_", ahn_upstream_dn_thres$conn),
	paste(ahn_upstream_dn_thres$type, ahn_upstream_dn_thres$side, sep="_"))
ahn_upstream_dn_thres$vertex_label = sapply(ahn_upstream_dn_thres$type_conn, FUN = function(x) sum(ahn_upstream_dn_thres$type_conn == x))
ahn_upstream_dn_thres$vertex_label = paste0(ahn_upstream_dn_thres$type_conn, "(", ahn_upstream_dn_thres$vertex_label, ")")
	
rownames(ahn_upstream_dn_adj_sum) = ahn_upstream_dn_thres$type_conn[match(rownames(ahn_upstream_dn_adj_sum), ahn_upstream_dn_thres$bodyid)]

ahn_upstream_dn_adj_sum = apply(ahn_upstream_dn_adj_sum, 2, function(x) tapply(x, ahn_upstream_dn_thres$type_conn, sum, na.rm = TRUE))
ahn_upstream_dn_adj_sum = squarify(ahn_upstream_dn_adj_sum)

ahn_upstream_dn_graph_sum = graph_from_adjacency_matrix(adjmatrix = ahn_upstream_dn_adj_sum, mode = "directed", weighted = TRUE, diag = FALSE)
V(ahn_upstream_dn_graph_sum)$label = ahn_upstream_dn_thres$vertex_label[match(V(ahn_upstream_dn_graph_sum)$name, ahn_upstream_dn_thres$type_conn)]
V(ahn_upstream_dn_graph_sum)$label[match(c(msahn, mtahn), V(ahn_upstream_dn_graph_sum)$name)] = c("R", "L", "R", "L")
V(ahn_upstream_dn_graph_sum)$color = ahn_upstream_dn_thres$vertex_color[match(V(ahn_upstream_dn_graph_sum)$name, ahn_upstream_dn_thres$type_conn)]
V(ahn_upstream_dn_graph_sum)$color[match(c(msahn, mtahn), V(ahn_upstream_dn_graph_sum)$name)] = graph_lut[c("MsAHN", "MsAHN", "MtAHN", "MtAHN")]
tkplot_graph = tkplot(ahn_upstream_dn_graph_sum, canvas.width = 600, canvas.height = 600,
	edge.arrow.size = 0.8, edge.arrow.width = 1.5, curved = TRUE, edge.label = E(ahn_upstream_dn_graph_sum)$weight, edge.width = 2, vertex.color = V(ahn_upstream_dn_graph_sum)$color,
	margin = c(0,0,0,0), layout = layout_with_graphopt(ahn_upstream_dn_graph_sum, spring.constant = 5, charge = 0.01, spring.length = 50), vertex.label.family = "sans")
readline(prompt="Adjust graph layout by dragging nodes around in tkplot window.\nPress ENTER once done with adjusting node positions: ")
tkplot_layout_upstream_dn = tkplot.getcoords(tkplot_graph, norm = TRUE)
tkplot.close(tkplot_graph)
svglite(paste0("AHN_upstream_DNs_MANC_summed_by_side_", Sys.Date(),".svg"), width = 5, height = 5)
print(plot.igraph(ahn_upstream_dn_graph_sum, edge.arrow.size = 0.4,
	edge.arrow.width = 1.5, curved = TRUE, edge.label = E(ahn_upstream_dn_graph_sum)$weight, edge.width = log1p(E(ahn_upstream_dn_graph_sum)$weight)/2, vertex.color = V(ahn_upstream_dn_graph_sum)$color,
	margin = c(0,0,0,0), layout = tkplot_layout_upstream_dn, vertex.label.family = "sans", vertex.label.cex = 0.65, edge.label.family = "sans", edge.label.cex = 0.6))
dev.off()

#calculate an ipsi-contra score for direct dn connections to ahns
temp_ic_mat = ahn_upstream_dn_adj_sum[!(rownames(ahn_upstream_dn_adj_sum) %in% c(msahn,mtahn)),as.character(c(msahn,mtahn))]
colnames(temp_ic_mat) = gsub("-" , "_", names(c(msahn,mtahn)))
temp_group = gsub("_([LRM]|ND)$", "", rownames(temp_ic_mat))
#check if dn groups have both L and R sides, if not exclude from ipsi-contra analysis
temp_row_groups_of_interest = sapply(temp_group, FUN = function(x) all(paste0(x, c("_L","_R")) %in% rownames(temp_ic_mat)))
temp_row_groups_of_interest = unique(names(temp_row_groups_of_interest)[temp_row_groups_of_interest])
#fails if any one group does not have both L and R sides
temp_ipsi_contra_mat = sapply(temp_row_groups_of_interest, FUN = function(x) sapply(c("MsAHN","MtAHN"), function(y)
		(temp_ic_mat[paste0(x,"_L"), paste0(y,"_L")] - temp_ic_mat[paste0(x,"_L"), paste0(y,"_R")]
		+ temp_ic_mat[paste0(x,"_R"), paste0(y,"_R")] - temp_ic_mat[paste0(x,"_R"), paste0(y,"_L")])/
		(temp_ic_mat[paste0(x,"_L"), paste0(y,"_L")] + temp_ic_mat[paste0(x,"_L"), paste0(y,"_R")]
		+ temp_ic_mat[paste0(x,"_R"), paste0(y,"_R")] + temp_ic_mat[paste0(x,"_R"), paste0(y,"_L")])
	))
temp_ipsi_contra_mat = as.data.frame(t(temp_ipsi_contra_mat))
temp_ipsi_contra_mat$dn_type = rownames(temp_ipsi_contra_mat)
temp_ipsi_contra_mat = pivot_longer(temp_ipsi_contra_mat, cols=c("MsAHN","MtAHN"), values_to = "ic_index", names_to="AHN")

for (i in c("MsAHN","MtAHN")) {
	temp_plot_dat = temp_ipsi_contra_mat[temp_ipsi_contra_mat$AHN == i & !is.nan(temp_ipsi_contra_mat$ic_index),]
	temp_plot_dat$dn_type = factor(temp_plot_dat$dn_type, levels = names(dn_lut)[names(dn_lut) %in% temp_plot_dat$dn_type])
	pdf(paste0("MANC_", i, "_us_DN_ipsi_contra_index", Sys.Date(), ".pdf"), width = 300/72, height = 200/72)
	print(ggplot(data=temp_plot_dat, aes(x=dn_type, y=ic_index, fill=dn_type)) +
	  geom_bar(stat="identity") +
	  scale_fill_manual(values = dn_lut[levels(temp_plot_dat$dn_type)]) +
	  theme_classic() +
	  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
	  ylim(-0.25,0.25) +
	  geom_hline(yintercept = 0))
	dev.off()
}