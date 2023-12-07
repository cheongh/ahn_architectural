library(igraph)
library(ggplot2)
library(RColorBrewer)
library(svglite)
library(googlesheets4)
library(nat)
library(fancr)

syn_thres = 3

#dn_lut = c(brewer.pal(n = 8, name = "Dark2"), brewer.pal(n = length(levels(ahn_downstream_dn_summary_sum$type))-8, name = "Accent"))
#dn_lut = c(DNg02 = "#3A9E79", DNp54 = "#D15E13", DNp08 = "#7670B1", DNg32 = "#DE2A88", DNp38 = "#6DA62D", "26730" = "#E1AB24", DNg29 = "#A0711C", "17526" = "#666666", MANC_11286 = "#daa520", "DN (summed)" = "#CBCBCB")
dn_lut = c(DNg02 = "#3A9E79", DNp54 = "#D15E13", DNp08 = "#7670B1", DNg32 = "#DE2A88", DNp38 = "#6DA62D", "26730" = "#E1AB24", DNg29 = "#A0711C", "MANC_17526" = "#666666", "29987" = "#6eb2e6", "DN (summed)" = "#CBCBCB")

#color lookup table for vertices
graph_lut = c(MsAHN = "#ffe099", MtAHN = "#ffff64", AN = "#9632fa", DN = "#5ce7e1", IN = "#ff6666", MN = "#1910c7", SN = "#db58d7", ASN = "#ffc9b4", EN = "#80c8ff", AEN = "#b9b4ff", ND = "#aaaaaa")
	

#mtahn & msahn names in adj matrix, written like this to reuse code from manc
mtahn = c("MtAHN-L" = "MtAHN-L", "MtAHN-R" = "MtAHN-R")
msahn = c("MsAHN-L" = "MsAHN-L", "MsAHN-R" = "MsAHN-R")

google_write_target_ss = "https://docs.google.com/spreadsheets/d/11Xb6NJxcshw_Xf3sIUovP8NC62JXq0hrjO6EgQeCrIg/edit#gid=1108943567"
ss_tabs_to_read = c("MsAHN_L_summary", "MsAHN_R_summary", "MtAHN_L_summary", "MtAHN_R_summary")
names(ss_tabs_to_read) = c("MsAHN-L", "MsAHN-R", "MtAHN-L", "MtAHN-R")
fanc_dn_type_ss = "https://docs.google.com/spreadsheets/d/1k0RZwgouAdXKO2lR7R1YuDgDgNRxRCG-pDKCNZrUBAk/edit#gid=2091967020"

ahn_downstream_table_list = list()
all_partners = data.frame()
for (i in 1:length(ss_tabs_to_read)) {
	ahn_downstream_table_list[[i]] = read_sheet(google_write_target_ss, sheet = ss_tabs_to_read[i])
	ahn_downstream_table_list[[i]]$autoid = unlist(ahn_downstream_table_list[[i]]$autoid)
	ahn_downstream_table_list[[i]] = ahn_downstream_table_list[[i]][ahn_downstream_table_list[[i]]$autoid != "0",]
	#apply syn_thres
	ahn_downstream_table_list[[i]] = ahn_downstream_table_list[[i]][ahn_downstream_table_list[[i]]$syn_count >= syn_thres,]
	#recode PN & LN as IN
	ahn_downstream_table_list[[i]]$matched_celltype[ahn_downstream_table_list[[i]]$matched_celltype == "PN"] = "IN"
	ahn_downstream_table_list[[i]]$matched_celltype[ahn_downstream_table_list[[i]]$matched_celltype == "LN"] = "IN"
	#get rid of ND
	ahn_downstream_table_list[[i]] = ahn_downstream_table_list[[i]][!(ahn_downstream_table_list[[i]]$matched_celltype %in% c("ND")),]
	#don't count self-synapses
	ahn_downstream_table_list[[i]] = ahn_downstream_table_list[[i]][!(ahn_downstream_table_list[[i]]$matched_celltype %in% c("MsAHN-l", "MsAHN-r", "MtAHN-l", "MtAHN-r")[i]),]
	all_partners = rbind(all_partners, data.frame(autoid = ahn_downstream_table_list[[i]]$autoid,
		abbrv = ahn_downstream_table_list[[i]]$matched_celltype,
		cell_type = ahn_downstream_table_list[[i]]$dn_type,
		synapse_weight = ahn_downstream_table_list[[i]]$syn_count,
		ahn = names(ss_tabs_to_read)[i]))
	
}

#construct adj matrix
ahn_downstream_bodyid = unique(all_partners$autoid[!(all_partners$abbrv %in% c("MsAHN-r", "MsAHN-l", "MtAHN-r","MtAHN-l"))])
ahn_downstream_adj = matrix(data = 0, nrow = length(c(names(ss_tabs_to_read), ahn_downstream_bodyid)),
	ncol = length(c(names(ss_tabs_to_read), ahn_downstream_bodyid)),
	dimnames = list(c(names(ss_tabs_to_read), ahn_downstream_bodyid), c(names(ss_tabs_to_read), ahn_downstream_bodyid)))
for (i in 1:length(ss_tabs_to_read)) {
	ahn_downstream_adj[names(ss_tabs_to_read)[i], ahn_downstream_table_list[[i]]$autoid[!(ahn_downstream_table_list[[i]]$matched_celltype %in% c("MsAHN-r", "MsAHN-l", "MtAHN-r","MtAHN-l"))]
		] = ahn_downstream_table_list[[i]]$syn_count[!(ahn_downstream_table_list[[i]]$matched_celltype %in% c("MsAHN-r", "MsAHN-l", "MtAHN-r","MtAHN-l"))]
}


#make group names based on type and connectivity

ahn_downstream_bodyid = data.frame(bodyid = ahn_downstream_bodyid)
ahn_downstream_bodyid$abbrv = all_partners$abbrv[match(ahn_downstream_bodyid$bodyid, all_partners$autoid)]
ahn_downstream_bodyid$type_conn = apply(ahn_downstream_bodyid, MARGIN = 1, 
	FUN = function(x) paste(x["abbrv"], paste(names(c(msahn, mtahn))[ahn_downstream_adj[names(c(msahn, mtahn)), x["bodyid"]] > 0], collapse = "_"), sep = "_"))
ahn_downstream_bodyid$vertex_label = apply(ahn_downstream_bodyid, MARGIN = 1, 
	FUN = function(x) sum(ahn_downstream_bodyid$type_conn == x["type_conn"]))
ahn_downstream_bodyid$vertex_color = graph_lut[ahn_downstream_bodyid$abbrv]

#sum ahn upstream by type and connectivity with AHNs
ahn_downstream_adj_sum = ahn_downstream_adj
ahn_downstream_adj_sum = t(apply(t(ahn_downstream_adj_sum), 2, function(x) tapply(x, c(msahn, mtahn, ahn_downstream_bodyid$type_conn), sum, na.rm = TRUE)))
ahn_downstream_adj_sum = apply(ahn_downstream_adj_sum, 2, function(x) tapply(x, c(msahn, mtahn, ahn_downstream_bodyid$type_conn), sum, na.rm = TRUE))

ahn_downstream_graph_sum = graph_from_adjacency_matrix(adjmatrix = ahn_downstream_adj_sum, mode = "directed", weighted = TRUE, diag = FALSE)
V(ahn_downstream_graph_sum)$label = ahn_downstream_bodyid$vertex_label[match(V(ahn_downstream_graph_sum)$name, ahn_downstream_bodyid$type_conn)]
V(ahn_downstream_graph_sum)$label[match(c(msahn, mtahn), V(ahn_downstream_graph_sum)$name)] = c("L", "R", "L", "R")
V(ahn_downstream_graph_sum)$color = ahn_downstream_bodyid$vertex_color[match(V(ahn_downstream_graph_sum)$name, ahn_downstream_bodyid$type_conn)]
V(ahn_downstream_graph_sum)$color[match(c(msahn, mtahn), V(ahn_downstream_graph_sum)$name)] = graph_lut[c("MsAHN", "MsAHN", "MtAHN", "MtAHN")]

tkplot_graph = tkplot(ahn_downstream_graph_sum, canvas.width = 600, canvas.height = 600,
	edge.arrow.size = 0.8, edge.arrow.width = 1.5, curved = TRUE, edge.label = E(ahn_downstream_graph_sum)$weight, edge.width = 2,
	vertex.color = V(ahn_downstream_graph_sum)$color, vertex.label = V(ahn_downstream_graph_sum)$name,
	margin = c(0,0,0,0), layout = layout_with_graphopt(ahn_downstream_graph_sum, spring.constant = 5, charge = 0.01, spring.length = 50), vertex.label.family = "sans")
readline(prompt="Adjust graph layout by dragging nodes around in tkplot window.\nPress ENTER once done with adjusting node positions: ")
tkplot_layout_upstream = tkplot.getcoords(tkplot_graph)
tkplot.close(tkplot_graph)
#manually tweak layout and get coordinates
#tkplot_layout_upstream = tkplot.getcoords(3)
#after saving coords, go back and relabel vertices to proper names
svglite(paste0("ahn_downstream_FANC_", Sys.Date(),".svg"), width = 5, height = 5)
print(plot.igraph(ahn_downstream_graph_sum, edge.arrow.size = 0.4,
	edge.arrow.width = 1.5, curved = TRUE, edge.label = E(ahn_downstream_graph_sum)$weight, edge.width = log1p(E(ahn_downstream_graph_sum)$weight)/2, vertex.color = V(ahn_downstream_graph_sum)$color,
	margin = c(0,0,0,0), layout = tkplot_layout_upstream, vertex.label.family = "sans", vertex.label.cex = 0.65, edge.label.family = "sans", edge.label.cex = 0.6))
dev.off()
 


#calculate summary statistics for upstream partners of each AHN
all_partners$abbrv[all_partners$abbrv %in% c("MsAHN-l", "MsAHN-r")] = "MsAHN"
all_partners$abbrv[all_partners$abbrv %in% c("MtAHN-l", "MtAHN-r")] = "MtAHN"
ahn_downstream_summary = aggregate(all_partners$synapse_weight, by = list(all_partners$ahn, all_partners$abbrv), FUN = sum)
names(ahn_downstream_summary) = c("AHN", "type", "synapse_weight")
ahn_downstream_summary$type = factor(ahn_downstream_summary$type, levels = names(graph_lut)[names(graph_lut) %in% ahn_downstream_summary$type])


svglite(paste0("ahn_downstream_FANC_synapse_frac_", Sys.Date(),".svg"), width = 6, height = 5)
print(ggplot(ahn_downstream_summary, aes(fill=type, y=synapse_weight, x=AHN)) + 
    geom_bar(position="fill", stat="identity") +
	xlab("") +
	ylab("% synapses") +
	scale_x_discrete() +
	theme(text = element_text(size=15), axis.line.y.left = element_line(color = "black"), axis.line.x.bottom = element_line(color = "black"), panel.background = element_blank()) +
	scale_fill_manual(values = graph_lut[names(graph_lut) %in% ahn_downstream_summary$type]))
dev.off()

write.csv(ahn_downstream_summary, paste0("ahn_partner_FANC_downstream_syn_thres_", syn_thres, "_", Sys.Date(),".csv"))


#extract only DNs from AHN upstream list, and add identities
ahn_downstream_dns = all_partners[all_partners$abbrv == "DN",]
ahn_downstream_dns$cell_type[is.na(ahn_downstream_dns$cell_type)] = "ND"
ahn_downstream_dns$vertex_label = sapply(ahn_downstream_dns$cell_type, FUN = function(x) length(unique(ahn_downstream_dns$autoid[ahn_downstream_dns$cell_type == x])))

#construct adj matrix
ahn_downstream_dn_bodyid = unique(ahn_downstream_dns$autoid)
ahn_downstream_dn_adj = matrix(data = 0, nrow = length(c(names(ss_tabs_to_read), ahn_downstream_dn_bodyid)),
	ncol = length(c(names(ss_tabs_to_read), ahn_downstream_dn_bodyid)),
	dimnames = list(c(names(ss_tabs_to_read), ahn_downstream_dn_bodyid), c(names(ss_tabs_to_read), ahn_downstream_dn_bodyid)))

for (i in 1:length(ss_tabs_to_read)) {
	ahn_downstream_dn_adj[ahn_downstream_dns$autoid[ahn_downstream_dns$ahn == names(ss_tabs_to_read)[i]], names(ss_tabs_to_read)[i]] = ahn_downstream_dns$synapse_weight[ahn_downstream_dns$ahn == names(ss_tabs_to_read)[i]]
}

rownames(ahn_downstream_dn_adj) = ifelse(is.na(ahn_downstream_dns$cell_type[match(rownames(ahn_downstream_dn_adj), ahn_downstream_dns$autoid)]),
	rownames(ahn_downstream_dn_adj),
	ahn_downstream_dns$cell_type[match(rownames(ahn_downstream_dn_adj), ahn_downstream_dns$autoid)])
colnames(ahn_downstream_dn_adj) = rownames(ahn_downstream_dn_adj)

ahn_downstream_dn_adj = apply(ahn_downstream_dn_adj, 2, function(x) tapply(x, rownames(ahn_downstream_dn_adj), sum, na.rm = TRUE))
ahn_downstream_dn_adj = t(apply(t(ahn_downstream_dn_adj), 2, function(x) tapply(x, colnames(ahn_downstream_dn_adj), sum, na.rm = TRUE)))

ahn_downstream_dn_graph = graph_from_adjacency_matrix(adjmatrix = ahn_downstream_dn_adj, mode = "directed", weighted = TRUE, diag = FALSE)
V(ahn_downstream_dn_graph)$label = ahn_downstream_dns$vertex_label[match(V(ahn_downstream_dn_graph)$name, ahn_downstream_dns$cell_type)]
V(ahn_downstream_dn_graph)$label = paste0(V(ahn_downstream_dn_graph)$name, "(", V(ahn_downstream_dn_graph)$label, ")")
V(ahn_downstream_dn_graph)$label[match(c(msahn, mtahn), V(ahn_downstream_dn_graph)$name)] = c("L", "R", "L", "R")
V(ahn_downstream_dn_graph)$color = graph_lut["DN"]
V(ahn_downstream_dn_graph)$color[match(c(msahn, mtahn), V(ahn_downstream_dn_graph)$name)] = graph_lut[c("MsAHN", "MsAHN", "MtAHN", "MtAHN")]

tkplot_graph = tkplot(ahn_downstream_dn_graph, canvas.width = 600, canvas.height = 600,
	edge.arrow.size = 0.8, edge.arrow.width = 1.5, curved = TRUE, edge.label = E(ahn_downstream_dn_graph)$weight, edge.width = 2, vertex.color = V(ahn_downstream_dn_graph)$color,
	margin = c(0,0,0,0), layout = layout_with_graphopt(ahn_downstream_dn_graph, spring.constant = 5, charge = 0.01, spring.length = 50), vertex.label.family = "sans")
readline(prompt="Adjust graph layout by dragging nodes around in tkplot window.\nPress ENTER once done with adjusting node positions: ")
tkplot_layout_upstream_dn = tkplot.getcoords(tkplot_graph)
tkplot.close(tkplot_graph)
#manually tweak layout and get coordinates
#tkplot_layout_upstream_dn = tkplot.getcoords(2)
#after saving coords, go back and relabel vertices to proper names
svglite(paste0("ahn_downstream_DNs_FANC_", Sys.Date(),".svg"), width = 5, height = 5)
print(plot.igraph(ahn_downstream_dn_graph, edge.arrow.size = 0.4,
	edge.arrow.width = 1.5, curved = TRUE, edge.label = E(ahn_downstream_dn_graph)$weight, edge.width = log1p(E(ahn_downstream_dn_graph)$weight)/2, vertex.color = V(ahn_downstream_dn_graph)$color,
	margin = c(0,0,0,0), layout = tkplot_layout_upstream_dn, vertex.label.family = "sans", vertex.label.cex = 0.4, edge.label.family = "sans", edge.label.cex = 0.4))
dev.off()


#calculate summary statistics for individual upstream DN types of each AHN
ahn_downstream_dn_summary = aggregate(ahn_downstream_dns$synapse_weight, by = list(ahn_downstream_dns$ahn, ahn_downstream_dns$cell_type), FUN = sum)
names(ahn_downstream_dn_summary) = c("AHN", "dn_identity", "synapse_weight")
#made this group-to-group connectivity instead of thresholding by individual AHN connectivity
ahn_downstream_dn_summary$AHN_type = ifelse(ahn_downstream_dn_summary$AHN %in% c("MsAHN-L", "MsAHN-R"), "MsAHN", "MtAHN")
ahn_downstream_dn_summary$less_than_5_percent = apply(ahn_downstream_dn_summary, MARGIN = 1, FUN = function(x)
	ifelse(sum(ahn_downstream_dn_summary$synapse_weight[ahn_downstream_dn_summary$dn_identity == x["dn_identity"] & ahn_downstream_dn_summary$AHN_type == x["AHN_type"]]) < 0.05*sum(ahn_downstream_dn_summary$synapse_weight[ahn_downstream_dn_summary$AHN_type == x["AHN_type"]]),
	"DN (summed)", x["dn_identity"]
	))
#new line 3-21-22: if connectivity of a DN with either AHN is above 5%, show DN for both
ahn_downstream_dn_summary$less_than_5_percent = sapply(ahn_downstream_dn_summary$dn_identity, FUN = function(x) {
	temp_eval = ahn_downstream_dn_summary$less_than_5_percent[ahn_downstream_dn_summary$dn_identity == x];
	if (any(temp_eval != "DN (summed)")) return(x) else return("DN (summed)")
	})
#temporarily set unidentified_DN to DN(summed)--remove this line once all DNs are properly identified
#ahn_downstream_dn_summary$less_than_5_percent[ahn_downstream_dn_summary$less_than_5_percent == "unidentified_DN"] = "DN (summed)"
ahn_downstream_dn_summary_sum = aggregate(ahn_downstream_dn_summary$synapse_weight, by = list(ahn_downstream_dn_summary$AHN, ahn_downstream_dn_summary$less_than_5_percent), FUN = sum)
colnames(ahn_downstream_dn_summary_sum) = c("AHN", "type", "synapse_weight")
ahn_downstream_dn_summary_sum$AHN_type = ifelse(ahn_downstream_dn_summary_sum$AHN %in% c("MsAHN-L", "MsAHN-R"), "MsAHN", "MtAHN")
dn_type_order = sapply(unique(ahn_downstream_dn_summary_sum$type), FUN = function(x) sum(ahn_downstream_dn_summary_sum$synapse_weight[ahn_downstream_dn_summary_sum$type == x]))
dn_type_order = sort(dn_type_order[names(dn_type_order) != "DN (summed)"], decreasing = TRUE)
ahn_downstream_dn_summary_sum$type = factor(ahn_downstream_dn_summary_sum$type, levels = c(names(dn_type_order), "DN (summed)"))

#plot summed upstream DNs
svglite(paste0("ahn_downstream_FANC_DN_synapse_frac_", Sys.Date(),".svg"), width = 6, height = 5)
print(ggplot(ahn_downstream_dn_summary_sum, aes(fill=type, y=synapse_weight, x=AHN)) + 
    geom_bar(position="fill", stat="identity") +
	xlab("") +
	ylab("% synapses") +
	scale_x_discrete() +
	theme(text = element_text(size=15), axis.line.y.left = element_line(color = "black"), axis.line.x.bottom = element_line(color = "black"), panel.background = element_blank()) +
	scale_fill_manual(values = dn_lut[names(dn_lut) %in% ahn_downstream_dn_summary_sum$type]))
dev.off()

#plot DNs
for (i in list(msahn, mtahn)) {
	ahn_downstream_dn_unique = unique(ahn_downstream_dns$autoid[(ahn_downstream_dns$ahn %in% i)
		& ahn_downstream_dns$cell_type %in% ahn_downstream_dn_summary$less_than_5_percent[ahn_downstream_dn_summary$AHN %in% i]])
	#get rid of "DN (summed)" category for now
	dn_mesh_col = dn_lut[ahn_downstream_dns$cell_type[match(ahn_downstream_dn_unique, ahn_downstream_dns$autoid)]]
	neuron_meshes = read_fanc_meshes(ahn_downstream_dn_unique)

	FANC.hxsurf = as.hxsurf(FANC.surf)

	#open3d(antialias=4)
	open3d()
	userMatrix = matrix(c(0.1058311, 0.99110168, -0.08062597, 0, -0.9838080, 0.09257191, -0.15341219, 0, -0.1445848, 0.09555651,  0.98485911, 0, 0,0,0,1), nrow = 4, ncol = 4, byrow = TRUE)
	par3d(windowRect=c(0,0,1920,800),zoom=0.5,userMatrix=userMatrix)
	plot3d(FANC.hxsurf,alpha=0.1)
	for (j in 1:length(neuron_meshes)) wire3d(neuron_meshes[[j]], col = dn_mesh_col[j], lit = FALSE)
	snapshot3d(paste0(i[1], "_dns_fanc_", Sys.Date(), ".png"),fmt="png",webshot = FALSE)
	#rgl.postscript("msahn_dns_fanc.svg",fmt="svg")
	rgl.close()
}






#make DN graph with summed categories
ahn_downstream_dn_adj = matrix(data = 0, nrow = length(c(names(ss_tabs_to_read), ahn_downstream_dn_bodyid)),
	ncol = length(c(names(ss_tabs_to_read), ahn_downstream_dn_bodyid)),
	dimnames = list(c(names(ss_tabs_to_read), ahn_downstream_dn_bodyid), c(names(ss_tabs_to_read), ahn_downstream_dn_bodyid)))

for (i in 1:length(ss_tabs_to_read)) {
	ahn_downstream_dn_adj[ahn_downstream_dns$autoid[ahn_downstream_dns$ahn == names(ss_tabs_to_read)[i]], names(ss_tabs_to_read)[i]] = ahn_downstream_dns$synapse_weight[ahn_downstream_dns$ahn == names(ss_tabs_to_read)[i]]
}

ahn_downstream_dn_thres = data.frame(bodyid = ahn_downstream_dn_bodyid)
ahn_downstream_dn_thres$type = ahn_downstream_dns$cell_type[match(ahn_downstream_dn_bodyid, ahn_downstream_dns$autoid)]
#determine categories by connectivity with AHNs for entire cell type rather than individual cell
ahn_downstream_dn_thres$conn = apply(ahn_downstream_dn_thres, MARGIN = 1,
	FUN = function(x) paste(names(c(msahn, mtahn))[colSums(ahn_downstream_dn_adj[ ahn_downstream_dn_thres$bodyid[ahn_downstream_dn_thres$type==x["type"]], names(c(msahn, mtahn)), drop = FALSE]) > 0], collapse = "_"))
ahn_downstream_dn_thres$type = ifelse(ahn_downstream_dn_thres$type %in% ahn_downstream_dn_summary$less_than_5_percent,
	ahn_downstream_dn_thres$type,
	"DN (summed)")
	
#only separate out type_conn for DN summed category
ahn_downstream_dn_thres$type_conn = ifelse(ahn_downstream_dn_thres$type == "DN (summed)",
	paste(ahn_downstream_dn_thres$type, ahn_downstream_dn_thres$conn, sep = "_"),
	ahn_downstream_dn_thres$type)
ahn_downstream_dn_thres$vertex_label = sapply(ahn_downstream_dn_thres$type_conn, FUN = function(x) sum(ahn_downstream_dn_thres$type_conn == x))
ahn_downstream_dn_thres$vertex_label = paste0(ahn_downstream_dn_thres$type, "(", ahn_downstream_dn_thres$vertex_label, ")")
ahn_downstream_dn_thres$vertex_color = dn_lut[ahn_downstream_dn_thres$type]
rownames(ahn_downstream_dn_adj) = c(names(ss_tabs_to_read), ahn_downstream_dn_thres$type_conn)
colnames(ahn_downstream_dn_adj) = rownames(ahn_downstream_dn_adj)

#sum ahn upstream by type and connectivity with AHNs
ahn_downstream_dn_adj = t(apply(t(ahn_downstream_dn_adj), 2, function(x) tapply(x, colnames(ahn_downstream_dn_adj), sum, na.rm = TRUE)))
ahn_downstream_dn_adj = apply(ahn_downstream_dn_adj, 2, function(x) tapply(x, rownames(ahn_downstream_dn_adj), sum, na.rm = TRUE))

write.csv(ahn_downstream_dn_adj, paste0("fanc_ahn_downstream_dn_adj_sum_syn_thres_", syn_thres, "_", Sys.Date(),".csv"))

ahn_downstream_dn_graph_sum = graph_from_adjacency_matrix(adjmatrix = ahn_downstream_dn_adj, mode = "directed", weighted = TRUE, diag = FALSE)
V(ahn_downstream_dn_graph_sum)$label = ahn_downstream_dn_thres$vertex_label[match(V(ahn_downstream_dn_graph_sum)$name, ahn_downstream_dn_thres$type_conn)]
V(ahn_downstream_dn_graph_sum)$label[match(c(msahn, mtahn), V(ahn_downstream_dn_graph_sum)$name)] = c("L", "R", "L", "R")
V(ahn_downstream_dn_graph_sum)$color = ahn_downstream_dn_thres$vertex_color[match(V(ahn_downstream_dn_graph_sum)$name, ahn_downstream_dn_thres$type_conn)]
V(ahn_downstream_dn_graph_sum)$color[match(c(msahn, mtahn), V(ahn_downstream_dn_graph_sum)$name)] = graph_lut[c("MsAHN", "MsAHN", "MtAHN", "MtAHN")]

tkplot_graph = tkplot(ahn_downstream_dn_graph_sum, canvas.width = 600, canvas.height = 600,
	edge.arrow.size = 0.8, edge.arrow.width = 1.5, curved = TRUE, edge.label = E(ahn_downstream_dn_graph_sum)$weight, edge.width = 2,
	vertex.color = V(ahn_downstream_dn_graph_sum)$color, vertex.label = V(ahn_downstream_dn_graph_sum)$name,
	margin = c(0,0,0,0), layout = layout_with_graphopt(ahn_downstream_dn_graph_sum, spring.constant = 5, charge = 0.01, spring.length = 50), vertex.label.family = "sans")
readline(prompt="Adjust graph layout by dragging nodes around in tkplot window.\nPress ENTER once done with adjusting node positions: ")
tkplot_layout_upstream = tkplot.getcoords(tkplot_graph)
tkplot.close(tkplot_graph)
#manually tweak layout and get coordinates
#tkplot_layout_upstream = tkplot.getcoords(3)
#after saving coords, go back and relabel vertices to proper names
svglite(paste0("ahn_downstream_DNs_FANC_summed", Sys.Date(),".svg"), width = 5, height = 5)
print(plot.igraph(ahn_downstream_dn_graph_sum, edge.arrow.size = 0.4,
	edge.arrow.width = 1.5, curved = TRUE, edge.label = E(ahn_downstream_dn_graph_sum)$weight, edge.width = log1p(E(ahn_downstream_dn_graph_sum)$weight)/2, vertex.color = V(ahn_downstream_dn_graph_sum)$color,
	margin = c(0,0,0,0), layout = tkplot_layout_upstream, vertex.label.family = "sans", vertex.label.cex = 0.65, edge.label.family = "sans", edge.label.cex = 0.6))
dev.off()





#specifically export matrix of DN upstream connectivity with the same DN categories as MANC
#this is hardcoded, need to update this if MANC DN categories change
manc_matched_dn_categories = c('DNg02','DNp54','DNp08','DNg32','DNp38','DNg29','MANC_29987','MANC_17526')
#make DN graph with summed categories
ahn_downstream_dn_adj = matrix(data = 0, nrow = length(c(names(ss_tabs_to_read), ahn_downstream_dn_bodyid)),
	ncol = length(c(names(ss_tabs_to_read), ahn_downstream_dn_bodyid)),
	dimnames = list(c(names(ss_tabs_to_read), ahn_downstream_dn_bodyid), c(names(ss_tabs_to_read), ahn_downstream_dn_bodyid)))

for (i in 1:length(ss_tabs_to_read)) {
	ahn_downstream_dn_adj[ahn_downstream_dns$autoid[ahn_downstream_dns$ahn == names(ss_tabs_to_read)[i]], names(ss_tabs_to_read)[i]] = ahn_downstream_dns$synapse_weight[ahn_downstream_dns$ahn == names(ss_tabs_to_read)[i]]
}

ahn_downstream_dn_thres = data.frame(bodyid = ahn_downstream_dn_bodyid)
ahn_downstream_dn_thres$type = ahn_downstream_dns$cell_type[match(ahn_downstream_dn_bodyid, ahn_downstream_dns$autoid)]
#determine categories by connectivity with AHNs for entire cell type rather than individual cell
ahn_downstream_dn_thres$conn = apply(ahn_downstream_dn_thres, MARGIN = 1,
	FUN = function(x) paste(names(c(msahn, mtahn))[colSums(ahn_downstream_dn_adj[ ahn_downstream_dn_thres$bodyid[ahn_downstream_dn_thres$type==x["type"]], names(c(msahn, mtahn)), drop = FALSE]) > 0], collapse = "_"))
ahn_downstream_dn_thres$type = ifelse(ahn_downstream_dn_thres$type %in% manc_matched_dn_categories,
	ahn_downstream_dn_thres$type,
	"DN (summed)")
	
#only separate out type_conn for DN summed category
ahn_downstream_dn_thres$type_conn = ifelse(ahn_downstream_dn_thres$type == "DN (summed)",
	paste(ahn_downstream_dn_thres$type, ahn_downstream_dn_thres$conn, sep = "_"),
	ahn_downstream_dn_thres$type)
ahn_downstream_dn_thres$vertex_label = sapply(ahn_downstream_dn_thres$type_conn, FUN = function(x) sum(ahn_downstream_dn_thres$type_conn == x))
ahn_downstream_dn_thres$vertex_label = paste0(ahn_downstream_dn_thres$type, "(", ahn_downstream_dn_thres$vertex_label, ")")
ahn_downstream_dn_thres$vertex_color = dn_lut[ahn_downstream_dn_thres$type]
rownames(ahn_downstream_dn_adj) = c(names(ss_tabs_to_read), ahn_downstream_dn_thres$type_conn)
colnames(ahn_downstream_dn_adj) = rownames(ahn_downstream_dn_adj)

#sum ahn upstream by type and connectivity with AHNs
ahn_downstream_dn_adj = t(apply(t(ahn_downstream_dn_adj), 2, function(x) tapply(x, colnames(ahn_downstream_dn_adj), sum, na.rm = TRUE)))
ahn_downstream_dn_adj = apply(ahn_downstream_dn_adj, 2, function(x) tapply(x, rownames(ahn_downstream_dn_adj), sum, na.rm = TRUE))

write.csv(ahn_downstream_dn_adj, paste0("ahn_downstream_DN_FANC_using_MANC_categories_syn_thres_", syn_thres, "_", Sys.Date(),".csv"))
