library(catmaid)
library(ggplot2)
library(plyr)
library(svglite)
library(igraph)

#note: Mt and MsAHN-r are actually -l using the FANC definition (soma side rather than projection side)--relabel the graphs and plots once made

msahn_neurons_of_interest = read.csv("FAFB MsAHN downstream neurons of interest 3.csv")
syn_thres = 3

mtahn_partners = catmaid_query_connected(3385431)
mtahn_partners$outgoing$cell_type = "ND"
mtahn_partners$outgoing$input_to_JON = FALSE
mtahn_partners$outgoing$output_from_JON = FALSE

for (i in 1:length(mtahn_partners$outgoing$partner)) {
	mtahn_partner_annotations = catmaid_get_annotations_for_skeletons(skid = mtahn_partners$outgoing$partner[i])
	#classify cell type by serially going through annotation types
	if ("Downstream of MtAHN-r (PN)" %in% mtahn_partner_annotations$annotation ||
		"Downstream of MtAHN-r (auditory PN)" %in% mtahn_partner_annotations$annotation ||
		"Downstream of MtAHN-r (multimodal PN)" %in% mtahn_partner_annotations$annotation) {
		mtahn_partners$outgoing$cell_type[i] = "IN"
		} else if ("Downstream of MtAHN-r (LN)" %in% mtahn_partner_annotations$annotation ||
		"Downstream of MtAHN-r (auditory LN)" %in% mtahn_partner_annotations$annotation ||
		"Downstream of MtAHN-r (multimodal LN)" %in% mtahn_partner_annotations$annotation) {
		mtahn_partners$outgoing$cell_type[i] = "IN"
		} else if ("Downstream of MtAHN-r (DN)" %in% mtahn_partner_annotations$annotation) {
		mtahn_partners$outgoing$cell_type[i] = "DN"
		} else if ("Downstream of MtAHN-r (AN)" %in% mtahn_partner_annotations$annotation) {
		mtahn_partners$outgoing$cell_type[i] = "AN"
		} else if ("Downstream of MtAHN-r (JO-A)" %in% mtahn_partner_annotations$annotation) {
		mtahn_partners$outgoing$cell_type[i] = "JON"
		} else if ("Downstream of MtAHN-r (JO-B)" %in% mtahn_partner_annotations$annotation) {
		mtahn_partners$outgoing$cell_type[i] = "JON"
		} else if ("Downstream of MtAHN-r (JON)" %in% mtahn_partner_annotations$annotation) {
		mtahn_partners$outgoing$cell_type[i] = "JON"
		} else if ("Downstream of MtAHN-r (efferent)" %in% mtahn_partner_annotations$annotation) {
		mtahn_partners$outgoing$cell_type[i] = "EN"
		} else if ("Downstream of MtAHN-r (unclassified)" %in% mtahn_partner_annotations$annotation) {
		mtahn_partners$outgoing$cell_type[i] = "unclassified"
	}
		
	#is neuron an output from JONs?
	if ("Downstream of MtAHN-r (auditory PN)" %in% mtahn_partner_annotations$annotation ||
		"Downstream of MtAHN-r (multimodal PN)" %in% mtahn_partner_annotations$annotation ||
		"Downstream of MtAHN-r (auditory LN)" %in% mtahn_partner_annotations$annotation ||
		"Downstream of MtAHN-r (multimodal LN)" %in% mtahn_partner_annotations$annotation ||
		"Downstream of MtAHN-r (output from JON)" %in% mtahn_partner_annotations$annotation) {
		mtahn_partners$outgoing$output_from_JON[i] = TRUE
		}
	
	#is neuron an input to JONs
	if ("Downstream of MtAHN-r (input to JON)" %in% mtahn_partner_annotations$annotation) {
		mtahn_partners$outgoing$input_to_JON[i] = TRUE
		}
	}



msahn_partners = catmaid_query_connected(2455455)
msahn_partners$outgoing$cell_type = "ND"
msahn_partners$outgoing$input_to_JON = FALSE
msahn_partners$outgoing$output_from_JON = FALSE

for (i in 1:length(msahn_partners$outgoing$partner)) {
	msahn_partner_annotations = catmaid_get_annotations_for_skeletons(skid = msahn_partners$outgoing$partner[i])
	#classify cell type by serially going through annotation types
	if ("Downstream of MsAHN-r (PN)" %in% msahn_partner_annotations$annotation) {
		msahn_partners$outgoing$cell_type[i] = "IN"
		} else if ("Downstream of MsAHN-r (LN)" %in% msahn_partner_annotations$annotation) {
		msahn_partners$outgoing$cell_type[i] = "IN"
		} else if ("Downstream of MsAHN-r (DN)" %in% msahn_partner_annotations$annotation) {
		msahn_partners$outgoing$cell_type[i] = "DN"
		} else if ("Downstream of MsAHN-r (AN)" %in% msahn_partner_annotations$annotation) {
		msahn_partners$outgoing$cell_type[i] = "AN"
		} else if ("Downstream of MsAHN-r (JO-A)" %in% msahn_partner_annotations$annotation) {
		msahn_partners$outgoing$cell_type[i] = "JON"
		} else if ("Downstream of MsAHN-r (JO-B)" %in% msahn_partner_annotations$annotation) {
		msahn_partners$outgoing$cell_type[i] = "JON"
		} else if ("Downstream of MsAHN-r (JON)" %in% msahn_partner_annotations$annotation) {
		msahn_partners$outgoing$cell_type[i] = "JON"
		} else if ("Downstream of MsAHN-r (Efferent)" %in% msahn_partner_annotations$annotation) { #note spelling of Efferent with big E (compared to MtAHN's with little e)
		msahn_partners$outgoing$cell_type[i] = "EN"
		} else if ("Downstream of MsAHN-r (unclassified)" %in% msahn_partner_annotations$annotation) {
		msahn_partners$outgoing$cell_type[i] = "unclassified"
		} else if ("Downstream of MsAHN-r (LPTC)" %in% msahn_partner_annotations$annotation) {
		msahn_partners$outgoing$cell_type[i] = "LPTC"
	}
		
	#is neuron downstream of JONs?
	if ("Downstream of MsAHN-r (output from JON)" %in% msahn_partner_annotations$annotation) {
		msahn_partners$outgoing$output_from_JON[i] = TRUE
		}
	
	#is neuron an input to JONs?
	if ("Downstream of MsAHN-r (input to JON)" %in% msahn_partner_annotations$annotation) {
		msahn_partners$outgoing$input_to_JON[i] = TRUE
		}
	}
	#not dealt with: (input and output to/from DNa02, input/output with GF)

#make adjacency matrix
bodyid_list = c("MsAHN-r", "MtAHN-r", unique(c(mtahn_partners$outgoing$partner[mtahn_partners$outgoing$cell_type != "ND"], msahn_partners$outgoing$partner[msahn_partners$outgoing$cell_type != "ND"])))
ahn_adj = matrix(data = 0, nrow = length(bodyid_list), ncol = length(bodyid_list), dimnames = list(bodyid_list, bodyid_list))
ahn_adj["MsAHN-r", as.character(msahn_partners$outgoing$partner[msahn_partners$outgoing$cell_type != "ND"])] = msahn_partners$outgoing$syn.count[msahn_partners$outgoing$cell_type != "ND"]
ahn_adj["MtAHN-r", as.character(mtahn_partners$outgoing$partner[mtahn_partners$outgoing$cell_type != "ND"])] = mtahn_partners$outgoing$syn.count[mtahn_partners$outgoing$cell_type != "ND"]

#start building graphs
graph_color_lut = c("#5ce7e1", "#9632fa", "#ff6666", "#1910c7", "#db58d7", "#ffe099", "#ffff64", "#b9b4ff", "#ffc9b4", "#ffA0BE", "#aaaaaa")
names(graph_color_lut) = c("DN", "AN", "IN", "MN", "JON", "MsAHN-r", "MtAHN-r", "EN", "ASN", "LPTC", "ND")

ahn_downstream_bodyid = data.frame(bodyid = bodyid_list, cell_type = NA)
ahn_downstream_bodyid$cell_type = ifelse(ahn_downstream_bodyid$bodyid %in% mtahn_partners$outgoing$partner,
	mtahn_partners$outgoing$cell_type[match(ahn_downstream_bodyid$bodyid, mtahn_partners$outgoing$partner)],
	msahn_partners$outgoing$cell_type[match(ahn_downstream_bodyid$bodyid, msahn_partners$outgoing$partner)])
ahn_downstream_bodyid$cell_type[1:2] = ahn_downstream_bodyid$bodyid[1:2]

ahn_downstream_bodyid$type_conn = ifelse(colSums(ahn_adj[c("MsAHN-r","MtAHN-r"),ahn_downstream_bodyid$bodyid] > 0) == 2, 
	paste(ahn_downstream_bodyid$cell_type,"AHN",sep = "_"),
	ifelse(ahn_adj["MsAHN-r", ahn_downstream_bodyid$bodyid] > 0, paste(ahn_downstream_bodyid$cell_type,"MsAHN-r",sep = "_"),
	ifelse(ahn_adj["MtAHN-r", ahn_downstream_bodyid$bodyid] > 0, paste(ahn_downstream_bodyid$cell_type,"MtAHN-r",sep = "_"),
	"ERROR"
	)))
ahn_downstream_bodyid$type_conn[1:2] = ahn_downstream_bodyid$bodyid[1:2]

#vertex label is cell count
ahn_downstream_bodyid$vertex_label = apply(ahn_downstream_bodyid, MARGIN = 1, 
	FUN = function(x) sum(ahn_downstream_bodyid$type_conn == x[3]))
ahn_downstream_bodyid$type_conn[1:2] = ahn_downstream_bodyid$bodyid[1:2]
ahn_downstream_bodyid$vertex_color = graph_color_lut[ahn_downstream_bodyid$cell_type]

#sum ahn downstream by type and connectivity with AHNs
#note this is not thresholded
ahn_downstream_adj_sum = ahn_adj
ahn_downstream_adj_sum = t(apply(t(ahn_downstream_adj_sum), 2, function(x) tapply(x, ahn_downstream_bodyid$type_conn, sum, na.rm = TRUE)))
ahn_downstream_adj_sum = apply(ahn_downstream_adj_sum, 2, function(x) tapply(x, ahn_downstream_bodyid$type_conn, sum, na.rm = TRUE))

ahn_downstream_graph_sum = graph_from_adjacency_matrix(adjmatrix = ahn_downstream_adj_sum, mode = "directed", weighted = TRUE, diag = FALSE)
V(ahn_downstream_graph_sum)$label = ahn_downstream_bodyid$vertex_label[match(V(ahn_downstream_graph_sum)$name, ahn_downstream_bodyid$type_conn)]
V(ahn_downstream_graph_sum)$label[V(ahn_downstream_graph_sum)$name %in% c("MsAHN-r", "MtAHN-r")] = c("R","R")
V(ahn_downstream_graph_sum)$color = ahn_downstream_bodyid$vertex_color[match(V(ahn_downstream_graph_sum)$name, ahn_downstream_bodyid$type_conn)]

tkplot_index = tkplot(ahn_downstream_graph_sum, canvas.width = 650, canvas.height = 650, width = 650, height = 650,
	edge.arrow.size = 0.8, edge.arrow.width = 1.5, curved = TRUE, edge.label = E(ahn_downstream_graph_sum)$weight, edge.width = 2,
	vertex.label = V(ahn_downstream_graph_sum)$name, vertex.color = V(ahn_downstream_graph_sum)$color,
	margin = c(0,0,0,0), layout = layout_with_graphopt(ahn_downstream_graph_sum, spring.constant = 5, charge = 0.01, spring.length = 50), vertex.label.family = "sans")
readline(prompt="Adjust graph layout by dragging nodes around in tkplot window. Press ENTER once done with adjusting node positions: ")
#manually tweak layout and get coordinates
tkplot_layout = tkplot.getcoords(tkplot_index)
#after saving coords, go back and relabel vertices to proper names
svglite(paste0("AHN_downstream_FAFB_", Sys.Date(),".svg"), width = 4, height = 4)
print(plot.igraph(ahn_downstream_graph_sum, vertex.size = 20, edge.arrow.size = 0.4,
	edge.arrow.width = 1.5, curved = TRUE, edge.label = E(ahn_downstream_graph_sum)$weight, edge.width = log1p(E(ahn_downstream_graph_sum)$weight)/2, vertex.color = V(ahn_downstream_graph_sum)$color,
	margin = c(0,0,0,0), layout = tkplot_layout, vertex.label.family = "sans", vertex.label.cex = 0.65, edge.label.family = "sans", edge.label.cex = 0.6))
dev.off()
tkplot.close(tkplot_index)

#also find further downstream connectivity to pathways such as GF, DNa02, LPTCins
#msahn_neurons_of_interest = c(GF-R = 4947529, DNa02-R = 12526673)
#msahn_neurons_of_interest$skid = as.character(msahn_neurons_of_interest$skid)
#add neurons of interest to adjacency matrix
#syn threshold applies here
#msahn
bodyid_list = c("MsAHN-r", msahn_partners$outgoing$partner[msahn_partners$outgoing$cell_type != "ND" & msahn_partners$outgoing$syn.count >= syn_thres])
bodyid_list = c(bodyid_list, msahn_neurons_of_interest$skid[!(msahn_neurons_of_interest$skid %in% bodyid_list)])
msahn_adj = matrix(data = 0, nrow = length(bodyid_list), ncol = length(bodyid_list), dimnames = list(bodyid_list, bodyid_list))
msahn_adj["MsAHN-r", as.character(msahn_partners$outgoing$partner[msahn_partners$outgoing$cell_type != "ND"  & msahn_partners$outgoing$syn.count >= syn_thres])] = msahn_partners$outgoing$syn.count[msahn_partners$outgoing$cell_type != "ND"  & msahn_partners$outgoing$syn.count >= syn_thres]

#for each neuron in msahn_neurons_of_interest, check if any have upstream or downstream connectivity with MsAHN downstream partners
for (i in msahn_neurons_of_interest$skid) {
	temp_neuron_connections = catmaid_query_connected(i, minimum_synapses = syn_thres)
	temp_partner_overlap = intersect(temp_neuron_connections$outgoing$partner, msahn_partners$outgoing$partner[msahn_partners$outgoing$cell_type != "ND" & msahn_partners$outgoing$syn.count >= syn_thres])
	if (length(temp_partner_overlap) > 0) {
		msahn_adj[as.character(i), as.character(temp_partner_overlap)] = temp_neuron_connections$outgoing$syn.count[match(temp_partner_overlap, temp_neuron_connections$outgoing$partner)]
	}
	temp_partner_overlap = intersect(temp_neuron_connections$incoming$partner, msahn_partners$outgoing$partner[msahn_partners$outgoing$cell_type != "ND" & msahn_partners$outgoing$syn.count >= syn_thres])
	if (length(temp_partner_overlap) > 0) {
		msahn_adj[as.character(temp_partner_overlap), as.character(i)] = temp_neuron_connections$incoming$syn.count[match(temp_partner_overlap, temp_neuron_connections$incoming$partner)]
	}
}
#get rid of neuron in msahn_neurons_of_interest with no connections
temp_no_conn = sapply(1:dim(msahn_adj)[1], FUN = function(x) all(msahn_adj[x,] == 0) & all(msahn_adj[,x] == 0))
msahn_adj = msahn_adj[!temp_no_conn, !temp_no_conn]
rm(temp_no_conn)

#make labels
# temp_msahn_labels = ifelse(colnames(msahn_adj) %in% msahn_partners$outgoing$partner,
		# paste(msahn_partners$outgoing$cell_type[match(colnames(msahn_adj), msahn_partners$outgoing$partner)], "MsAHN", sep = "_"),
		# ifelse(colnames(msahn_adj) %in% msahn_neurons_of_interest$skid,
			# paste(msahn_neurons_of_interest$cell_type[match(colnames(msahn_adj), msahn_neurons_of_interest$skid)], "ds", sep = "_"),
			# colnames(msahn_adj)))
			
temp_msahn_labels = ifelse(colnames(msahn_adj) %in% msahn_partners$outgoing$partner & msahn_adj[1,] >= syn_thres,
		paste(msahn_partners$outgoing$cell_type[match(colnames(msahn_adj), msahn_partners$outgoing$partner)], "MsAHN", sep = "_"),
		ifelse(colnames(msahn_adj) %in% msahn_neurons_of_interest$skid,
			paste(msahn_neurons_of_interest$cell_type[match(colnames(msahn_adj), msahn_neurons_of_interest$skid)], "ds", sep = "_"),
			colnames(msahn_adj)))

#sum ahn downstream by type and connectivity with AHNs
msahn_adj_sum = msahn_adj
msahn_adj_sum = t(apply(t(msahn_adj_sum), 2, function(x) tapply(x, temp_msahn_labels, sum, na.rm = TRUE)))
msahn_adj_sum = apply(msahn_adj_sum, 2, function(x) tapply(x, temp_msahn_labels, sum, na.rm = TRUE))

msahn_downstream_graph_sum = graph_from_adjacency_matrix(adjmatrix = msahn_adj_sum, mode = "directed", weighted = TRUE, diag = FALSE)
V(msahn_downstream_graph_sum)$label = sapply(V(msahn_downstream_graph_sum)$name, FUN = function(x) sum(temp_msahn_labels == x))
V(msahn_downstream_graph_sum)$label[V(msahn_downstream_graph_sum)$name == "MsAHN-r"] = "R"
V(msahn_downstream_graph_sum)$color = graph_color_lut[match(sapply(V(msahn_downstream_graph_sum)$name, FUN = function(x) strsplit(x, split = "_")[[1]][1]), names(graph_color_lut))]
V(msahn_downstream_graph_sum)$label[is.na(V(msahn_downstream_graph_sum)$color)] = paste0(V(msahn_downstream_graph_sum)$name[is.na(V(msahn_downstream_graph_sum)$color)], "(",
	V(msahn_downstream_graph_sum)$label[is.na(V(msahn_downstream_graph_sum)$color)], ")")
V(msahn_downstream_graph_sum)$color[is.na(V(msahn_downstream_graph_sum)$color)] = "#FFFFFF"

tkplot_index = tkplot(msahn_downstream_graph_sum, canvas.width = 650, canvas.height = 650, width = 650, height = 650,
	edge.arrow.size = 0.8, edge.arrow.width = 1.5, edge.curved = TRUE, edge.label = E(msahn_downstream_graph_sum)$weight, edge.width = 2,
	vertex.label = V(msahn_downstream_graph_sum)$name, vertex.color = V(msahn_downstream_graph_sum)$color,
	margin = c(0,0,0,0), layout = layout_with_graphopt(msahn_downstream_graph_sum, spring.constant = 5, charge = 0.01, spring.length = 50), vertex.label.family = "sans")
readline(prompt="Adjust graph layout by dragging nodes around in tkplot window. Press ENTER once done with adjusting node positions: ")
#manually tweak layout and get coordinates
tkplot_layout = tkplot.getcoords(tkplot_index)
#after saving coords, go back and relabel vertices to proper names
svglite(paste0("MsAHN_downstream_FAFB_", Sys.Date(),".svg"), width = 4, height = 4)
print(plot.igraph(msahn_downstream_graph_sum, vertex.size = 20, edge.arrow.size = 0.4,
	edge.arrow.width = 1.5, edge.curved = 0.1, edge.label = E(msahn_downstream_graph_sum)$weight, edge.width = log1p(E(msahn_downstream_graph_sum)$weight)/2, vertex.color = V(msahn_downstream_graph_sum)$color,
	margin = c(0,0,0,0), layout = tkplot_layout, vertex.label.family = "sans", vertex.label.cex = 0.65, edge.label.family = "sans", edge.label.cex = 0.6))
dev.off()
tkplot.close(tkplot_index)




#also find further downstream connectivity to pathways such as GF, DNa02, LPTCins
#msahn_neurons_of_interest = c(GF-R = 4947529, DNa02-R = 12526673)
#mtahn to GF
gf_bodyid = 4947529
bodyid_list = c("MtAHN-r", mtahn_partners$outgoing$partner[mtahn_partners$outgoing$cell_type != "ND" & mtahn_partners$outgoing$syn.count >= syn_thres])
#gf is already in bodyid_list

if (!(gf_bodyid %in% bodyid_list)) bodyid_list = c(bodyid_list, gf_bodyid)
mtahn_adj = matrix(data = 0, nrow = length(bodyid_list), ncol = length(bodyid_list), dimnames = list(bodyid_list, bodyid_list))
mtahn_adj["MtAHN-r", as.character(mtahn_partners$outgoing$partner[mtahn_partners$outgoing$cell_type != "ND"  & mtahn_partners$outgoing$syn.count >= syn_thres])] = mtahn_partners$outgoing$syn.count[mtahn_partners$outgoing$cell_type != "ND"  & mtahn_partners$outgoing$syn.count >= syn_thres]

#for each neuron in msahn_neurons_of_interest, check if any have upstream or downstream connectivity with MsAHN downstream partners
temp_neuron_connections = catmaid_query_connected(gf_bodyid, minimum_synapses = syn_thres)
temp_partner_overlap = intersect(temp_neuron_connections$incoming$partner, mtahn_partners$outgoing$partner[mtahn_partners$outgoing$cell_type != "ND" & mtahn_partners$outgoing$syn.count >= syn_thres])
if (length(temp_partner_overlap) > 0) {
	mtahn_adj[as.character(temp_partner_overlap), as.character(gf_bodyid)] = temp_neuron_connections$incoming$syn.count[match(temp_partner_overlap, temp_neuron_connections$incoming$partner)]
}

#make labels
# temp_msahn_labels = ifelse(colnames(msahn_adj) %in% msahn_partners$outgoing$partner,
		# paste(msahn_partners$outgoing$cell_type[match(colnames(msahn_adj), msahn_partners$outgoing$partner)], "MsAHN", sep = "_"),
		# ifelse(colnames(msahn_adj) %in% msahn_neurons_of_interest$skid,
			# paste(msahn_neurons_of_interest$cell_type[match(colnames(msahn_adj), msahn_neurons_of_interest$skid)], "ds", sep = "_"),
			# colnames(msahn_adj)))


temp_mtahn_labels = ifelse(colnames(mtahn_adj) %in% mtahn_partners$outgoing$partner & mtahn_adj[1,] >= syn_thres,
		paste(mtahn_partners$outgoing$cell_type[match(colnames(mtahn_adj), mtahn_partners$outgoing$partner)], "MtAHN", sep = "_"),
		colnames(mtahn_adj))
temp_mtahn_labels[names(temp_mtahn_labels) == gf_bodyid] = "GF-R"
#make copy of temp labels to make variant of adj matrix that separates GF-R upstream cells
temp_mtahn_labels_gf = ifelse(mtahn_adj[,as.character(gf_bodyid)] >= syn_thres,
	paste(temp_mtahn_labels, "GF", sep = "_"),
	temp_mtahn_labels)

#sum ahn downstream by type and connectivity with AHNs
mtahn_adj_sum = mtahn_adj
mtahn_adj_sum = t(apply(t(mtahn_adj_sum), 2, function(x) tapply(x, temp_mtahn_labels, sum, na.rm = TRUE)))
mtahn_adj_sum = apply(mtahn_adj_sum, 2, function(x) tapply(x, temp_mtahn_labels, sum, na.rm = TRUE))

#variant that separates GF us cells
#mtahn_adj_sum = mtahn_adj
#mtahn_adj_sum = t(apply(t(mtahn_adj_sum), 2, function(x) tapply(x, temp_mtahn_labels_gf, sum, na.rm = TRUE)))
#mtahn_adj_sum = apply(mtahn_adj_sum, 2, function(x) tapply(x, temp_mtahn_labels_gf, sum, na.rm = TRUE))

mtahn_downstream_graph_sum = graph_from_adjacency_matrix(adjmatrix = mtahn_adj_sum, mode = "directed", weighted = TRUE, diag = FALSE)
V(mtahn_downstream_graph_sum)$label = sapply(V(mtahn_downstream_graph_sum)$name, FUN = function(x) sum(temp_mtahn_labels == x))
V(mtahn_downstream_graph_sum)$label[V(mtahn_downstream_graph_sum)$name == "MtAHN-r"] = "R"
V(mtahn_downstream_graph_sum)$color = graph_color_lut[match(sapply(V(mtahn_downstream_graph_sum)$name, FUN = function(x) strsplit(x, split = "_")[[1]][1]), names(graph_color_lut))]
V(mtahn_downstream_graph_sum)$label[is.na(V(mtahn_downstream_graph_sum)$color)] = paste0(V(mtahn_downstream_graph_sum)$name[is.na(V(mtahn_downstream_graph_sum)$color)], "(",
	V(mtahn_downstream_graph_sum)$label[is.na(V(mtahn_downstream_graph_sum)$color)], ")")
V(mtahn_downstream_graph_sum)$color[is.na(V(mtahn_downstream_graph_sum)$color)] = "#FFFFFF"

tkplot_index = tkplot(mtahn_downstream_graph_sum, canvas.width = 650, canvas.height = 650, width = 650, height = 650,
	edge.arrow.size = 0.8, edge.arrow.width = 1.5, edge.curved = FALSE, edge.label = E(mtahn_downstream_graph_sum)$weight, edge.width = 2,
	vertex.label = V(mtahn_downstream_graph_sum)$name, vertex.color = V(mtahn_downstream_graph_sum)$color,
	margin = c(0,0,0,0), layout = layout_with_graphopt(mtahn_downstream_graph_sum, spring.constant = 5, charge = 0.01, spring.length = 50), vertex.label.family = "sans")
readline(prompt="Adjust graph layout by dragging nodes around in tkplot window. Press ENTER once done with adjusting node positions: ")
#manually tweak layout and get coordinates
tkplot_layout = tkplot.getcoords(tkplot_index)
#after saving coords, go back and relabel vertices to proper names
svglite(paste0("MtAHN_downstream_FAFB_", Sys.Date(),".svg"), width = 4, height = 4)
print(plot.igraph(mtahn_downstream_graph_sum, vertex.size = 20, edge.arrow.size = 0.4,
	edge.arrow.width = 1.5, edge.curved = FALSE, edge.label = E(mtahn_downstream_graph_sum)$weight, edge.width = log1p(E(mtahn_downstream_graph_sum)$weight)/2, vertex.color = V(mtahn_downstream_graph_sum)$color,
	margin = c(0,0,0,0), layout = tkplot_layout, vertex.label.family = "sans", vertex.label.cex = 0.65, edge.label.family = "sans", edge.label.cex = 0.6))
dev.off()
tkplot.close(tkplot_index)



#sum_cell_type = aggregate(mtahn_partners$outgoing$syn.count, by=list(mtahn_partners$outgoing$cell_type), FUN=sum)
#plot cell type for all neurons (tracing progress)
#ggplot(sum_cell_type, aes(fill=Group.1, y=x, x=1)) + 
#	geom_bar(position="fill", stat="identity")
#this is for mtahn
sum_cell_type = ddply(mtahn_partners$outgoing, "cell_type", summarise,synapses= sum(syn.count))

#msahn
sum_cell_type_msahn = ddply(msahn_partners$outgoing, "cell_type", summarise,synapses= sum(syn.count))

#remove orphans and untraced postsynaptic sites (mtahn)
mtahn_partners_traced = mtahn_partners$outgoing[ which(mtahn_partners$outgoing$cell_type!="ND"),]
sum_cell_type_traced = ddply(mtahn_partners_traced, "cell_type", summarise,synapses= sum(syn.count))
for (i in 1:nrow(mtahn_partners_traced)) {
	if (mtahn_partners_traced$input_to_JON[i] && mtahn_partners_traced$output_from_JON[i]) {
		mtahn_partners_traced$JON_pre_post[i] = "both"
		} else if (mtahn_partners_traced$input_to_JON[i]) {
		mtahn_partners_traced$JON_pre_post[i] = "presynaptic"
		} else if (mtahn_partners_traced$output_from_JON[i]) {
		mtahn_partners_traced$JON_pre_post[i] = "postsynaptic"
		} else {
		mtahn_partners_traced$JON_pre_post[i] = "neither"
		}
	}
sum_JON_pre_post_traced = ddply(mtahn_partners_traced, "JON_pre_post", summarise,synapses= sum(syn.count))
	
#remove below syn thres
mtahn_partners_certain = mtahn_partners_traced[mtahn_partners_traced$syn.count>=syn_thres,]
sum_cell_type_certain = ddply(mtahn_partners_certain, "cell_type", summarise,synapses= sum(syn.count))
sum_JON_pre_post_certain = ddply(mtahn_partners_certain, "JON_pre_post", summarise,synapses= sum(syn.count))


#do the same for msahn
msahn_partners_traced = msahn_partners$outgoing[ which(msahn_partners$outgoing$cell_type!="ND"),]
msahn_sum_cell_type_traced = ddply(msahn_partners_traced, "cell_type", summarise,synapses= sum(syn.count))
for (i in 1:nrow(msahn_partners_traced)) {
	if (msahn_partners_traced$input_to_JON[i] && msahn_partners_traced$output_from_JON[i]) {
		msahn_partners_traced$JON_pre_post[i] = "both"
		} else if (msahn_partners_traced$input_to_JON[i]) {
		msahn_partners_traced$JON_pre_post[i] = "presynaptic"
		} else if (msahn_partners_traced$output_from_JON[i]) {
		msahn_partners_traced$JON_pre_post[i] = "postsynaptic"
		} else {
		msahn_partners_traced$JON_pre_post[i] = "neither"
		}
	}
msahn_sum_JON_pre_post_traced = ddply(msahn_partners_traced, "JON_pre_post", summarise,synapses= sum(syn.count))
#remove cells below syn thres
msahn_partners_certain = msahn_partners_traced[msahn_partners_traced$syn.count>=syn_thres,]
msahn_sum_cell_type_certain = ddply(msahn_partners_certain, "cell_type", summarise,synapses= sum(syn.count))
#sum_JON_pre_post_certain = ddply(msahn_partners_certain, "JON_pre_post", summarise,synapses= sum(syn.count))
	
#save plots
#combine Ms and Mt tables
sum_cell_type$ahn = "MtAHN-r"
sum_cell_type_msahn$ahn = "MsAHN-r"
ahn_sum_cell_type = rbind(sum_cell_type_msahn, sum_cell_type)
ahn_sum_cell_type$cell_type = factor(ahn_sum_cell_type$cell_type, levels=names(graph_color_lut)[names(graph_color_lut) %in% ahn_sum_cell_type$cell_type])
sum_cell_type_traced$ahn = "MtAHN-r"
msahn_sum_cell_type_traced$ahn = "MsAHN-r"
ahn_sum_cell_type_traced = rbind(msahn_sum_cell_type_traced, sum_cell_type_traced)
ahn_sum_cell_type_traced$cell_type = factor(ahn_sum_cell_type_traced$cell_type, levels=names(graph_color_lut)[names(graph_color_lut) %in% ahn_sum_cell_type_traced$cell_type])
sum_cell_type_certain$ahn = "MtAHN-r"
msahn_sum_cell_type_certain$ahn = "MsAHN-r"
ahn_sum_cell_type_certain = rbind(msahn_sum_cell_type_certain, sum_cell_type_certain)
ahn_sum_cell_type_certain$cell_type = factor(ahn_sum_cell_type_certain$cell_type, levels=names(graph_color_lut)[names(graph_color_lut) %in% ahn_sum_cell_type_certain$cell_type])

#tracing progress
#cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#009292" , "#222222")

png(paste0("FAFB AHN tracing progress_", Sys.Date(), ".png"), height=500,width=300)
print(ggplot(ahn_sum_cell_type, aes(fill=cell_type, y=synapses, x=ahn)) + 
	geom_bar(position="fill", stat="identity") +
	xlab("") +
	ylab("% synapses") +
	theme(text = element_text(size=20), axis.line.y.left = element_line(color = "black"), axis.line.x.bottom = element_line(color = "black"), panel.background = element_blank()) +
	scale_fill_manual(values=graph_color_lut[names(graph_color_lut) %in% ahn_sum_cell_type$cell_type]))
ggsave(file=paste0("FAFB AHN tracing progress_", Sys.Date(), ".svg"))
dev.off()

# png(paste0("MsAHN-r tracing progress_", Sys.Date(), ".png"),height=500,width=200)
# print(ggplot(sum_cell_type_msahn, aes(fill=cell_type, y=synapses, x=1)) + 
	# geom_bar(position="fill", stat="identity", width = 1) +
	# xlab("MsAHN-r") +
	# ylab("% synapses") +
	# scale_x_discrete() +
	# theme(text = element_text(size=20)) +
	# scale_fill_manual(values=cbPalette))
# ggsave(file=paste0("MsAHN-r tracing progress_", Sys.Date(), ".svg"))
# dev.off()

#cell type with no orphans/untraced
png(paste0("FAFB AHN downstream cell types_", Sys.Date(), ".png"),height=500,width=300)
print(ggplot(ahn_sum_cell_type_traced, aes(fill=cell_type, y=synapses, x=ahn)) + 
	geom_bar(position="fill", stat="identity") +
	xlab("") +
	ylab("% synapses") +
	theme(text = element_text(size=20), axis.line.y.left = element_line(color = "black"), axis.line.x.bottom = element_line(color = "black"), panel.background = element_blank()) +
	scale_fill_manual(values=graph_color_lut[names(graph_color_lut) %in% ahn_sum_cell_type_traced$cell_type]))
ggsave(file=paste0("FAFB AHN downstream cell types_", Sys.Date(), ".svg"))
dev.off()

# png(paste0("MsAHN-r downstream cell types_", Sys.Date(), ".png"),height=500,width=200)
# print(ggplot(msahn_sum_cell_type_traced, aes(fill=cell_type, y=synapses, x=1)) + 
	# geom_bar(position="fill", stat="identity", width = 1) +
	# xlab("MsAHN-r") +
	# ylab("% synapses") +
	# scale_x_discrete()) +
	# theme(text = element_text(size=20)) +
	# scale_fill_manual(values=cbPalette)
# ggsave(file=paste0("MsAHN-r downstream cell types_", Sys.Date(), ".svg"))
# dev.off()

#cell type at syn_thres or above
png(paste("FAFB MsAHN-r MtAHN-r downstream cell types syn_thres ", syn_thres, ".png"),height=500,width=300)
print(ggplot(ahn_sum_cell_type_certain, aes(fill=cell_type, y=synapses, x=ahn)) + 
	geom_bar(position="fill", stat="identity") +
	xlab("") +
	ylab("% synapses") +
	theme(text = element_text(size=20), axis.line.y.left = element_line(color = "black"), axis.line.x.bottom = element_line(color = "black"), panel.background = element_blank()) +
	scale_fill_manual(values=graph_color_lut[names(graph_color_lut) %in% ahn_sum_cell_type_certain$cell_type]))
ggsave(file=paste("FAFB MtAHN-r MtAHN-r downstream cell types syn_thres ", syn_thres, ".svg"))
dev.off()


#JON pre post > 1 synapse
# png("MtAHN-r downstream JON pre post more than 1 synapse.png",height=500,width=300)
# print(ggplot(sum_JON_pre_post_certain, aes(fill=JON_pre_post, y=synapses, x=1)) + 
	# geom_bar(position="fill", stat="identity") +
	# xlab("") +
	# ylab("% synapses") +
	# scale_x_discrete()) +
	# theme(text = element_text(size=20)) +
	# scale_fill_manual(values=cbPalette)
# ggsave(file="MtAHN-r downstream JON pre post more than 1 synapse.svg")
# dev.off()


#plot the neurons
# library(elmr)
# userMatrix = matrix(c(1,0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1), nrow = 4, ncol = 4)

# neurons = read.neurons.catmaid(mtahn_partners_certain[ which(mtahn_partners_certain$cell_type == "DN"), 'partner'])
# plot3d(neurons, soma = TRUE, lwd = 1, line_antialias = TRUE)
# plot3d(FAFB14, depth_mask = FALSE)
# par3d(windowRect=c(85,85,800,800),zoom=0.78,userMatrix=userMatrix)
# rgl.snapshot("MtAHN-r DNs.png",fmt="png")
# rgl.close()

# neurons = read.neurons.catmaid(mtahn_partners_certain[ which(mtahn_partners_certain$cell_type == "AN"), 'partner'])
# plot3d(neurons, soma = TRUE, lwd = 1, line_antialias = TRUE)
# plot3d(FAFB14, depth_mask = FALSE)
# par3d(windowRect=c(85,85,800,800),zoom=0.78,userMatrix=userMatrix)
# rgl.snapshot("MtAHN-r ANs.png",fmt="png")
# rgl.close()

# neurons = read.neurons.catmaid(mtahn_partners_certain[ which(mtahn_partners_certain$cell_type == "LN"), 'partner'])
# plot3d(neurons, soma = TRUE, lwd = 1, line_antialias = TRUE)
# plot3d(FAFB14, depth_mask = FALSE)
# par3d(windowRect=c(85,85,800,800),zoom=0.78,userMatrix=userMatrix)
# rgl.snapshot("MtAHN-r LNs.png",fmt="png")
# rgl.close()

# neurons = read.neurons.catmaid(mtahn_partners_certain[ which(mtahn_partners_certain$cell_type == "PN"), 'partner'])
# plot3d(neurons, soma = TRUE, lwd = 1, line_antialias = TRUE)
# plot3d(FAFB14, depth_mask = FALSE)
# par3d(windowRect=c(85,85,800,800),zoom=0.78,userMatrix=userMatrix)
# rgl.snapshot("MtAHN-r PNs.png",fmt="png")
# rgl.close()

# neurons = read.neurons.catmaid(mtahn_partners_certain[ which(mtahn_partners_certain$cell_type == "JO-A" | mtahn_partners_certain$cell_type == "JON (other)" | mtahn_partners_certain$cell_type == "JO-B"), 'partner'])
# plot3d(neurons, soma = TRUE, lwd = 1, line_antialias = TRUE)
# plot3d(FAFB14, depth_mask = FALSE)
# par3d(windowRect=c(85,85,800,800),zoom=0.78,userMatrix=userMatrix)
# rgl.snapshot("MtAHN-r JONs.png",fmt="png")
# rgl.close()

# neurons = read.neurons.catmaid(mtahn_partners_certain[ which(mtahn_partners_certain$cell_type == "efferent"), 'partner'])
# plot3d(neurons, soma = TRUE, lwd = 1, line_antialias = TRUE)
# plot3d(FAFB14, depth_mask = FALSE)
# par3d(windowRect=c(85,85,800,800),zoom=0.78,userMatrix=userMatrix)
# rgl.snapshot("MtAHN-r efferents.png",fmt="png")
# rgl.close()
