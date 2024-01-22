library(neuprintr)
library(igraph)
library(svglite)
library(malevnc)
library(gplots)
library(viridis)
#library(plyr)
library(dplyr)
library(googlesheets4)

us_syn_thres = 10
ds_syn_thres = 3
msahns = c(13926, 12536)

tect_in_groups = c(22289, 21307, 18519, 14502, 14625, 15337, 22521, 21381, 13060, 15788, 15113)

temp_list = neuprintr::neuprint_list2df(neuprintr::neuprint_fetch_custom(cypher=paste0("MATCH (a:Neuron)-[w:ConnectsTo]->(b:Neuron) WHERE b.bodyId IN [", paste0(msahns, collapse = ","), "] AND a.type IN ['DNg02','DNp54'] AND w.weight >=", us_syn_thres," RETURN DISTINCT a.bodyId AS bodyid, a.class AS class, a.group AS group, a.type AS temp_label"), timeout=2000))
temp_list_2 = neuprintr::neuprint_list2df(neuprintr::neuprint_fetch_custom(cypher=paste0("MATCH (a:Neuron)-[w:ConnectsTo]->(b:Neuron) WHERE a.bodyId IN [", paste0(msahns, collapse = ","), "] AND w.weight >=", ds_syn_thres," AND b.group IN [", paste0(tect_in_groups, collapse = ","), "] RETURN DISTINCT b.bodyId AS bodyid, b.class AS class, b.group AS group"), timeout=2000))
temp_list_2$temp_label = "Tect IN"
temp_list_3 = neuprintr::neuprint_list2df(neuprintr::neuprint_fetch_custom(cypher=paste0("MATCH (a:Neuron) WHERE a.type IN ['DLMn a b', 'DLMn c-f', 'DVMn 1a-c','DVMn 2a b','DVMn 3a b'] RETURN a.bodyId AS bodyid, a.class AS class, a.group AS group, a.type AS temp_label"), timeout=2000))
temp_list = rbind(data.frame(bodyid=msahns, class="Ascending Interneuron", group="MsAHN", temp_label="MsAHN"), temp_list, temp_list_2, temp_list_3)

dn_adj = neuprint_get_adjacency_matrix(temp_list$bodyid, timeout = 2000)

dn_adj[dn_adj<ds_syn_thres] = 0
dn_adj_combined = t(apply(t(dn_adj), 2, function(x) tapply(x, temp_list$temp_label, sum, na.rm = TRUE)))
dn_adj_combined = apply(dn_adj_combined, 2, function(x) tapply(x, temp_list$temp_label, sum, na.rm = TRUE))
#clean up connectivity to only show important connections, others are really low
dn_adj_combined["MsAHN", colnames(dn_adj_combined) == "DNg02"] = 0
dn_adj_combined[rownames(dn_adj_combined) == "Tect IN", "MsAHN"] = 0
dn_adj_combined[c("DLMn a b", "DLMn c-f", "DVMn 1a-c", "DVMn 2a b", "DVMn 3a b"),] = 0
dn_adj_combined["Tect IN", "DNg02"] = 0

colnames(dn_adj_combined) = sapply(colnames(dn_adj_combined), FUN = function(x) paste0(x, "(", sum(temp_list$temp_label == x), ")"))
rownames(dn_adj_combined) = sapply(rownames(dn_adj_combined), FUN = function(x) paste0(x, "(", sum(temp_list$temp_label == x), ")"))

dn_graph = graph_from_adjacency_matrix(adjmatrix = dn_adj_combined, mode = "directed", weighted = TRUE, diag = FALSE)

	tkplot_graph = tkplot(canvas.height = 800, canvas.width = 800, dn_graph, vertex.size = 20, vertex.label.cex = 1, vertex.color = "grey",
		edge.arrow.size = 0.8, edge.arrow.width = 1.5, curved = TRUE, edge.color = "grey", edge.label = E(dn_graph)$weight, edge.width = 0.1,
		margin = c(0,0,0,0), layout = layout_with_graphopt(dn_graph, spring.constant = 5, charge = 0.01, spring.length = 50), vertex.label.family = "sans")
	readline(prompt="Adjust graph layout by dragging nodes around in tkplot window.\nPress ENTER once done with adjusting node positions: ")
	tkplot_layout = tkplot.getcoords(tkplot_graph)/800
	tk_close(tkplot_graph)

	svglite(paste0("manc_msahn_dng02_tect_in_mn_graph_", Sys.Date(),".svg"), width = 15, height = 15)
	plot.igraph(dn_graph, vertex.size = 7, vertex.label.cex = 1, vertex.color = "grey",
		edge.arrow.size = 1.5, edge.arrow.width = 1, edge.curved = 0.2, edge.label = E(dn_graph)$weight, edge.width = log(E(dn_graph)$weight),
		edge.label.cex = 1, margin = c(0,0,0,0), layout = tkplot_layout, vertex.label.family = "sans", edge.label.family = "sans", rescale = FALSE)
	dev.off()
