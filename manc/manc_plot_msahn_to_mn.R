library(malevnc)
library(neuprintr)
library(igraph)
library(ggplot2)
library(RColorBrewer)
library(svglite)
library(googlesheets4)
library(nat)

syn_thres = 10
neurons_of_interest = read.csv("MANC_MsAHN_partners_to_plot.csv")

#type_1_ds = neuprint_connection_table(neurons_of_interest$bodyid[neurons_of_interest$type == "Tect IN"], prepost = "POST", threshold = syn_thres)
type_1_ds = manc_connection_table(neurons_of_interest$bodyid[neurons_of_interest$type == "Tect IN"], partners = 'outputs', threshold = syn_thres, moredetails = T)
type_1_ds = unique(type_1_ds$partner[type_1_ds$class == "Motor neuron"])
type_1_ds = type_1_ds[!is.na(type_1_ds)]

google_mn_subclass_ss = "https://docs.google.com/spreadsheets/d/18xOAxfptGLjnMhrVQKcXOIffSfhjTQCQRoWovTjKdd4/edit#gid=0"
google_mn_subclass_ss = read_sheet(google_mn_subclass_ss, sheet = 1)

type_1_ds_type = ifelse(is.na(google_mn_subclass_ss$synonyms[match(type_1_ds, google_mn_subclass_ss$bodyid)]),
	paste(google_mn_subclass_ss$group[match(type_1_ds, google_mn_subclass_ss$bodyid)], "MN"),
	google_mn_subclass_ss$synonyms[match(type_1_ds, google_mn_subclass_ss$bodyid)])

msahn_adj = neuprint_get_adjacency_matrix(bodyids = c(neurons_of_interest$bodyid, type_1_ds))
msahn_adj[msahn_adj<syn_thres] = 0

colnames(msahn_adj) = ifelse(colnames(msahn_adj) %in% neurons_of_interest$bodyid,
	neurons_of_interest$type[match(colnames(msahn_adj), neurons_of_interest$bodyid)],
	type_1_ds_type[match(colnames(msahn_adj), type_1_ds)])
colnames(msahn_adj) = sapply(colnames(msahn_adj), FUN = function(x) paste0(x,"(", sum(colnames(msahn_adj)==x), ")"))
rownames(msahn_adj) = colnames(msahn_adj)

msahn_adj = t(apply(t(msahn_adj), 2, function(x) tapply(x, colnames(msahn_adj), sum, na.rm = TRUE)))
msahn_adj = apply(msahn_adj, 2, function(x) tapply(x, rownames(msahn_adj), sum, na.rm = TRUE))
	
msahn_graph = graph_from_adjacency_matrix(adjmatrix = msahn_adj, mode = "directed", weighted = TRUE, diag = FALSE)

tkplot_graph = tkplot(msahn_graph, canvas.width = 600, canvas.height = 600,
	edge.arrow.size = 0.8, edge.arrow.width = 1.5, curved = TRUE, edge.label = E(msahn_graph)$weight, edge.width = 2, vertex.color = "cyan",
	margin = c(0,0,0,0), layout = layout_with_graphopt(msahn_graph, spring.constant = 5, charge = 0.01, spring.length = 50), vertex.label.family = "sans")
readline(prompt="Adjust graph layout by dragging nodes around in tkplot window.\nPress ENTER once done with adjusting node positions: ")
tkplot_layout_upstream_dn = tkplot.getcoords(tkplot_graph)
tkplot.close(tkplot_graph)
svglite(paste0("MsAHN_graph_with_wing_IN_MN_", Sys.Date(),".svg"), width = 5, height = 5)
print(plot.igraph(msahn_graph, edge.arrow.size = 0.4,
	edge.arrow.width = 1.5, curved = TRUE, edge.label = E(msahn_graph)$weight, edge.width = log1p(E(msahn_graph)$weight)/2, vertex.color = "cyan",
	margin = c(0,0,0,0), layout = tkplot_layout_upstream_dn, vertex.label.family = "sans", vertex.label.cex = 0.65, edge.label.family = "sans", edge.label.cex = 0.6))
dev.off()

dn_mesh = read_manc_meshes(neurons_of_interest$bodyid[neurons_of_interest$type == "Tect IN"])
#MANC.tissue.surf.scale = MANC.tissue.surf
#MANC.tissue.surf.scale$Vertices = MANC.tissue.surf.scale$Vertices*1000
MANC.surf.scale = MANC.surf
MANC.surf.scale$Vertices = MANC.surf.scale$Vertices*1000
open3d()
userMatrix = matrix(c(0.03377458, 0.98989213, 0.13771242, 0, -0.4494019, 0.1381172, -0.8825836, 0, -0.89268696, -0.03207915, 0.44952551, 0, 0, 0, 0, 1), nrow = 4, ncol = 4)
par3d(windowRect=c(0,0,1920,500),zoom=0.45,userMatrix=userMatrix)
for (i in 1:length(dn_mesh)) wire3d(dn_mesh[[i]],col = i, lit = FALSE)
#plot3d(MANC.tissue.surf.scale, alpha = 0.1)
wire3d(MANC.surf.scale, alpha = 0.02, col = "grey")
rgl.snapshot(paste0("MANC_Tect_IN_", Sys.Date(), ".png"),fmt="png")
rgl.close()


dn_mesh = read_manc_meshes(type_1_ds)
#MANC.tissue.surf.scale = MANC.tissue.surf
#MANC.tissue.surf.scale$Vertices = MANC.tissue.surf.scale$Vertices*1000
MANC.surf.scale = MANC.surf
MANC.surf.scale$Vertices = MANC.surf.scale$Vertices*1000
open3d()
userMatrix = matrix(c(0.03377458, 0.98989213, 0.13771242, 0, -0.4494019, 0.1381172, -0.8825836, 0, -0.89268696, -0.03207915, 0.44952551, 0, 0, 0, 0, 1), nrow = 4, ncol = 4)
par3d(windowRect=c(0,0,1920,500),zoom=0.45,userMatrix=userMatrix)
for (i in 1:length(dn_mesh)) wire3d(dn_mesh[[i]],col = i, lit = FALSE)
#plot3d(MANC.tissue.surf.scale, alpha = 0.1)
wire3d(MANC.surf.scale, alpha = 0.02, col = "grey")
rgl.snapshot(paste0("MANC_Tect_IN_downstream_wing_MN_", Sys.Date(), ".png"),fmt="png")
rgl.close()