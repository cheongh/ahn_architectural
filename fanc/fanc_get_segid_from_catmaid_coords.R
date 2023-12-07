library(googlesheets4)
library(fancr)
library(reticulate)
library(catmaid)
source_python("fanc_realignment.py")

catmaid_url = "https://radagast.hms.harvard.edu/catmaidvnc/"
c_token = "f057a5e08562abdd827f1b342d38733a454e8dc1"
fanc_pid = 13
fanc_voxel_dim_xy = 4.3
fanc_voxel_dim_z = 45

google_write_target_ss = "https://docs.google.com/spreadsheets/d/11Xb6NJxcshw_Xf3sIUovP8NC62JXq0hrjO6EgQeCrIg/edit#gid=0"

fanc_conn = catmaid_login(server = catmaid_url, token = c_token)
ahn_skid = c(MsAHN_L = 237078 , MsAHN_R = 402598, MtAHN_L = 313368, MtAHN_R = 250373)

for (i in 1:length(ahn_skid)) {
	ahn_connectors = catmaid_get_connector_table(ahn_skid[i], direction = "outgoing", get_partner_names = TRUE, conn = fanc_conn, pid = fanc_pid)[c("connector_id", "skid", "treenode_id", "partner_skid", "partner_name")]
	ahn_connectors = ahn_connectors[!is.na(ahn_connectors$partner_skid),]
	ahn_connectors$partner_treenode_id = apply(ahn_connectors[c("partner_skid", "connector_id", "skid")], MARGIN = 1, FUN = function(x) {
		temp = catmaid_get_connector_table(x["partner_skid"], direction = "incoming", conn = fanc_conn, pid = fanc_pid);
		temp[temp$connector_id == x["connector_id"] & temp$partner_skid == x["skid"], "treenode_id"][1]
		})
	ahn_connectors[,c("x","y","z")] = catmaid_get_treenodes_detail(tnids = ahn_connectors$partner_treenode_id, conn = fanc_conn, pid = fanc_pid)[,c("x","y","z")]
	ahn_connectors[,c("x","y")] = ahn_connectors[,c("x","y")]/fanc_voxel_dim_xy
	ahn_connectors["z"] = round(ahn_connectors["z"]/fanc_voxel_dim_z)
	temp_xform_coords = fanc3_to_4(ahn_connectors[,c("x","y","z")], precision = 0.001)
	ahn_connectors$autoid = fanc_xyz2id(temp_xform_coords, rawcoords = TRUE)
	ahn_connectors$fanc4_coords = apply(temp_xform_coords, MARGIN = 1, FUN = function(x) paste(x, collapse = ", "))
	ahn_connectors$ng_link = sapply(ahn_connectors$autoid, FUN = fanc_scene, open = FALSE)
	
	write_sheet(ahn_connectors, ss = google_write_target_ss, sheet = names(ahn_skid)[i])
}

summarize_ahn_downstream = function(ahn_skid) {
	for (i in 1:length(ahn_skid)) {
		temp_ss = read_sheet(google_write_target_ss, sheet = names(ahn_skid)[i])
		summarize_ahn_partners = aggregate(temp_ss$autoid, by = list(temp_ss$autoid), FUN = length)
		colnames(summarize_ahn_partners) = c("autoid", "syn_count")
		summarize_ahn_partners = summarize_ahn_partners[order(summarize_ahn_partners$syn_count, decreasing = TRUE),]
		summarize_ahn_partners$ng_link = sapply(summarize_ahn_partners$autoid, FUN = fanc_scene, open = FALSE)
		write_sheet(summarize_ahn_partners, ss = google_write_target_ss, sheet = paste(names(ahn_skid)[i], "summary", sep = "_"))
	}
}

fanc_update_autoid = function(ss, sheet, colname, write_range) {
	input_ss = read_sheet(ss, sheet)
	if (!(colname %in% colnames(input_ss))) stop("fanc_update_autoid: target sheet not found in input spreadsheet")
	input_ss[colname][is.na(input_ss[colname])] = "0"
	input_ss[colname] = unlist(input_ss[colname])
	#autoid = sapply(unlist(input_ss[colname][!is.na(input_ss[colname])]), FUN = function(x) fanc_latestid(x))
	autoid = vector(mode = "character", length = length(input_ss[[colname]]))
	autoid[1] = fanc_latestid(input_ss[[1, colname]])
	j = 1
	for (i in 2:length(input_ss[[colname]])) {
		if (input_ss[[i, colname]] %in% input_ss[[colname]][1:(i-1)]) {
			autoid[i] = autoid[1:(i-1)][match(input_ss[[i, colname]], input_ss[[colname]][1:(i-1)])]
		} else {
			autoid[i] = fanc_latestid(input_ss[[i, colname]])
		}
		if (j %% 100 == 0) range_write(ss, data.frame(autoID = autoid), sheet, write_range, reformat = FALSE)
		j = j + 1
	}
	range_write(ss, data.frame(autoID = autoid), sheet, write_range, reformat = FALSE)
}

autoid_target_ss = "https://docs.google.com/spreadsheets/d/1jd54V7H1f7RAQr9o03nf_Foxv5TxTIYrEj2V13PWBn0/edit#gid=1765223037"
autoid_target_sheet = c("MsAHN-l (KB)", "MsAHN-r (JR)", "MtAHN-l (AP)", "MtAHN-r (AC)")
#this is downstream
autoid_target_range = "P1"
for (i in 1:4) fanc_update_autoid(autoid_target_ss, autoid_target_sheet[i], "Final NG ID", autoid_target_range)

#this is upstream
autoid_target_ss = "https://docs.google.com/spreadsheets/d/19adhRzBUiiERoLX2PWxv-to0TcLIpCJZkQoKV-DAbQc/edit#gid=887186462"
autoid_target_sheet = c("250374 KCD MtAHN-r", "313369 KCD MtAHN-l", "205793 KCD MsAHN-r", "237079 KCD MsAHN-l")
autoid_target_range = "P1"
for (i in 1:4) fanc_update_autoid(autoid_target_ss, autoid_target_sheet[i], "NG Segment ID", autoid_target_range)

autosynapse_read_ss = "https://docs.google.com/spreadsheets/d/1jd54V7H1f7RAQr9o03nf_Foxv5TxTIYrEj2V13PWBn0/edit#gid=655526454"
autosynapse_sheets = c("MsAHN-l (KB)", "MsAHN-r (JR)", "MtAHN-l (AP)", "MtAHN-r (AC)")
autoseg_from_catmaid_ss = "https://docs.google.com/spreadsheets/d/1EYwtLjiZzBq0hae7aYlNu4_V5f7BqeQWy6EbikophR4/edit#gid=1436073928"
autoseg_from_catmaid_sheets = c("MsAHN_L_summary", "MsAHN_R_summary", "MtAHN_L_summary", "MtAHN_R_summary")

fanc_match_old_data = function(read_ss, read_sheets, write_ss, write_sheet) {
	data = data.frame()
	for (i in read_sheets) {
		temp_ss = read_sheet(read_ss, i)
		if (!("Neuron Class abb" %in% colnames(temp_ss))) stop("fanc_match_old_data: target column not found in sheet")
		data = rbind(data, temp_ss[c("Neuron Class abb", "Initials", "autoID")])
	}
	data = data[!is.na(data$autoID),]
	for (i in write_sheet) {
		write_target_data = read_sheet(write_ss, i)
		write_target_data$matched_celltype = sapply(write_target_data$autoid, function(x) {
			num_unique_types = length(unique(data$"Neuron Class abb"[data$autoID == x & !(is.na(data$"Neuron Class abb"))]));
			if (num_unique_types == 1) {
				data$"Neuron Class abb"[match(x, data$autoID)]
			} else if (num_unique_types > 1) {
				return("inconsistent_label")
			} else {return("")}})
		write_target_data$matched_initials = sapply(write_target_data$autoid, function(x) {
			if(length(data$Initials[data$autoID == x]) > 0) {
				paste(unique(data$Initials[data$autoID == x]), collapse = "/")
			} else {return("")}})
		write_sheet(write_target_data, write_ss, i)
	}
}

fanc_match_old_data(autosynapse_read_ss, autosynapse_sheets, autoseg_from_catmaid_ss, autoseg_from_catmaid_sheets)

sheets_to_sample = c("MsAHN_L", "MsAHN_R", "MtAHN_L", "MtAHN_R")
base_catmaid_url = "https://radagast.hms.harvard.edu/catmaidvnc/?pid=13"
replace_ng_coords = function(ng_scene, coords) {
	gsub('voxelCoordinates%22:\\[\\d*\\.?\\d+,\\d*\\.?\\d+,\\d*\\.?\\d+\\]',
		paste0("voxelCoordinates%22:[", coords, "]"),
		ng_scene)
}

construct_catmaid_url = function(skid, treenode, x, y, z, url) {
	x = x*fanc_voxel_dim_xy
	y = y*fanc_voxel_dim_xy
	z = z*fanc_voxel_dim_z
	return(paste0(url, "&zp=", z, "&yp=", y, "&xp=", x, "&tool=tracingtool&active_skeleton_id=", skid, "&active_node_id=", treenode, "&sid0=10&s0=0"))
}

sample_fanc_downstream_ss = function(ss, sheets) {
	for (i in sheets) {
		temp_sheet = read_sheet(ss, i)
		temp_output = temp_sheet[sample(1:dim(temp_sheet)[1], 50, replace=FALSE),]
		temp_output$ahn_partner = i
		temp_output$ng_link = mapply(temp_output$ng_link, temp_output$fanc4_coords, FUN = replace_ng_coords)
		temp_output$catmaid_link = construct_catmaid_url(skid = temp_output$partner_skid, treenode = temp_output$partner_treenode_id, x = temp_output$x, y = temp_output$y, z = temp_output$z, url = base_catmaid_url)
		write_sheet(temp_output, ss, paste(i, "sample"))
	}
}

sample_fanc_downstream_ss(autoseg_from_catmaid_ss, sheets_to_sample)


#https://radagast.hms.harvard.edu/catmaidvnc/?pid=13&zp=49185&yp=536349.75&xp=225294.19999999998&tool=tracingtool&active_skeleton_id=342349&active_node_id=8108590&sid0=10&s0=0
#voxelCoordinates%22:[48848.1719,114737.2109,2690]
#%22zoomFactor%22:8%7D,%22perspectiveZoom%22:2230.6094


#fanc_ahn_downstream_ss = read_sheet("https://docs.google.com/spreadsheets/d/19adhRzBUiiERoLX2PWxv-to0TcLIpCJZkQoKV-DAbQc/edit#gid=829557905", sheet = "AHNs Wing-Raven")
#temp_coords = sapply(fanc_ahn_downstream_ss$"Coord (Paste in NG)", FUN = function(x) strsplit(x, split = ", "))
#temp_coords = lapply(temp_coords, FUN = function(x) if (all(is.na(x))) c(0,0,0) else x)
#temp_coords = t(matrix(unlist(temp_coords), nrow = 3, ncol = length(temp_coords)))
#temp_xform_coords = fanc3_to_4(temp_coords, precision = 0.01)
#fanc_ahn_downstream_ss$autoid = fanc_xyz2id(temp_xform_coords, rawcoords = TRUE)
#fanc_ahn_downstream_ss$autoid_manual_match = (fanc_ahn_downstream_ss$"NG Segment ID" == fanc_ahn_downstream_ss$autoid)
#fanc_ahn_downstream_ss$fanc3_to_4_coords = apply(temp_xform_coords, MARGIN = 1, FUN = function(x) paste(x, collapse = ", "))

