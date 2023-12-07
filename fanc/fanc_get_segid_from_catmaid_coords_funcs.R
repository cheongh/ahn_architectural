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

fanc_conn = catmaid_login(server = catmaid_url, token = c_token)
ahn_skid = c(MsAHN_L = 237078 , MsAHN_R = 402598, MtAHN_L = 313368, MtAHN_R = 250373)

base_catmaid_url = "https://radagast.hms.harvard.edu/catmaidvnc/?pid=13" #for constructing a catmaid link to a neuron
fanc_dn_type_ss = "https://docs.google.com/spreadsheets/d/1k0RZwgouAdXKO2lR7R1YuDgDgNRxRCG-pDKCNZrUBAk/edit#gid=2091967020"

#gets autoids for downstream partners of a skid, and writes it to a google ss with sheet name same as ahn_skid name
fanc_get_autoid_from_catmaid_connectors = function(ahn_skid, fanc_conn, google_write_target_ss, prepost = "POST", connectors = NULL) {
	return_list = list()
	for (i in 1:length(ahn_skid)) {
		ahn_connectors = catmaid_get_connector_table(ahn_skid[i], direction = ifelse(prepost == "PRE", "incoming", "outgoing"), get_partner_names = TRUE, conn = fanc_conn, pid = fanc_pid)[c("connector_id", "skid", "treenode_id", "partner_skid", "partner_name")]
		if (!is.null(connectors)) ahn_connectors = ahn_connectors[ahn_connectors$connector_id %in% unlist(connectors),]
		ahn_connectors = ahn_connectors[!is.na(ahn_connectors$partner_skid),]
		ahn_connectors$partner_treenode_id = apply(ahn_connectors[c("partner_skid", "connector_id", "skid")], MARGIN = 1, FUN = function(x) {
			temp = catmaid_get_connector_table(x["partner_skid"], direction = ifelse(prepost == "PRE", "outgoing", "incoming"), conn = fanc_conn, pid = fanc_pid);
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
}

summarize_ahn_downstream = function(ahn_skid, google_write_target_ss) {
	dn_type_ss = read_sheet(fanc_dn_type_ss, sheet = "AHNs")
	for (i in 1:length(ahn_skid)) {
		temp_ss = read_sheet(google_write_target_ss, sheet = names(ahn_skid)[i])
		summarize_ahn_partners = aggregate(temp_ss$autoid, by = list(temp_ss$autoid), FUN = length)
		colnames(summarize_ahn_partners) = c("autoid", "syn_count")
		summarize_ahn_partners$dn_type = dn_type_ss$nblast_dn_type[match(summarize_ahn_partners$autoid, dn_type_ss$autoID)]
		summarize_ahn_partners$dn_type[is.na(summarize_ahn_partners$dn_type)] = ""
		summarize_ahn_partners = summarize_ahn_partners[order(summarize_ahn_partners$syn_count, decreasing = TRUE),]
		summarize_ahn_partners$ng_link = sapply(summarize_ahn_partners$autoid, FUN = fanc_scene, open = FALSE)
		write_sheet(summarize_ahn_partners, ss = google_write_target_ss, sheet = paste(names(ahn_skid)[i], "summary", sep = "_"))
	}
}

fanc_update_autoid = function(ss, sheet, colname, write_range) {
	input_ss = read_sheet(ss, sheet)
	if (!(colname %in% colnames(input_ss))) stop("fanc_update_autoid: target sheet not found in input spreadsheet")
	autoid = sapply(unlist(input_ss[colname][!is.na(input_ss[colname])]), FUN = function(x) fanc_latestid(x))
	range_write(ss, data.frame(autoID = autoid), sheet, write_range, reformat = FALSE)
}

fanc_match_old_data = function(read_ss, read_sheets, write_ss, write_sheet) {
	data = data.frame()
	for (i in read_sheets) {
		temp_ss = read_sheet(read_ss, i)
		if (!("matched_celltype" %in% colnames(temp_ss))) stop("fanc_match_old_data: target column not found in sheet")
		data = rbind(data, temp_ss[c("matched_celltype", "matched_initials", "autoid", "notes")])
	}
	data = data[!is.na(data$autoid),]
	for (i in write_sheet) {
		write_target_data = read_sheet(write_ss, i)
		write_target_data$matched_celltype = sapply(write_target_data$autoid, function(x) {
			num_unique_types = length(unique(data$"matched_celltype"[data$autoid == x & !(is.na(data$"matched_celltype"))]));
			if (num_unique_types == 1) {
				unique(data$"matched_celltype"[data$autoid == x & !(is.na(data$"matched_celltype"))])
			} else if (num_unique_types > 1) {
				return("inconsistent_label")
			} else {return("")}})
		write_target_data$matched_initials = sapply(write_target_data$autoid, function(x) {
			if(length(data$matched_initials[data$autoid == x & !(is.na(data$matched_initials))]) > 0) {
				paste(unique(data$matched_initials[data$autoid == x & !(is.na(data$matched_initials))]), collapse = "/")
			} else {return("")}})
		write_target_data$matched_initials[write_target_data$matched_initials == "NA"] = ""
		write_target_data$notes = sapply(write_target_data$autoid, function(x) {
			if(length(data$notes[data$autoid == x & !sapply(data$notes, FUN = is.null)]) > 0) {
				paste(unique(data$notes[data$autoid == x & !sapply(data$notes, FUN = is.null)]), collapse = ";")
			} else {return("")}})
		write_sheet(write_target_data, write_ss, i)
	}
}

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

#sheets_to_sample = c("MsAHN_L", "MsAHN_R", "MtAHN_L", "MtAHN_R")
#sample_fanc_downstream_ss(autoseg_from_catmaid_ss, sheets_to_sample)

