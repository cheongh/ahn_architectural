library(googlesheets4)
library(fancr)

google_write_target_ss = "https://docs.google.com/spreadsheets/d/19adhRzBUiiERoLX2PWxv-to0TcLIpCJZkQoKV-DAbQc/edit#gid=887186462"
ss_tabs_to_read = c("250374 KCD MtAHN-r", "313369 KCD MtAHN-l", "205793 KCD MsAHN-r", "237079 KCD MsAHN-l")
fanc_dn_type_ss = "https://docs.google.com/spreadsheets/d/1k0RZwgouAdXKO2lR7R1YuDgDgNRxRCG-pDKCNZrUBAk/edit#gid=2091967020"
fanc_dn_type_ss = read_sheet(fanc_dn_type_ss, sheet = "AHNs")
fanc_dn_type_ss$autoID = unlist(fanc_dn_type_ss$autoID)

summarize_ahn_upstream = function(ss, ss_tab) {
	for (i in ss_tab) {
		temp_ss = read_sheet(google_write_target_ss, sheet = i)
		summarize_ahn_partners = aggregate(unlist(temp_ss$autoID), by = list(unlist(temp_ss$autoID)), FUN = length)
		colnames(summarize_ahn_partners) = c("autoid", "syn_count")
		summarize_ahn_partners$class = sapply(summarize_ahn_partners$autoid, function(x) {if (length(unique(temp_ss$"Cell type"[temp_ss$autoID == x & !(is.na(temp_ss$"Cell type"))])) == 1) 
			{temp_ss$"Cell type"[match(x, temp_ss$autoID)]} else {
			return("inconsistent_label")}})
		summarize_ahn_partners$cell_type = fanc_dn_type_ss$nblast_dn_type[match(summarize_ahn_partners$autoid, fanc_dn_type_ss$autoID)]
		summarize_ahn_partners = summarize_ahn_partners[order(summarize_ahn_partners$syn_count, decreasing = TRUE),]
		summarize_ahn_partners$ng_link = sapply(summarize_ahn_partners$autoid, FUN = fanc_scene, open = FALSE)
		write_sheet(summarize_ahn_partners, ss = ss, sheet = paste("summary", i, "HC", sep = " "))
	}
}

summarize_ahn_upstream(ss = google_write_target_ss, ss_tab = ss_tabs_to_read)

inconsistent_label_review_sheet = "inconsistent_label"

correct_ahn_upstream_labels = function(ss, ss_tab, label_review_sheet) {
	label_review_dat = read_sheet(google_write_target_ss, sheet = label_review_sheet)
	for (i in ss_tab) {
		temp_ss = read_sheet(google_write_target_ss, sheet = i)
		temp_ss$autoID = unlist(temp_ss$autoID)
		temp_ss$"Cell type" = ifelse(temp_ss$autoID %in% label_review_dat$autoID,
			label_review_dat$"cell_type (abbreviated)"[match(temp_ss$autoID, label_review_dat$autoID)],
			temp_ss$"Cell type")
		#write_sheet(temp_ss, ss = ss, sheet = i)
		range_write(temp_ss["Cell type"], ss = ss, sheet = i, range = paste0(intToUtf8(64 + which(colnames(temp_ss) == "Cell type")), 1), reformat = FALSE)
	}
}

correct_ahn_upstream_labels(ss = google_write_target_ss, ss_tab = ss_tabs_to_read, label_review_sheet = inconsistent_label_review_sheet)