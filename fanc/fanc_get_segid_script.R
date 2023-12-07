source("fanc_get_segid_from_catmaid_coords_funcs.R")

ahn_skid = c(MsAHN_L = 237078 , MsAHN_R = 402598, MtAHN_L = 313368, MtAHN_R = 250373)
google_write_target_ss = "https://docs.google.com/spreadsheets/d/1lQbUnT08_OBMdBSKeFw9rhw7xL4raGmzT683Ybqn-eo/edit#gid=0"

wing_connectors = read_sheet(google_write_target_ss, sheet = "Wing")["Connector ID"]
#get autoids for ahn_skid, match connectors only with those from connectors variable, and write to ss with sheet name from ahn_skid names
fanc_get_autoid_from_catmaid_connectors(ahn_skid = ahn_skid[c("MsAHN_L", "MsAHN_R", "MtAHN_L")], fanc_conn = fanc_conn, google_write_target_ss = google_write_target_ss, connectors = wing_connectors)

#spreadsheet you want to match old data from
autosynapse_read_ss = "https://docs.google.com/spreadsheets/d/1EYwtLjiZzBq0hae7aYlNu4_V5f7BqeQWy6EbikophR4/edit#gid=1994178847"
autosynapse_sheets = c("MsAHN_L_summary", "MsAHN_R_summary", "MtAHN_L_summary", "MtAHN_R_summary")
#write target sheets
autoseg_from_catmaid_sheets = c("MsAHN_L", "MsAHN_R", "MtAHN_L")
fanc_match_old_data(autosynapse_read_ss, autosynapse_sheets, google_write_target_ss, autoseg_from_catmaid_sheets)




