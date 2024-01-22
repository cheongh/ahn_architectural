library(neuprintr)
library(malevnc)
box::use(mancfuns/han)

syn_thres=3

msahn = c(12536, 13926)
names(msahn) = c("MsAHN","MsAHN")#c("MsAHN-R", "MsAHN-L")

#complete list of tect ins, check which have input from msahns
type_1_sorted = c(23510,22289,22472,101206,21307,23628,20774,18519,157981,20590,18048,14502,16970,19315,18370,14625,32328,15337,23485,18657,18517,20537,15586,36805,22521,21381,23462,24945,153911,23191,23638,13060,22220,24442,155196,20586,16638,15788,17055,15113,101682,17945)
tect_ins = type_1_sorted[type_1_sorted %in% neuprint_connection_table(msahn, prepost = "POST", threshold = syn_thres)$partner]
names(tect_ins)=rep('tect_in',length(tect_ins))
dnp54_dng02=neuprintr::neuprint_list2df(neuprintr::neuprint_fetch_custom(cypher=paste0("MATCH (a:Neuron) WHERE a.type IN ['DNg02','DNp54','DNxn085'] RETURN a.bodyId AS bodyid, a.class AS class, a.group AS group, a.type AS type, a.subclass AS subclass, a.synonyms AS synonyms, a.rootSide AS root_side"), timeout=2000))
dns=dnp54_dng02$bodyid
names(dns)=dnp54_dng02$type
names(dns)[names(dns)=="DNxn085"] = "DNp54"
ifmn=neuprintr::neuprint_list2df(neuprintr::neuprint_fetch_custom(cypher=paste0("MATCH (a:Neuron) WHERE a.type IN ['DLMn a, b','DLMn c-f','DVMn 1a-c', 'DVMn 2a, b', 'DVMn 3a, b'] RETURN a.bodyId AS bodyid, a.class AS class, a.type AS type, a.group AS group, a.subclass AS subclass, a.synonyms AS synonyms, a.rootSide AS root_side"), timeout=2000))
ifmns=ifmn$bodyid
names(ifmns)=ifmn$type

all_bodyids = c(msahn, tect_ins, dns, ifmns)

graph = han$graph_from_bodyid(all_bodyids, save_name = "msahn_dng02_dnp54_tect_in_ifmn_graph.svg", threshold = syn_thres, by_group = TRUE, custom_groups = names(all_bodyids), append_group_name = FALSE, nt_cutoff = 0.7)
