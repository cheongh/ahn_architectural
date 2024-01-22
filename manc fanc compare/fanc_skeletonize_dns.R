library(fancr)
library(fafbseg)
library(reticulate)
library(googlesheets4)
library(malevnc)
library(nat)
#sk = import("skeletor")
#tm = import("trimesh")

#use below to install required python packages
#reticulate::py_install('meshparty', pip = T, pip_options='--upgrade --prefer-binary --user --ignore-installed certifi')

fanc_dn_ss = "https://docs.google.com/spreadsheets/d/1k0RZwgouAdXKO2lR7R1YuDgDgNRxRCG-pDKCNZrUBAk/edit#gid=2091967020"
fanc_dn_data = read_sheet(fanc_dn_ss, "AHNs")
fanc_skel_path = paste(getwd(), "fanc_skels", sep = "\\\\")
if (!dir.exists(fanc_skel_path)) dir.create(fanc_skel_path)
fanc_skel_subpath = paste(getwd(), "fanc_skels", "dn", sep = "\\\\")
if (!dir.exists(fanc_skel_subpath)) dir.create(fanc_skel_subpath)

#tweak these settings
with_fanc(download_neuron_obj(unlist(fanc_dn_data$autoID), save.obj = fanc_skel_subpath))

#fanc_dn_skels = with_fanc(meshparty_skeletonize(unlist(fanc_dn_data$autoID[1])))
#skeletor function bugged
#fanc_dn_skels = with_fanc(skeletor(obj = fanc_skel_subpath, method = "wavefront", waves = 1, radius = FALSE))
skeletor_wavefront = function(obj, waves, step_size, heal = TRUE) {
	#obj can be file name or directory name
	if (!file_test("-f", obj)) #check if file
		filelist = list.files(path = obj)
	else {
		filelist = unlist(strsplit(obj, "[\\/]"))
		filelist = filelist[length(filelist)]
	}
	swclist = neuronlist()
	#check that all have valid extensions by strsplit()
	filelist = filelist[sapply(filelist, FUN = function(x) {temp = unlist(strsplit(x, "\\."));
		if (length(temp) < 2 | temp[length(temp)] != "obj") return(FALSE) else return(TRUE)})]
	reticulate::py_run_string("import skeletor as sk")
	reticulate::py_run_string("import trimesh as tm")
	for (i in filelist) {
		reticulate::py_run_string(sprintf('mesh = tm.load_mesh("%s", process = False)',
			paste(obj, i, sep = '\\\\')))
		reticulate::py_run_string(sprintf("mesh = sk.pre.fix_mesh(mesh = mesh, remove_disconnected = 5, inplace = True)"))
		reticulate::py_run_string(sprintf("swc = sk.skeletonize.by_wavefront(mesh = mesh, waves = %s, step_size = %s, progress=False)",
            waves, step_size))
		swc = reticulate::py$swc$swc
		#swc = skel$swc
		colnames(swc) = c("PointNo","Parent","X","Y","Z","W")
		neuron = nat::as.neuron(swc)
		#do we need some way to add segmentID to neuron
		if(heal){
			neuron = suppressMessages(nat::stitch_neurons_mst(x = neuron, threshold = Inf, k = 10L))
		}
		nl = neuronlist(neuron)
		names(nl) = strsplit(i, "\\.")[[1]][1]
		swclist = append(swclist, nl)
	}
	return(swclist)
}

swcs = skeletor_wavefront(fanc_skel_subpath, waves = 1, step_size = 1)
write.neurons(swcs, dir = fanc_skel_subpath, format = "swc")
