library(ggplot2)
library(signal)
#library(coin)
library(exactRankTests)

fly_sex = c(M = 'male', F = 'female')

calcium_data_folder = "E:/calcium imaging ahns/msahn epi flight exact sync/calcium measurements"
video_sync_info = read.csv("E:/calcium imaging ahns/msahn epi flight exact sync/video sync info.csv")#, fileEncoding="UTF-8-BOM")

t_after = 0.5 #time after flight start used for finding slope maxima
t_before = 0.5
camera_frame_rate = 100 #Hz #in this experiment, behavior camera is triggered by epifluorescence imaging camera frame-by-frame
epi_frame_rate = camera_frame_rate #Hz

lowpass_butter = butter(4, 5/camera_frame_rate, "low") #5 Hz cutoff
lowpass_filter = function(signal) filtfilt(lowpass_butter, x=signal)

df_f = function(dat, baseline) (dat-baseline)/baseline*100

#convert to takeoff time
video_sync_info$takeoff_frame = (video_sync_info$takeoff_frame-1)/camera_frame_rate



msahn_summary_dat = data.frame()
all_calcium_dat = data.frame()

for (i in unique(video_sync_info$Fly)) {
	if (all(!is.na(video_sync_info$takeoff_frame[video_sync_info$Fly==i]))) {
		for (j in video_sync_info$video[video_sync_info$Fly == i]) {
			takeoff_frame = video_sync_info$takeoff_frame[video_sync_info$Fly == i & video_sync_info$video == j]
			temp_calcium = read.csv(paste0(calcium_data_folder, "/", fly_sex[video_sync_info$sex[video_sync_info$Fly==i][1]], " msahn ", i, "-", j, ".csv"))
			temp_calcium$t = 0:(dim(temp_calcium)[1]-1)/epi_frame_rate
			#reduce data to period 5s before takeoff until 20s after
			temp_calcium$t = temp_calcium$t-takeoff_frame
			temp_calcium = temp_calcium[temp_calcium$t >= -5 & temp_calcium$t <= 20,]
			
			#process ahn, 1 or 2 cells #not the smartest way to do this
			temp_ahn_col = c('Mean1','Mean3')
			temp_ahn_col_2 = c('Mean2','Mean4')
			temp_ahn_mean_col = c('mean1.2','mean3.4')
			temp_filt = c('mean1.2_filt','mean3.4_filt')
			for (k in 1:2) {
				if (temp_ahn_col[k] %in% colnames(temp_calcium)) {
					temp_calcium[temp_ahn_mean_col[k]] = temp_calcium[[temp_ahn_col[k]]] - temp_calcium[[temp_ahn_col_2[k]]]
					#lowpass filter
					temp_calcium_mean = mean(temp_calcium[temp_calcium$t<=-0.5 & temp_calcium$t>=-5, temp_ahn_mean_col[k]]) #in 5s period before takeoff
					temp_calcium[temp_ahn_mean_col[k]] = df_f(temp_calcium[temp_ahn_mean_col[k]], temp_calcium_mean)
					temp_calcium[temp_filt[k]] = lowpass_filter(temp_calcium[[temp_ahn_mean_col[k]]])

					#test fit polynomial curve #fit to 0.5s after first significant rise
					diffdiff = diff(temp_calcium[[temp_filt[k]]], differences=2)
					t_rise = which(diff(sign(diff(temp_calcium[[temp_filt[k]]], differences=3)))==-2)+3 #find local maxima
					t_rise = temp_calcium$t[t_rise] #convert to t
					t_rise = t_rise[t_rise >= -t_before & t_rise <= t_after] #find maxima of diffdiff within a 1 s window around takeoff
					t_rise = t_rise[which.max(diffdiff[match(t_rise,temp_calcium$t)+2])]
					#to avoid blips due to z-motion, move t_rise to the right until we find non-negative fluorescence
					temp_positive_calcium = (temp_calcium$t>=t_rise & temp_calcium[[temp_filt[k]]]>=0)
					t_rise = temp_calcium$t[min(which(temp_positive_calcium))]
	
					
					msahn_summary_dat = rbind(msahn_summary_dat, data.frame(fly=i, rep=j, cell=k, t_diff = t_rise))
					temp_calcium$temp_plot = temp_calcium[[temp_ahn_mean_col[k]]]
					temp_calcium$temp_filt = temp_calcium[[temp_filt[k]]]
					
					all_calcium_dat = rbind(all_calcium_dat, data.frame(fly=i, rep=j, cell=k, t = temp_calcium$t, calcium=temp_calcium[[temp_ahn_mean_col[k]]]))
					
					pdf(paste0("msahn_calcium_vs_takeoff_fly_", i, "_", j, "-", k, ".pdf"), width = 400/72, height = 400/72)
					print(ggplot(temp_calcium, aes(x=t,y=temp_plot)) +
						geom_line() +
						geom_line(aes(x=t,y=temp_filt), color='magenta') +
						xlim(-1, 2) + ylim(-50,400)+
						facet_grid(rows=, space="free_y") +
						geom_vline(xintercept = 0, color="red") + #takeoff
						geom_vline(xintercept = t_rise, color="blue") + #calcium rise
						ylab("dF/F (%)") +
						xlab("t (s)"))
					dev.off()
					png(paste0("msahn_calcium_vs_takeoff_fly_", i, "_", j, "-", k, ".png"), width = 400, height = 400)
					print(ggplot(temp_calcium, aes(x=t,y=temp_plot)) +
						geom_line() +
						geom_line(aes(x=t,y=temp_filt), color='magenta') +
						xlim(-1, 2) + ylim(-50,400)+
						facet_grid(rows=, space="free_y") +
						geom_vline(xintercept = 0, color="red") + #takeoff
						geom_vline(xintercept = t_rise, color="blue") + #calcium rise
						ylab("dF/F (%)") +
						xlab("t (s)"))
					dev.off()
				}
			}
		}
	}
}

all_calcium_dat$t = round(all_calcium_dat$t, digits = 2) #why was this not exact?
temp_mean = aggregate(calcium~t, data=all_calcium_dat, FUN = mean)
pdf(paste0("msahn_calcium_vs_takeoff_all_fly_reps_plot.pdf"), width = 400/72, height = 400/72)
print(ggplot(all_calcium_dat, aes(x=t,y=calcium,group=interaction(fly,cell,rep))) +
	geom_line(size = 0.25, color='grey') +
	geom_line(aes(x=t, y=calcium, group='none'), size = 0.5, color='black', data=temp_mean) +
	xlim(-5, 15) + ylim(-100,500)+
	facet_grid(rows=, space="free_y") +
	geom_vline(xintercept = 0, color="red") + #takeoff
	ylab("dF/F (%)") +
	xlab("t (s)") +
	theme_classic())
dev.off()

pdf(paste0("msahn_calcium_vs_takeoff_all_fly_reps_plot_zoom.pdf"), width = 200/72, height = 200/72)
print(ggplot(all_calcium_dat, aes(x=t,y=calcium,group=interaction(fly,cell,rep))) +
	geom_line(size = 0.25, color='grey') +
	geom_line(aes(x=t, y=calcium, group='none'), size = 0.5, color='black', data=temp_mean) +
	xlim(-0.5, 0.5) + ylim(-50,200)+
	facet_grid(rows=, space="free_y") +
	geom_vline(xintercept = 0, color="red") + #takeoff
	ylab("dF/F (%)") +
	xlab("t (s)") +
	theme_classic())
dev.off()
					
					
msahn_summary_dat$fly = factor(msahn_summary_dat$fly)

png(paste0("msahn_calcium_vs_takeoff_diff.png"), width = 200, height = 400)
print(ggplot(msahn_summary_dat, aes(x="MsAHN",y=t_diff, color=fly)) +
	geom_boxplot(aes(x="MsAHN",y=t_diff, color=NA), outlier.shape = NA) +
	geom_jitter(width=0.4, height = 0) +
	ylab("calcium rise start - first wing movement timing (s)") +
	xlab("AHN"))
dev.off()


msahn_summary_dat$ahn = "MsAHN"
all_summary_dat = msahn_summary_dat

png(paste0("ahn_calcium_vs_takeoff_diff.png"), width = 300, height = 600)
print(ggplot(all_summary_dat, aes(x=ahn,y=t_diff, color=fly)) +
	geom_boxplot(aes(x=ahn,y=t_diff, color=NA), outlier.shape = NA) +
	geom_jitter(width=0.4, height = 0) +
	ylab("calcium rise start - first wing movement timing (s)") +
	xlab("AHN")) +
	theme(text = element_text(size = 20)) 
dev.off()

#summarize by fly means, then do t-test to test if significantly below zero
per_fly_summary = aggregate(t_diff~fly+ahn, data = all_summary_dat, FUN = mean, na.action = na.omit)
pdf(paste0("per_fly_ahn_calcium_vs_takeoff_diff_curve_fit.pdf"), width = 250/72, height = 600/72)
print(ggplot(per_fly_summary, aes(x=ahn,y=t_diff)) +
	geom_boxplot(aes(x=ahn,y=t_diff), outlier.shape = NA) +
	geom_jitter(width=0.4, height = 0, size = 3) +
	ylab("calcium rise start - first wing movement timing (s)") +
	xlab("AHN") +
	theme(text = element_text(size = 20)) )
dev.off()

#t.test(per_fly_summary$t_diff[per_fly_summary$ahn=="MsAHN"], y = NULL, alternative = "less") #possibly non-normal

wt = wilcox.exact(per_fly_summary$t_diff[per_fly_summary$ahn=="MsAHN"], mu = 0, alternative = "less", exact=TRUE, conf.int = TRUE)
