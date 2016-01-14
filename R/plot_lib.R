volcanoPlot = function(	dat,
						key_columns = c(pvalue="dropout_rate_diff_pvalue", fdr="dropout_rate_diff_fdr", effect_size="dropout_rate_diff_estimate"),
						fdr_cutoff=0.1,
						relative_effect_size = TRUE,
						key_targets_minus = NA,
						key_targets_plus = NA,
            show_key_target_symbols = TRUE,
						plottitle="",
						xlabel="",
						ylabel="",
						fdr_label="",
						xlimits=NA,
						ylimits=NA,
						y_scale_by=2,
						plot_colours = c(insignificant="#DDDDDD", significant="#888888", key_target_minus="royalblue", key_target_plus="darkgreen")
						) {

	library(ggplot2)

	dat$pvalue = dat[,key_columns["pvalue"]]
	dat$fdr = dat[,key_columns["fdr"]]
	dat$effect_size = dat[,key_columns["effect_size"]]

	dat$log_pvalue = -log10(dat$pvalue)
	dat$log_fdr = -log10(dat$fdr)
	dat$fdr_status = ifelse(dat$fdr <= fdr_cutoff, "significant", "insignificant")

	highlight_minus = ifelse(length(key_targets_minus) < 2 & is.na(key_targets_minus), FALSE, TRUE)[1]
	highlight_plus = ifelse(length(key_targets_plus) < 2 & is.na(key_targets_plus), FALSE, TRUE)[1]

	minVal = min(dat$effect_size)
	maxVal = max(dat$effect_size)

	# Adjust x axis for relative dropout rate
	if(relative_effect_size) {
		min_effect = floor(min(dat$effect_size))
		max_effect = ceiling(max(dat$effect_size))
		xlimits = c(min_effect+1, max_effect-1)
		effect_breaks = (min_effect+1):(max_effect-1)
		effect_labels = c(min_effect:-2, "same", 2:max_effect)
		dat$effect_size = ifelse(dat$effect_size < 0, dat$effect_size+1, dat$effect_size-1)
#		print(effect_breaks)
#		print(effect_labels)
	} else {
		# For example, in expression vs. essentiality analyses, the effect size
		# is substantially smaller than (-1;1), so the x-axis ticks need to be
		# set differently.
		if(minVal > -1 & maxVal < 1) {
			min_effect = floor(minVal*10)/10
			max_effect = ceiling(maxVal*10)/10
			xlimits = c(min_effect, max_effect)
			effect_breaks_minus = rev(-(0:(-minVal*10)/10))
			effect_breaks_plus = (0:(maxVal*10))/10
			effect_breaks = c(effect_breaks_minus, effect_breaks_plus[-1])
#			effect_breaks = ((minVal*10):(maxVal*10))/10
			effect_labels = as.character(round(effect_breaks, digits=1))
		} else {
			min_effect = floor(minVal)
			max_effect = ceiling(maxVal)
			xlimits = c(min_effect, max_effect)
			effect_breaks = min_effect:max_effect
			effect_labels = as.character(effect_breaks)			
		}
	}

	if(highlight_minus) {
		dat$key_target_minus = ifelse(dat$symbol %in% key_targets_minus, "key_target_minus", "unknown")
		dat$fdr_status = ifelse(dat$key_target_minus == "key_target_minus", "key_target_minus", dat$fdr_status)
		dat$key_target_effect_minus = ifelse(dat$key_target_minus == "key_target_minus", dat$effect_size, NA)
		dat$key_target_log_pvalue_minus = ifelse(dat$key_target_minus == "key_target_minus", dat$log_pvalue, NA)
	}

	if(highlight_plus){
		dat$key_target_plus = ifelse(dat$symbol %in% key_targets_plus, "key_target_plus", "unknown")
		dat$fdr_status = ifelse(dat$key_target_plus == "key_target_plus", "key_target_plus", dat$fdr_status)
		dat$key_target_effect_plus = ifelse(dat$key_target_plus == "key_target_plus", dat$effect_size, NA)
		dat$key_target_log_pvalue_plus = ifelse(dat$key_target_plus == "key_target_plus", dat$log_pvalue, NA)
	}
	
	###################################################################
	# Generate Gene level volcano plot
	###################################################################
	
	temp = table(dat$fdr_status)

	color_levels = c()

	if(!is.na(temp["insignificant"])) color_levels = c(color_levels, 1)
	if(!is.na(temp["significant"])) color_levels = c(color_levels, 2)
	if(highlight_minus) color_levels = c(color_levels, 3)
	if(highlight_plus) color_levels = c(color_levels, 4)

	color_labels = c("insignificant", "significant", "key_target_minus", "key_target_plus")
	color_labels = color_labels[color_levels]

	colors_fdr = plot_colours[color_levels]
	dat$fdr_status = factor(dat$fdr_status, levels=color_labels)

	maxY = min(floor(max(dat$log_pvalue)), 60)

	ylimits = c(0,ceiling(max(dat$log_pvalue)))
	ybreaks = seq(0,maxY,by=y_scale_by)
	yticksAll = c(1,
				expression(10^{-1}), expression(10^{-2}), expression(10^{-3}), expression(10^{-4}), expression(10^{-5}),
				expression(10^{-6}), expression(10^{-7}), expression(10^{-8}), expression(10^{-9}), expression(10^{-10}),
				expression(10^{-11}), expression(10^{-12}), expression(10^{-13}), expression(10^{-14}), expression(10^{-15}),
				expression(10^{-16}), expression(10^{-17}), expression(10^{-18}), expression(10^{-19}), expression(10^{-20}),
				expression(10^{-21}), expression(10^{-22}), expression(10^{-23}), expression(10^{-24}), expression(10^{-25}),
				expression(10^{-26}), expression(10^{-27}), expression(10^{-28}), expression(10^{-29}), expression(10^{-30}),
				expression(10^{-31}), expression(10^{-32}), expression(10^{-33}), expression(10^{-34}), expression(10^{-35}),
				expression(10^{-36}), expression(10^{-37}), expression(10^{-38}), expression(10^{-39}), expression(10^{-40}), 
				expression(10^{-41}), expression(10^{-42}), expression(10^{-43}), expression(10^{-44}), expression(10^{-45}), 
				expression(10^{-46}), expression(10^{-47}), expression(10^{-48}), expression(10^{-49}), expression(10^{-50}),
				expression(10^{-51}), expression(10^{-52}), expression(10^{-53}), expression(10^{-54}), expression(10^{-55}),
				expression(10^{-56}), expression(10^{-57}), expression(10^{-58}), expression(10^{-59}), expression(10^{-60}))
	yticks = yticksAll[ybreaks+1]


  maxValue = max(dat[dat$fdr <= fdr_cutoff, "pvalue"], na.rm=T)
	fdr_line = ifelse(!is.na(maxValue), -log10(maxValue), NA)
	fdr_label = fdr_label
	fdr_label_x = min(dat$effect_size) + 0.1*(max(dat$effect_size)-min(dat$effect_size))
	fdr_label_y = fdr_line-0.5
	#label_minus = "<= More Essential\nWith Increased Expression"

	range_x = range(dat$effect_size)
	range_y = range(dat$pvalue)

	if(minVal > -1 & maxVal < 1) {
		jitter_x = abs(0.03*(range_x[2]-range_x[1]))
	} else {
		jitter_x = abs(0.15*(range_x[2]-range_x[1]))		
	}
	
	jitter_y = abs(0.15*(range_y[2]-range_y[1]))

	pl = ggplot(data=dat, aes(x=effect_size, y=log_pvalue))
	pl = pl + geom_point(aes(color=fdr_status), size=3)
	pl = pl + theme_bw()
	pl = pl + xlab(xlabel) + ylab(ylabel) + labs(title=plottitle)
  pl = pl + scale_color_manual(values=colors_fdr)
  pl = pl + theme(
#    plot.title = element_text(colour="black", size=18, face="bold"),
    axis.line=element_line(size = 0.2, linetype = 'solid', colour="grey"),
#     axis.text.x=element_text(size=14, colour="black"),
#     axis.text.y=element_text(size=14, colour="black"),
#     axis.title.x=element_text(size=16, colour="black", face="bold"),
#     axis.title.y=element_text(size=16, colour="black", face="bold", angle=90),
    strip.text.x=element_text(colour="black", size=14, face="bold"),
    strip.text.y=element_text(colour="black", size=14, face="bold"),
    panel.background=element_blank(),
    #  legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

  if(!is.na(fdr_line)) {
    pl = pl + geom_hline(yintercept=fdr_line, linetype="dashed")
    pl = pl + annotate("text", label=fdr_label, x=fdr_label_x, y=fdr_label_y, size=6, fontface="italic")
  }

	if(!is.na(xlimits[1])) pl = pl + scale_x_continuous(limits=xlimits, breaks=effect_breaks, labels=effect_labels)
	if(!is.na(ylimits[1])) pl = pl + scale_y_continuous(limits=ylimits, breaks=ybreaks, labels=yticks)

	if(highlight_minus) {
		pl = pl + geom_point(aes(x=key_target_effect_minus, y=key_target_log_pvalue_minus, color=fdr_status), size=5)

    if(show_key_target_symbols) {
      pl = pl + geom_text(aes(label=symbol, x=key_target_effect_minus, y=key_target_log_pvalue_minus, color=fdr_status), position=position_jitter(h=jitter_x, w=jitter_y), size=7, fontface="bold")      
    }
	}

	if(highlight_plus) {
		pl = pl + geom_point(aes(x=key_target_effect_plus, y=key_target_log_pvalue_plus, color=fdr_status), size=5)

		if(show_key_target_symbols) {
      pl = pl + geom_text(aes(label=symbol, x=key_target_effect_plus, y=key_target_log_pvalue_plus, color=fdr_status), position=position_jitter(h=jitter_x, w=jitter_y), size=7, fontface="bold")
		}
	}

	return(pl)
}
