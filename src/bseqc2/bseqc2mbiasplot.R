# vim: set noexpandtab tabstop=2:

pdf(NULL)
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
plotstrand=function() {
	f$methCG=f$numC_CG/(f$numC_CG+f$numT_CG)
	f$methCH=f$numC_CH/(f$numC_CH+f$numT_CH)
	f$tag=factor(f$tag, levels=c('++', '-+', '+-', '--'))
	N=nrow(f)
	pcg=ggplot(data=f, aes(x=position, y=methCG, color=tag)) +
	geom_line() +
	scale_x_continuous(breaks=scales::pretty_breaks(n=min(10, N))) +
	scale_y_continuous(breaks=scales::pretty_breaks(n=min(10, N)), limits=c(0, 1)) +
	xlab(xlab) +
	ylab(ylab) +
	ggtitle('Methylation levels of CG') +
	theme_bw() +
	theme(
		plot.background = element_blank()
		, panel.grid.major = element_blank()
		, panel.grid.minor = element_blank()
		, panel.border = element_blank()
		, plot.title = element_text(hjust=0.5)
		, axis.line = element_line(color='black')
		, axis.text = element_text(color='black')
		, legend.position = 'top'
		, legend.title = element_blank()
		)

	pch=ggplot(data=f, aes(x=position, y=methCH, color=tag)) +
	geom_line() +
	scale_x_continuous(breaks=scales::pretty_breaks(n=min(10, N))) +
	scale_y_continuous(breaks=scales::pretty_breaks(n=min(10, N)), limits=c(0, 1)) +
	xlab(xlab) +
	ylab(ylab) +
	ggtitle('Methylation levels of CH') +
	theme_bw() +
	theme(
		plot.background = element_blank()
		, panel.grid.major = element_blank()
		, panel.grid.minor = element_blank()
		, panel.border = element_blank()
		, plot.title = element_text(hjust=0.5)
		, axis.line = element_line(color='black')
		, axis.text = element_text(color='black')
		, legend.position = 'top'
		, legend.title = element_blank()
		)

	ml=marrangeGrob(list(pcg, pch), nrow=1, ncol=2, top='')
	ggsave(ml, file=outfile, width=width, height=height, useDingbats=F)
}

plotpairs=function() {
	f$end=sprintf('End%s', f$end+1)
	f$methCG=f$numC_CG/(f$numC_CG+f$numT_CG)
	f$methCH=f$numC_CH/(f$numC_CH+f$numT_CH)
	if (pico) {
		f$tag=factor(f$tag
			, levels=c(
				'++,+-', '+-,++', '-+,--', '--,-+', '++,N', '+-,N', 'N,++', 'N,+-', '-+,N', '--,N', 'N,-+', 'N,--'
				, '++,-+', '-+,++', '++,--', '--,++', '+-,--', '--,+-', '+-,-+', '-+,+-', '++,++', '+-,+-', '-+,-+', '--,--'
				)
			)
	} else {
		f$tag=factor(f$tag
			, levels=c(
				'++,+-', '-+,--', '++,N', '-+,N', 'N,+-', 'N,--'
				, '++,--', '--,++', '++,-+', '+-,++', '+-,N', 'N,++', '--,N', 'N,-+', '--,-+', '-+,++', '+-,--', '--,+-', '+-,-+', '-+,+-', '++,++', '+-,+-', '-+,-+', '--,--'
				)
			)
	}
	f=f[order(f$tag),]

	plots=lapply(unique(f$tag)
		, function (tag) {
			subf=f[f$tag==tag,]
			N=nrow(subf)
			pcg=ggplot(data=subf, aes(x=position, y=methCG, color=end)) +
			geom_line() +
			scale_x_continuous(breaks=scales::pretty_breaks(n=min(10, N))) +
			scale_y_continuous(breaks=scales::pretty_breaks(n=min(10, N)), limits=c(0, 1)) +
			xlab(xlab) +
			ylab(ylab) +
			ggtitle('Methylation levels of CG') +
			theme_bw() +
			theme(
				plot.background = element_blank()
				, panel.grid.major = element_blank()
				, panel.grid.minor = element_blank()
				, panel.border = element_blank()
				, plot.title = element_text(hjust=0.5)
				, axis.line = element_line(color='black')
				, axis.text = element_text(color='black')
				, legend.position = 'top'
				, legend.title = element_blank()
				)

			pch=ggplot(data=subf, aes(x=position, y=methCH, color=end)) +
			geom_line() +
			scale_x_continuous(breaks=scales::pretty_breaks(n=min(10, N))) +
			scale_y_continuous(breaks=scales::pretty_breaks(n=min(10, N)), limits=c(0, 1)) +
			xlab(xlab) +
			ylab(ylab) +
			ggtitle('Methylation levels of CH') +
			theme_bw() +
			theme(
				plot.background = element_blank()
				, panel.grid.major = element_blank()
				, panel.grid.minor = element_blank()
				, panel.border = element_blank()
				, plot.title = element_text(hjust=0.5)
				, axis.line = element_line(color='black')
				, axis.text = element_text(color='black')
				, legend.position = 'top'
				, legend.title = element_blank()
				)

			ml=marrangeGrob(list(pcg, pch), nrow=1, ncol=2, top=sprintf('Mapping pair: %s', tag))
			ml
		})
	pdf(file=outfile, width=width, height=height, useDingbats=F)
	print(plots)
	dev.off()
}

f=read.table(infile, header=T, sep='\t', quote="", check.names=F, comment.char='#', stringsAsFactors=F)
if ('end' %in% names(f)) {
	plotpairs()
} else {
	plotstrand()
}
