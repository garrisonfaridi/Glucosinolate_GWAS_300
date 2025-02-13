# Manhatton Function

manhattan_single = function(rs, chr, ps, pval, max.y, outfile, snps2highlight, sig.thresh, sig.thresh.multi, point.scale = 1, gap.size = 1e6){
  
  ## create object storing positions to plot
  
  ps.plot = ps
  
  # loop over chromosomes, assuming they have consecutive numeric names starting at 1
  for (i in 2:length(unique(chr))){
    
    # adjust positions of each SNP on x-axis to add "gap" between each chromosome
    ps.plot[chr >= i] = ps[chr >= i] + max(ps.plot[chr == i-1]) + gap.size
    
  }
  
  # # qq Plot (not used; ignore)
  # 
  # CMplot(data.frame(rs, chr, ps, pval), plot.type = "q",
  #        main = file_path_sans_ext(infilename), memo = file_path_sans_ext(infilename))
  
  ## Manhattan Plot
  
  # convert between -log10 and unscaled p values
  pval[-log10(pval) > max.y] = 10^-max.y
  
  # set filename for output Manhattan plot
  png(filename = outfile, width = 4800, height = 2400)
  
  # plot points for all SNPs
  print( plot(-log10(pval) ~ ps.plot, ylim = c(0,max.y+1), axes = F, xlab = "", ylab = "",  
              pch = 16, col = "gray90", cex = 2.2 * point.scale) )
  
  # plot points for significant SNPs  
  print( points(-log10(pval[pval<sig.thresh]) ~ ps.plot[pval<sig.thresh], 
                pch = 16, col = "gray75", cex = 2.2 * point.scale) )
  
  # plot points for non-significant SNPs in candidate genes
  print( points(-log10(pval[rs %in% snps2highlight & pval > sig.thresh.multi]) ~ 
                  ps.plot[rs %in% snps2highlight & pval > sig.thresh.multi], 
                pch = 16, col = "dodgerblue", cex = 3.5 * point.scale) ) 
  
  # plot points for significant SNPs in candidate genes
  print( points(-log10(pval[rs %in% snps2highlight & pval < sig.thresh.multi & pval != 10^-max.y]) ~ 
                  ps.plot[rs %in% snps2highlight & pval < sig.thresh.multi & pval != 10^-max.y], 
                pch = 16, col = "dodgerblue", cex = 5.0 * point.scale) )
  
  # optionally, if p-values are equal to the y axis limit, make these points larger
  # (can be useful if setting p-values above the limit to equal the limit to avoid stretched axis)
  print( points(-log10(pval[rs %in% snps2highlight & pval == 10^-max.y]) ~ 
                  ps.plot[rs %in% snps2highlight & pval == 10^-max.y], 
                pch = 17, col = "dodgerblue", cex = 6.5 * point.scale) ) 
  
  # draw lines for single GWAS and and multi-GWAS significance thresholds
  print( abline(h = -log10(sig.thresh), col = "red", lty = 2, lwd = 4) )
  print( abline(h = -log10(sig.thresh.multi), col = "black", lty = 2, lwd = 4) )
  
  # plot axes
  print( axis(side = 2, cex.axis = 5, cex.lab = 7, lwd = 5))
  
  dev.off()
  
}