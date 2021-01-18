# Author: Ning Liu
from_densematrix_to_hmplotting <- function(matrix, chr = NULL){
  require(tidyverse)
  require(data.table)
  data <- fread(matrix, sep = '\t') %>%
    as_tibble()
  if(nrow(data) + 3 != ncol(data)){
    stop("You need input with proper format.")
    return(NULL)
  }
  if(length(unique(data$V1)) > 1){
    if(is.null(chr) == TRUE){
      stop("Please indicate one chromosome.")
      return(NULL)
    }
    data <- filter(data, V1 == chr)
  }
  
  name <- data$V3
  data <- data[4:length(data)]
  colnames(data) <- name
  matrix <- as.matrix(data)
  rownames(matrix) <- name
  
  return(matrix)
}



hicHeatmap <- function (hicdata, chrom, chromstart, chromend, max_y = 30, zrange = NULL, colpalette = colorRampPalette(inferno(256)), flip = FALSE) 
{
  require(viridisLite)
  rows = as.numeric(rownames(hicdata))
  cols = as.numeric(colnames(hicdata))
  
  # get the subset of matrix from requested start and end
  hicregion = hicdata[which(rows >= chromstart & rows <= chromend), 
                      which(cols >= chromstart & cols <= chromend)]
  
  # how many bins are in the target matrix
  nbins = nrow(hicregion)
  halfbin = abs(chromstart - chromend)/(2 * nbins)
  hicmx = as.matrix(hicregion)
  if (is.null(zrange) == TRUE) {
    max_z = max(hicmx)
    min_z = min(hicmx)
  }
  else {
    min_z = zrange[1]
    max_z = zrange[2]
  }
  
  # map to color
  
  if (is.null(zrange) == TRUE) {
    breaks <- seq(min(hicmx), max(hicmx), length.out = 100)
  }
  else {
    hicmx[which(hicmx < zrange[1])] = zrange[1]
    hicmx[which(hicmx > zrange[2])] = zrange[2]
    breaks <- seq(zrange[1], zrange[2], length.out = 100)
  }
  
  palettes <- colpalette(length(breaks) - 1)
  coloredhicmx <- as.character(cut(hicmx, breaks, labels = palettes))
  coloredhicmx[is.na(coloredhicmx)] <- colpalette(1)
  colormx <- matrix(coloredhicmx, nrow = nrow(hicmx))
  
  # make canvas
  nf <- c(0, max_y)
  f <- c(-max_y, 0)
  plot(1, 1, xlim = c(chromstart, chromend), xlab = "", 
       ylab = "", type = "n", bty = "n", xaxt = "n", yaxt = "n", 
       ylim = get(ifelse(flip == FALSE, "nf", "f")))
  
  for (rownum in (1:nrow(hicmx))) {
    # by making diamond shape for each bin
    y <- ifelse(flip == FALSE, -0.5, 0.5)
    x <- chromstart + (rownum - 1) * (2 * halfbin)
    for (colnum in (rownum:ncol(hicmx))) {
      x <- x + halfbin
      y <- ifelse(flip == FALSE, y + 0.5, y - 0.5)
      x2plot = c(x - halfbin, x, x + halfbin, x, x - halfbin)
      y2plot = c(y, y + 0.5, y, y - 0.5, y)
      polygon(x2plot, y2plot, col = colormx[colnum, rownum], border = NA)
    }
  }
}

plot_sub_triangle <- function(start,end,res){
  max_y <- ((end - start)/2)/res
  mid <- (end - start)/2 + start
  height <- ((end - start)/(end - start))*max_y
  
  segments(start, 0, end, 0, col = "green", lwd = 2, lty=2)
  segments(start, 0, mid, height, col = "green", lwd = 2, lty = 2) 
  segments(mid, height, end, 0, col = "green", lwd = 2, lty = 2) 
}

plotTADs2HM <- function(tadfile,chr,chromstart,chromend,res,tad_col = "red"){
  max_y <- ((chromend - chromstart)/2)/res
  
  tad <- fread(tadfile) %>%
    as_tibble() %>%
    filter(name == "domain") %>%
    dplyr::select(-name)
  
  tads2plot <- tad %>% 
    filter(chrom == chr) %>% 
    filter(chromEnd >= chromstart & chromEnd <= chromend | 
             chromStart >= chromstart & chromStart <= chromend) %>%
    magrittr::set_colnames(c("chr","start","end")) %>%
    mutate(mid = (end - start)/2 + start, 
           height = ((end - start)/(chromend - chromstart))*max_y) %>%
    dplyr::select(-chr) %>%
    as.matrix()
  
  for(i in seq_len(nrow(tads2plot))){
    if(tads2plot[i,1] < chromstart){
      segments(chromstart,0,tads2plot[i,2], 0, col = tad_col, lwd = 2)
      h = ((tads2plot[i,2] - chromstart)/(chromend - chromstart))*max_y
      segments((tads2plot[i,2]-chromstart)/2+chromstart,h,tads2plot[i,2], 0, col = tad_col, lwd = 2)
    } else if(tads2plot[i,2] > chromend){
      segments(tads2plot[i,1], 0, chromend, 0, col = tad_col, lwd = 2)
      h = ((chromend - tads2plot[i,1])/(chromend - chromstart))*max_y
      segments((chromend-tads2plot[i,1])/2+tads2plot[i,1],h,tads2plot[i,1], 0, col = tad_col, lwd = 2)
    } else {
      segments(tads2plot[i,1], 0, tads2plot[i,2], 0, col = tad_col, lwd = 2)
      segments(tads2plot[i,1], 0, tads2plot[i,3], tads2plot[i,4], col = tad_col, lwd = 2)
      segments(tads2plot[i,3], tads2plot[i,4], tads2plot[i,2], 0, col = tad_col, lwd = 2)
    }
  }
}

heatmapLegend <- function(col = inferno(256), min, max=-min, ruler = FALSE, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(col)-1)/(max-min)
  
  dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  if(isTRUE(ruler)){
    axis(2, ticks, las=1)
  }
  for (i in 1:(length(col)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=col[i], border=NA)
  }
}
