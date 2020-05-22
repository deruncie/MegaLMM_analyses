# code taken and adapted from the phytools package: https://cran.r-project.org/web/packages/phytools/index.html
my_phylo2map = function (x, type = c("phylogram", "direct"), ...)
{
  type <- type[1]
  if (class(x) == "phylo.to.map") {
    tree <- x$tree
    map <- x$map
    coords <- x$coords
  }
  else stop("x should be an object of class \"phylo.to.map\"")
  if (hasArg(xlim))
    xlim <- list(...)$xlim
  else xlim <- map$range[1:2]
  if (hasArg(ylim))
    ylim <- list(...)$ylim
  else ylim <- map$range[3:4]
  if (hasArg(fsize))
    fsize <- list(...)$fsize
  else fsize <- 1
  if (hasArg(split))
    split <- list(...)$split
  else split <- c(0.4, 0.6)
  if (hasArg(psize))
    psize <- list(...)$psize
  else psize <- 1
  if (hasArg(cex.points)) {
    cex.points <- list(...)$cex.points
    if (length(cex.points) == 1)
      cex.points <- c(0.6 * cex.points, cex.points)
  }
  else cex.points <- c(0.6 * psize, psize)
  if (hasArg(mar))
    mar <- list(...)$mar
  else mar <- rep(0, 4)
  if (hasArg(asp))
    asp <- list(...)$asp
  else asp <- 1
  if (hasArg(ftype))
    ftype <- list(...)$ftype
  else ftype <- "reg"
  ftype <- which(c("off", "reg", "b", "i", "bi") == ftype) -
    1
  if (!ftype)
    fsize = 0
  if (hasArg(from.tip))
    from.tip <- list(...)$from.tip
  else from.tip <- FALSE
  if (hasArg(colors))
    colors <- list(...)$colors
  else colors <- "red"
  if (length(colors) == 1)
    colors <- rep(colors[1], 2)
  if (length(colors) == 2 && type == "phylogram") {
    colors <- matrix(rep(colors, nrow(coords)), nrow(coords),
                     2, byrow = TRUE)
    rownames(colors) <- rownames(coords)
  }
  else if (is.vector(colors) && (length(colors) == Ntip(tree))) {
    # recover()
    COLS <- matrix("red", nrow(coords), 2, dimnames = list(rownames(coords)))
    for (i in 1:length(colors)) COLS[which(rownames(COLS) ==
                                             names(colors)[i]), 1:2] <- colors[i]
    colors <- COLS
  }
  if (hasArg(direction))
    direction <- list(...)$direction
  else direction <- "downwards"
  if (hasArg(pch))
    pch <- list(...)$pch
  else pch <- 21
  if (length(pch) == 1)
    pch <- rep(pch, 2)
  if (hasArg(lwd))
    lwd <- list(...)$lwd
  else lwd <- c(2, 1)
  if (length(lwd) == 1)
    lwd <- rep(lwd, 2)
  if (hasArg(lty))
    lty <- list(...)$lty
  else lty <- "dashed"
  if (hasArg(pts))
    pts <- list(...)$pts
  else pts <- TRUE
  if (type == "phylogram") {
    if (x$direction == "downwards" && direction == "rightwards") {
      cat("\"phylo.to.map\" direction is \"downwards\" but plot direction has been given as \"rightwards\".\n")
      cat("Re-optimizing object....\n")
      cc <- aggregate(coords, by = list(rownames(coords)),
                      mean)
      cc <- matrix(as.matrix(cc[, 2:3]), nrow(cc), 2, dimnames = list(cc[,
                                                                         1], colnames(cc)[2:3]))
      tree <- minRotate(tree, cc[, 1])
    }
    else if (x$direction == "rightwards" && direction ==
             "downwards") {
      cat("\"phylo.to.map\" direction is \"rightwards\" but plot direction has been given as \"downwards\".\n")
      cat("Re-optimizing object....\n")
      cc <- aggregate(coords, by = list(rownames(coords)),
                      mean)
      cc <- matrix(as.matrix(cc[, 2:3]), nrow(cc), 2, dimnames = list(cc[,
                                                                         1], colnames(cc)[2:3]))
      tree <- minRotate(tree, cc[, 2])
    }
  }
  if (type == "phylogram") {
    if (direction == "downwards") {
      if (!ftype)
        ylim <- c(ylim[1], ylim[2] + 0.03 * diff(ylim))
      ylim <- c(ylim[1], ylim[2] + split[1]/split[2] *
                  (ylim[2] - ylim[1]))
    }
    else if (direction == "rightwards") {
      if (!ftype)
        xlim <- c(xlim[1] - 0.03 * diff(xlim), xlim[2])
      xlim <- c(xlim[1] - split[1]/split[2] * (xlim[2] -
                                                 xlim[1]), xlim[2])
    }
  }
  if (all(mar == 0))
    mar <- mar + 0.01
  plot.new()
  par(mar = mar)
  plot.window(xlim = xlim, ylim = ylim, asp = asp)
  map(map, add = TRUE, fill = TRUE, col = "gray95", mar = rep(0,
                                                              4))
  if (type == "phylogram") {
    # recover()
    cw = tree
    # cw <- reorder(tree, "cladewise")
    if (!is.binary(cw))
      cw <- multi2di(cw)
    n <- Ntip(cw)
    if (direction == "downwards") {
      dx <- abs(diff(xlim))
      rect(xlim[1] - 1.04 * dx, ylim[2] - split[1] * (ylim[2] -
                                                        ylim[1]), xlim[2] + 1.04 * dx, ylim[2], col = "white",
           border = "white")
      pdin <- par()$din[2]
      dy = diff(ylim)
      segments(xlim[1]-.03*dx,ylim[2]-split[1]*dy,xlim[1]-.03*dx,ylim[2])
      breaks = seq(0,1,by=.2)
      segments(rep(xlim[1]-.03*dx,length(breaks)),ylim[2]-split[1]*dy + breaks*dy*split[1],rep(xlim[1]-.04*dx,length(breaks)),ylim[2]-split[1]*dy + breaks*dy*split[1])
      # recover()
      text(rep(xlim[1]-.04*dx,length(breaks)),ylim[2]-split[1]*dy + breaks*dy*split[1],labels = sprintf('%0.1f',breaks),adj=1.5,cex=.5)
      text(xlim + dx/2,ylim[2],labels = list(...)$main,adj=c(.5,1))
      # print('asdf')
      # text(rep(xlim[1]-.04*dx,length(breaks)),ylim[2]-split[1]*dy + breaks*dy*split[1],labels = breaks,adj=1.5,cex=.5)

      sh <- (fsize * strwidth(paste(" ", cw$tip.label,
                                    sep = "")) + 0.3 * fsize * strwidth("W")) * (par()$din[1]/par()$din[2]) *
        (diff(par()$usr[3:4])/diff(par()$usr[1:2]))
      cw$edge.length <- cw$edge.length/max(nodeHeights(cw)) *
        (split[1] * (ylim[2] - ylim[1]) - max(sh))
      pw <- reorder(cw, "postorder")
      x <- vector(length = n + cw$Nnode)
      x[cw$edge[cw$edge[, 2] <= n, 2]] <- 0:(n - 1)/(n -
                                                       1) * (xlim[2] - xlim[1]) + xlim[1]
      nn <- unique(pw$edge[, 1])
      for (i in 1:length(nn)) {
        xx <- x[pw$edge[which(pw$edge[, 1] == nn[i]),
                        2]]
        x[nn[i]] <- mean(range(xx))
      }
      Y <- ylim[2] - nodeHeights(cw)
      coords <- coords[, 2:1]
      for (i in 1:nrow(coords)) {
        tip.i <- which(cw$tip.label == rownames(coords)[i])
        lines(c(x[tip.i], coords[i, 1]), c(Y[which(cw$edge[,
                                                           2] == tip.i), 2] - if (from.tip) 0 else sh[tip.i],
                                           coords[i, 2]), col = colors[i, 1], lty = lty,
              lwd = lwd[2])
      }
      points(coords, pch = pch, cex = cex.points[2], bg = colors[,2])
      for (i in 1:nrow(Y)) lines(rep(x[cw$edge[i, 2]],
                                     2), Y[i, ], lwd = lwd[1], lend = 2)
      for (i in 1:cw$Nnode + n) lines(range(x[cw$edge[which(cw$edge[,
                                                                    1] == i), 2]]), Y[which(cw$edge[, 1] == i), 1],
                                      lwd = lwd[1], lend = 2)
      for (i in 1:n) {
        if (ftype)
          text(x[i], Y[which(cw$edge[, 2] == i), 2],
               paste(" ", sub("_", " ", cw$tip.label[i]),
                     sep = ""), pos = 4, offset = c(0, 1), srt = -90,
               cex = fsize, font = ftype)
        if (pts)
          points(x[i], Y[which(cw$edge[, 2] == i), 2],pch = 21, bg = colors[cw$tip.label, ][i,2], cex = cex.points[1],col=NA)
      }
      PP <- list(type = "phylogram", use.edge.length = TRUE,
                 node.pos = 1, show.tip.label = if (ftype) TRUE else FALSE,
                 show.node.label = FALSE, font = ftype, cex = fsize,
                 adj = 0, srt = 0, no.margin = FALSE, label.offset = fsize *
                   strwidth(" ")/(par()$usr[2] - par()$usr[1]) *
                   (par()$usr[4] - par()$usr[3]), x.lim = par()$usr[1:2],
                 y.lim = par()$usr[3:4], direction = direction,
                 tip.color = "black", Ntip = Ntip(cw), Nnode = cw$Nnode,
                 edge = cw$edge, xx = x, yy = sapply(1:(Ntip(cw) +
                                                          cw$Nnode), function(x, y, z) y[match(x, z)],
                                                     y = Y, z = cw$edge))
    }
    else {
      dy <- abs(diff(ylim))
      rect(xlim[1], ylim[1], xlim[1] + split[1] * (xlim[2] -
                                                     xlim[1]), ylim[2], col = "white", border = "red")
      sh <- fsize * strwidth(paste(" ", cw$tip.label, sep = "")) +
        0.2 * fsize * strwidth("W")
      cw$edge.length <- cw$edge.length/max(nodeHeights(cw)) *
        (split[1] * (xlim[2] - xlim[1]) - max(sh))
      pw <- reorder(cw, "postorder")
      y <- vector(length = n + cw$Nnode)
      y[cw$edge[cw$edge[, 2] <= n, 2]] <- 0:(n - 1)/(n -
                                                       1) * (ylim[2] - ylim[1]) + ylim[1]
      nn <- unique(pw$edge[, 1])
      for (i in 1:length(nn)) {
        yy <- y[pw$edge[which(pw$edge[, 1] == nn[i]),
                        2]]
        y[nn[i]] <- mean(range(yy))
      }
      H <- nodeHeights(cw)
      X <- xlim[1] + H
      coords <- coords[, 2:1]
      for (i in 1:nrow(coords)) {
        tip.i <- which(cw$tip.label == rownames(coords)[i])
        lines(c(X[which(cw$edge[, 2] == tip.i), 2] +
                  if (from.tip) 0 else sh[tip.i], coords[i, 1]),
              c(y[tip.i], coords[i, 2]), col = colors[i,
                                                      1], lty = lty, lwd = lwd[2])
      }
      points(coords, pch = pch, cex = cex.points[2], bg = colors[,
                                                                 2])
      for (i in 1:nrow(X)) lines(X[i, ], rep(y[cw$edge[i,
                                                       2]], 2), lwd = lwd[1], lend = 2)
      for (i in 1:cw$Nnode + n) lines(X[which(cw$edge[,
                                                      1] == i), 1], range(y[cw$edge[which(cw$edge[,
                                                                                                  1] == i), 2]]), lwd = lwd[1], lend = 2)

      for (i in 1:n) {
        if (ftype)
          text(X[which(cw$edge[, 2] == i), 2], y[i],
               paste(" ", sub("_", " ", cw$tip.label[i]),
                     sep = ""), pos = 4, offset = 0.1, cex = fsize,
               font = ftype)
        if (pts)
          points(X[which(cw$edge[, 2] == i), 2], y[i],
                 pch = 21, bg = colors[cw$tip.label, ][i,
                                                       2], cex = cex.points[1])
      }
      PP <- list(type = "phylogram", use.edge.length = TRUE,
                 node.pos = 1, show.tip.label = if (ftype) TRUE else FALSE,
                 show.node.label = FALSE, font = ftype, cex = fsize,
                 adj = 0, srt = 0, no.margin = FALSE, label.offset = 0.1,
                 x.lim = par()$usr[1:2], y.lim = par()$usr[3:4],
                 direction = direction, tip.color = "black", Ntip = Ntip(cw),
                 Nnode = cw$Nnode, edge = cw$edge, xx = sapply(1:(Ntip(cw) +
                                                                    cw$Nnode), function(x, y, z) y[match(x, z)],
                                                               y = X, z = cw$edge), yy = y)
    }
    assign("last_plot.phylo", PP, envir = .PlotPhyloEnv)
  }
  else if (type == "direct") {
    phylomorphospace(tree, coords[, 2:1], add = TRUE, label = "horizontal",
                     node.size = c(0, psize), lwd = lwd[2], control = list(col.node = setNames(rep(colors[2],
                                                                                                   max(tree$edge)), 1:max(tree$edge)), col.edge = setNames(rep(colors[1],
                                                                                                                                                               nrow(tree$edge)), tree$edge[, 2])), ftype = c("off",
                                                                                                                                                                                                             "reg", "b", "i", "bi")[ftype + 1], fsize = fsize)
  }
}
