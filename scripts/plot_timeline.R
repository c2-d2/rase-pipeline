#!/usr/bin/env Rscript

library(optparse)


###############
# CLI PARSING #
###############

kIsRStudio <- Sys.getenv("RSTUDIO") == "1"

if (kIsRStudio) {
  src.exp <- "08_norwich_6766"
  src.exp <- "01_sparc_01"
  src.exp <- "06_phili_04"
  src.file <- "tests/01_sparc_01.predict.tsv"
} else {
  parser <- OptionParser(usage = "%prog [options] timeline.tsv plot.pdf")
  arguments <- parse_args(parser, positional_arguments = 2)

  opt <- arguments$options

  src.file <- arguments$args[1]
  out.file <- arguments$args[2]

  kWidth <- 4
  kHeight <- 10

  pdf(out.file,
      width = kWidth,
      height = kHeight)
}


#################
# CONFIGURATION #
#################

set.seed(42)

# first and second panels
kFirstMinutes <- 15
kLastHours <- 2

# time of the first snapshot
kFirstSnapshotTime <- as.difftime(c(1), units = "mins")

# remove endpoints more than ... far away
kEndpointFilter=10

kLWD <- 2

kIndicSize <- 0.75
kIndicPos <- -1.2

kYLabDist <- 2.3

par(cex = 0.6)
par(oma = c(4, 4, 0.5, 0.5))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))



#############
# FUNCTIONS #
#############

LoadTimelineData <- function(src.file) {
  df <- read.delim(src.file, header = TRUE)
  df$datetime <-
    as.POSIXct(strptime(df$datetime, "%Y-%m-%d %H:%M:%S"))
  first.datetime <- df$datetime[1] - kFirstSnapshotTime

  df$time.mins <-
    difftime(df$datetime, first.datetime, units = "mins")

  # is it really sorted? should be..
  stopifnot(sort(df$datetime) == df$datetime)

  # remove too distant end points
  while(diff(tail(df, 2)$time) >= kEndpointFilter){
    df=head(df,-1)
  }

  df$inv.time.mins <- max(df$time.mins) - df$time.mins

  df
}

DfToAnts <- function(dataframe) {
  cols <- colnames(dataframe)
  antcols <- cols[grepl("_cat", cols)]
  ants <- gsub("_cat", "", antcols)
  ants <- ants[ants != "CHL"]
  ants <- ants[ants != "chl"]
  ants
}

#
# plots vertical ablines
#
TimeAblines <- function(x) {
  abline(
    v = x,
    lty = 2,
    col = "grey",
    lwd = 2
  )
}

#
# plots horizontal ablines
#
ThresholdAbline <- function(y) {
  abline(h = c(y),
         lty = 1,
         col = "grey")
}


#
# set margins for a subfigure
#
margin <- function(i) {
  if (i == 1) {
    par(mar = c(0, 0, 0, 0))
  }
  else{
    par(mar = c(0, 8, 0, 0))
  }
}

#
# plots curve for nb of reads
#

PlotReads <- function(i) {
  margin(i)
  if (i == 1) {
    par(bty = "[")
    plot(
      df1$time.mins,
      df1$read.count / 1000,
      xlim = l.xlim,
      ylim = c(0, max(df$read.count) / 1000),
      type = 'l',
      las = 1,
      xaxt = 'n',
      ylab = NA,
      xlab = NA,
      lwd = kLWD
    )
    mtext(
      "#reads (thousands)",
      side = 2,
      line = kYLabDist,
      cex.lab = 1,
      cex = 0.7,
      las = 3
    )
  }
  else {
    par(bty = "]")
    plot(
      df2$time.mins / 60,
      df2$read.count / 1000,
      xlim = r.xlim,
      ylim = c(0, max(df$read.count) / 1000),
      type = 'l',
      las = 1,
      xaxt = 'n',
      yaxt = 'n',
      ylab = NA,
      xlab = NA,
      lwd = kLWD
    )
  }
  TimeAblines(kVerticalAblines)
}


#
# plots curve for PG score
#

PlotPG <- function(i) {
  margin(i)
  if (i == 1) {
    par(bty = "[")
    plot(
      df1$time.mins,
      2 * df1$PG_score - 1,
      type = 'l',
      xlim = l.xlim,
      ylim = c(0, 1),
      las = 1,
      ylab = NA,
      xlab = NA,
      xaxt = 'n',
      lwd = kLWD
    )


    mtext(
      "PG score",
      side = 2,
      line = kYLabDist,
      cex.lab = 1,
      cex = 0.7,
      las = 3
    )

    mtext(
      "fail",
      side = 2,
      line = kIndicPos,
      cex = kIndicSize,
      at = 0.3
    )

    mtext(
      "pass",
      side = 2,
      line = kIndicPos,
      cex = kIndicSize,
      at = 0.85
    )
  }
  else {
    par(bty = "]")
    plot(
      df2$time.mins / 60,
      2 * df2$PG_score - 1,
      type = 'l',
      xlim = r.xlim,
      ylim = c(0, 1),
      xlab = NA,
      ylab = NA,
      las = 1,
      xaxt = 'n',
      yaxt = 'n',
      lwd = kLWD
    )
  }
  ThresholdAbline(0.6)
  TimeAblines(kVerticalAblines)
}




#
# plots curve for an antibiotic
#
PlotAntibiotic <- function(ant, i, is.last) {
  antcol = paste(ant, "_susc_score", sep = "")
  print(paste(ant, antcol))

  par(bty = "l")
  margin(i)
  if (i == 1) {
    plot(
      df1$time.mins,
      df1[, antcol],
      type = 'l',
      xlim = l.xlim,
      ylim = c(0.0, 1),
      ylab = NA,
      xlab = NA,
      yaxt = "s",
      las = 1,
      xaxt = 'n',
      bty = "[",
      lwd = kLWD
    )

    mtext(
      paste(ant, "susc score"),
      side = 2,
      line = kYLabDist,
      cex.lab = 1,
      cex = 0.7,
      las = 3
    )

    mtext(
      "non-susc",
      side = 2,
      line = kIndicPos,
      cex = kIndicSize,
      at = 0.6 / 2
    )
    mtext(
      "susc",
      side = 2,
      line = kIndicPos,
      cex = kIndicSize,
      at = 0.8
    )
  }
  else{
    plot(
      df2$time.mins / 60,
      df2[, antcol],
      type = 'l',
      xlim = r.xlim,
      ylim = c(0.0, 1),
      yaxt = "n",
      xlab = NA,
      ylab = NA,
      las = 1,
      xaxt = 'n',
      bty = "]",
      lwd = kLWD
    )
  }


  ThresholdAbline(0.6)
  TimeAblines(kVerticalAblines)


  # last row => plot labels
  if (is.last) {
    if (i == 1) {
      axis(1, lwd = 0.5)

      mtext("minutes",
            side = 1,
            line = 2)

    }
    else {
      axis(1, lwd = 0.5)

      mtext("hours",
            side = 1,
            line = 2)

    }
  }
}


############
# PLOTTING #
############

#src.file <- paste("./input/snapshots/", src.exp, ".tsv", sep = "")

df <- LoadTimelineData(src.file)

df1 <- df[df$time.min <= kFirstMinutes,]
df2 <- df[df$inv.time.min <= 60 * kLastHours,]

last.min <- max(df$time.mins)

l.xlim <- c(0, kFirstMinutes)
r.xlim <- c(last.min / 60 - kLastHours, last.min / 60)
kVerticalAblines <- c(1, 5, last.min / 60)

ants <- DfToAnts(df)

par(mfrow = c(length(ants) + 2, 2), tcl = -0.5)


################
# READS
################

for (i in c(1, 2)) {
  PlotReads(i)
}

################
# PG SCORE
################

for (i in c(1, 2)) {
  PlotPG(i)
}

################
# RESISTANCE
################

last.ant <- tail(ants, 1)
for (ant in ants) {
  is.last <- ant == last.ant
  for (i in c(1, 2)) {
    PlotAntibiotic(ant, i, is.last)
  }
}

if (!kIsRStudio) {
  dev.off()
}
