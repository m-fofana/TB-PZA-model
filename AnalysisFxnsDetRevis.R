# Mariam Fofana
# November 2015
# Multi-strain TB model("DST" project)
# Functions for analysis of model output from deterministic model

Convert3D <- function(x, n, yrs) {
  # Function to convert trajectory data from 2D (when read from csv) to 3D array
  # When converted to csv, data in 3rd dimension are appended as add'l columns
  # This function selects appropriate columns and moves them back to 3rd dim
  # Args are: x--array to be converted
  #           n--length of data vector for each year
  #           yrs--number of years of data (= # of years simulated + 1)
  # Returns: 3D array w/ 1 row for each trajectory, 1 col for each year
  #           # of slices = # of output variables collected each year
  
  # Error handling
  if (is.numeric(x) != 1) {
    stop ("Array x must be numeric")
  }
  if (is.numeric(n) != 1) {
    stop ("Argument n must be numeric")
  }
  if (is.numeric(yrs) != 1) {
    stop ("Argument yrs must be numeric")
  }
  if ((yrs * n) != ncol(x)) {
    stop ("Check that arguments match desired dimensions")
  }
  if (require(abind) != 1) {
    stop ("Must install and load package 'abind'")
  }
  
  # Convert!
  out3D <- x[, 1:yrs]
  for (i in 2:n) {
    print(paste("i=", i))  # Counter to track progress
    min <- (i - 1) * yrs + 1  # col at which data for yr begin
    max <- i * yrs  # col at which data for yr end
    temp <- x [,min:max]  # all cols containing data for yr
    out3D <- abind(out3D, temp, along=3)
    remove(temp)
  }
  return(out3D)
}  # End of function


PlotDet <- function(outvar, outpop, data1, inityr1, endyr1, 
					          data2=NULL, inityr2=NULL, endyr2=NULL,  
					          plotname=NULL, filename=NULL, savepath=NULL, save=F, 
					          axisyr=NULL, ymax=NULL, col1="gray", col2="red") {
  # Function to plot trajectory data for deterministic output
  # Requires that data1 & data2 be 3D arrays w/ 1 row per sim
  # Assumes that both sets of sims starting from year 0
	# Range of plotting years for data2 should be subset of range for data1
	
  # Args are: data1--data for 1st set of sims
  #           outvar--outcome to be plotted
  #             1 = Prevalence per 100,000
  #             2 = Proportion of active TB cases
  #             3 = Incidence
  #           outpop--Population for which outcomes to be plotted
  #             1 = All active TB cases
  #             2 = MDR TB cases
  #             3 = Pre-XDR TB cases  
  #           inityr1--year to start plotting data1 (sim yr, not calendar yr!)
  #           endyr1--year to end plotting data1 (sim yr, not calendar yr!)
  #           data2--data for 2nd set of sims
  #           inityr2--year to start plotting data2, if =/= inityr1 
  #           endyr2--year to end plotting data2, if =/= endyr1
  #           axisyr--year of upper end of x-axis
  #           plotname--title for plot
  #           filename--file name if saved
  #           savepath--directory path to save plots  
  #           save--Should plots be saved?
  #           ymax--upper limit of y-axis
  #           col1--plot color for data1
  #           col2--plot color for data2

  # Produces plot with x-axis = calendar year, y-axis = outcome
  # Overlapping sets of trajectories (data2 plotted over data1) 
  
  # Error handling
  if (outvar %in% c(1, 2, 3) != 1) {
    stop ("Argument 'outvar' must be '1=Prev', or '2=Prop'")
  }
  if (outpop %in% c(1, 2, 3) != 1) {
    stop ("Argument 'outpop' must be '1=All', '2=MDR', or '3=Pre-XDR'")
  }  
  if (is.numeric(data1) != 1) {
    stop ("Data arrays must be numeric")
  }
  if (is.numeric(inityr1) != 1 || is.numeric(endyr1) != 1) {
    stop ("Arguments 'inityr1' and 'endyr1' must be numeric")
  }
  if (is.null(data2) != 1 && is.numeric(data2) != 1) {
    stop ("Data arrays must be numeric")
  }
  if ((!is.null(inityr2) && !is.numeric(inityr2)) ||
      (!is.null(endyr2) && !is.numeric(endyr2))) {
    stop ("Arguments 'inityr2' and 'endyr2' must be numeric")
  }
  
  if (is.null(inityr2) && !is.null(data2)) {
    inityr2 <- inityr1
  }
  if (is.null(endyr2) && !is.null(data2)) {
    endyr2 <- endyr1
  }

  # Time
  plotyr1 <- data1[, inityr1:endyr1, 1] + 1954
  calyr1 <- inityr1 - 1 + 1954 
  calyr2 <- endyr1 - 1 + 1954   
  
  # Column #s for total active TB cases
  indtot <- c(4, 10, 16, 22, 28, 34, 40, 46, 7, 13, 19, 25, 31, 37, 43, 49)

  # Column numbers for compartments to tally   
  if (outvar < 3) {
    if (outpop == 1) {
      indnum <- c(4, 10, 16, 22, 28, 34, 40, 46, 7, 13, 19, 25, 31, 37, 43, 49)
    } else if (outpop == 2) {
      indnum <- c(10, 28, 34, 46, 13, 31, 37, 49)
    } else if (outpop == 3) {
      indnum <- c(28, 31, 46, 49)
    }
  }

  if (outvar == 3) {
    if (outpop == 1) {
      indnum <- 50
    } else if (outpop == 2) {
      indnum <- 51
    } else if (outpop == 3) {
      indnum <- 52
    }
  }
  
	# Compute outcomes
  if (outvar < 3) {
    out1 <- t(apply(data1[, inityr1:endyr1, indnum], 1, rowSums))
    denom1 <- t(apply(data1[, inityr1:endyr1, indtot], 1, rowSums)) / 100
  }
  
  if (outvar == 3) {
    out1 <- data1[, inityr1:endyr1, indnum]
  }
  
  # Prevalence vs. proportion
  if (outvar != 2) {denom1 <- 1} 
  out1 <- out1 / denom1

  
  # If plotting a 2nd set of trajectories
  if (!is.null(data2)) {
    plotyr2 <- data2[, inityr2:endyr2, 1] + 1954 
    if (outvar < 3) {
      out2 <- t(apply(data2[, inityr2:endyr2, indnum], 1, rowSums))
      denom2 <- t(apply(data2[, inityr2:endyr2, indtot], 1, rowSums)) / 100
    }
    
    if (outvar == 3) {
      out2 <- data2[, inityr2:endyr2, indnum]
    }
    
    if (outvar != 2) {denom2 <- 1} 
    out2 <- out2 / denom2
  }  
  
  # Plot parameters
  ylab1 <- c("", "%", "")
  ylab2 <- c("TB", "MDR", "Pre-XDR")
  ylab3 <- c("Prevalence per 100K", "", "Incidence per 100K")
  if (is.null(ymax)) {
    ymax <- ceiling(max(out1[out1 < Inf], out2[out2 < Inf], na.rm=T) * 1.2)
  }
  if (is.null(axisyr)) {
    axisyr <- calyr2
  }
  par(mai=rep(1.4, 4), xpd=F, mgp=c(3.5, 1, 0))
  plot(0, 0, xlab="Year", cex.sub=0.8, main=plotname, type="n",
       ylab=paste(ylab1[outvar], ylab2[outpop], ylab3[outvar], sep=" "), 
       xlim=c(calyr1, axisyr), ylim=c(0, ymax), cex.axis=1, las=1,
       xaxp=c(calyr1, axisyr, (axisyr - calyr1) / 5)) 

  # Plot
  for (i in 1:nrow(out1)) {
    lines(plotyr1[i, ], out1[i, ], col=col1, lwd=1)
  }
  if (!is.null(data2)) {
    for (i in 1:nrow(out2)) {
      lines(plotyr2[i, ], out2[i, ], col=col2, lwd=1)
    }
  }
  
  # Save
  if (save == T) {
    dev.copy(png, paste(savepath, filename, ".png", sep=""), width=600, height=350)
    dev.off()
  }
}  # End of function
  

vif_func <- function(in_frame, thresh=10, trace=T) {  
  # Creates function for stepwise regression based on
  #   variance inflation factor (VIF)
  # Adapted from: https://gist.github.com/fawda123/4717702#file-vif_fun-r
  
  if (class(in_frame) != 'data.frame') {
    in_frame<-data.frame(in_frame)
  }
  
  # Get initial VIF value for all comparisons of variables
  vif_init <- NULL
  for (val in names(in_frame)) {
    form_in <- formula(paste(val, ' ~ .'))
    vif_init <- rbind(vif_init, c(val, VIF(lm(form_in, data=in_frame))))
  }
  vif_max <- max(as.numeric(vif_init[, 2]))
  
  if (vif_max < thresh) {
    if (trace==T) { 
      # Print output of each iteration
      prmatrix(vif_init, collab=c('var', 'vif'), 
               rowlab=rep('', nrow(vif_init)), quote=F)
      cat('\n')
      cat(paste('All VIFs < ', thresh, ', max VIF ', round(vif_max, 2), sep=''),'\n')
    }
    return(names(in_frame))
  }
  else {
    in_dat<-in_frame
    # Backwards selection of explanatory variables
    # Stops when all VIF values are below 'thresh'
    while (vif_max >= thresh) {
      vif_vals <- NULL
      for (val in names(in_dat)) {
        form_in <- formula(paste(val, ' ~ .'))
        vif_add <- VIF(lm(form_in, data=in_dat))
        vif_vals <- rbind(vif_vals, c(val, vif_add))
      }
      max_row <- which(vif_vals[, 2] == max(as.numeric(vif_vals[, 2])))[1]
      vif_max <- as.numeric(vif_vals[max_row, 2])
      if (vif_max < thresh) {break}
      if (trace==T) { 
        # Print output of each iteration
        prmatrix(vif_vals, collab=c('var', 'vif'), rowlab=rep('', nrow(vif_vals)), quote=F)
        cat('\n')
        cat('removed: ', vif_vals[max_row, 1], vif_max, '\n')
        flush.console()
      }
      
      in_dat<-in_dat[, !names(in_dat) %in% vif_vals[max_row, 1]]
    }
    return(names(in_dat))
  }
}