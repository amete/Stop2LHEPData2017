# First prepare the data

truth.dir <- "./TruthJobs"
truth.out <- "./CleanedTruthJobs"
comb.out  <- "./CombinedTruthJobs"
untar.files <- FALSE
merge.files <- FALSE
make.df     <- FALSE 
make.plots  <- TRUE

#
# Don't touch below!!!
#

###########
# Loop over all subfolders under truth.dir to find the job outputs.
# Then untar all jobs and move the cleaned files to truth.out
###########
manipulate_files <- function(x, y, n) {
    # Construct the file name
    file.name <- paste(x, y, sep = '/')
    # Check if file exists
    if(!file.exists(file.name)) {
        next        
    }
    # Untar the file - this goes to the workdir!
    untar(file.name)
    # Find the masses
    dsid   <- gsub('.*user.amete.(.*).MGPy8EG.*','\\1',file.name)
    mass.x <- gsub('.*bWN_(.*)_(.*)_MadSpin.*','\\1',file.name)
    mass.y <- gsub('.*bWN_(.*)_(.*)_MadSpin.*','\\2',file.name)
    # Rename
    out.name <- paste(dsid, mass.x, mass.y, n, 'txt', sep = '.')
    file.rename(from = 'StopTwoLepton2016.txt',
                to   = paste(truth.out, out.name, sep = '/'))
                
}

extract_files <- function() {
    folders <- list.dirs(path = truth.dir, recursive = FALSE)
    
    lapply(seq_along(folders),
           function(x) {
               lapply( seq_along(list.files(folders[[x]])),
                       function(y) {
                           manipulate_files(folders[[x]],
                                            list.files(folders[[x]])[[y]],
                                            y)
                       } )
           } )
}
###########

###########
# Loop over the extracted files and merge the subjob outputs
###########
add_files <- function(x, y) {
    # x : DSID
    # y : File
    if( grepl(x,y) ) {
        #print(paste(x,y))
        mx <- gsub('(.*)\\.(.*)\\.(.*)\\.(.*).txt','\\2',y)
        my <- gsub('(.*)\\.(.*)\\.(.*)\\.(.*)\\.txt','\\3',y)
        combined.file <- paste(comb.out,paste(x, mx, my, 'txt', sep = '.'), sep = '/')
        if(!file.exists(combined.file)) {
            #file.create(combined.file)
            my.df <- read.csv(paste(truth.out, y, sep = '/'))
            my.df <- my.df [ ,(names(my.df) %in% c('SR','events'))]
            write.csv(my.df, file = combined.file, row.names=FALSE)
        } else {
            my.df <- read.csv(combined.file)
            my.df_temp <- read.csv(paste(truth.out, y, sep = '/'))
            my.df$events <- my.df$events + my.df_temp$events
            write.csv(my.df, file = combined.file, row.names=FALSE)
        }
    }
}

combine_files <- function() {
    dsids <- 389940:389981
    files <- list.files(path = truth.out, recursive = FALSE)
    
    # Loop over disds and files and add sub-outputs    
    lapply( dsids, function(x) { lapply( files, function(y) { add_files(x,y) } ) } )
}
###########

###########
# Read the combined files to make data frames
###########
make_df <- function(x, y, xsections) {
    # x : DSID
    # y : File
    if( grepl(x,y) ) {
        # Read the cross-section and eff
        xsec <- xsections[xsections$DSID == x,]$CROSSECTION.pb
        filt <- xsections[xsections$DSID == x,]$FILTEREFF
                
        # Read masses
        mx <- gsub('(.*)\\.(.*)\\.(.*)\\.txt','\\2',y)
        my <- gsub('(.*)\\.(.*)\\.(.*)\\.txt','\\3',y)

        # Read the file
        my.df <- read.csv(paste(comb.out, y, sep = '/'))
        total <- as.numeric(my.df[my.df$SR == 'All',]$events)
        # Keep interesting rows
        my.df <- filter(my.df, grepl('SR3Body',SR))
        # Add Extra information
        my.df <- mutate(my.df, DSID = x, 
                        MX = as.numeric(mx), 
                        MY = as.numeric(my),
                        Total      = total,
                        Acceptance = as.numeric(events)/total,
                        Error      = Acceptance*sqrt(1./as.numeric(events) + 
                                                     1./total),
                        CrossSecion = xsec,
                        FilterEff   = filt,
                        NormEvents  = as.numeric(events)*36100.0*CrossSecion*FilterEff/Total)
        # Reorder columns (not super necessary)
        my.df <- my.df[,c(3,4,5,6,9,10,1,2,11,7,8)]

        # Write into file
        combined.file <- 'Fancy_Data.csv'
        if(!file.exists(combined.file)) {
            write.csv(my.df, file = combined.file, row.names=FALSE)
        }
        else {
            current.df <- read.csv(combined.file)
            total.df <- rbind(current.df,my.df)
            write.csv(total.df, file = combined.file, row.names=FALSE)
        }
    }
}

build_dataframes <- function() {
    # Use dplyr and tidyr
    if(!require(dplyr)) {
        install.packages("dplyr")
        require(dplyr)
    }
    if(!require(tidyr)) {
        install.packages("tidyr")
        require(tidyr)
    }
    
    dsids <- 389940:389981
    files <- list.files(path = comb.out, recursive = FALSE)

    # Delete if file exists
    if(file.exists('Fancy_Data.csv')) {
        file.remove('Fancy_Data.csv')
    }
    
    # Load cross-sections
    xsections <- read.csv('MGPy8EG_A14N23LO_TT_bWN.txt')
        
    # Loop over disds and files and make the dataframes
    lapply( dsids, function(x) { lapply( files, function(y) { make_df(x,y,xsections) } ) } )
    
}

###########

# Main
if(untar.files) {
    extract_files()
}

if(merge.files) {
    combine_files()    
}

if(make.df) {
    build_dataframes()
}

# Clean-up the workspace
clean_list <- ls()
clean_list <- clean_list[!grepl("make.plots",clean_list)]
rm(list = clean_list)


make_plot <- function(df,val,region,file) {
    # Acceptance
    if(grepl(val,"Acceptance")) {
        p <- qplot(MX, MY, data = df, color = Acceptance*100)
        title <- "Acceptance [%]"
        p <- p + geom_text(aes(label=round(Acceptance*100,2)), vjust = -1)
        p <- p + scale_colour_gradientn(colours = c('orange','red','black'),
                                        limits = c(0,1))
    } else if (grepl(val,"Reco. Eff.")) {
        p <- qplot(MX, MY, data = df, color = YIELD/NormEvents*100)
        title <- "Reco. Eff. [%]"
        p <- p + geom_text(aes(label=round(YIELD/NormEvents*100,0)), vjust = -1)
        p <- p + scale_colour_gradientn(colours = c('orange','red','black'),
                                        limits = c(0,100))
    }
    p <- p + geom_point(size=5)
    p <- p + labs(title = paste(title,"for",region),
                  x = "Stop Mass [GeV]", y = "LSP Mass [GeV]",
                  color = title)
    p <- p + theme_light()
    p <- p + theme(plot.title = element_text(hjust = 0.5,size=14,face="bold"),
                   axis.title = element_text(size=14,face="bold"),
                   axis.text  = element_text(size=14,face="bold"))
    pdf(file, width = 10, height = 7)
    print(p)
    dev.off()
}

if(make.plots) {
    # ggplot2
    if(!require(ggplot2)) {
        install.packages("ggplot2")
        require(ggplot2)
    }
    # Use dplyr and tidyr
    if(!require(dplyr)) {
        install.packages("dplyr")
        require(dplyr)
    }
    # Use gdata
    if(!require(gdata)) {
        install.packages("gdata")
        require(gdata)
    }   
    
    # First read in the data then subset 
    my.df <- read.csv('Fancy_Data.csv')

    ###########
    # SR3BodyWSF
    ###########
    my.df.SR3BodyWSF <- my.df[grepl('SR3BodyWSF',my.df$SR),]
    # Get the reco
    reco.df.SR3BodyWSF <- read.csv('RecoJobs/raw_reco_yields_srwSF_weighted.txt',sep = "\t")
    reco.df.SR3BodyWSF <- mutate(reco.df.SR3BodyWSF, SR = "SR3BodyWSF")
    # Merge
    my.df.SR3BodyWSF <- merge(my.df.SR3BodyWSF,
                              reco.df.SR3BodyWSF,
                              by = c("DSID","MX","MY","SR"))
    # Write
    # write.fwf(my.df.SR3BodyWSF,file = 'SR3BodyWSF_Summary.txt')
    # Acceptance
    make_plot(my.df.SR3BodyWSF,"Acceptance","SR3BodyWSF","Plots/SR3BodyWSF_Acc.pdf")
    make_plot(my.df.SR3BodyWSF,"Reco. Eff.","SR3BodyWSF","Plots/SR3BodyWSF_RecoEff.pdf")

    ###########
    # SR3BodyWDF
    ###########
    my.df.SR3BodyWDF <- my.df[grepl('SR3BodyWDF',my.df$SR),]
    # Get the reco
    reco.df.SR3BodyWDF <- read.csv('RecoJobs/raw_reco_yields_srwDF_weighted.txt',sep = "\t")
    reco.df.SR3BodyWDF <- mutate(reco.df.SR3BodyWDF, SR = "SR3BodyWDF")
    # Merge
    my.df.SR3BodyWDF <- merge(my.df.SR3BodyWDF,
                              reco.df.SR3BodyWDF,
                              by = c("DSID","MX","MY","SR"))
    # Write
    # write.fwf(my.df.SR3BodyWDF,file = 'SR3BodyWDF_Summary.txt')
    # Acceptance
    make_plot(my.df.SR3BodyWDF,"Acceptance","SR3BodyWDF","Plots/SR3BodyWDF_Acc.pdf")
    make_plot(my.df.SR3BodyWDF,"Reco. Eff.","SR3BodyWDF","Plots/SR3BodyWDF_RecoEff.pdf")

    ###########
    # SR3BodyTSF
    ###########
    my.df.SR3BodyTSF <- my.df[grepl('SR3BodyTSF',my.df$SR),]
    # Get the reco
    reco.df.SR3BodyTSF <- read.csv('RecoJobs/raw_reco_yields_srtSF_weighted.txt',sep = "\t")
    reco.df.SR3BodyTSF <- mutate(reco.df.SR3BodyTSF, SR = "SR3BodyTSF")
    # Merge
    my.df.SR3BodyTSF <- merge(my.df.SR3BodyTSF,
                              reco.df.SR3BodyTSF,
                              by = c("DSID","MX","MY","SR"))
    # Write
    # write.fwf(my.df.SR3BodyTSF,file = 'SR3BodyTSF_Summary.txt')
    # Acceptance
    make_plot(my.df.SR3BodyTSF,"Acceptance","SR3BodyTSF","Plots/SR3BodyTSF_Acc.pdf")
    make_plot(my.df.SR3BodyTSF,"Reco. Eff.","SR3BodyTSF","Plots/SR3BodyTSF_RecoEff.pdf")

    ###########
    # SR3BodyTDF
    ###########
    my.df.SR3BodyTDF <- my.df[grepl('SR3BodyTDF',my.df$SR),]
    # Get the reco
    reco.df.SR3BodyTDF <- read.csv('RecoJobs/raw_reco_yields_srtDF_weighted.txt',sep = "\t")
    reco.df.SR3BodyTDF <- mutate(reco.df.SR3BodyTDF, SR = "SR3BodyTDF")
    # Merge
    my.df.SR3BodyTDF <- merge(my.df.SR3BodyTDF,
                              reco.df.SR3BodyTDF,
                              by = c("DSID","MX","MY","SR"))
    # Write
    # write.fwf(my.df.SR3BodyTDF,file = 'SR3BodyTDF_Summary.txt')    
    # Acceptance
    make_plot(my.df.SR3BodyTDF,"Acceptance","SR3BodyTDF","Plots/SR3BodyTDF_Acc.pdf")
    make_plot(my.df.SR3BodyTDF,"Reco. Eff.","SR3BodyTDF","Plots/SR3BodyTDF_RecoEff.pdf")
    
}