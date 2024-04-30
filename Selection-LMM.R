#Purpose: Differential selection through Linear Mixed Models 
#Author: Daniela Fernandez Salinas
#Version: 7
#Date: October 17 2022
#Last modification: 
#Usage Examples:
#Dependencies: R:
	#lme4,dyplr,optparse
#Developer Notes:

# 	Packages
require(lme4)
require(dplyr)
require(optparse)
require(tidyverse)

# 	Functions
read_time<-function(name,file){
    t0<-Sys.time()
    temp<-read.table(file, sep = '\t', header = TRUE, row.names=1)
    tf<-Sys.time()-t0
    banner<-paste(file, "...stored in: ",name,". Reading time: ", round(tf,3),"minutes.", sep="")
    print(banner)
    assign(name,temp, envir=parent.frame())
    rm(temp)   
}

selection_model<-function(name,covType,expressed,covF=NULL,covR=NULL){
#Defining regression variables (predictor and covariates)
	#Predictor variable
	x<-as.factor(meta_matrix[,name])
	
	#Adding covariates
	if(!is.null(covType)){
		#Notation for covariates with random effect
		if(covType=="1" || covType=="0"){
		covRandom <- paste(paste("(1|meta_matrix$",covR,")",sep=""),collapse=" + ")
		}
		if(covType=="1"){
		covs <- paste(" ", covRandom, sep="")
		}
		#Notation for covariates with fixed effect
		if(covType=="-1" || covType=="0"){
	    covs<- paste(paste("meta_matrix$", covF, sep=""), collapse=" + ")
	  	}

	  	#Using covariates with multiple effects
	  	if(covType=="0"){
	  	covs <- paste(covs, covRandom, sep=" + ")
	  	}
	}
#Perform test and control regression per peak
 	results <- apply(expressed,1,function(y){
 		if(!is.null(covType)){
	 	#Differential Model with Covariates
	 	  	form1 <- as.formula(paste("y ~ x +", covs))
	 	  	print(form1)
	 	#Control Model
	 	  	form0 <- as.formula(paste("y ~", covs))
	 	  	if(covType=="-1"){
	 	  		lm1 <- lm(form1)
	 			lm0 <- lm(form0)
	 	  	}else{
	 	  		lm1 <- lmer(form1, REML=FALSE)
	 			lm0 <- lmer(form0, REML=FALSE)
	 	  	}	  	
 	  	#Perform likelihood ratio test using anova function
 			anv <- anova(lm0, lm1)
 		#Output Stats
 			c(summary(lm1)$coefficients[2,], anv$P[2])
 		}
 		else{
 		#Simple Differential Model with No Covariates
 			form1<-as.formula("y ~ x")
 			lm1<-lm(form1)
 		#Output Stats
 			c(summary(lm1)$coefficients[2,])
 		}  	
 	})           
 	results <- t(results)
 	if(covType=="-1"){
 		colnames(results)<-c("x1","Std","x1_t","Pr(>|t|)","p_value")
 	}else{
 		colnames(results)<-c("x1","Std","x1_t","p_value")
 	}
	if(nrow(results)!=nrow(expressed)){
	 	print("WARNING: Input and output lines differ!")
	}
return(results) 	
}


################### 	MAIN 	########################
#Parameter Definition
parser=list(
	make_option(c("--data","-d"),help="Main input containing counts matrix.",metavar="FILE"),
	make_option("--simulation",action="store_true",help="Use if you don't intend to run the entire data set. (Only first 100 lines will be used)",default=FALSE),
	make_option("--oDir",help="Directory where to write the output.",default=".",metavar="/PATH/", type="character"),
	make_option("--oFile",help="Name for the output file.",default="Differential-Output", metavar="FILENAME", type="character"),
	make_option(c("--name","-n"),help="Use to run single analysis mode. Indicate the name of the (column) in your metamatrix for which to select.",type="character",metavar="STRING"),
	make_option(c("--meta","-m"),help="Binary meta_matrix.",type="character",metavar="FILE"),
	make_option(c("--random","-r"),help="Add column names used for covariates with random effect."),
	make_option(c("--fixed","-f"),help="Add column names used for covariates with random effect."),
	make_option("--no_covars",action="store_true",help="Use to run a model without covariates.",default=FALSE),
	make_option(c("--all","-a"),action="store_true",help="Use to test all rows in the counts matrix without filtering."),
	make_option("--number",help="Use to test all instances that have a value greater than 1 in at least N samples."),
	make_option("--sorting_column",help="Provide the name of the column in your metamatrix that matches the column names of the counts matrix for sorting.", default="Library", metavar="STRING"),
	make_option(c("--columns","-c"),type="character",help="Use to indicate columns to keep from your counts matrix.", metavar="FILE",default=NULL)
)
args=parse_args(OptionParser(option_list=parser))

#Required Parameters Check
if(is.character(args$data)==FALSE){
	print("ERROR: NO counts matrix provided.")
	stop()
}
if(is.character(args$meta)==FALSE){
	print(args$meta)
	print("ERROR: No meta matrix provided.")
	stop()
}
if(is.character(args$name)==FALSE){
	print("ERROR: No column indicated for selection (name).")
	stop()
}

#Printing Parameters
print("-----------DIFFERENTIAL SELECTION------------")
print(paste0("Counts Matrix:",args$data))
print(paste0("Writting output to:",args$oDir,"/",args$oFile))
print(paste0("File to meta_matrix: ",args$m))
print(paste0("Sorting meta_matrix by: ",args$sorting_column))
print(paste0("Performing selection for: ",args$name))

if(args$simulation){
	print("Running simulation MODE")
}

fixed<-NULL
random<-NULL
covars<-0

if(args$no_covars){
	covars<-NULL
}else{
	if(is.character(args$r)){
		covars<-covars+1
		random<-args$r
	}
	if(is.character(args$f)){
		covars<-covars-1
		if(grepl(",",args$f)){
			fixed<-str_split(args$f,",",simplify=T)
			#print(paste(paste("meta_matrix$", fixed, sep=""), collapse=" + "))
		}else{
			fixed<-args$f
		}
	}
}



#Reading Input
print("Reading input...")
meta_matrix<-read.table(args$meta,header=T,sep='\t')
print("Meta matrix successfully read!")
print("Reading counts...")
read_time("logcounts",args$data)

#Subseting Counts Matrix
if(!is.null(args$columns)){
  print(paste("Subsetting columns from",args$columns,sep=" "))
  selected_columns<-read.table(args$columns)
  logcounts <- logcounts %>% select(selected_columns$V1)
  #meta_matrix<-meta_matrix[which(meta_matrix$sorting_column %in% selected_columns$V1),]
}

#Sorting Metamatrix
if(ncol(logcounts)!=sum(colnames(logcounts) %in% meta_matrix[,args$sorting_column])){
	print("ERROR: Counts column names don't match the meta_matrix.")
	print("Levels in sorting_column:\n")
	print(meta_matrix[,args$sorting_column])
	print("Columns in counts:")
	print(colnames(logcounts))
	stop()
}else{
	row.names(meta_matrix)<-meta_matrix[,args$sorting_column]
	meta_matrix<-meta_matrix[colnames(logcounts),]
	#print(meta_matrix)
}

#Simulation Mode
if(args$simulation){
	logcounts<-logcounts[1:100,]
	print("Subsetting first 100 peaks...")
}

#Defining Output Path
output_path<-args$oDir
if(!dir.exists(output_path)){
 dir.create(output_path)
 print(paste("Created output directory:",output_path,sep=""))
}

if(isTRUE(args$a)){
	log_filtered<-logcounts
	print(paste0("Testing all instances in the counts matrix ",nrow(log_filtered)))
}else{

	minimum<-as.integer(args$number)
	print(paste0("Filtering instances with a value greater than 1 in at least ",minimum," columns..."))
	log_filtered<- logcounts[which(rowSums(logcounts>1)>=minimum),]
	print(paste0("Testing ",nrow(log_filtered)," rows."))
	
}
print(meta_matrix)
print(head(log_filtered))
#Running Regression
if(is.null(covars)){
	print("Running regression without covariates...")
	model_results<-selection_model(name=args$name,expressed=log_filtered,covType=NULL)
	print(model_results)
}
if(!is.null(covars)){
	model_results<-selection_model(name=args$name,expressed=log_filtered,covR=random,covF=fixed,covType=as.character(covars))
}

#Writing output file
filename_regression<-paste(output_path,"/",args$oFile,"_mixedmodel.txt", sep="")
print(filename_regression)
write.table(model_results, filename_regression, row.names=T,col.names=T,quote=F, sep = "\t")
 
print("Done!")
