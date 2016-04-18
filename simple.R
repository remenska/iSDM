#!/usr/bin/Rscript
# Fetch command line arguments
myArgs <- commandArgs(trailingOnly = TRUE)
#cat(myArgs)
# Convert to numerics
nums = as.numeric(myArgs)
# cat will write the result to the stdout stream

 newDef <- function(a,b)
 {
     x = runif(10,a,b)
     mean(x)
 }

sayHello <- function(){
   print('hello')
}

sayHello()
cat("maximum is:", max(nums),"\n")

