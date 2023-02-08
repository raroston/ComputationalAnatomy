# Multichannel PCA, multiple
# Author: R.A. Roston, Ph.D.
# Date: 2023-01-30

# This script runs a multichannelPCA on a list of transforms and saves the results. 
# Transforms & TotalTransforms were generated with "1_registration.R"

runID <- # set runID
wd <- # set working directory
omit <-  # pattern(s) used to select specific transforms to omit from PCA

transform.type <- c("Total", "Syn1")
masks <- c("head", "limbs", "torso", "wholebody")

setwd(wd)

runPCA <- function(wd = getwd(), transform.type, mask, runID, omit){
  
  library(ANTsR)
  
  createtransformlist <- function(wd, folder, omit){
    transformlist <- paste0(wd, folder, dir(paste0(wd, folder)))
    
    for (i in 1:length(omit)){
      transformlist <- transformlist[-grep(patt = omit[i], transformlist)]
    }
    return(transformlist)
  }
  
  doPCA <- function(dir.PCAresults, transformlist, pca_mask){
    
    tlist <- list()
    for (i in 1:length(transformlist)) {
      tlist[[i]] <- resampleImage(antsImageRead(transformlist[i]), c(0.054,0.054,0.054))
    }
    
    pca <- multichannelPCA(x = tlist, mask = pca_mask, pcaOption = "randPCA")
    
    # SAVE PCA RESULTS
    # save() does not work for ANTs images, so pcaWarps has to be saved separately
    dir.create(dir.PCAresults)
    dir.create(paste0(dir.PCAresults, "/pcaWarps"))
    
    save(objects = transformlist, file = paste0(dir.PCAresults, "/transformlist"))
    save(objects = pca, file = paste0(dir.PCAresults, "/PCA"))
    for(i in 1:length(pca$pcaWarps)){
      antsImageWrite(image = pca$pcaWarps[[i]], filename = paste0(dir.PCAresults, "/pcaWarps/pcaWarp_", i, ".nrrd"))
    }
  }
  
  maskfile <- paste0(wd, "/Masks/Embryo_Atlas_mask_", mask, ".nrrd")
  if (file.exists(maskfile)) {
      
      pca_mask <- resampleImage(antsImageRead(maskfile), c(0.054,0.054,0.054), interpType = 1)
      
      if (transform.type == "Total"){
        folder <- "/TotalTransforms/"
        transformlist = createtransformlist(wd = wd, folder = folder, omit = omit)
        dir.create(paste0(wd, "/PCA/", runID))
        dir.PCAresults <- paste0(wd, "/PCA/", runID, "/", transform.type, "-", mask)
        doPCA(dir.PCAresults = dir.PCAresults, transformlist = transformlist, pca_mask = pca_mask)
      } else if (transform.type == "Syn1") {
        folder <- "/Transforms/"
        transformlist = createtransformlist(wd = wd, folder = folder, omit = omit)
        dir.PCAresults <- paste0(wd, "/PCA/", runID, "/", transform.type, "-", mask)
        doPCA(dir.PCAresults = dir.PCAresults, transformlist = transformlist, pca_mask = pca_mask)
      } else { 
        stop(paste0( "Transform.type '", transform.type, "' is invalid."))
      }
    
  } else { 
    stop(paste0("File for mask '", mask, "' does not exist."))
  }
  
}

for (i in 1:length(transform.type)) {
  for (j in 1:length(masks)) {
    print(paste0("Starting ", transform.type[i], "-", masks[j], " : ", Sys.time()))
    runPCA(transform.type = transform.type[i], 
           mask = masks[j],
           runID = runID)
    print(paste0("Completed ", transform.type[i], "-", masks[j], " : ", Sys.time()))
  }
}




