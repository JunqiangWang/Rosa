#' 
#'@description 
#'@parameter regul A list of regulatory networks
#'@author Junqiang Wang
#'@export 
NetworkIntegration<-function(regul,
                               network.integration.method="maxLikelihood"
                               ){
  
  require(dplyr)
  
  # get all the regulators 
  
  regulators <- unique(unlist(lapply(regul, function(x) names(x)), use.names = FALSE))
  
  unlist(regulators) %>% unname() %>% sort() %>% unique() -> regulators
  
  # merge the regulons, and generate the new regul
  # regulators<-regulators[1:2]
  
 l<-lapply(regulators, FUN=function(regulator, regul){
    
    regul<-unname(regul)
    
   tfmode<-lapply(regul, FUN=function(x, regulator){
      
      x<-x[[regulator]]$tfmode

    }, regulator=regulator)
   
      tfmode<-unlist(tfmode)
   
      
      likelihood<-lapply(regul, FUN=function(x, regulator){
        
        x<-x[[regulator]]$likelihood
        
      }, regulator=regulator)
      
 likelihood<-unlist(likelihood)
 
l<-list()

l$tfmode=tfmode
l$likelihood=likelihood

return(l)
      
  }, regul=regul)
 
 names(l)<-regulators
 
 return(l)
 
}


  #----
  #
  #----

  
  
  
  