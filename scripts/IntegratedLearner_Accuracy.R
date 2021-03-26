
#' Title Accuracy Plot for Integrated learner
#'
#' @param SL_fit Nested list of Super Learner fits, 
#' first layer for SL.library, second layer for omics layer
#' @param method_vector Labels of SL.library 
#' @param layers_vector Labels of omics layers
#'
IntegratedLearner_Accuracy<-function(SL_fit,method_vector,layers_vector)
{
  Accuracy<-vector("numeric")
  Methods<- vector()
  for(i in 1:length(method_vector)) {
    Methods<-c(Methods, rep(method_vector[i], length(layers_vector)))
  }
Layers <- rep(layers_vector , length(method_vector))
for(i in 1:length(SL_fit)){
  for(j in 1:length(SL_fit[[i]])) {
    Accuracy[(i-1)*length(SL_fit[[i]]) + j]= 1-SL_fit[[i]][[j]]$cvRisk
  }
}
data <- data.frame(Methods,Layers,Accuracy)
ggplot(data, aes(fill=Layers, y=Accuracy, x=Methods)) + 
  geom_bar(position="dodge", stat="identity")
}


