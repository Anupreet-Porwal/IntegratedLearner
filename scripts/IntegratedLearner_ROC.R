#' @param SL_fit_layers SuperLearner object List with list names
#' @param y Outcome vector if not already included in the SL object.
ILWrapper_SL_ROC<-function(SL_fit_layers_list,y )
{
  require(pROC)
  require(ggplot2)
  require(ck37r)
  require(ggpubr)
  #The loop to calculate ROC's and add them as new layers
  cols <- palette()
  p <- ggplot() + ggplot2::labs(
                             x = "False positive % (1 - specificity)",
                             y = "True positive % (sensitivity)")
  q<- ggplot() + ggplot2::labs(
                            x = "Recall",
                            y = "Precision")
  auroc<-vector(mode="numeric",length = length(SL_fit_layers_list))
  for(i in 1:length(SL_fit_layers_list)){
    preds<-SL_fit_layers_list[[i]]$Z
    pred = ROCR::prediction(preds, y)
    perf1 = ROCR::performance(pred, "sens", "spec") 
    perf2 = ROCR::performance(pred, "prec", "rec") 
    auroc = round(ck37r::auc_table(SL_fit_layers_list[[i]],y)$auc, 3)
    auprc = round(ck37r::prauc_table(SL_fit_layers_list[[i]],y)$prauc, 3)
    sens_spec <- data.frame(spec=methods::slot(perf1, "y.values")[[1]],
                            sens=1 - methods::slot(perf1, "x.values")[[1]],
                            data=paste(names(SL_fit_layers_list)[i], " AUROC = ", auroc, sep=""))
    prec_rec <- data.frame(prec=methods::slot(perf1, "x.values")[[1]],
                            rec = methods::slot(perf1, "y.values")[[1]],
                            data=paste(names(SL_fit_layers_list)[i], " AUPR = ", auprc, sep=""))
    
    p = p + 
      geom_line(data=sens_spec, aes(x=sens, y=spec, color=data))
    q = q + geom_line(data=prec_rec, aes(x=rec, y=prec, color=data))
  }
  p +  ggplot2::annotate("segment", x = 0, xend = 1, y = 0, yend = 1)
  figure <- ggpubr::ggarrange(p, q, 
                      labels = c("AUROC", "AUPR"),
                      ncol = 1, nrow = 2)
  figure
  
}
  
    

