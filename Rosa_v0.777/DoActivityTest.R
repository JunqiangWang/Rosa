#'@description This function do activity marker test using T test
#'@export
DoActivityMarkerTest<-function(seu,
                               group.by=NULL,
                               assay="Activity",
                               slot="data",
                               test="t.test"
){

  require(dplyr)

  if(is.null(group_by)){stop("please assign the group.by parameter")}


  tmp<-match(group.by, colnames(seu@meta.data))
  labels<-seu@meta.data[,tmp]

  groups<-unique(labels)

  mat<-GetAssayData(object = seu, assay = assay, slot = slot)


  out.list<-list()

  for(i in 1:length(groups)){

    message("performing t test using activity metric on group ", groups[i], "...")

    group.i<-groups[i]

    test.id<-which(labels==group.i)
    control.id<-which(labels!=group.i)

    mat.test<-mat[, test.id]
    mat.control<-mat[, control.id]

    out<-matrix(0, nrow=nrow(mat.test), ncol=5)
    rownames(out)<-rownames(mat.test)
    colnames(out)<-c("p_value", "mean_x", "mean_y", "mean_diff", "t")

    #

    if(test!="t.test"){message("only T test is performed in current version")}

    for(j in 1:nrow(mat.test)){

      t_stats<-t.test(mat.test[j,], mat.control[j,], alternative="greater")
      t_stats$p.value
      t_stats$estimate
      mean_diff<-t_stats$estimate[1]-t_stats$estimate[2]
      mean_diff
     # relative_diff<-abs((t_stats$estimate[1]-t_stats$estimate[2])/t_stats$estimate[1])
     # relative_diff

      #
      out[j,]<-unname(c(t_stats$p.value,t_stats$estimate[1], t_stats$estimate[2], mean_diff, t_stats$statistic))

    }

    out<-as.data.frame(out)

    out$cluster<-rep(group.i, times=nrow(out))

    out$gene<-rownames(out)

    out<-out %>% arrange(desc(mean_diff))

    out

    out.list[[i]]<-out

  }

  names(out.list)<-groups

  df <- do.call(rbind, out.list)

  df$p_value_fdr<-p.adjust(df$p_value, method="fdr")

  return(df)

}


# example
# markers<-DoActivityMarkerTest(seu,
#                              group.by="Class",
#                              assay="Activity",
#                              slot="data",
#                              test="t.test"
# )

# markers

