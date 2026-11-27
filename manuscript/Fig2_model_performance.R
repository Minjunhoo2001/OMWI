#!/usr/bin/env Rscript
# Available panels a, b, c and f. No obsolete study-exclusion or LOOCV inputs are used.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pROC)
  library(patchwork)
  library(scales)
})
options(stringsAsFactors = FALSE)

# 只需要在这里改路径
ROOT_DIR <- Sys.getenv("OMWI_ROOT", getwd())
IN_DIR <- Sys.getenv("OMWI_MODEL_RESULTS", file.path(ROOT_DIR, "data/model_outputs"))
OUT_DIR <- Sys.getenv("OMWI_FIG2_OUT", file.path(ROOT_DIR, "fig2"))
if (!dir.exists(IN_DIR)) stop("Frozen model folder missing: ", IN_DIR)
TABLE_DIR <- file.path(OUT_DIR, "tables")
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)

# 这里全部是绘图代码，不会启动稳定性模型重训。
# ROC: LOOCV & 10-fold CV 相近，采用透明度 + 虚线区分。
INCLUDE_ROC_CI <- TRUE
ROC_BOOT_N <- 250L

COL_HEALTH  <- "#4DBBD5"  # 原版青蓝
COL_DISEASE <- "#E64B35"  # 原版砖红
COL_TRAIN   <- "#DE491E"  # 原版10-fold CV橙色
COL_TEST    <- "#51999F"  # 06.fig2原版外部验证青绿
COL_GREY    <- "#7E7E7E"
COL_RET     <- "#3C5488"  # 右轴、样本保留率：NPG深蓝
COL_BA      <- "#00A087"  # balanced accuracy：原版绿色调
COL_ROC <- c(
  "LOOCV OOF" = COL_GREY,
  "10-fold CV OOF" = COL_TRAIN,
  "External validation" = COL_TEST
)
COL_GROUP <- c("Healthy" = COL_HEALTH, "Non-healthy" = COL_DISEASE)

# 直接复用原始 06.fig2(2).R 中的灰框、白色分面条、刻度线设计。
theme_cell <- function(base_size=10) {
  theme_bw(base_size=base_size, base_family="sans") +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(fill=NA,color="grey50",linewidth=.45),
          axis.line=element_blank(),
          axis.ticks=element_line(color="grey40",linewidth=.4),
          axis.text=element_text(color="black"),
          axis.title=element_text(color="black"),
          plot.title=element_blank(),
          legend.title=element_blank(),
          legend.background=element_blank(),
          legend.key=element_blank(),
          strip.background=element_rect(fill="white",color="grey50",linewidth=.45),
          strip.text=element_text(face="bold",color="black"),
          plot.tag=element_text(face="bold",size=13))
}

save_both <- function(p, base, width, height, dpi=600) {
  ggsave(file.path(OUT_DIR,paste0(base,".pdf")), p,
         width=width,height=height,bg="white",limitsize=FALSE)
  ggsave(file.path(OUT_DIR,paste0(base,".png")), p,
         width=width,height=height,dpi=dpi,bg="white",limitsize=FALSE)
}
write_csv <- function(x, file) {
  fwrite(x,file.path(TABLE_DIR,paste0(file,".csv")),bom=TRUE)
}
read_csv <- function(path) {
  if (!file.exists(path)) stop("缺少所需输入文件：",path)
  # 大部分 CSV 为 UTF8；原始训练表部分文本为 GB18030，优先保留 ASCII 核心字段。
  # fread 的 Latin-1 兼容读取可避免含个别非 UTF8 的非核心元数据导致报错。
  tryCatch(fread(path, data.table=FALSE, check.names=FALSE),
           error=function(e) fread(path, data.table=FALSE, check.names=FALSE,
                                     encoding="Latin-1"))
}
assert_cols <- function(x, names_req, file_label) {
  missing <- setdiff(names_req,names(x))
  if (length(missing)) stop(file_label," 缺少列：",paste(missing,collapse=", "))
}
resolve <- function(rel) read_csv(file.path(IN_DIR, rel))

coef <- resolve("final/species_nonzero_coefficients.csv")
train <- resolve("final/species_train_in_sample_scores.csv")
ext <- resolve("final/external_independent_predictions.csv")
oof <- resolve("summary/nested_repeated10fold_OOF.csv")
selected_params <- resolve("final/species_selected_parameters.csv")
external_metrics <- resolve("final/external_independent_metrics.csv")
cv_each <- resolve("summary/nested_repeated10fold_each.csv")
assert_cols(coef,c("feature","coefficient"),"coefficients")
assert_cols(train,c("sample_id","y","OMHI","prob_healthy"),"training")
assert_cols(ext,c("sample_id","disease","OMHI","prob_healthy"),"external")
assert_cols(oof,c("sample_id","y_true","repeat_id","OMHI","prob_healthy"),"nested CV")
train$y_true <- as.integer(train$y == 1)
ext$y_true <- as.integer(tolower(trimws(ext$disease)) == "healthy")
stopifnot(nrow(train)==1077L,nrow(ext)==716L,nrow(oof)==10770L, all(oof$alpha==0.1),
          !anyDuplicated(train$sample_id),!anyDuplicated(ext$sample_id),
          !anyDuplicated(oof[,c("sample_id","repeat_id")]),
          setequal(unique(oof$sample_id),train$sample_id),!length(intersect(train$sample_id,ext$sample_id)))
for (d in list(train,ext,oof)) stopifnot(all(is.finite(d$OMHI)),all(is.finite(d$prob_healthy)),max(abs(plogis(d$OMHI)-d$prob_healthy))<1e-6)
cv_df <- oof %>% group_by(sample_id) %>% summarise(y_true=first(y_true),n_repeats=n_distinct(repeat_id),prob_healthy=mean(prob_healthy),.groups="drop") %>% mutate(OMHI=qlogis(pmin(pmax(prob_healthy,1e-12),1-1e-12)))
stopifnot(all(cv_df$n_repeats==10L),all(cv_df$y_true==train$y_true[match(cv_df$sample_id,train$sample_id)]))
# a. 80 species: TWO lollipop columns in the SAME panel a
#    Left: 39 health-associated taxa; Right: 41 disease-associated taxa.
#    Follow original 06.fig2(2).R: cyan/red, dashed zero reference, italic species.
# ==========================================================================
coef_plot <- coef %>%
  mutate(feature=as.character(feature), coefficient=as.numeric(coefficient)) %>%
  filter(feature!="(Intercept)", is.finite(coefficient), coefficient!=0) %>%
  mutate(taxon=gsub("_"," ",sub("^s__","",feature)),
         direction=ifelse(coefficient>0,"Health-associated","Disease-associated"))
if (anyDuplicated(coef_plot$taxon)) stop("物种名称显示重复，请核对原始feature")
if (nrow(coef_plot)!=80L) stop("Expected 80 nonzero species")
write_csv(coef_plot,"Fig2a_coefficients_80_taxa")

make_coef_panel <- function(dat, direction, color, tag=NULL) {
  d <- dat %>% filter(.data$direction == .env$direction) %>%
    arrange(desc(abs(coefficient)), taxon)
  d$taxon <- factor(d$taxon, levels=rev(d$taxon))
  p <- ggplot(d,aes(x=coefficient,y=taxon)) +
    geom_vline(xintercept=0,linetype="dashed",color="grey60",linewidth=.55) +
    geom_segment(aes(x=0,xend=coefficient,yend=taxon),color=color,linewidth=.65) +
    geom_point(color=color,size=1.9) +
    labs(x="Coefficient",y=NULL,tag=tag) +
    theme_cell(15) +#a大小
    theme(axis.text.y=element_text(size=15,face="italic"),
          axis.text.x=element_text(size=14),
          axis.title.x=element_text(size=14),
          axis.ticks.y=element_blank(),
          plot.margin=margin(5,6,5,5))
  if (direction=="Health-associated") {
    p <- p + scale_x_continuous(limits=c(0,max(d$coefficient)*1.10),
                                 expand=expansion(mult=c(.03,.04)))
  } else {
    p <- p + scale_x_continuous(limits=c(min(d$coefficient)*1.10,0),
                                 expand=expansion(mult=c(.04,.03)))
  }
  # White strip header equivalent to original facet strip.
  p + facet_wrap(~direction,ncol=1) +
    theme(strip.text = element_text(face="bold", size=15),
    plot.tag = element_text(face="bold", size=18))
}
p_a_health <- make_coef_panel(coef_plot,"Health-associated",COL_HEALTH,tag="a")
p_a_disease <- make_coef_panel(coef_plot,"Disease-associated",COL_DISEASE)
p_a <- (p_a_health | p_a_disease) + plot_layout(widths=c(1,1.06))
save_both(p_a_health,"Fig2a_health_associated_39_taxa",7.3,13.0)
save_both(p_a_disease,"Fig2a_disease_associated_41_taxa",7.3,13.0)
save_both(p_a,"Fig2a_two_columns_health_vs_disease",16,13.0)

# ==========================================================================
# c. Two original-style violins: Training | External validation
#    Same palette and original Cell-like panel border; italic P and d only.
# ==========================================================================
violin_dat <- bind_rows(
  train %>% transmute(set="Training", OMHI, y_true),
  ext   %>% transmute(set="External validation", OMHI, y_true)
) %>%
  mutate(set=factor(set, levels=c("Training","External validation")),
         health_status=factor(ifelse(y_true==1,"Healthy","Non-healthy"),
                              levels=c("Healthy","Non-healthy")))
cliff_delta <- function(d) {
  x <- d$OMHI[d$y_true == 1]
  y <- d$OMHI[d$y_true == 0]

  n1 <- length(x)
  n2 <- length(y)

  U <- as.numeric(wilcox.test(x, y, exact = FALSE)$statistic)

  delta <- 2 * U / (n1 * n2) - 1

  return(delta)
}
calc_ba <- function(d,thr=0) {
  yy <- as.integer(d$OMHI>thr)
  mean(c(mean(yy[d$y_true==1]==1),mean(yy[d$y_true==0]==0)))
}
violin_stats <- violin_dat %>% group_by(set) %>%
  group_modify(~data.frame(n_healthy=sum(.x$y_true==1),
                           n_nonhealthy=sum(.x$y_true==0),
                           P=wilcox.test(OMHI~health_status,data=.x,exact=FALSE)$p.value,
                           delta=cliff_delta(.x),BA=calc_ba(.x))) %>% ungroup()
write_csv(violin_stats,"Fig2c_violin_statistics")
format_p_plotmath <- function(x) {
  if (x<0.001) {
    parts <- strsplit(sprintf("%.2e",x),"e",fixed=TRUE)[[1]]
    paste0(parts[1]," %*% 10^(",as.integer(parts[2]),")")
  } else sprintf("%.3f",x)
}
# Plotmath expressions: P and d italic; BA roman. Use three lines for reliable parsing.
annotation_df <- bind_rows(lapply(seq_len(nrow(violin_stats)),function(i) {
  st <- violin_stats[i,,drop=FALSE]
  data.frame(set=st$set,
             annotation=c(paste0("italic(P) == ",format_p_plotmath(st$P)),
                          sprintf("delta == %.2f",st$delta),
                          sprintf("BA == %.1f*'%%'",100*st$BA)),
             line_id=1:3)
}))
rng <- range(violin_dat$OMHI)
yspan <- diff(rng); if (yspan<=0) yspan<-1
annotation_top <- rng[2]+.23*yspan
annotation_df$y <- annotation_top-(annotation_df$line_id-1)*.064*yspan
annotation_df$x <- 1.5
counts <- violin_dat %>% count(set,health_status,name="n") %>%
  mutate(axis_label=paste0(as.character(health_status),"\n(n=",n,")"))
violin_dat <- violin_dat %>% left_join(counts,by=c("set","health_status"))
violin_dat$axis_label <- factor(violin_dat$axis_label,levels=counts$axis_label[order(counts$set,counts$health_status)])
set.seed(704)
p_c <- ggplot(violin_dat,aes(x=axis_label,y=OMHI,fill=health_status)) +
  geom_hline(yintercept=0,linetype="dashed",color=COL_GREY,linewidth=.65) +
  geom_violin(aes(color = health_status),width=.82,trim=FALSE,alpha=.65,linewidth = .6) +
  geom_jitter(aes(color=health_status),width=.095,size=.9,alpha=.84) +
  geom_boxplot(width=.17,outlier.shape=NA,fill=NA,color="black",linewidth=.55) +
  geom_text(data=annotation_df,aes(x=x,y=y,label=annotation),
            inherit.aes=FALSE,parse=TRUE,size=5) +
  facet_grid(.~set,scales="free_x",space="free_x") +
  scale_fill_manual(values=COL_GROUP,guide="none") +
  scale_color_manual(values=COL_GROUP,guide="none") +
  coord_cartesian(ylim=c(rng[1]-.12*yspan,annotation_top+.06*yspan)) +
  labs(x=NULL,y="OMWI",tag="c") +
  theme_cell(10) +
  theme(
        axis.text.x=element_text(size=13),
        axis.text.y=element_text(size=13),
        axis.title.y = element_text(size = 15, face = "plain"),
        strip.text = element_text(size = 14, face = "bold"),
        plot.tag = element_text(size = 18, face = "bold"),
        panel.spacing=grid::unit(.15,"cm"))
save_both(p_c,"Fig2c_violin_two_sets_original_style",6.3,5.5)
save_both(p_c,"Fig2C_Training_vs_External_validation",6.3,5.5)

#    Green solid left axis (BA) and blue dashed right axis (retention).
#    Both series have HOLLOW circles; BA numbers only, no retention labels.
# ==========================================================================
cutoffs <- c(0, .5, 1)
cutoff_ba <- function(data, dataset, cutoff) {
  retained <- data[abs(data$OMHI)>=cutoff,,drop=FALSE]
  ypred <- as.integer(retained$OMHI>0)
  healthy <- retained$y_true==1
  disease <- retained$y_true==0
  sen <- if(any(healthy)) mean(ypred[healthy]==1) else NA_real_
  spe <- if(any(disease)) mean(ypred[disease]==0) else NA_real_
  data.frame(set=dataset,cutoff=cutoff,
             BA=if(all(is.finite(c(sen,spe)))) mean(c(sen,spe)) else NA_real_,
             n_total=nrow(data),n_retained=nrow(retained),
             coverage=if(nrow(data)) nrow(retained)/nrow(data) else NA_real_)
}
cut_input <- list("Training"=train,"External validation"=ext,"10-fold CV"=cv_df)
e_data <- bind_rows(lapply(names(cut_input),function(name)
  bind_rows(lapply(cutoffs,function(x)cutoff_ba(cut_input[[name]],name,x))))) %>%
  mutate(set=factor(set,levels=names(cut_input)),
         BA_pct=100*BA, coverage_pct=100*coverage)
write_csv(e_data,"Fig2f_cutoff_balanced_accuracy_and_coverage")
if (any(!is.finite(e_data$BA_pct)))
  warning("部分cutoff导致单类全部被排除：该处BA为NA；不会捏造数值")
# Label only the green BA points. At high BA values place the text below
# the point to avoid clipping near 100%; otherwise place it above.
# Samples retained remain unlabeled, as requested.
e_data <- e_data %>%
  mutate(BA_label = ifelse(is.finite(BA_pct), sprintf("%.1f%%", BA_pct), NA_character_),
         BA_label_y = ifelse(BA_pct > 92, BA_pct - 7, BA_pct + 5.7))
p_e <- ggplot(e_data,aes(x=cutoff)) +
  geom_line(aes(y=BA_pct,group=set),color=COL_BA,linewidth=1.0) +
  geom_point(aes(y=BA_pct),shape=21,size=3,stroke=.86,
             fill="white",color=COL_BA) +
  geom_line(aes(y=coverage_pct,group=set),color=COL_RET,
            linewidth=.93,linetype="dashed") +
  geom_point(aes(y=coverage_pct),shape=21,size=3,stroke=.86,
             fill="white",color=COL_RET) +
  geom_text(aes(y=BA_label_y,label=BA_label),
            color="black",size=4,na.rm=TRUE,
            fontface="plain",show.legend=FALSE) +
  facet_grid(.~set) +
  scale_x_continuous(breaks=cutoffs,expand=expansion(mult=c(.09,.09))) +
  scale_y_continuous(limits=c(0,100),breaks=seq(0,100,20),
                     expand=expansion(add=c(1.0,2.0)),
                     name="Balanced accuracy (%)",
                     sec.axis=sec_axis(~.,name="Samples retained (%)",
                                       breaks=seq(0,100,20))) +
  labs(x="|OMWI| cutoff",tag="f") +
  theme_cell(9.4) +
  theme(panel.grid.major.y=element_line(color="grey90",linewidth=.28),
        axis.title.y.left=element_text(size=13,color=COL_BA,margin=margin(r=2)),
        axis.text.y.left=element_text(size=13,color="black"),
        axis.title.y.right=element_text(size=13,color=COL_RET,margin=margin(l=4)),
        axis.text.y.right=element_text(size=13,color=COL_RET),
        axis.ticks.y.right=element_line(color=COL_RET),
        axis.line.y.right=element_line(color=COL_RET,linewidth=.45),
        axis.text.x=element_text(size=13),
        axis.title.x=element_text(size=13),
        strip.text=element_text(size=15),
         plot.tag = element_text(face="bold", size=18),
        panel.spacing=grid::unit(.12,"cm"),
        plot.margin=margin(4,6,5,4))
save_both(p_e,"Fig2f_cutoff_dual_y_axis_original_style",9.2,3.5)


# Fig. 2b: pooled repeated out-of-fold scores; corrected LOOCV is optional until available.
roc_inputs <- list(`10-fold CV OOF`=oof,`External validation`=ext)
loo_path <- Sys.getenv("OMWI_LOOCV_PREDICTIONS",file.path(ROOT_DIR,"results/loocv/summary/nested_LOOCV_OOF.csv"))
if (file.exists(loo_path)) {
  loo <- read_csv(loo_path)
  assert_cols(loo,c("sample_id","y_true","prob_healthy","alpha"),"LOOCV")
  stopifnot(nrow(loo)==1077L,!anyDuplicated(loo$sample_id),all(loo$alpha==0.1),setequal(loo$sample_id,train$sample_id),all(loo$y_true==train$y_true[match(loo$sample_id,train$sample_id)]))
  roc_inputs <- c(list(`LOOCV OOF`=loo),roc_inputs)
} else message("Corrected LOOCV predictions unavailable; drawing the two available ROC curves only.")
roc_results <- lapply(names(roc_inputs),function(nm) {
  d <- roc_inputs[[nm]]; r <- pROC::roc(d$y_true,d$prob_healthy,levels=c(0,1),direction="<",quiet=TRUE)
  data.frame(dataset=nm,fpr=1-r$specificities,tpr=r$sensitivities,AUC=as.numeric(pROC::auc(r)))
})
roc_plot <- bind_rows(roc_results)
roc_labels <- unique(roc_plot[,c("dataset","AUC")]); write_csv(roc_labels,"Fig2b_AUC_summary")
labels <- setNames(sprintf("%s (AUC = %.2f)",roc_labels$dataset,roc_labels$AUC),roc_labels$dataset)
p_b <- ggplot(roc_plot,aes(fpr,tpr,color=dataset))+geom_abline(slope=1,intercept=0,linetype="dashed",color="grey70")+geom_path(linewidth=1)+scale_color_manual(values=COL_ROC,labels=labels)+coord_equal()+labs(x="False positive rate",y="True positive rate",tag="b")+theme_cell(14)+theme(legend.position="bottom")
save_both(p_b,"Fig2b_ROC",6.3,5.4)
