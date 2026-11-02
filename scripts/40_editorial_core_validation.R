#!/usr/bin/env Rscript
# Fixed-set external gene validation and paired deconvolution effects.
source("scripts/00_config.R")
suppressPackageStartupMessages({library(limma); library(dplyr); library(readr); library(tidyr); library(ggplot2)})
set.seed(20261031)
x <- readRDS(out_path("results", "analysis_inputs.rds"))$gse31348
e <- as.matrix(x$expression)
m <- x$metadata[match(colnames(e), x$metadata$geo_accession), ]
stopifnot(identical(colnames(e), m$geo_accession), ncol(e) == 135L, length(unique(m$SubjectID)) == 27L)
m$TimePoint <- factor(m$TimePoint, levels = c("Week 0", "Week 1", "Week 2", "Week 4", "Week 26"))
d <- model.matrix(~0+TimePoint, data=m); colnames(d) <- make.names(colnames(d))
dc <- duplicateCorrelation(e, d, block=m$SubjectID)
fit <- lmFit(e, d, block=m$SubjectID, correlation=dc$consensus.correlation)
fit <- eBayes(contrasts.fit(fit, makeContrasts(Week26=TimePointWeek.26-TimePointWeek.0, levels=d)))
stats <- topTable(fit, coef="Week26", number=Inf, sort.by="none") %>% tibble::rownames_to_column("Gene")
sets <- list(restoration_core=list(file="restoration_core_genes.txt", direction=1), disease_driver=list(file="disease_driver_genes.txt", direction=-1))
gene_rows <- lapply(names(sets), function(nm) {
  s <- sets[[nm]]; g <- unique(readLines(out_path("results", s$file), warn=FALSE))
  tibble(Set=nm, Gene=g, ExpectedDirection=if(s$direction==1) "increase" else "decrease") %>%
    left_join(stats %>% select(Gene, logFC, AveExpr, P.Value, adj.P.Val), by="Gene") %>%
    mutate(Measured=!is.na(logFC), DirectionConsistent=Measured & sign(logFC)==s$direction)
}) %>% bind_rows()
write_tsv(gene_rows, out_path("results", "gse31348_fixed_gene_direction.tsv"))
# Expression-matched empirical reference; selected sets remain fixed, never reselected in GSE31348.
stats$bin <- ntile(stats$AveExpr, 10)
summary_rows <- lapply(names(sets), function(nm) {
  z <- gene_rows %>% filter(Set==nm, Measured)
  bins <- stats$bin[match(z$Gene, stats$Gene)]
  wanted <- sets[[nm]]$direction
  observed <- sum(sign(z$logFC)==wanted)
  null <- replicate(5000, {
    sampled <- unlist(lapply(unique(bins), function(b) sample(which(stats$bin==b), sum(bins==b))))
    sum(sign(stats$logFC[sampled])==wanted)
  })
  tibble(Set=nm, N_defined=sum(gene_rows$Set==nm), N_measured=nrow(z), N_consistent=observed,
         Proportion=observed/nrow(z), N_FDR05_expected=sum(z$DirectionConsistent & z$adj.P.Val<.05),
         Matched_random_mean=mean(null), Empirical_P=(1+sum(null>=observed))/(length(null)+1),
         Validation="Week 26 versus Week 0; 27 paired participants; expression-matched random gene sets")
}) %>% bind_rows()
write_tsv(summary_rows, out_path("results", "gse31348_fixed_gene_direction_summary.tsv"))

cib <- read_tsv(out_path("results", "gse31348_cibersort_lm22_legacy_settings.tsv"), show_col_types=FALSE)
q <- read_tsv(out_path("results", "gse31348_quantiseq_explicit_nonlog.tsv"), show_col_types=FALSE)
cat("CIBERSORT columns: ", paste(names(cib),collapse=", "), "\n")
neut <- grep("^Neutrophils?_CIBERSORT$", names(cib), value=TRUE, ignore.case=TRUE)
if(length(neut)!=1L) stop("Cannot find CIBERSORT neutrophil column")
cell <- tibble(SampleID=cib$SampleID, SubjectID=cib$SubjectID, TimePoint=cib$TimePoint,
  CIBERSORT_B=cib$CIBERSORT_total_B, CIBERSORT_Neutrophil=cib[[neut]]) %>%
  inner_join(q %>% transmute(SampleID, quanTIseq_B=`B cell`, quanTIseq_Neutrophil=Neutrophil), by="SampleID")
stopifnot(nrow(cell)==135L)
long <- cell %>% pivot_longer(matches("^(CIBERSORT|quanTIseq)_(B|Neutrophil)$"), names_to=c("Method","Lineage"), names_sep="_", values_to="Fraction")
paired <- long %>% filter(TimePoint %in% c("Week 0", "Week 26")) %>%
  select(SubjectID, TimePoint, Method, Lineage, Fraction) %>% pivot_wider(names_from=TimePoint, values_from=Fraction) %>%
  mutate(Change=`Week 26`-`Week 0`)
stopifnot(nrow(paired)==108L, all(is.finite(paired$Change)))
effects <- paired %>% group_by(Method,Lineage) %>% summarise(N_participants=n(), Mean_change=mean(Change),
    SD_change=sd(Change), SE=sd(Change)/sqrt(n()), CI_low=Mean_change-qt(.975,n()-1)*SE,
    CI_high=Mean_change+qt(.975,n()-1)*SE, P_value=t.test(Change)$p.value, .groups="drop")
write_tsv(paired, out_path("results", "gse31348_paired_cell_fraction_changes.tsv"))
write_tsv(effects, out_path("results", "gse31348_paired_cell_fraction_effects.tsv"))
p <- ggplot(effects, aes(Mean_change*100, Method, colour=Method)) +
  geom_vline(xintercept=0, linetype=2, colour="grey55") +
  geom_errorbar(aes(xmin=CI_low*100,xmax=CI_high*100),orientation="y",width=.12,linewidth=.9) +
  geom_point(size=3) + facet_wrap(~Lineage,ncol=2,scales="free_x") +
  scale_colour_manual(values=c(CIBERSORT="#3B6D8E",quanTIseq="#C46743"),guide="none") +
  labs(title="Estimated cell fractions after 26 weeks of treatment", subtitle="Mean paired change and 95% CI; 27 patients",
       x="Change from baseline (percentage points)",y=NULL) +
  theme_classic(base_size=14) + theme(strip.text=element_text(face="bold",size=14),plot.title=element_text(face="bold"))
ggsave(out_path("figures","FigureS4_GSE31348_CIBERSORT_LM22_Bcell_Trajectory.pdf"),p,width=8,height=4.4,device=cairo_pdf)
print(summary_rows); print(effects)
