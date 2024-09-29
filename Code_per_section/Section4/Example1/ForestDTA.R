if (!require(pacman)) {install.packages("pacman"); library(pacman)}
p_load(tidyverse,
       magrittr,
       glue,
       mada,
       lme4,
       broom,
       broom.mixed,
       ggnewscale)

fnDiamond <- function(m, l, h, y, w){
  tibble(xx = c(m, l, m, h),
         yy = c(y+w, y, y-w, y))
}

fnMetanDTA <- function(tblIn, UID, TP, FN, FP, TN,
                       level              = 0.95,
                       correction         = 0.5, 
                       correction.control = "all") {
  if (missing(tblIn) & (missing(TP) | missing(FN) | missing(FP) | missing(FN) | missing(UID))) stop("no valid input")
  if (missing(tblIn)) {
    tblIn <- tibble(TP, FN, FP, TN, UID)
  } else {
    if (!is.data.frame(tblIn)) {
      stop("no valid input")
    } else {
      tblIn %<>% select(all_of(c("UID", "TP", "FN", "FP", "TN")))}
  }
  
  if (correction.control != "none") {
    tblIn %<>%
      {if (correction.control != "all") rowwise(.) else .} %>%
      mutate(tmp = ifelse(any(c_across(TP:TN) == 0), correction, 0),
             across(TP:TN, ~ .x + tmp)) %>%
      select(-tmp)
  }
  
  tbl2x2Long <- tblIn %>%
    mutate(refP = TP+FN,
           refN = TN+FP) %>%
    pivot_longer(cols = c("TP", "FP"),
                 values_to = "indP") %>%
    mutate(sens = as.numeric(name == "TP"),
           fpr  = 1-sens,
           n    = ifelse(sens, refP, refN))
  
  metan <- tbl2x2Long %>%
    glmer(data    = .,
          formula = cbind(indP, n - indP) ~ 0 + sens + fpr + (0 + sens + fpr|UID),
          family  = binomial)
  
  tidyCoef <- metan %>%
    tidy() %>%
    bind_rows({.} %>%
                filter(term == "fpr") %>%
                mutate(term     = "spec",
                       estimate = -estimate)) %>%
    mutate(term = str_replace_all(term, "(sens|spec|fpr)", "ln\\1"),
           CIlo = estimate - std.error * qnorm(0.975),
           CIhi = estimate + std.error * qnorm(0.975)) %>%
    bind_rows({.} %>%
                filter(str_starts(term, "lns")) %>%
                mutate(term = str_remove_all(term, "ln"),
                       across(c(estimate, CIlo, CIhi), plogis),
                       across(std.error:p.value, ~NA)))
  
  vcm <- vcov(metan) 
  rownames(vcm) <- colnames(vcm) <- c("lnsens", "lnfpr")
  
  list(data = tblIn, metan = metan, tidyCoef = tidyCoef, vcov = vcm)
}

fnForestDTA <- function(dfIn, TP, FN, FP, TN, UID, w_sens, w_spec,
                        summ_sens_est, summ_sens_CI, summ_spec_est, summ_spec_CI,
                        accuracy = 0.001,
                        correction.control = "none",
                        correction = 0.5,
                        wMethod = c("sample_size", "inverse_variance", "custom", "none"),
                        ci.method = "clopper-pearson",
                        lw = 0.5,
                        ts = 10,
                        blnStripes = TRUE,
                        blnSavePlot = FALSE,
                        strSavePath) {
  
  wMethod <- match.arg(wMethod)
  
  if (missing(dfIn) & (missing(TP) | missing(FN) | missing(FP) | missing(FN) | missing(UID))) stop("no valid input")
  if (missing(dfIn)) {
    tblIn <- tibble(TP, FN, FP, TN, UID)
  } else {
    if (!is.data.frame(dfIn)) stop("no valid input: dfIn is not a dataframe")
    if (!"TP" %in% colnames(dfIn)) stop("no valid input: dfIn does not have a column called 'TP'")
    if (!"FN" %in% colnames(dfIn)) stop("no valid input: dfIn does not have a column called 'FN'")
    if (!"FP" %in% colnames(dfIn)) stop("no valid input: dfIn does not have a column called 'FP'")
    if (!"TN" %in% colnames(dfIn)) stop("no valid input: dfIn does not have a column called 'TN'")
    if (!"UID" %in% colnames(dfIn)) stop("no valid input: dfIn does not have a column called 'UID'")
    tblIn <- dfIn %>%
      select(all_of(c("UID", "TP", "FN", "FP", "TN")), any_of(c("w_sens", "w_spec")))
  }
  
  if (wMethod=="custom") {
    if (missing(w_sens) | missing(w_spec)) {
      if (!"w_sens" %in% colnames(dfIn) | !"w_spec" %in% colnames(dfIn)) {
        stop("custom weights specified but none provided")
      } else {
        tblIn %<>% 
          rename(sns_w = w_sens, spc_w = w_spec)}
    } else {
      tblIn %<>% 
        add_column(sns_w = w_sens) %>% 
        add_column(spc_w = w_spec)
    }
  } else {
    tblIn %<>% 
      add_column(sns_w = 0) %>% 
      add_column(spc_w = 0)
  }
  
  nStuds <- nrow(tblIn)

  mdd <- tblIn %>% 
    rename(names = UID) %>% 
    madad(data,
          method = ci.method,
          correction.control = correction.control,
          correction = correction)
  
  tblHet <- tibble(term = c("sns", "spc"),
                   hh = c(mdd$sens.htest %$%
                            paste0("\u03c7\u00b2 = ", signif(statistic, 5), ", df = ", parameter, ", p ", ifelse(p.value < 0.001, "< 0.001", paste("=", scales::label_number(0.001)(p.value)))),
                          mdd$spec.htest %$%
                            paste0("\u03c7\u00b2 = ", signif(statistic, 5), ", df = ", parameter, ", p ", ifelse(p.value < 0.001, "< 0.001", paste("=", scales::label_number(0.001)(p.value)))))) %>% 
    rowwise() %>% 
    mutate(hh = ifelse(nStuds==1, "", hh))
  
  tblMA <- mdd %$%
    tibble(UID = names,
           dat = data,
           sns_estimate = sens$sens,
           sns_CIlo = sens$sens.ci[,1],
           sns_CIhi = sens$sens.ci[,2],
           spc_estimate = spec$spec,
           spc_CIlo = spec$spec.ci[,1],
           spc_CIhi = spec$spec.ci[,2]) %>% 
    unnest(dat) %>% 
    left_join(tblIn,
              by = join_by(UID, TP, FN, FP, TN))
  
  lev <- tblIn %>% 
    arrange(UID) %>% 
    pull(UID) %>% 
    c("", ., "Total")
  
  tblMA %<>% 
    {if (correction.control != "all") rowwise(.) else .} %>%
    mutate(tmp = ifelse(any(c_across(TP:TN) == 0), correction, 0),
           across(TP:TN, ~ .x + tmp, .names = "{.col}c")) %>% 
    pivot_longer(cols = matches("_"),
                 names_to = c("term", "v"),
                 names_sep = "_") %>% 
    mutate(cT = ifelse(term=="sns", TPc, TNc),
           cF = ifelse(term=="sns", FNc, FPc),
           ff = "plain") %>% 
    select(-(TP:TNc)) %>% 
    pivot_wider(values_from = value,
                names_from = v) %>% 
    mutate(across(c(estimate, CIlo, CIhi),
                  .fns = ~scales::label_number(accuracy = accuracy)(.x),
                  .names = "{.col}_ps"),
           strEM = glue("{estimate_ps} ({CIlo_ps} to {CIhi_ps})")) %>% 
    group_by(UID, term) %>% 
    mutate(n = sum(c_across(c(cT, cF))),
           UID = factor(UID, levels = lev)) %>% 
    group_by(term) %>% 
    mutate(se = sqrt(1/cT + 1/cF),
           iv = se^-2,
           nw = n / sum(n),
           ivw = iv / sum(iv),
           cw = w / sum(w)) %>%
    rowwise() %>% 
    mutate(ww = case_match(wMethod,
                           "sample_size"      ~ nw,
                           "inverse_variance" ~ ivw,
                           "custom"           ~ cw,
                           "none"             ~ 0.1),
           across(c(cT, cF), as.character)) %>% 
    ungroup() %>% 
    add_row(term = "sns", UID = factor(""), cT = "TP", cF = "FN", strEM = "Sensitivity (95%CI)", ff = "bold") %>% 
    add_row(term = "spc", UID = factor(""), cT = "TN", cF = "FP", strEM = "Specificity (95%CI)", ff = "bold") 
  
  MA <- tblIn %>%
    fnMetanDTA(correction.control = correction.control,
               correction = correction)
  
  tblSumm <- MA$tidyCoef %>% 
    filter(term %in% c("sens", "spec")) %>% 
    select(term, estimate, CIlo, CIhi)
  
  if (!missing(summ_sens_est) & !missing(summ_sens_CI)) {
    tblSumm %<>% 
      rows_update(tibble_row(term = "sens", estimate = summ_sens_est, CIlo = summ_sens_CI[1], CIhi = summ_sens_CI[2]), by = "term")
  }
  if (!missing(summ_spec_est) & !missing(summ_spec_CI)) {
    tblSumm %<>% 
      rows_update(tibble_row(term = "spec", estimate = summ_spec_est, CIlo = summ_spec_CI[1], CIhi = summ_spec_CI[2]), by = "term")
  }
  
  tblMA %<>%
    bind_rows(tblSumm %>%
                mutate(across(-term,
                              .fns = ~scales::label_number(accuracy = accuracy)(.x),
                              .names = "{.col}_s"),
                       strEM = glue("{estimate_s} ({CIlo_s} to {CIhi_s})"),
                       ff = "bold",
                       UID = factor("Summary estimate"),
                       term = str_replace_all(term, c("sens" = "sns", "spec" = "spc"))) %>%
                select(UID, term, ff, estimate, CIlo, CIhi, strEM)) %>%
    arrange(UID)
  
  f <- tblMA %>% 
    ggplot(aes(x = estimate, y = UID)) +
    {if (blnStripes) {geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
              data = tibble(xmin = rep(0, nStuds),
                            xmax = rep(2, nStuds),
                            ymin = 1:nStuds+0.5,
                            ymax = 1:nStuds+1.5) %>%
                mutate(fill = as.character(row_number() %% 2)),
              inherit.aes = FALSE)}} +
    scale_fill_manual(values = c("white", "grey90")) +
    new_scale_fill() +
    geom_path(aes(x = estimate, y = yy, group = term), colour = "gray50", linetype = 2, linewidth = lw,
              data = tblSumm %>% 
                mutate(term = str_replace_all(term, c("sens" = "sns", "spec" = "spc"))) %>% 
                select(term, estimate) %>% 
                slice(rep(1:2, each = 2)) %>% 
                bind_cols(yy = rep(c(nStuds+1.5, 0.5), 2))) +
    geom_point(aes(size = ww), fill = "white", shape = 22, stroke = 0) +
    geom_errorbar(aes(xmin = CIlo, xmax = CIhi), width = 0, linewidth = lw, data = . %>% filter(UID != "Summary estimate")) +
    geom_point(aes(size = ww, fill = term, colour = term), shape = 22, stroke = NA, alpha = 0.5) +
    scale_fill_manual(values = c("red", "blue")) +
    geom_text(aes(x = 1.10, label = cT, fontface = ff), hjust = 0.5, size = ts, size.unit = "pt" ) +
    geom_text(aes(x = 1.25, label = cF, fontface = ff), hjust = 0.5, size = ts, size.unit = "pt" ) +
    geom_text(aes(x = 1.35, label = strEM, fontface = ff), hjust = 0, size = ts, size.unit = "pt" ) +
    scale_x_continuous(name = "",
                       limits = c(0, 2),
                       breaks = 0:5/5,
                       expand = c(0,0)) +
    scale_y_discrete(name = "",
                     limits = rev,
                     expand = expansion(add = c(0, 0.5))) +
    geom_path(aes(x = xx, y = yy, group = xx), linewidth = lw, data = tibble(xx = rep(0:1, each=2),
                                                                             yy = rep(c(nStuds+1.5, 0.5), 2))) +
    geom_path(aes(x = xx, y = yy, group = yy), linewidth = lw, data = tibble(xx = rep(0:1, 2),
                                                                             yy = rep(c(nStuds+1.5, 0.5), each = 2))) +
    geom_polygon(aes(x = xx, y = yy, fill = term),
                 data = tblSumm %>% 
                   group_by(term) %>% 
                   mutate(xx = list(fnDiamond(estimate, CIlo, CIhi, 1, 0.4)),
                          term = str_replace_all(term, c("sens" = "sns", "spec" = "spc"))) %>% 
                   unnest(xx) %>% 
                   select(term, xx, yy)) +
    coord_cartesian(ylim = c(0.5,NA), expand = FALSE, clip = 'off') +
    geom_text(aes(x = 0.5, y = -1, label = ll), size = ts/.pt, fontface = "bold",
              data = tibble(term = c("sns", "spc"), ll = c("Sensitivity", "Specificity"))) +
    geom_text(aes(x = 1.1, y = -0, label = hh),
              fontface = "italic",
              size = ts/.pt,
              vjust = 0.5,
              hjust = 0, 
              data = tblHet) +
    facet_wrap(vars(term), labeller = as_labeller(
      c(tsens = "Sensitivity", tspec = "Specificity")),
      strip.position = "bottom") +
    theme_minimal() +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          strip.background = element_blank(),
          axis.text = element_text(size = ts),
          axis.text.y = element_text(hjust = 0, face = c("bold", rep("plain", nStuds+1))),
          strip.placement = "none",
          panel.spacing.x = unit(0, "cm"),
          strip.text.x = element_blank(),
          axis.ticks.x = element_line(colour = "black", linewidth = lw))
  
  f
  if (blnSavePlot) {
    fnSaveGG(plt = f,
             strName = strSavePath,
             cmW = 30,
             cmH = 2 + (nStuds+1)*0.5,
             typ = "svg")
  }
  
  return(f)
}

fnSaveGG <- function(plt, strName, cmW, cmH, dpi=500, typ = c("svg", "pdf", "jpg")) {
  if (missing(plt)) {plt <- last_plot()}
  for (t in typ) {
    ggsave(filename = glue("{strName}.{t}"),
           plot     = plt,
           device   = t,
           scale    = 1,
           width    = cmW,
           height   = cmH,
           units    = "cm",
           dpi      = dpi,
           bg       = "transparent")
  }
}

tblBNP100 <- tribble(
  ~Study,  ~TP, ~FN,  ~FP,  ~TN,
  "Alibay, 2005",  59L,  1L,  53L,  47L,
  "Arques, 2005",  31L,  1L,  14L,  24L,
  "Barcase, 2004",  55L,  2L,   4L,  37L,
  "Blonde-Cynober, 2011",  23L,  3L,  12L,  26L,
  "Chenevier-Gobeaux, 2010", 114L,  1L, 155L, 108L,
  "Chung, 2006",  72L,  0L,  42L,  29L,
  "Dao, 2001",  91L,  6L,   9L, 144L,
  "Davis, 1994",  26L,  6L,   2L,  18L,
  "Lainchbury, 2003",  68L,  2L,  69L,  66L,
  "Logeart, 2002", 110L,  5L,  33L,  15L,
  "Lokuge, 2010", 252L, 22L, 166L, 172L,
  "Maisel, 2002", 670L, 74L, 202L, 640L,
  "Maisel, 2010", 543L, 25L, 409L, 664L,
  "Mueller, 2005", 132L,  5L,  44L,  70L,
  "Parab, 2005",  45L,  2L,  17L,   6L,
  "Ray, 2004", 127L, 14L,  68L,  99L,
  "Rogers, 2009", 353L, 15L, 115L, 257L,
  "Sanz, 2006",  43L,  2L,   3L,  27L,
  "Wang, 2010",  46L,  3L,  23L,  12L) %>%
  rename(UID = Study)

tblBNP100 %>% 
  fnForestDTA(summ_sens_est = 0.95, summ_sens_CI = c(0.93, 0.97),
              summ_spec_est = 0.63, summ_spec_CI = c(0.50, 0.75),
              accuracy = 0.01,
              wMethod = "samp",
              blnStripes = T,
              blnSavePlot = T,
              strSavePath = here::here("@@TstDelete"))

tblBNP100 %$%
  fnForestDTA(TP = TP, FN = FN, FP = FP, TN = TN, UID = UID,
              summ_sens_est = 0.95, summ_sens_CI = c(0.93, 0.97),
              summ_spec_est = 0.63, summ_spec_CI = c(0.50, 0.75),
              accuracy = 0.01,
              blnSavePlot = F,
              strSavePath = here::here("@@TstDelete"))
