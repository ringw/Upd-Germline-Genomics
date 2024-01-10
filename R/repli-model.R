#source('C:/Users/dringwa1/OneDrive - Johns Hopkins/Documents/scTJ/Rollutils.R')

repli_fit_hmm = function(repli, ident='Tj') {
  chrs = rowData(repli)$chr %>% recode(
    `2L`='2', `2R`='2', `3L`='3', `3R`='3'
  )
  # chrs are continuous
  chr_lens = table(chrs)
  repli = assay(repli, 'pct')[, repli$ident == ident]
  repli_kmeans = kmeans(repli[is.finite(repli[,1]), ], 3)
  repli = na.approx(repli, na.rm=F) %>% replace(is.na(.), 0)
  
  hmmlearn = reticulate::import('hmmlearn')
  np = reticulate::import('numpy')
  np$random$seed(0L)
  hmm_dirichlet = 1000*diag(nrow=3) + matrix(0.0001, 3, 3)
  hmm = hmmlearn$hmm$GaussianHMM(n_components=3L, init_params='stc', transmat_prior = hmm_dirichlet, 
                                 means_prior=repli_kmeans$centers,
                                 means_weight=matrix(10000, nrow=nrow(repli_kmeans$centers), ncol=ncol(repli_kmeans$centers)))
  # hmm$means_ = np$asarray(repli_kmeans$centers)
  hmm = hmm$fit(repli, chr_lens)
  hmm_proba = hmm$predict_proba(repli, chr_lens)
  list(
    means=matrix(
      hmm$means_,
      nrow=3,
      ncol=4,
      dimnames=list(c('late','indeterminate','early'), colnames(repli))
    ),
    covars=list(
      hmm$covars_[1,,],
      hmm$covars_[2,,],
      hmm$covars_[3,,]
    ),
    predict=matrix(
      hmm_proba,
      nrow=nrow(hmm_proba),
      ncol=ncol(hmm_proba),
      dimnames=list(rownames(repli), c('late','indeterminate','early'))
    )
  )
}