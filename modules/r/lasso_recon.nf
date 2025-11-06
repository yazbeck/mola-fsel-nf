process LASSO_RECON {
  tag { "${project}-${gene}" }

  input:
    tuple val(project), val(gene), val(mofa_dir)

  script:
  """
  WD=\$PWD; cd ${projectDir}
  Rscript scripts/04_lasso_recon.R \
    --project ${project} \
    --genes   ${gene} \
    --results ${params.results_base} \
    --in_data ${params.data_base} \
    --nfolds  ${params.nfolds} \
    --seed    ${params.seed} \
    --topn_expr ${params.topn_expr} \
    --topn_mut  ${params.topn_mut} \
    --topn_cnv  ${params.topn_cnv} 
  cd "\$WD"
  """
}
