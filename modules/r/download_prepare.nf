process DL_PREP {
  tag { "${project}-${gene}" }

  input:
    tuple val(project), val(gene)

  // hand off the PROJECT dir as a VALUE
  output:
    tuple val(project), val(gene), val("${params.data_base}/${project}")

  script:
  """
  WD=\$PWD; cd ${projectDir}
  Rscript scripts/legacy/DownloadPrepareMatrics.R \
    --project ${project} --gene ${gene} --outdir ${params.data_base}
  cd "\$WD"

  # ensure MOFA expects these at <data_base>/<PROJECT>/
  mkdir -p ${params.data_base}/${project}
  for f in expression_se.rds mutations_maf.rds cnv_thresholded_by_genes.rds; do
    [ -e "${params.data_base}/${project}/\$f" ] || \
      ln -sf \$(find "${params.data_base}/${project}" -maxdepth 3 -name "\$f" | head -n1) \
              "${params.data_base}/${project}/\$f" || true
  done
  """
}
