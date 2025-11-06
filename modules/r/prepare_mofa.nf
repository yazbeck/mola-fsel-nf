process PREPARE_MOFA {
  tag { "${project}" }

  input:
    tuple val(project), val(project_data_dir)

  // emit MOFA dir as a VALUE (project-level)
  output:
    tuple val(project), val("${params.results_base}/${project}/mofa")

  script:
  """
  WD=\$PWD; cd ${projectDir}
  Rscript scripts/legacy/PrepareMOFA.R \
    --project ${project} \
    --gene ${params.genes?.split(',')?.first() ?: 'NA'} \
    --in_data ${params.data_base} \
    --out_mofa ${params.results_base}/${project}/mofa \
    --topn ${params.mofa_topn} \
    ${ params.predef      ? "--predef '${params.predef}'" : "" } \
    ${ params.predef_file ? "--predef_file '${params.predef_file}'" : "" }
  cd "\$WD"
  """
}
