process GET_TOP_WEIGHTS {
  tag { "${project}-${gene}" }

  input:
    tuple val(project), val(gene), val(mofa_dir)

  output:
    // pass-through so downstream LASSO has the same (p,g,mofa_dir)
    tuple val(project), val(gene), val(mofa_dir)	

  script:

  def outdir = "${params.results_base}/${project}/weights/${gene}"
  """
  set -euo pipefail
  mkdir -p "${outdir}"	
  WD=\$PWD; cd ${projectDir} 
  Rscript ${workflow.projectDir}/scripts/legacy/GetTopWeights.R \
    --project ${project} --gene ${gene} \
    --mofa_dir ${mofa_dir} \
    --outdir   ${outdir} \
    --n_top    ${params.n_top} \
    ${ params.predef      ? "--predef '${params.predef}'" : "" } \
    ${ params.predef_file ? "--predef_file '${params.predef_file}'" : "" }
  cd "\$WD"
  """
}
