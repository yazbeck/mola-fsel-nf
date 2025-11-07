nextflow.enable.dsl=2

include { DL_PREP         } from "./modules/r/download_prepare.nf"
include { PREPARE_MOFA    } from "./modules/r/prepare_mofa.nf"
include { GET_TOP_WEIGHTS } from "./modules/r/get_top_weights.nf"
include { LASSO_RECON     } from "./modules/r/lasso_recon.nf"

workflow {
  // params → lists
  def projects = (params.projects ?: 'TCGA-BRCA').split(',').collect{ it.trim() }.findAll{ it }
  def genes    = (params.genes    ?: 'TP53').split(',').collect{ it.trim() }.findAll{ it }

  // (project,gene) pairs
  Channel.fromList(projects.collectMany { p -> genes.collect { g -> tuple(p,g) } })
         .set { pg_pairs }

  // 1) download/prepare (per project,gene is fine)
  DL_PREP(pg_pairs)    // emits: (project, gene, "<data_base>/<project>")

  // 2) unique projects → MOFA (run once per project)
  DL_PREP.out
    .map { p, g, projdir -> tuple(p, projdir) }
    .distinct()
    .set { perProject }                       // (project, "<data_base>/<project>")

  PREPARE_MOFA(perProject)                    // emits: (project, "<results_base>/<project>/mofa")

  // 3) join MOFA back to (project,gene) → (project, gene, mofa_dir)
  pg_pairs
    .combine( PREPARE_MOFA.out )                 // join by 'project' key... now use combine instead of join
    .map { p1, g, p2, mofa_dir -> tuple(p1, g, mofa_dir) }  // Extract only the needed elements
    //.map { p, g, mofa_dir -> tuple(p, g, mofa_dir) }
    .set { perGeneWithMofa }


  GET_TOP_WEIGHTS(perGeneWithMofa)
  LASSO_RECON(GET_TOP_WEIGHTS.out)   // feed LASSO from weights' output (no .into needed)
}