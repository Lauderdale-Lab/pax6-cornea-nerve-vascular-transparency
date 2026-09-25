###############################################################################
## run_all.R -- runs the pipeline end to end and records provenance
##
## Analysis pipeline for:
##   Nerve and vascular abnormalities precede loss of corneal transparency
##   in Pax6-haploinsufficient mice
##   Sneha K. Mohan, James D. Lauderdale
##
## James D. Lauderdale, PhD  (ORCID 0000-0001-7503-0528)
## Department of Cellular Biology, University of Georgia
## Athens, GA 30602, USA
##
## Repository : https://github.com/Lauderdale-Lab/pax6-cornea-nerve-vascular-transparency
## Archived   : <ZENODO_DOI>
## Licence    : MIT (see LICENSE)
## Contact    : James D. Lauderdale, jdlauder@uga.edu
##
## Run the pipeline with run_all.R. Scripts are numbered in execution order and
## share one R session by design; see run_all.R for why.
###############################################################################

###############################################################################
## run_all.R
##
## Runs the analysis end to end and records what was run.
##
## EXECUTION ORDER IS THE FILE ORDER. Reading the numbers top to bottom gives
## the order the steps actually execute in, which is the only ordering a reader
## should have to keep in their head.
##
##   00_config.R                 every parameter; computes nothing
##   00b_panel_universes.R       builds the three curated gene panels
##   01_load.R                   metadata, counts, pool sex composition, panels
##   01b_qc_library_depth.R      library complexity at matched depth
##   01c_contamination_indices.R off-target tissue and stromal share per
##                               library (Supp. Table 6); also defines the
##                               rule 09 applies to GSE183742
##   02_shared_gene_assignment.R primary panel for genes on more than one
##   03_expression_status.R      expressed / not expressed / indeterminate
##   04_differential.R           the two models, plus omnibus and interactions
##   05_program_level.R        abundance-matched permutation tests
##   06_developmental.R          trajectories within genotype
##   07_concordance.R            the status calls against a replicate-aware
##                               analysis of the same data; its verdicts are
##                               drawn on Figures 7 and 8B
##   08_figures.R                Figures 7-10 (optional; needs patchwork)
##   09_cross_dataset_replication.R  the calls against GSE183742 (optional;
##                               runs only when PAX6_DUNCAN_COUNTS and
##                               PAX6_DUNCAN_META are set)
##
## 07 precedes the figures because they read it: a call the replicate-aware
## model does not confirm is drawn as an open symbol. 09 is a check that
## nothing consumes, and it needs external data; a step whose inputs are not
## set is recorded as skipped, not failed, and the chain continues.
##
## EVERY STEP RUNS IN ONE ENVIRONMENT, deliberately. 05 needs transcriptome-wide
## expression status and 06 needs the primary-panel assignment, neither of which
## is written to disk -- one is too large to be worth it and the other would
## then exist in two places that could disagree. Sourcing into separate
## environments would break both.
##
## THE CHAIN STOPS AT THE FIRST FAILURE. Later steps consume earlier outputs, so
## continuing past an error would produce results built on a partial run, which
## is worse than no results.
##
## Provenance is written to the output directory on every run: the resolved
## configuration, the input files with their sizes and modification times, and
## the full session information. A result that cannot be traced to the code and
## inputs that produced it is not reproducible, whatever else is true of it.
##
## THE DRIVER IS A FUNCTION, and that is not stylistic. Steps are sourced into
## the global environment because they must share objects, which means any step
## can assign over a name the runner is using. A step with a top-level loop over
## `i`, or one that defines its own `STEPS`, would silently corrupt the runner's
## own bookkeeping and misreport what ran. Keeping the runner's state in a
## function frame puts it out of reach: source() still evaluates each script in
## the global environment, but the loop index, the step table and the results
## cannot be touched from there.
###############################################################################

run_started <- Sys.time()

SCRIPT_DIR <- if (nzchar(Sys.getenv("PAX6_SCRIPT_DIR"))) {
  Sys.getenv("PAX6_SCRIPT_DIR")
} else getwd()

## Steps in execution order. `required` marks a step whose absence is fatal;
## the two front-end steps are optional so that a working tree missing them
## still runs, with the gap reported rather than silently tolerated.
## `needs_env` lists environment variables a step cannot run without; when any
## is unset the step is recorded as skipped and the chain continues. Only the
## external-data check uses it.
STEPS <- data.frame(
  file = c("00_config.R", "00b_panel_universes.R", "01_load.R",
           "01b_qc_library_depth.R", "01c_contamination_indices.R",
           "02_shared_gene_assignment.R",
           "03_expression_status.R", "04_differential.R",
           "05_program_level.R", "06_developmental.R", "07_concordance.R",
           "08_figures.R", "09_cross_dataset_replication.R"),
  label = c("configuration", "panel universes", "load inputs",
            "library depth QC", "off-target tissue indices",
            "shared-gene assignment",
            "expression status", "differential expression",
            "program level", "developmental trajectories",
            "concordance with a replicate-aware analysis", "figures",
            "cross-dataset replication"),
  ## 07 is required: the figures read its verdicts. 01c is required: it writes
  ## Supp. Table 6 and defines the rule 09 uses.
  required = c(TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
               FALSE, FALSE),
  needs_env = c("", "", "", "", "", "", "", "", "", "", "", "",
                "PAX6_DUNCAN_COUNTS,PAX6_DUNCAN_META"),
  stringsAsFactors = FALSE)

message(strrep("=", 74))
message("PAX6 cornea RNA-seq -- three-state analysis")
message("Scripts: ", SCRIPT_DIR)
message(strrep("=", 74))

missing <- STEPS$file[!file.exists(file.path(SCRIPT_DIR, STEPS$file))]
if (length(missing)) {
  fatal <- intersect(missing, STEPS$file[STEPS$required])
  if (length(fatal)) {
    stop("Required script(s) not found in ", SCRIPT_DIR, ": ",
         paste(fatal, collapse = ", "))
  }
  message("NOT PRESENT, will be skipped: ", paste(missing, collapse = ", "))
}

## Any other numbered script in the directory is a version this run is not
## using. A stale copy left beside the current one is how a pipeline quietly
## runs the wrong code.
## Matched on any leading digit, not just 0: a stray 10_*.R sat beside the
## pipeline unreported because the pattern was "0*.R".
present <- sort(basename(Sys.glob(file.path(SCRIPT_DIR, "[0-9]*.R"))))
stale <- setdiff(present, STEPS$file)
if (length(stale)) {
  message("PRESENT BUT NOT RUN: ", paste(stale, collapse = ", "))
}

## ---------------------------------------------------------------------------

run_pipeline <- function(steps, script_dir) {
  results <- vector("list", nrow(steps))

  for (i in seq_len(nrow(steps))) {
    ## Everything the loop needs is read from the function's own arguments and
    ## frame, never from a global a sourced script could have overwritten.
    this_file  <- steps$file[i]
    this_label <- steps$label[i]
    this_path  <- file.path(script_dir, this_file)

    if (!file.exists(this_path)) {
      results[[i]] <- data.frame(step = i, script = this_file, ok = NA,
                                 seconds = NA, stringsAsFactors = FALSE)
      next
    }

    ## A step whose inputs are not set is skipped, and says so, rather than
    ## stopping the chain with an error about a missing file.
    env_needed <- if (nzchar(steps$needs_env[i])) {
      trimws(strsplit(steps$needs_env[i], ",")[[1]])
    } else character(0)
    env_unset <- env_needed[!nzchar(Sys.getenv(env_needed))]
    if (length(env_unset)) {
      message("\n[", i, "/", nrow(steps), "] ", this_label, " -- ", this_file,
              "\n  SKIPPED: not set: ", paste(env_unset, collapse = ", "))
      results[[i]] <- data.frame(step = i, script = this_file, ok = NA,
                                 seconds = NA, stringsAsFactors = FALSE)
      next
    }

    message("\n", strrep("-", 74))
    message("[", i, "/", nrow(steps), "] ", this_label, " -- ", this_file)
    message(strrep("-", 74))

    started <- Sys.time()
    ok <- tryCatch({
      source(this_path, local = FALSE, echo = FALSE)
      TRUE
    }, error = function(e) {
      message("\n*** FAILED: ", conditionMessage(e))
      FALSE
    })
    elapsed <- round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1)

    results[[i]] <- data.frame(step = i, script = this_file, ok = ok,
                               seconds = elapsed, stringsAsFactors = FALSE)
    if (!ok) {
      message("\nStopping. Later steps consume this step's output, so continuing ",
              "would build results on a partial run.")
      break
    }
  }

  out <- do.call(rbind, results)
  ## A step that neither ran nor was recorded means the loop lost track of
  ## itself, which is exactly the failure this function exists to prevent.
  stopifnot(nrow(out) == length(unique(out$step)))
  out
}

status <- run_pipeline(STEPS, SCRIPT_DIR)

## ---------------------------------------------------------------------------
## Provenance
## ---------------------------------------------------------------------------

if (exists("OUT_ROOT")) {
  input_files <- c(metadata = if (exists("FILE_METADATA")) FILE_METADATA else NA,
                   counts   = if (exists("FILE_COUNTS")) FILE_COUNTS else NA)
  input_rows <- vapply(names(input_files), function(nm) {
    p <- input_files[[nm]]
    if (is.na(p) || !file.exists(p)) return(paste0("  ", nm, ": NOT FOUND"))
    inf <- file.info(p)
    ## Content hash, because a modification time is not a reliable identity
    ## on a synced volume (iCloud rewrites mtime on sync). md5sum is in base
    ## R (tools) so this adds no dependency.
    md5 <- unname(tools::md5sum(p))
    sprintf("  %s: %s\n    %.1f MB, modified %s\n    md5 %s", nm, p,
            inf$size / 1024^2, format(inf$mtime, "%Y-%m-%d %H:%M:%S"), md5)
  }, character(1))

  provenance <- c(
    "PAX6 cornea RNA-seq -- run provenance",
    strrep("=", 70), "",
    paste("Run started :", format(run_started, "%Y-%m-%d %H:%M:%S %Z")),
    paste("Run finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("Script dir  :", SCRIPT_DIR),
    paste("Output dir  :", OUT_ROOT), "",
    "INPUT FILES", input_rows, "",
    "SCRIPTS RUN",
    sprintf("  %-32s %-6s %8s", "script", "ok", "seconds"),
    sprintf("  %-32s %-6s %8s", status$script,
            ifelse(is.na(status$ok), "skip", ifelse(status$ok, "TRUE", "FALSE")),
            ifelse(is.na(status$seconds), "-", format(status$seconds))), "",
    "PARAMETERS",
    sprintf("  detection exclusions : %s",
            if (length(DETECTION_EXCLUDE)) paste(DETECTION_EXCLUDE, collapse = ", ") else "(none)"),
    sprintf("  DE exclusions        : %s",
            if (length(DE_EXCLUDE)) paste(DE_EXCLUDE, collapse = ", ") else "(none)"),
    sprintf("  expression band      : ON at lower limit >= %s CPM, OFF at upper limit < %s CPM",
            T_HI, T_LO),
    sprintf("  Poisson alpha        : %s", POISSON_ALPHA),
    sprintf("  stability gate       : %s", STABILITY_GATE),
    sprintf("  categorical min fold : %s (conservative bound)", CATEGORICAL_MIN_FOLD),
    sprintf("  FDR                  : %s", FDR),
    sprintf("  equivalence bound    : +/-%s log2", EQUIV_BOUND),
    sprintf("  sex covariate        : %s", USE_SEX_COVARIATE),
    sprintf("  shared-gene rule     : %s", paste(PANEL_PRECEDENCE, collapse = " > ")),
    sprintf("  shared-gene overrides: %s",
            if (nrow(PANEL_PRIMARY_OVERRIDE))
              paste(PANEL_PRIMARY_OVERRIDE$symbol, "->",
                    PANEL_PRIMARY_OVERRIDE$panel, collapse = "; ") else "(none)"),
    sprintf("  permutations / bins  : %s / %s", N_PERM, N_BINS),
    sprintf("  seed                 : %s", SEED), "",
    "SESSION INFO", strrep("-", 70),
    utils::capture.output(utils::sessionInfo()))

  writeLines(provenance, file.path(OUT_ROOT, "RunProvenance.txt"))
}

## ---------------------------------------------------------------------------

message("\n", strrep("=", 74))
message("SUMMARY")
print(status, row.names = FALSE)
message(strrep("=", 74))

if (all(status$ok[!is.na(status$ok)])) {
  message("All steps completed in ", round(sum(status$seconds, na.rm = TRUE)), " seconds.")
  if (exists("OUT_ROOT")) {
    message("\nOutputs and RunProvenance.txt in:\n  ", OUT_ROOT)
    message("\nManuscript numbers come from:")
    message("  Contamination_Indices.csv    off-target tissue per library (Supp. Table 6)")
    message("  ExpressionCalls.csv          status, pooled counts, Poisson limits")
    message("  Assessability.csv            what statement each gene supports")
    message("  CategoricalGenes.csv         on in one group, off in the other")
    message("  DE_Results.csv               tested genes, with bounded vs undetermined")
    message("  InteractionResults.csv       does the genotype effect change with age")
    message("  AdultOmnibusLRT.csv          any difference among the three adult groups")
    message("  ProgramLevel.csv           panel-level, abundance-matched")
    message("  Developmental_*.csv          trajectories and status crossings")
    message("  Concordance/                 which status calls a replicate-aware model")
    message("                               confirms; drawn on Figures 7 and 8B")
    message("\nCheck on the analysis (nothing above depends on it):")
    message("  CrossDataset/                the calls against GSE183742")
    message("\nEvery negative in these files is either BOUNDED, with the excluded")
    message("effect stated, or UNDETERMINED. Undetermined is not evidence of no")
    message("effect and must not be reported as one.")
  }
} else {
  message("RUN INCOMPLETE. Nothing from this run should be quoted.")
}
