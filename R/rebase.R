

load_rebase <- function(file="data-raw/itype2.txt") {
  rebase <- readr::read_delim(file, delim="\t", skip=10,
                              col_names = c("enzyme", "prototype", "sequence",
                                            "methylation_site", "commercial_source", "references"),
                              col_types = "cccccc") %>%
    dplyr::mutate(base_sequence = stringr::str_replace_all(sequence, "[^A-Z]", ""),
                  is_palindrome = reverse_complement(base_sequence) == base_sequence,
                  is_ACTG = !stringr::str_detect(base_sequence, "[^ACTG]"),
                  length = nchar(base_sequence))

  return(rebase)
}

simulate_length <- function(rebase, nsites=1, nreps=10, lengths=c(20,40,60,80,100), timer="elapsed") {
  results <- dplyr::tibble(
    n_sites = rep(nsites, each=nreps*length(lengths)),
    oligo_lengths = rep(lengths, each=nreps),
    rep = rep(1:nreps, times=length(lengths)),
    sites = lapply(n_sites, function(n) sample(rebase$base_sequence, n)),
    time = 0.0
  )

  for (i in 1:nrow(results)) {
    solution <- cutfree(m=results$oligo_lengths[i],
                        sites=results$sites[[i]],
                        min_blocks=1,
                        re_randomize=FALSE,
                        quiet=TRUE)
    results$time[i] <- solution$time[timer]
    print(results$oligo_lengths[i])
  }

  return(results)
}

simulate_sites <- function(rebase, nsites=1:5, nreps=20, oligo_length=20, timer="elapsed") {
  results <- dplyr::tibble(
    n_sites = rep(nsites, each=nreps),
    rep = rep(1:nreps, times=length(nsites)),
    sites = lapply(n_sites, function(n) sample(rebase$base_sequence, n)),
    time = 0.0,
    degeneracy = NA
  )

  for (i in 1:nrow(results)) {
    solution <- cutfree(m=oligo_length,
                        sites=results$sites[[i]],
                        min_blocks=1,
                        re_randomize=FALSE,
                        quiet=TRUE)
    results$time[i] <- solution$time[timer]
    if (!is.na(solution$code))
      results$degeneracy[i] <- degeneracy(solution$code)
  }

  return(results)
}

run_cutree_paper <- function() {
  rebase <- load_rebase()

  rebase_pal_com <- rebase %>%
    dplyr::filter(!is.na(commercial_source),
                  is_palindrome,
                  length==6,
                  is_ACTG) %>%
    dplyr::filter(!duplicated(base_sequence))

  results <- simulate_sites(rebase_pal_com)
  results_lengths <- simulate_length(rebase_pal_com)

  theme_white <- theme_bw() + theme(panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank())

  plot1 <- ggplot2::ggplot(results, aes(x=n_sites, y=time)) +
    ggplot2::scale_y_log10() +
    ggplot2::geom_point() +
    ggplot2::xlab("Blocked Sites") +
    ggplot2::ylab("Runtime [s]") +
    ggplot2::geom_smooth(method="lm", se=FALSE)

  plot2 <- ggplot2::ggplot(results, aes(x=n_sites, y=degeneracy)) +
    ggplot2::scale_y_log10() +
    ggplot2::geom_point() +
    ggplot2::xlab("Blocked Sites") +
    ggplot2::ylab("Possible Sequences") +
    ggplot2::geom_smooth(method="lm", se=FALSE)

  plot0 <- ggplot2::ggplot(results_lengths, aes(x=oligo_lengths, y=time)) +
    ggplot2::scale_y_log10() +
    ggplot2::geom_point() +
    ggplot2::xlab("Randomer Length [bp]") +
    ggplot2::ylab("Runtime [s]") +
    ggplot2::geom_smooth(method="lm", se=FALSE)

  results %>%
    dplyr::group_by(n_sites) %>%
    dplyr::select(n_sites, time, degeneracy) %>%
    dplyr::summarise_all(dplyr::funs(mean, sd))

  cowplot::plot_grid(plot1, plot2, plot0, labels=c("A","B","C"), align="h")
}



