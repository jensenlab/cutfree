
IUB_CODES <- list(
  A = c("A"),
  C = c("C"),
  G = c("G"),
  T = c("T"),
  R = c("A", "G"),
  Y = c("C", "T"),
  M = c("A", "C"),
  K = c("G", "T"),
  S = c("C", "G"),
  W = c("A", "T"),
  H = c("A", "C", "T"),
  B = c("C", "G", "T"),
  V = c("A", "C", "G"),
  D = c("A", "G", "T"),
  N = c("A", "C", "G", "T")
)

printf <- function(fmt, ...) {
  cat(sprintf(fmt, ...))
}

# Find codes that do not contain a code.
# Example: if code is A, then blocking codes are
#   C,G,T,Y,K,S,B (all codes without A)
# If names==TRUE, return the codes (character vector of bases);
# otherwise return the code names.
get_blocking_codes <- function(code, codes=IUB_CODES, names=FALSE) {
  to_block <- codes[code]
  blocking <- Filter(Negate(function(bases) any(to_block %in% bases)), codes)
  if (names) {
    return(names(blocking))
  } else {
    return(blocking)
  }
}

# Turn a string into a vector of single characters.
str_to_vector <- function(str) {
  if (length(str) == 1) {
    strsplit(str, "")[[1]]
  } else {
    str
  }
}

# Create a text block showing which bases are allowed
# by the oligo. Returns four strings denoting the presence
# of each base.
# Example: make_oligo_block("ACGTRYMKSWHBVDN")
#                 A                 C                 G                 T
# "A---A-A--AA-AAA" "-C---CC-C-CCC-C" "--G-G--GG--GGGG" "---T-T-T-TTT-TT"
make_oligo_block <- function(oligo, codes=IUB_CODES) {
  oligo <- str_to_vector(oligo)

  bases <- c("A", "C", "G", "T")
  strs <- c(A="", C="", G="", T="")
  for (base in bases) {
    strs[base] <- paste(sapply(oligo, function(code) if (base %in% codes[[code]]) base else "-"),
                        collapse="")
  }
  return(strs)
}

# Show an oligo followed by the allowed bases
# (from make_oligo_block).
# Example: print_oligo_block("ACGTRYMKSWHBVDN")
#   ACGTRYMKSWHBVDN
#   A---A-A--AA-AAA
#   -C---CC-C-CCC-C
#   --G-G--GG--GGGG
#   ---T-T-T-TTT-TT
print_oligo_block <- function(oligo, indent="") {
  cat(paste0(indent, oligo, collapse=""), "\n")
  cat(paste0(indent, make_oligo_block(oligo), collapse="\n"))
}

`%<in>%` <- function(x,y) all(x %in% y) & all(y %in% x)
which_name <- function(tf) names(tf)[which(tf)]

subcodes <- function(code, codes=IUB_CODES) {
  which_name(sapply(codes, function(x) all(x %in% codes[[code]])))
}

complement <- function(strings, codes=IUB_CODES) {
  base_comps <- c(A="T", C="G", G="C", T="A")
  complemented <- lapply(codes, function(x) unname(base_comps[x]))
  comp_codes <- sapply(codes, function(x) "")
  for (code in names(comp_codes)) {
    comp_codes[code] <- which_name(sapply(codes, `%<in>%`, complemented[[code]]))
  }
  # comp_codes now contains the complement code for each code in codes
  # For IUB_CODES:
  #   A   C   G   T   R   Y   M   K   S   W   H   B   V   D   N
  #  "T" "G" "C" "A" "Y" "R" "K" "M" "S" "W" "D" "V" "B" "H" "N"

  chartr(paste(names(comp_codes), collapse=""),
         paste(comp_codes, collapse=""),
         strings)
}

reverse <- function(strings) {
  sapply(lapply(strsplit(strings, NULL), rev), paste, collapse="")
}

reverse_complement <- function(strings) reverse(complement(strings))

expand_asymmetric <- function(sites) {
  unique(c(sites, reverse_complement(sites)))
}

str2char <- function(x) strsplit(x, NULL)[[1]]
char2str <- function(x) paste(x, collapse="")

degeneracy <- function(sequence, codes=IUB_CODES) {
  dups <- vapply(codes, length, 0)
  sapply(sequence, function(x) prod(dups[str2char(x)]))
}

cutfree <- function(m=20, sites=c(), codes=IUB_CODES,
                    starting_oligo=strrep("N",m),
                    min_blocks=1,
                    obj_weights=log(1:4),
                    re_randomize=TRUE,
                    seed=NULL,
                    obj_frac=1.0,
                    quiet=FALSE,
                    maxtime=30) {
  sites <- expand_asymmetric(sites)  # include both strands if asymmetric
  ks <- sapply(sites, nchar)  # length of each site
  nB <- length(codes)  # number of variables for position in randomer
  starting_oligo <- str2char(starting_oligo)
  m <- length(starting_oligo)  # number of positions in randomer

  ncons <- m + sum(m - ks + 1)
  nvars <- nB * m
  idx <- function(b,i) nB*(i-1) + which(names(codes) == b)

  A <- matrix(0, nrow=ncons, ncol=nvars,
              dimnames=list(NULL,
                            paste0(names(codes), rep(1:m, each=nB))))
  rhs = rep(0, ncons)
  sense <- rep(">", ncons)
  sense[1:m] <- "="
  lb = numeric(nvars)
  ub = numeric(nvars)

  # open up only codes allowed by the starting oligo
  for (i in 1:m) {
    subs <- subcodes(starting_oligo[i])
    for (b in subs) {
      ub[idx(b,i)] <- 1
    }
  }

  # allow only a single code per location
  for (i in 1:m) {
    A[i,(nB*(i-1)+1):(nB*i)] <- 1
    rhs[i] <- 1
  }

  # add constraints to remove the RS
  curr_row <- m + 1
  for (s in seq_along(sites)) {
    rs <- str2char(sites[s])
    for (j in 1:(m-ks[s]+1)) {
      for (i in 1:ks[s]) {
        blocked <- get_blocking_codes(rs[i], names=TRUE)
        for (b in blocked) {
          A[curr_row,idx(b,i+j-1)] <- 1
        }
      }
      rhs[curr_row] <- min_blocks
      curr_row <- curr_row + 1
    }
  }

  # define the objective
  codelens <- obj_weights[vapply(codes, length, 0)]
  obj <- rep(codelens, m)

  model <- list(
    A = A,
    obj = obj,
    modelsense = "max",
    rhs = rhs,
    lb = lb,
    ub = ub,
    sense = sense,
    vtype = "B"
  )
  params <- list(OutputFlag=as.integer(!quiet))

  times <- system.time(
    orig_result <- gurobi::gurobi(model, params)
  )

  if (re_randomize && orig_result$status == "OPTIMAL") {
    model2 <- add_objective_constraint(model, orig_result, frac=obj_frac)
    set.seed(seed)
    model2$obj <- runif(length(model2$obj))
    result <- gurobi::gurobi(model2, params)
  } else {
    result <- orig_result
  }

  if (result$status == "OPTIMAL")
    code <- paste(rep(names(codes),m)[result$x[1:nvars] == 1], collapse="")
  else
    code <- NA

  return(list(code = code,
              model = model,
              time = times,
              sol_orig = orig_result,
              sol_final = result))
}

create_barcodes <- function(randomer, n=100, min_distance=5,
                            file=NULL, codes=IUB_CODES,
                            seed=NULL) {
  randomer <- str2char(randomer)
  m <- length(randomer)

  nvars <- 4*m
  ncons <- m + n
  nuc_names <- c("A", "C", "G", "T")
  idx <- function(b,i) 4*(i-1) + which(nuc_names == b)

  model = list(
    modelsense = "max",
    vtype = "B"
  )
  model$A <- matrix(0, nrow=ncons, ncol=nvars,
                    dimnames=list(NULL,
                                  paste0(nuc_names, rep(1:m, each=4))))
  model$rhs = rep(0, ncons)
  model$rhs[1:m] <- 1
  model$sense <- rep("=", ncons)
  model$lb = numeric(nvars)
  model$ub = numeric(nvars)

  # open up only codes allowed by the starting oligo
  for (i in 1:m) {
    allowed <- codes[randomer[i]]
    for (b in allowed) {
      model$ub[idx(b,i)] <- 1
    }
  }

  # allow only a single code per location
  for (i in 1:m) {
    model$A[i,(4*(i-1)+1):(4*i)] <- 1
  }

  barcodes <- as.character(rep(NA, n))

  params <- list(OutputFlag=0)
  set.seed(seed)
  for (i in 1:n) {
    model$obj <- runif(nvars) - 0.5
    elapsed <- system.time(
      result <- gurobi::gurobi(model, params),
      gcFirst = FALSE
    )
    if (result$status != "OPTIMAL") {
      printf("   No additional barcodes possible.\n")
      return(barcodes)
    }
    barcodes[i] <- paste(rep(nuc_names,m)[result$x == 1], collapse="")

    printf("   Found barcode %i (%s) in %f seconds.\n", i, barcodes[i], elapsed[3])

    if (i < n) {
      # add existing barcode as constraint
      model$A[m+i, ] <- result$x
      model$sense[m+i] <- "<"
      model$rhs[m+i] <- m - min_distance
    }
  }
  return(barcodes)
}


# result <- cutfree(m=20, sites=c("GGTCTC", "GATATC"), min_blocks=2)
# bcs <- create_barcodes("NNNNNNN", min_distance=3, n=100)
