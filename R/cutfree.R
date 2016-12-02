
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
print_oligo_block <- function(oligo) {
  cat(paste(oligo, collapse=""), "\n")
  cat(paste(make_oligo_block(oligo), collapse="\n"))
}

make_random_oligo <- function(m=20, rs="GGTCTC", codes=IUB_CODES) {
  k <- nchar(rs)
  nB <- length(codes)
  ncons <- 2*m - k + 1
  A <- matrix(0, nrow=ncons, ncol=nB*m, dimnames=list(NULL, paste0(names(codes), rep(1:m, each=nB))))
  rhs <- rep(0, ncons)

  sense <- c(rep("=", m), rep(">", m-k+1))

  # allow only a single code per location
  for (i in 1:m) {
    A[i,(nB*(i-1)+1):(nB*i)] <- 1
    rhs[i] <- 1
  }

  # add constraints to remove the RS
  idx <- function(b,i) nB*(i-1) + which(names(codes) == b)
  rs <- str_to_vector(rs)
  for (j in 1:(m-k+1)) {
    for (i in 1:k) {
      blocked <- get_blocking_codes(rs[i], names=TRUE)
      for (b in blocked) {
        A[m+j,idx(b,i+j-1)] <- 1
      }
    }
    rhs[m+j] <- 1
  }

  # define the objective
  codelens <- vapply(codes, length, 0)
  obj <- rep(codelens, m)

  model <- list(
    A = A,
    obj = obj,
    modelsense = "max",
    rhs = rhs,
    sense = sense,
    vtype = "B"
  )
  params <- list(OutputFlag=1)

  result <- gurobi::gurobi(model, params)
  code <- rep(names(codes),m)[result$x == 1]
  return(list(code = code,
              model = model))
}

result <- make_random_oligo()
