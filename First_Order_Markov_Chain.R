#File Import Methods for Seqinr and Biostrings
library(seqinr)
s <- read.fasta(file = "C://path//to//sequence.fasta")
s <- s[[1]]
s

#Sequence Length
length(s)
#Nucleotide Frequencies:
table(s)
#Determine GC Content
GC(s)

require(Biostrings)
x <-readDNAStringSet(filepath = "C://path//to//sequence.fasta", nrec=-1L, skip=0L, use.names=TRUE)
x <- x[[1]]
x

#Nucleotide Frequencies:
alphabetFrequency(x, baseOnly = TRUE)
#ZeroOrder Markov:
alphabetFrequency(x, as.prob = TRUE, baseOnly = TRUE)
A <- 1115
C <- 647
G <- 773
T <- 1081

#Dinucleotide Frequencies:
dinucleotideFrequency(x)
dinucleotideFrequency(x, as.matrix = TRUE)

#FirstOrder Markov for "AN" = AN/A
A_first = matrix(c(341, 210, 237, 326)/A)          # Given nucleotide A in the 1st position
C_first = matrix(c(280, 136, 26, 205)/C)          # Given nucleotide C in the 1st position
G_first = matrix(c(220, 119, 198, 236)/G)      # Given nucleotide G in the 1st position
T_first = matrix(c(274, 182, 312, 313)/T)      # Given nucleotide T in the 1st position
L = list(A_first, C_first, G_first, T_first)
L

#Assign names to matrix for First Order Markov Transition Dinucleotide Probablities: 
transitionmatrix <- matrix(c(A_first, C_first, G_first, T_first), 4, 4, byrow = TRUE) # Create a 4 x 4 matrix
nucleotides <- c("A", "C", "G", "T")
rownames(transitionmatrix) <- nucleotides
colnames(transitionmatrix) <- nucleotides
transitionmatrix
