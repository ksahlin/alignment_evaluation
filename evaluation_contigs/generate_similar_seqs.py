import random
import sys

if len(sys.argv) != 4:
    print("Usage: python generate_similar_seqs.py <nr_seqs> <length> <SNV_rate> <SV_rate> <outfile> ") 

n = int(sys.argv[1])
L = int(sys.argv[2])
mut_freq = float(sys.argv[3])
sv_freq = float(sys.argv[4])
outfile = sys.argv[5]

f = open(outfile, "w")
# L = 100000
# n = 20
seq = "".join([random.choice("ACGT") for i in range(L)])
for i in range(n):
  f.write(">contig{0}\n".format(i))
  muts = set(random.sample(range(len(seq)),int(L*mut_freq)))
  seq_mut = "".join([seq[i] if i not in muts else random.choice(['', random.choice("ACGT"), seq[i] + random.choice("ACGT")]) for i in range(len(seq))]) 
  
  # large indels (exons)
  muts = set(random.sample(range(len(seq)),int(L*sv_freq)))
  seq_mut2 = seq_mut
  for j in sorted(muts):
    if j < len(seq_mut2) - 1001:
      seq_mut2 = seq_mut2[:j] + seq_mut2[j+ random.randint(300,1000):] 
  f.write("{0}\n".format(seq_mut2))
#  for chunk in range(0,len(seq_mut)-1, 60):
#    f.write("{0}\n".format(seq_mut[chunk:chunk+60]))



