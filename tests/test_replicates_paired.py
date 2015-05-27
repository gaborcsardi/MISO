
import setup_miso
import misopy
import misopy.gff_utils
import pysplicing
import random as plainrandom

from scipy import *
from numpy import *

gene = pysplicing.createGene( ((1, 100), (201, 300), (401, 500)),
                              ((0, 1, 2), (0, 2)) )

def do_replicates_paired(reads):
    result = pysplicing.doMISOPaired(
        GFF = gene, gene = 0L, reads = reads, read_len = 33L,
        mean_frag_len = 60, frag_variance = 5, num_sds = 4,
        num_iters = 5000L, burn_in = 1000L, lag = 10L,
        prior = pysplicing.MISO_PRIOR_LOGISTIC,
        replicate_mean_prior_mean = 0.0, replicate_mean_prior_var = 100000,
        replicate_var_prior_numobs = 0.0001, replicate_var_prior_var = 0.0001)

    psi = transpose(array(result[0]))
    meanpsi = mean(psi, 0)
    varpsi = var(psi, 0)
    print "rep mean", meanpsi, "var", varpsi
    assert(abs(varpsi[0] - 0.01172719) < 1e-8 and
           abs(varpsi[1] - 0.01172719) < 1e-8)
    return meanpsi

def do_one_paired(reads):
    results = pysplicing.doMISOPaired(
        gene, 0L, (reads,), read_len = 33L, num_iters = 5000L,
        burn_in = 1000L, lag = 10L, mean_frag_len = 60, frag_variance = 5,
        num_sds = 4)
    psi = transpose(array(results[0]))
    return psi

def do_onebyone_paired(reads):
    psi = [ do_one_paired(r) for r in reads ]
    meanpsi = matrix([ mean(p, 0) for p in psi ])
    varpsi = var(meanpsi, 0).ravel().tolist()[0]
    print "one mean", meanpsi, "var", varpsi
    assert(abs(varpsi[0] - 0.023567944694320754) < 1e-8 and
           abs(varpsi[1] - 0.023567944694320754) < 1e-8)
    return mean(meanpsi, 0).ravel().tolist()[0]

def test_replicates():
    plainrandom.seed(42)
    reads1 = pysplicing.simulatePairedReads(gene, 0L, (0.1, 0.9), 100L, 33L, 60, 5, 4)
    reads2 = pysplicing.simulatePairedReads(gene, 0L, (0.2, 0.8), 100L, 33L, 60, 5, 4)
    reads3 = pysplicing.simulatePairedReads(gene, 0L, (0.3, 0.7), 100L, 33L, 60, 5, 4)
    reads4 = pysplicing.simulatePairedReads(gene, 0L, (0.4, 0.6), 100L, 33L, 60, 5, 4)
    reads5 = pysplicing.simulatePairedReads(gene, 0L, (0.5, 0.5), 100L, 33L, 60, 5, 4)

    reads = (tuple(reads1[1:3]), tuple(reads2[1:3]), tuple(reads3[1:3]),
             tuple(reads4[1:3]), tuple(reads5[1:3]))
    psi_rep = do_replicates_paired(reads)
    psi_one = do_onebyone_paired(reads)

    print "psi_rep", psi_rep
    assert(abs(psi_rep[0] - 0.28372519) < 1e-8 and
           abs(psi_rep[1] - 0.71627481) < 1e-8)

    print "psi_one", psi_one
    assert(abs(psi_one[0] - 0.29992007896967365) < 1e-8 and
           abs(psi_one[1] - 0.7000799210303265) < 1e-8)
