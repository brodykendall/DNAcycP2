from inst.python.dnacycp_python import pred

def main():
    # out = pred.cycle_fasta("inst/data/ex1.fasta", "inst/python/irlstm_smooth", 1000, 1)
    # out = pred.cycle_fasta("inst/data/ex1.fasta", "inst/python/irlstm_smooth", 1000, 2)
    # out = pred.cycle_fasta("inst/data/ex1.fasta", "inst/python/irlstm_smooth", 1000, 4)
    # out = pred.cycle_fasta("inst/data/ex1.fasta", "inst/python/irlstm_smooth", 1000, 8)
    # out = pred.cycle_fasta("/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/saccer3/ncbi_dataset/data/GCF_000146045.2/GCF_000146045.2_R64_genomic.fna", "inst/python/irlstm_smooth", 100000, 1)
    # out = pred.cycle_fasta("/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/saccer3/ncbi_dataset/data/GCF_000146045.2/GCF_000146045.2_R64_genomic.fna", "inst/python/irlstm_smooth", 100000, 8)
    # out = pred.cycle_fasta("/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/human/chr21.fa", "inst/python/irlstm_smooth", 100000, 1)
    out = pred.cycle_fasta("/Users/Brody1/Dropbox/Northwestern/DNA_Cyclizability/human/chr21.fa", "inst/python/irlstm_smooth", 100000, 8)
    print(out)

if __name__ == "__main__":
    main()