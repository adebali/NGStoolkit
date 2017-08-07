import os

class referenceGenomePath(object):
    def __init__(self):
        mm10_ucsc_root = "/proj/seq/data/MM10_UCSC/Sequence"
        mm10_custom_root = "/proj/sancarlb/users/ogun/seq/mm10"

        self.mm10 = {
            "name": "",
            "bowtie": os.path.join(mm10_ucsc_root, "BowtieIndex/genome"),
            "fasta": os.path.join(mm10_ucsc_root, "WholeGenomeFasta/genome.fa"),
            "limits": os.path.join(mm10_ucsc_root, "WholeGenomeFasta/genome.fa.fai"),
            "genes": os.path.join(mm10_custom_root, "geneList_yy_6.bed"),
            "chmm": {
                "liver": "NA",
                "spleen": "NA",
                "kidney": "NA",
                "testes": "NA"
            }
        }

        mm9_ucsc_root = "/proj/seq/data/MM9_UCSC/Sequence"
        mm9_custom_root = "/proj/sancarlb/users/ogun/seq/mm9"
        mm9_chromatinStates = "/nas/longleaf/home/adebali/ogun/seq/mm9/chromatin_states_chromHMM_mm9"
        self.mm9 = {
            "name": "_mm9",
            "bowtie": os.path.join(mm9_ucsc_root, "BowtieIndex/genome"),
            "fasta": os.path.join(mm9_ucsc_root, "WholeGenomeFasta/genome.fa"),
            "limits": os.path.join(mm9_ucsc_root, "WholeGenomeFasta/genome.fa.fai"),
            "genes": "NA",
            "chmm": {
                "liver": os.path.join(mm9_chromatinStates, "liver_cStates_HMM"),
                "spleen": os.path.join(mm9_chromatinStates, "spleen_cStates_HMM"),
                "kidney": os.path.join(mm9_chromatinStates, "kidney_cStates_HMM"),
                "testes": os.path.join(mm9_chromatinStates, "testes_cStates_HMM")
            }
        }
        
        TAIR10_root = "/proj/seq/data/TAIR10_NCBI"
        TAIR10_custom_root = "/proj/sancarlb/users/ogun/seq/TAIR10"
        self.TAIR10 = {
            "name": "_TAIR10",
            "bowtie": os.path.join(TAIR10_root, "Sequence/BowtieIndex/genome"),
            "fasta": os.path.join(TAIR10_root, "Sequence/WholeGenomeFasta/genome.fa"),
            "limits": os.path.join(TAIR10_root, "Sequence/WholeGenomeFasta/genome.fa.fai"),
            "genes": "NA",
            "chmm": "NA"
        }

        TAIR9_root = "/proj/seq/data/TAIR9/Arabidopsis_thaliana/Ensembl/TAIR9"
        TAIR9_custom_root = "/proj/sancarlb/users/ogun/seq/TAIR9"
        self.TAIR9 = {
            "name": "_TAIR9",
            "bowtie": os.path.join(TAIR9_root, "Sequence/BowtieIndex/genome"),
            "fasta": os.path.join(TAIR9_root, "Sequence/WholeGenomeFasta/genome.fa"),
            "limits": os.path.join(TAIR9_root, "Sequence/WholeGenomeFasta/genome.fa.fai"),
            "transcripts": os.path.join(TAIR9_custom_root, "transcripts.bed"),
            "TSS": os.path.join(TAIR9_custom_root, "transcripts_gt4K_sorted_TSS_w100.bed"),
            "TES": os.path.join(TAIR9_custom_root, "transcripts_gt4K_sorted_TES_w100.bed"),
            "genes": "NA",
            "chmm": "NA"
        }

    def get(self, key): return self.__getattribute__(key)

if __name__ == "__main__":
    path = referenceGenomePath()
    # print(path.get("mm9"))
    