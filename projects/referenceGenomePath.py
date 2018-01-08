import os

class referenceGenomePath(object):
    def __init__(self):

        hg19_root = "/proj/seq/data/HG19_UCSC"
        hg19_custom_root = "/proj/sancarlb/users/ogun/seq/hg19"
        self.hg19 = {
            "name": "_hg19",
            "bowtie": os.path.join(hg19_root, "Sequence/BowtieIndex/genome"),
            "salmon_quasi": os.path.join(hg19_custom_root, 'transcriptome', 'Homo_sapiens.GRCh38.all.quasi_k31_index'),
            "fasta": os.path.join(hg19_root, "Sequence/WholeGenomeFasta/genome.fa"),
            "limits": os.path.join(hg19_root, "Sequence/WholeGenomeFasta/genome.fa.fai"),
            "transcripts": os.path.join(hg19_custom_root, "geneTranscripts_sorted.bed"),
            "TSS": "NA",
            "TES": "NA",
            # "genes": "NA",
            "chromatinStates": "NA",
            "LAD": os.path.join(hg19_custom_root, "LADcenter_w100bp_noY_bed6.bed") # http://compbio.med.harvard.edu/modencode/webpage/lad/human.fibroblast.DamID.hg19.bed

        }

        self.hg19nuc = self.hg19

        mm10_ucsc_root = "/proj/seq/data/MM10_UCSC/Sequence"
        mm10_custom_root = "/proj/sancarlb/users/ogun/seq/mm10"
        self.mm10 = {
            "name": "",
            "bowtie": os.path.join(mm10_ucsc_root, "BowtieIndex/genome"),
            "fasta": os.path.join(mm10_ucsc_root, "WholeGenomeFasta/genome.fa"),
            "limits": os.path.join(mm10_ucsc_root, "WholeGenomeFasta/genome.fa.fai"),
            "genes": os.path.join(mm10_custom_root, "geneList_yy_6.bed"),
            "genesNR": os.path.join(mm10_custom_root, "mm10_genesNR.bed"),
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
            "limits_bed": os.path.join(TAIR10_custom_root, "genome.bed"),
            "transcripts": os.path.join(TAIR10_custom_root, "rna_singleIsoform_chr.bed"),       
            "genes": os.path.join(TAIR10_custom_root, "protein_coding_genes_chr.bed"),
            "chromatinStates": os.path.join(TAIR10_custom_root, "chromatinStatesSorted.bed"),
            # "chmm36": os.path.join(TAIR10_custom_root, "Liu2017", "chrStates_noEdge_short.bed"),
            "chmm36": os.path.join(TAIR10_custom_root, "Liu2017", "chrStates_noEdge.bed"),
            "chmm36_shuffled": os.path.join(TAIR10_custom_root, "Liu2017", "chrStates_noEdge_shuffled.bed"),
            "epigeneticMarkers": os.path.join(TAIR10_custom_root, "Liu2017", "epigeneticMarks", "epigeneticMarkers.bed"),
            "epigeneticMarkers_shuffled": os.path.join(TAIR10_custom_root, "Liu2017", "epigeneticMarks", "epigeneticMarkers_shuffled.bed"),
            "CCA1potentialGenes": os.path.join(TAIR10_custom_root, "Liu2017", "epigeneticMarks", "CCA1-SRX997425_potentialGenes.bed"),
            # "dnase": os.path.join(TAIR10_custom_root, "DHSs_chr_center1kb_w10.bed"),
            "dnase": os.path.join(TAIR10_custom_root, "DHSs_chr.bed"),
            "TSS": os.path.join(TAIR10_custom_root, "Hetzel2016_noHeader_sorted_NascentRNA_noOv2K_gt2K_w100_TSS.bed"),
            "TES": os.path.join(TAIR10_custom_root, "Hetzel2016_noHeader_sorted_NascentRNA_noOv2K_gt2K_w100_TES.bed"),
            "TSS_Q1": os.path.join(TAIR10_custom_root, "Hetzel2016_noHeader_sorted_NascentRNA_noOv2K_gt2K_sortByScore.bed_Q1.bed_w100_TSS.bed"),
            "TES_Q1": os.path.join(TAIR10_custom_root, "Hetzel2016_noHeader_sorted_NascentRNA_noOv2K_gt2K_sortByScore.bed_Q1.bed_w100_TES.bed"),
            "TSS_Q2": os.path.join(TAIR10_custom_root, "Hetzel2016_noHeader_sorted_NascentRNA_noOv2K_gt2K_sortByScore.bed_Q2.bed_w100_TSS.bed"),
            "TES_Q2": os.path.join(TAIR10_custom_root, "Hetzel2016_noHeader_sorted_NascentRNA_noOv2K_gt2K_sortByScore.bed_Q2.bed_w100_TES.bed"),
            "TSS_Q3": os.path.join(TAIR10_custom_root, "Hetzel2016_noHeader_sorted_NascentRNA_noOv2K_gt2K_sortByScore.bed_Q3.bed_w100_TSS.bed"),
            "TES_Q3": os.path.join(TAIR10_custom_root, "Hetzel2016_noHeader_sorted_NascentRNA_noOv2K_gt2K_sortByScore.bed_Q3.bed_w100_TES.bed"),
            "TSS_Q4": os.path.join(TAIR10_custom_root, "Hetzel2016_noHeader_sorted_NascentRNA_noOv2K_gt2K_sortByScore.bed_Q4.bed_w100_TSS.bed"),
            "TES_Q4": os.path.join(TAIR10_custom_root, "Hetzel2016_noHeader_sorted_NascentRNA_noOv2K_gt2K_sortByScore.bed_Q4.bed_w100_TES.bed"),
            "transcripts2": os.path.join(TAIR10_custom_root, "Hetzel_sorted_noEdge10K.bed"),
            "transcriptsNoNeighbor5K": os.path.join(TAIR10_custom_root, "Hetzel_sorted_noEdge10K_noNeighborIn5K.bed"),
            # "transcripts2": os.path.join(TAIR10_custom_root, "Hetzel2016_noHeader_sorted_NascentRNA_noOv2K_gt2K.bed"),
            "originOfReplication": os.path.join(TAIR10_custom_root, "Gutierrez_DNA_replication_origin_TAIR10_GBROWSE.bed"),
            "originOfReplication_shuffled": os.path.join(TAIR10_custom_root, "Gutierrez_DNA_replication_origin_TAIR10_GBROWSE_shuffled.bed"),
            "geo": "/pine/scr/a/d/adebali/GEO/oadebali@gmail.com/"
        }

        TAIR9_root = "/proj/seq/data/TAIR9/Arabidopsis_thaliana/Ensembl/TAIR9"
        TAIR9_custom_root = "/proj/sancarlb/users/ogun/seq/TAIR9"
        self.TAIR9 = {
            "name": "_TAIR9",
            "bowtie": os.path.join(TAIR9_root, "Sequence/BowtieIndex/genome"),
            "fasta": os.path.join(TAIR9_root, "Sequence/WholeGenomeFasta/genome.fa"),
            "limits": os.path.join(TAIR9_root, "Sequence/WholeGenomeFasta/genome.fa.fai"),
            "transcripts": os.path.join(TAIR9_custom_root, "transcripts.bed"),
            "TSS": os.path.join(TAIR9_custom_root, "transcripts_noOv4K_gt4K_w100_TSS.bed"),
            "TES": os.path.join(TAIR9_custom_root, "transcripts_noOv4K_gt4K_w100_TES.bed"),
            # "genes": "NA",
            "chromatinStates": os.path.join(TAIR9_custom_root, "chromatinStatesSorted.bed")
        }

        NC_000913_2_root = "/nas02/home/a/d/adebali/ncbi/ecoli"
        self.NC_000913_2 = {
            "name": "NC000913.2",
            "bowtie": os.path.join(NC_000913_2_root, "NC_000913.2", "NC_000913.2"),
            "fasta": os.path.join(NC_000913_2_root, "NC_000913.2", "NC_000913.2.chr.fa"),
            "geneTT": os.path.join(NC_000913_2_root, "NC_000913.2", "NC_000913.2.genes.AA_TT.bed"),
            "geneNames": os.path.join(NC_000913_2_root, "NC_000913.2", "NC_000913.2.genes.txt"),
            "limits": os.path.join(NC_000913_2_root, "NC_000913.2", "NC_000913.2.chr.sizes"),
            "genes": os.path.join(NC_000913_2_root, "NC_000913.2", "NC_000913.2.chr.genes.bed6.bed"),
            "lacOperon": os.path.join(NC_000913_2_root, "NC_000913.2", "NC_000913.2.chr.lacOperon.bed")
        }

        S288C_R64_2_1_custom_root = "/proj/sancarlb/users/ogun/seq/S288C/S288C_reference_genome_R64-2-1_20150113"
        self.S288C_R64_2_1 = {
            "name": "S288C_R64_2_1",
            "bowtie": "/proj/sancarlb/users/wentao/S288C_reference_genome/Bt1_SacCer3/yeast_bt1.SacCer3",
            "genes": os.path.join(S288C_R64_2_1_custom_root, "genes.bed"),
            "scoredGenes": "/proj/sancarlb/users/ogun/sancarlabutils/projects/yeast/RNA/dataDir/1801/SRR2005997.cu.bo.coToBa.coToBe.unSo.geMa.noCo.adTr.maScBe.bed",
            "genes_noNeighborIn500bp": os.path.join(S288C_R64_2_1_custom_root, "genes_noNeighborIn500bp.bed"),
            "fasta": os.path.join(S288C_R64_2_1_custom_root, "genome.fa"),
            "limits": os.path.join(S288C_R64_2_1_custom_root, "genome.fa.fai")
        }

    def get(self, key): return self.__getattribute__(key)

if __name__ == "__main__":
    path = referenceGenomePath()
    # print(path.get("mm9"))
    