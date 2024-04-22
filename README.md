
# sgrna_selector.py 

Selects sgRNA to generate a library, currently only supports Human. sgrna_selector 
takes sgRNAs from several genome-wide sgRNA libraries.

# Installation

None required, just use python sgrna_selector.py

NOTE: Currently it requires glbase3, I will remove that in future versions

# Usage

Usage: sgrna_selector.py genenames_file.txt libraryname num_sgrnas_per_gene

# Example

```
sgrna_selector % python sgrna_selector.py test_data.tsv test 3
INFO    : NumExpr defaulting to 4 threads.
INFO    : genelist: loaded 'test_data.tsv' found 803 items
INFO    : genelist: loaded 'datasets/CEG2.txt.gz' found 683 items
INFO    : genelist: loaded 'datasets/sgRNA_libraries/MinLibCas9.txt.gz' found 37,722 items
INFO    : addEmptyKey: Added a new key 'original_library'
INFO    : genelist: loaded 'datasets/sgRNA_libraries/Brunello.txt.gz' found 76,441 items
INFO    : addEmptyKey: Added a new key 'original_library'
INFO    : genelist: loaded 'datasets/sgRNA_libraries/tkov3_guide_sequence.txt.gz' found 71,090 items
INFO    : addEmptyKey: Added a new key 'original_library'
INFO    : genelist: loaded 'datasets/sgRNA_libraries/Source_sgRNA_table.txt.gz' found 354,715 items
INFO    : addEmptyKey: Added a new key 'original_library'
INFO    : genelist: loaded 'datasets/sgRNA_libraries/CRISPRab.packed.tsv.gz' found 146,154 items
INFO    : addEmptyKey: Added a new key 'original_library'
INFO    : addEmptyKey: Added a new key 'rank'
INFO    : addEmptyKey: Added a new key 'rank'
INFO    : addEmptyKey: Added a new key 'rank'
INFO    : addEmptyKey: Added a new key 'rank'
INFO    : addEmptyKey: Added a new key 'rank'
INFO    : map: 'test_data' vs 'MinLibCas9.txt', using 'name', found: 31733 items
Rescued ABRAXAS1/FAM175A
Rescued BABAM2/BRE
Rescued ELP1/IKBKAP
Rescued EMSY/C11orf30
Rescued KAT14/CSRP2BP
Rescued KMT5A/SETD8
Rescued KMT5B/SUV420H1
Rescued KMT5C/SUV420H2
Rescued NSD2/WHSC1
Rescued NSD3/WHSC1L1
Rescued OGA/MGEA5
Rescued PPP4R3A/SMEK1
Rescued PPP4R3B/SMEK2
Rescued RIOX2/MINA
Rescued RTL5/RGAG4
Rescued SGF29/CCDC101
Rescued SHTN1/KIAA1598
Rescued SLF1/ANKRD32
Rescued TOMM70/TOMM70A
Rescued PWWP3A/MUM1
Rescued AC106886.5/SRCAP
INFO    : removeDuplicates: 14119 duplicates, list now 17816 items long
INFO    : Saved 'sgrna_test-all.tsv'
INFO    : removeDuplicates: 17013 duplicates, list now 803 items long
INFO    : getColumns: got only the columns: name
INFO    : Saved 'sgrna_test-genes.tsv'
INFO    : removeDuplicates: 17013 duplicates, list now 803 items long
INFO    : map: 'CEG2.txt' vs 'MinLibCas9.txt', using 'name', found: 1166 items
INFO    : removeDuplicates: 1113 duplicates, list now 53 items long

CEG2 match: CDC73
CEG2 match: CHAF1A
CEG2 match: CHAF1B
CEG2 match: CRNKL1
CEG2 match: CTR9
CEG2 match: DDX10
CEG2 match: DDX18
CEG2 match: DDX20
CEG2 match: DDX21
CEG2 match: DDX27
CEG2 match: DDX41
CEG2 match: DDX47
CEG2 match: DDX49
CEG2 match: DDX55
CEG2 match: DDX56
CEG2 match: DMAP1
CEG2 match: DNMT1
CEG2 match: DYNC1I2
CEG2 match: EIF4A3
CEG2 match: GEMIN5
CEG2 match: GRWD1
CEG2 match: GTF3C2
CEG2 match: HCFC1
CEG2 match: HDAC3
CEG2 match: HJURP
CEG2 match: KANSL3
CEG2 match: KAT8
CEG2 match: LAS1L
CEG2 match: MED12
CEG2 match: MYBBP1A
CEG2 match: NAA50
CEG2 match: NAT10
CEG2 match: OGT
CEG2 match: PABPC1
CEG2 match: PCNA
CEG2 match: PRMT1
CEG2 match: PRMT5
CEG2 match: RAD51C
CEG2 match: RAD51D
CEG2 match: RPS12
CEG2 match: RUVBL2
CEG2 match: SAP18
CEG2 match: SART3
CEG2 match: SNAPC2
CEG2 match: SRRM1
CEG2 match: SRSF1
CEG2 match: SS18L2
CEG2 match: SUPT6H
CEG2 match: TONSL
CEG2 match: TOP2A
CEG2 match: UPF1
CEG2 match: WDR77
CEG2 match: XAB2



Full set:
  All genes                        : 806
  Genes with at least 1 sgRNA      : 803
  Number of sgRNAs                 : 17816
  Average sgRNAs/gene              : 22.2
  Number of genes with <= 2 sgRNAs : 0
  Number of genes with > 5 sgRNAs  : 792
  CEG2 hits                        : 53

INFO    : Saved 'sgrna_test-final.tsv'
INFO    : Renamed key 'sgrna' to 'seq'
INFO    : Saved FASTA file: sgrna_test-final.fasta

Final set:
  Number of sgRNAs                 : 2541
  Average sgRNAs/gene              : 3.2
  Total number of CTRL sgRNAs      : 141
```


