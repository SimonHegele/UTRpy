<p align="center">
  <img src="UTRpyLogo.png" alt="Meine Bildunterschrift" width="300"/>
</p>

<p align="center">
  <img src="UTRpy.drawio.png" alt="Meine Bildunterschrift" width="500"/>
  <br>
  <em>UTRpy workflow</em>
</p>

Extending flanking transcript exons in genome annotations from protein orthology based gene<br>
prediction using exons in genome annotations from reference based transcriptome assembly.<br>
These extensions are the Un-Translated Regions (UTRs) of the transcripts.

## Workflow

**Step 1:**<br>
Identifying pairs of exons where exons from the transcriptome assembly overlap exons from the<br>
gene prediction.

**Step 2:**<br>
Deciding if in a pair of overlapping exons the exon from the transcriptome assembly (exon_ta) is<br>
an UTR-extension for the exon from the gene prediction (exon_gp) based on the following criteria:
1. The length of exon_ta does not exceed a certain limit
2. Their strand information must match (More details in the Parameters-section)
3. Exon_ta and exon_gp have either shared start or end positions
4. Exon_ta actually extends tran_gp

**Step 3:**<br>
Merging the annotations

## Installation

Tested with Python 3.12.8, but any Python >= 3.10 should be fine too.

```
git clone https://github.com/SimonHegele/UTRpy
cd UTRpy
pip install .
```

## Usage

```
usage: utrpy [-h] [-mel MAXIMUM_EXON_LENGTH] [-meo MINIMUM_EXON_OVERLAP] [-t THREADS] [-c] gff_prediction gff_assembly gff_utrpy

Genome annotions from protein orthology based gene prediction tools lack UTRs. Exons in these kind of gene annotations can be extended toalso cover UTRs using exons generated with reference based assembly tools
such as StringTie.

positional arguments:
  gff_prediction        GFF-format genome annotation from gene prediction.
  gff_assembly          GFF-format genome annotation from transcriptome assembly.
  gff_utrpy             Output file path

options:
  -h, --help            show this help message and exit
  -mel MAXIMUM_EXON_LENGTH, --maximum_exon_length MAXIMUM_EXON_LENGTH
                        Maximum exon length prevents the use of unreasonably long exons. [Default:20000]
  -meo MINIMUM_EXON_OVERLAP, --minimum_exon_overlap MINIMUM_EXON_OVERLAP
                        Minimum overlap of exons to be considered [Default:3]
  -t THREADS, --threads THREADS
                        Number of threads to use [Default:4]
  -c, --conservative    Will not use exons if their strand is not known [Default:True]
```

### Parameters

1. **maximum_exon_length**<br>
To my best knowlege the longest confirmed exon is of the MUC16 gene and has a length of<br>
~21kb (https://doi.org/10.1093/nar/gks652). One might therefore pose a limit on the exon<br>
lengths to prevent the use of false exons.
2. **minimum_exon_overlap**<br>
Minimum overlap length of pairs of exons to be further investigated controls sensitivity.<br>
4. **conservative**<br>
If True (default), exons from the gene prediction are only replaced with exons from the<br>
transcriptome assembly if they are known to be of the same strand.<br>
If False it is sufficient if they are not known to be of different strands (The strand of<br>
one or both of them might not be known).

### Output

The gene prediction gff with extended Exons:

Exons:<br>
1	seqid	<seqid_transcriptome_assembly> (same as <seqid_gene_prediction>)<br>
2	source  <source_gene_prediction> + <source_transcriptome_assembly> (UTRpy)<br>
3	type	exon<br>
4	start	<start_transcriptome_assembly><br>
5	end	<end_transcriptome_assembly><br>
6	score	<score_transcriptome_assembly><br>
7	strand <strand_transcriptome_assembly> (same as <strand_transcriptome_assembly> if conservative=True)<br>
8	phase	.<br>
9	attributes  <attributes_gene_prediction><br>

Transcripts<br>
Start and end positions updated according to the exon extension

Genes<br>
Start and end positions updated according to the exon extension

**Example:**

Gene prediction input exon:     
`Chr_1 AUGUSTUS    exon    2207399   2208325   .   -   .   ID=agat-exon-51;Parent=ga_chond_ext_ncbi_g49.t1;gene_id=g_p_13;transcript_id=ga_chond_ext_ncbi_g49.t1`

Transcriptome assembly input exon:     
`Chr_1  StringTie   exon    2206428 2208325 1000.0  -   .   transcript_id "STRG.130.1"; gene_id "STRG.130";`

UTRpy output exon:    
`Chr_1 AUGUSTUS + StringTie (UTRpy)    exon    2206428   2208325 1000.0  -   .   ID=agat-exon-51;Parent=ga_chond_ext_ncbi_g49.t1;gene_id=g_p_13;transcript_id=ga_chond_ext_ncbi_g49.t1`
