# ENA submission
Simple python script to write and validate files for ENA upload from outputs of skim2mito and skim2rrna.

## Overview
Script takes the outputs of skim2mito/skim2ribo & prepares everything needed for submission to ENA. It generates EMBL-format flatfiles from FASTA + GFF3 annotations, manifest files, and chromosome list files, validates them with ena-webin-cli, and produces the submission commands.

Requies samples to be registered on ENA first.

## Set up conda environment
```
conda create -n ena_submission_env -c conda-forge -c bioconda emblmygff3=2.4 ena-webin-cli=9.0.1 openjdk=17.0.15 pyyaml pandas --yes
conda activate ena_submission_env
```

## Run code 
```
python3 ena_submission.py \
    --basedir .test/skim2mito/ \
    --samples .test/sample_list.txt \
    --output example_output_skim2mito \
    --config .test/skim2mito/config/config.yaml \
    --project <project_id> \
    --locus_tag <locus_tag> \
    --overwrite \
    --ena_username <username> \
    --ena_password <password> 
```


## Arguments
--basedir — Root directory of the pipeline output. The script expects the specific skim2mito subdirectories beneath this

--pipeline — Either skim2mito or skim2ribo. This controls branching logic throughout: which GFF features get counted, how assembly type and topology are determined, what goes into the manifest, and whether chromosome lists are generated.

--samples — A TSV with three columns: ID (your internal sample identifier used to find files), scientificName (the species name, passed to EMBLmyGFF3 and used in description lines), and SampleID (the ENA sample accession like an ERS number).

--config — Path to a Snakemake config YAML. It's loaded but interestingly never seemingly used anywhere in the current code.

--project — ENA study/project accession (e.g. PRJEB...), written into every manifest.

--locus_tag — Locus tag prefix registered with ENA, passed to EMBLmyGFF3 for feature naming. Seems like no specific rules, but ENA have some recommendations - https://ena-docs.readthedocs.io/en/latest/faq/locus_tags.html

--output — Where all generated files go. Gets wiped and recreated if --overwrite is set; otherwise the script exits if the directory already exists.

--ena_username / --ena_password — Credentials passed directly to ena-webin-cli for validation (and would be used for submission).

--platform — Sequencing platform string for the manifest, defaults to Illumina.

--assembly_type — Either isolate or clone, defaults to isolate. Written into the manifest's ASSEMBLYTYPE field. ENA mention here - https://ena-docs.readthedocs.io/en/latest/submit/assembly/genome.html

--trans_table - Specify the genetic code translation table to be used in EMBLmyGFF3. Default: 4

## Process
First verifies that EMBLmyGFF3 is on PATH (exits with an error if not), that the base directory exists, and that the samples TSV and config YAML can be read.

For each row in the samples file, the script then:
1. Locates and validates input files (in {basedir}/results/assembled_sequence/{sample_id}.fasta & finds associated .gff files).
2. Merges and filters GFF annotations. All GFF files concatenated into a temporary .gff3 file, stripping comment lines. Any annotation where start>end (indicating a feature that bridges the circular breakpoint) is excluded with a warning (comment in script suggesting this might be addressed later by rotating the sequence).
3. Counts relevant annotations. For skim2mito, it counts gene features. For skim2ribo, it counts rRNA features. The names are extracted from the Name= attribute.
4. Determines assembly type and topology. For skim2mito, a single sequence with "circular" in its ID → mitochondrial complete (circular topology); a single sequence without → mitochondrial partial (linear); multiple sequences → mitochondrial contig (linear). For skim2ribo, if there are exactly three rRNA annotations matching the set {5_8S_rRNA, 18S_rRNA, 28S_rRNA} → ribosomal RNA complete; otherwise → ribosomal RNA partial. The description line is also constructed here — for skim2ribo it's a formatted string with species name and gene names. For skim2mito it's currently hardcoded as 'XXX' (presumably TODO???).
5. Runs EMBLmyGFF3. Converts the merged GFF3 + FASTA into an EMBL-format flatfile. Translation table is selected. Output and logs go into a per-sample subdirectory under the output folder.
6. Post-processes the EMBL file. For ribo assemblies, lines starting with DT (date lines) are stripped out to avoid an ena-webin-cli validation error. The EMBL file is then gzipped.
7. Extracts coverage (mitochondrial). For complete or partial mitochondrial assemblies, it reads the blobtools table at {basedir}/results/blobtools/{sample_id}/table.tsv, takes coverage column, and computes the mean. If the blobtools file is missing, the script exits.
8. Writes a chromosome list (mitochondrial complete/partial only). A gzipped file mapping sequence IDs to chromosome types. For circular assemblies it writes a single MIT entry; for linear ones it writes one entry per contig with incrementing MIT0, MIT1, etc.
9. Writes the manifest file. The content varies by assembly type. Mitochondrial manifests include SAMPLE, ASSEMBLYNAME, COVERAGE, PROGRAM (hardcoded as GetOrganelle), PLATFORM, and ASSEMBLYTYPE. Ribosomal manifests include NAME instead. Both include STUDY, DESCRIPTION, and FLATFILE. Only complete/partial mitochondrial assemblies get a CHROMOSOME_LIST entry.
10. Validates with ena-webin-cli. Runs the CLI in -validate mode with -context genome for mitochondrial or -context sequence for ribosomal assemblies. Output goes to a log file.
11. Checks validation result. If the log contains "Submission(s) validated successfully.", the corresponding submit command is added to a list. The pass/fail status is written to the summary.

## Outputs
- Per sample: a subdirectory containing the merged .gff3, the .embl.gz flatfile, an EMBLmyGFF3 log, the manifest file, a validation log, and (for mitochondrial) a gzipped chromosome list.
- Overall: an assembly_summary.tsv in the output dir with columns for sample ID, assembly type, sequence count, total length, annotation count, and validation status. At the end, all validated submission commands are printed to stdout (so next, you'd pipe or copy them to actually run the submissions).

