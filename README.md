# ENA submission
Simple python script to write and validate files for ENA upload from outputs of skim2mito and skim2rrna.

#### Set up conda environment
```
conda create -n ena_submission_env -c conda-forge -c bioconda emblmygff3=2.4 ena-webin-cli=9.0.1 openjdk=17.0.15 pyyaml pandas --yes
conda activate ena_submission_env
```

#### Run code 
```
python3 ena_submission.py \
   -s .test/sample_list.txt \
   -c .test/skim2mito/config/config.yaml \
   -o example_output \
   -b .test/skim2mito \
   --project PRJEB76850 \
   --platform Illumina \
   --moleculetype "genomic DNA" \
   --overwrite \
   --ena_username Webin-67784 \
   --ena_password xyx

```



