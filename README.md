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



