#!/usr/bin/env python3
"""
Create ENA submission files for mitochondrial genome or ribosomaml assemblies from the skim2mito and skim2rrna pipeline outputs.
"""

# load modules
import argparse
from Bio import SeqIO
import yaml
import pandas as pd
import os
import subprocess
import sys
import gzip
import shutil

# argument parser
parser = argparse.ArgumentParser(
    description="""
    Create ENA submission files for mitochondrial genome or ribosomal gene assemblies from the skim2mito and skim2rrna pipeline outputs.
    """,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""
Example usage:
  %(prog)s \\
    --basedir .test/skim2mito/ \\
    --samples .test/sample_list.txt \\
    --output example_output/ \\
    """
    )

parser.add_argument('--basedir', required=True,
    help='Base directory for skim2mito or skim2rrna with input files')  

parser.add_argument('--samples', required=True,
    help='Path to TSV file with sample information (ID, scientificName, SampleID)')

parser.add_argument('--config', required=True,
    help='Path to snakemake config YAML file')

parser.add_argument('--project', required=True,
    help='ENA project accession')

parser.add_argument('--locus_tag', required=True,
    help='Locus tag prefix for ENA submission')

parser.add_argument('--output', required=True,
    help='Output directory for generated files')

parser.add_argument('--overwrite', action='store_true',
    help='Overwrite existing output directory')

parser.add_argument('--ena_username', required=True, 
    help='ENA username for webin-cli')

parser.add_argument('--ena_password', required=True,
    help='ENA password for webin-cli')

parser.add_argument('--platform', required=False,
    help='Sequencing platform',
    default='Illumina')

parser.add_argument('--assembly_type', default='isolate',
    help='Assembly type (default: isolate)',
    choices=['isolate', 'clone'])
    
args = parser.parse_args()


# check if EMBLmyGFF3 is available
try:
    subprocess.run(['EMBLmyGFF3', '--help'], check=True, 
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print("EMBLmyGFF3 is available in PATH.")
except FileNotFoundError:
    raise EnvironmentError(
        "EMBLmyGFF3 is not found in the system PATH."
        "Please install it and ensure it's accessible."
    )

# check if base directory exists
if not os.path.exists(args.basedir):
    print(f"Error: Base directory does not exist: {args.basedir}")
    sys.exit(1)
print(f"Base directory exists: {args.basedir}")

# get basename of directory which indicates if it's skim2mito or skim2rrna
pipeline_name = os.path.basename(os.path.normpath(args.basedir))
if pipeline_name not in ['skim2mito', 'skim2rrna']:
    print(f"Error: Base directory should be named 'skim2mito' or 'skim2rrna', but got '{pipeline_name}'")
    sys.exit(1)
print(f"Detected pipeline: {pipeline_name}")

# read samples file
try:
    samples = pd.read_csv(args.samples, sep='\t')
    print(f"Loaded {len(samples)} samples from {args.samples}")
except Exception as e:
    print(f"Error reading samples file: {e}")
    sys.exit(1)

# read config file
try:
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)
        print(f"Loaded config from {args.config}")
except Exception as e:
    print(f"Error reading config file: {e}")
    sys.exit(1)


# create output directory if it doesn't exist
if not os.path.exists(args.output):
    os.makedirs(args.output)
    print(f"Created output directory: {args.output}")
else:
    if args.overwrite:
        print(f"Overwriting existing output directory: {args.output}")
        shutil.rmtree(args.output)
        os.makedirs(args.output)
    else:
        print(f"Output directory already exists and overwrite is not enabled: {args.output}")
        sys.exit(1)



def run_emblmygff3(input_fasta, input_gff, topology, moleculetype, trans_code, locus_tag, scientific_name, project_id, description, output_embl, output_log):
    """Run EMBLmyGFF3 to convert GFF3 and FASTA files to EMBL format for ENA submission."""
    cmd = [
        'EMBLmyGFF3', input_gff, input_fasta,
        '--topology', topology,
        '--molecule_type', f'{moleculetype}',
        '--transl_table', str(trans_code),
        '--species', f'"{scientific_name}"',
        '--locus_tag', locus_tag,
        '--project_id', project_id,
        '--de', description,
        '--shame', '-o', output_embl
    ]
    with open(output_log, 'w') as log:
        result = subprocess.run(cmd, stdout=log, stderr=log, text=True)


cmds_to_submit = []

# create assembly summary file
with open(os.path.join(args.output, 'assembly_summary.tsv'), 'w') as summary:
    summary.write("SampleID\tAssemblyType\tFastaCount\tFastaSumLen\tAnnotations\tena_validation\n")

    # loop through samples and process each one
    for _, row in samples.iterrows():

        # get sample ID
        sample_id = row['ID']

        # get scientific name
        scientific_name = row['scientificName']

        # get ena sample ID
        sample_ena_id = row['SampleID']

        print(f"Processing sample {sample_id} with scientific name {scientific_name}...")

        # define fasta and annotation paths
        fasta_path = os.path.join(args.basedir, f"results/assembled_sequence/{sample_id}.fasta")
        annotation_dir = os.path.join(args.basedir, f"results/annotations/{sample_id}/")

        # check if fasta file exists
        if not os.path.isfile(fasta_path):
            print(f"   Warning: FASTA file not found for sample {sample_id} at {fasta_path}. Skipping.")
            continue
        
        # read fasta file as list so it can be used more than once without re-reading from disk
        fasta = list(SeqIO.parse(fasta_path, "fasta"))
        
        # count number of sequences in the fasta file
        fasta_count = sum(1 for _ in fasta)

        # get sequence IDs from fasta file
        sequence_ids = [record.id for record in fasta]

        # sum length of all sequences in the fasta file
        total_length = sum(len(record.seq) for record in fasta)

        # get paths to gff annotation files
        gff_files = []
        for root, dirs, files in os.walk(annotation_dir):
            for file in files:
                if file.endswith('.gff'):
                    gff_files.append(os.path.join(root, file))

        # number of sequence in the fasta file should equal the number of gff files
        if fasta_count != len(gff_files):
            print(f"   Warning: Number of sequences in FASTA ({fasta_count}) does not match number of GFF files ({len(gff_files)}) for sample {sample_id}. Skipping.")
            continue

        # make output subdirectory for sample
        sample_output_dir = os.path.join(args.output, sample_id)
        os.makedirs(sample_output_dir)

        # write temporary gff3 file for EMBLmyGFF3
        temp_gff3_path = os.path.join(sample_output_dir, f"{sample_id}.gff3")
        with open(temp_gff3_path, 'w') as gff3:
            for gff_file in gff_files:
                with open(gff_file, 'r') as gf:
                    for line in gf:
                        if not line.startswith('#'):
                            # get columns from gff3 file
                            seqname, source, feature, start, end, score, strand, frame, attribute = line.rstrip("\n").split("\t")
                            gene_name = attribute.split("Name=")[-1].split(";")[0]
                            if int(start) > int(end):
                                print(f"   Warning: Excluding {gene_name} annotation which bridges circular breakpoint.") ### TBC rotate feature to correct for reverse strand
                                continue
                            gff3.write(line)

        # count the number of annotations in the gff3 file for gene feature for skim2mito and rRNA feature for skim2rrna
        annotation_count = 0
        annotation_names = []
        with open(temp_gff3_path, 'r') as gff3:
            for line in gff3:
                line = line.rstrip("\n").split("\t")
                # get columns from gff3 file
                seqname, source, feature, start, end, score, strand, frame, attribute = line
                # get annotation name from attribute column
                name = attribute.split("Name=")[-1].split(";")[0]
                # count gene feature for skim2mito and rRNA feature for skim2rrna
                if pipeline_name == 'skim2mito' and feature == 'gene':
                    annotation_count += 1
                    annotation_names.append(name)
                elif pipeline_name == 'skim2rrna' and feature == 'rRNA':
                    annotation_count += 1
                    annotation_names.append(name)
        
        # define assembly type
   
        #### how to handle Arc_AOLE_QLD_5 skim2mito 
        #### single contig with 347 length

        if pipeline_name == 'skim2mito':
            molecule_type = 'genomic DNA'
            de_line = 'XXX'

            # mitochondrial complete
            if len(sequence_ids) == 1 and 'circular' in sequence_ids[0]:
                assembly_type = 'mitochondrial complete'
                topology = 'circular'

            # mitochondrial partial
            elif len(sequence_ids) == 1 and 'circular' not in sequence_ids[0]:
                assembly_type = 'mitochondrial partial'
                topology = 'linear'
            
            # mitochondrial with multiple contigs
            elif len(sequence_ids) > 1:
                    assembly_type = 'mitochondrial contig'
                    topology = 'linear'

        elif pipeline_name == 'skim2rrna':
            molecule_type = 'rRNA'
            topology = 'linear'
            
            # ribosomal complete
            if annotation_count == 3 and set(annotation_names) == {'5_8S_rRNA', '18S_rRNA', '28S_rRNA'}:
                assembly_type = 'ribosomal RNA complete'
            # ribosomal partial
            else:
                assembly_type = 'ribosomal RNA partial'
            
            # format names for description line
            annotation_names_description = list(set(annotation_names))
            # remove "_rRNA" and replace with " rRNA gene"
            annotation_names_description = [name.replace('_rRNA', ' rRNA gene') for name in annotation_names_description]
            # remove "5_8S" and replace with "5.8S"
            annotation_names_description = [name.replace('5_8S', '5.8S') for name in annotation_names_description]
            # repl
            de_line = f"{scientific_name} {' '.join(annotation_names_description)}, {sample_id} Isolate"

        # run EMBLmyGFF3 to create EMBL file for ENA submission
        run_emblmygff3(input_fasta = fasta_path, 
            input_gff = temp_gff3_path, 
            topology = topology, 
            moleculetype = molecule_type,
            trans_code = 4, 
            locus_tag = args.locus_tag, 
            scientific_name = scientific_name, 
            project_id = args.project, 
            description = de_line,
            output_embl = os.path.join(sample_output_dir, f"{sample_id}.embl"),
            output_log= os.path.join(sample_output_dir, f"{sample_id}_EMBLmyGFF3.log"))
        
        # remove the line starting with DT for ribsomal flatfiles from the EMBL file to avoid ena-webin-cli validation error
        if assembly_type == 'ribosomal RNA complete' or assembly_type == 'ribosomal RNA partial':
            with open(os.path.join(sample_output_dir, f"{sample_id}.embl"), 'r') as embl_file:
                lines = embl_file.readlines()
            with open(os.path.join(sample_output_dir, f"{sample_id}.embl"), 'w') as embl_file:
                for line in lines:
                    if not line.startswith('DT'):
                        embl_file.write(line)   

        # gzip EMBL file
        if os.path.exists(os.path.join(sample_output_dir, f"{sample_id}.embl.gz")):
            os.remove(os.path.join(sample_output_dir, f"{sample_id}.embl.gz"))
        subprocess.run(['gzip', os.path.join(sample_output_dir, f"{sample_id}.embl")])

        # write gzip copy fasta to output directory
        #SeqIO.write(fasta, os.path.join(sample_output_dir, f"{sample_id}.fasta"), "fasta")
        #subprocess.run(['gzip', os.path.join(sample_output_dir, f"{sample_id}.fasta")])
        
        if assembly_type == 'mitochondrial complete' or assembly_type == 'mitochondrial partial':

            # get sequence coverage
            coverage_values = []
            blobtools_table = os.path.join(args.basedir, f'results/blobtools/{sample_id}/table.tsv')
            if os.path.isfile(blobtools_table):
                with open(blobtools_table, 'r') as bt:
                    for line in bt:
                        if not line.startswith('index'): 
                            line = line.rstrip("\n").split("\t")
                            coverage_values.append(float(line[4]))              
                average_coverage = sum(coverage_values) / len(coverage_values)
            else:
                print(f"   Warning: Blobtools file not found for sample {sample_id} at {blobtools_table}. Skipping coverage extraction.")
                sys.exit(1)

            # write chromosome list
            with gzip.open(f'{sample_output_dir}/{sample_id}_chromosome_list.txt.gz', 'wt') as chrom_file:
                if topology == 'circular':
                    chrom_file.write(f"{sample_id}_circular\tMIT\tCircular-Chromosome\tMitochondrion\n")
                else:
                    for n, chr in enumerate(sequence_ids):
                        chrom_file.write(f"{chr}\tMIT{n}\tLinear-Chromosome\tMitochondrion\n")

        # write manifest file
        with open(f'{sample_output_dir}/{sample_id}_manifest.txt', 'w') as manifest:
            manifest.write(f'STUDY\t{args.project}\n')
            manifest.write(f'DESCRIPTION\t{assembly_type}\n')
            if assembly_type == 'mitochondrial complete' or assembly_type == 'mitochondrial partial' or assembly_type == 'mitochondrial contig':
                manifest.write(f'SAMPLE\t{sample_ena_id}\n')
                manifest.write(f'ASSEMBLYNAME\t{sample_id}\n')
                manifest.write(f'COVERAGE\t{average_coverage}\n')
                manifest.write(f'PROGRAM\tGetOrganelle\n')
                manifest.write(f'PLATFORM\t{args.platform}\n')
                manifest.write(f'ASSEMBLYTYPE\t{args.assembly_type}\n')
            if assembly_type == 'ribosomal RNA complete' or assembly_type == 'ribosomal RNA partial':
                 manifest.write(f'NAME\t{sample_ena_id}\n')
            manifest.write(f'FLATFILE\t{sample_id}.embl.gz\n')
            
            # manifest.write(f'FASTA\t{sample_id}.fasta.gz\n')
            if assembly_type == 'mitochondrial complete' or assembly_type == 'mitochondrial partial':
                manifest.write(f'CHROMOSOME_LIST\t{sample_id}_chromosome_list.txt.gz\n')

        # validate ena submission files using ena-validate
        if assembly_type == 'mitochondrial complete' or assembly_type == 'mitochondrial partial' or assembly_type == 'mitochondrial contig':
            cmd_ena = ['ena-webin-cli', '-context', 'genome', '-validate',
                      '-manifest', f'{sample_output_dir}/{sample_id}_manifest.txt', 
                      '-userName', f'{args.ena_username}', '-password', f'{args.ena_password}',
                      '-inputDir', f'{sample_output_dir}']
        if assembly_type == 'ribosomal RNA complete' or assembly_type == 'ribosomal RNA partial':
            cmd_ena = ['ena-webin-cli', '-context', 'sequence', '-validate',
                      '-manifest', f'{sample_output_dir}/{sample_id}_manifest.txt', 
                      '-userName', f'{args.ena_username}', '-password', f'{args.ena_password}',
                      '-inputDir', f'{sample_output_dir}']
        log_file_ena_webin = f'{sample_output_dir}/{sample_id}_ena_webin.log'
        with open(log_file_ena_webin, 'w') as log:
            result = subprocess.run(cmd_ena, stdout=log, stderr=log, text=True)
        
        # define cmd to upload to ENA using ena-webin-cli
        if assembly_type == 'mitochondrial complete' or assembly_type == 'mitochondrial partial' or assembly_type == 'mitochondrial contig':
            cmd_ena_submit = ['ena-webin-cli', '-context', 'genome', '-submit',
                             '-manifest', f'{sample_output_dir}/{sample_id}_manifest.txt', 
                             '-userName', f'{args.ena_username}', '-password', f'{args.ena_password}',
                             '-inputDir', f'{sample_output_dir}']
        if assembly_type == 'ribosomal RNA complete' or assembly_type == 'ribosomal RNA partial':
            cmd_ena_submit = ['ena-webin-cli', '-context', 'sequence', '-submit',
                             '-manifest', f'{sample_output_dir}/{sample_id}_manifest.txt', 
                             '-userName', f'{args.ena_username}', '-password', f'{args.ena_password}',
                             '-inputDir', f'{sample_output_dir}']

        # check if en_webin log contains "INFO : Submission(s) validated successfully."
        with open(log_file_ena_webin, 'r') as log:
            log_content = log.read()
            if "INFO : Submission(s) validated successfully." in log_content:
                pass_status = True
                cmds_to_submit.append(cmd_ena_submit)
            else:
                pass_status = False

        summary.write(f"{sample_id}\t{assembly_type}\t{fasta_count}\t{total_length}\t{annotation_count}\t{pass_status}\n")
    
for c in cmds_to_submit:
    print(' '.join(c)) 
