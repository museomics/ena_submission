#!/usr/bin/env python3
"""
Create ENA submission files for mitochondrial genome assemblies.

This script processes assembled mitochondrial genomes and generates EMBL files,
manifests, and other required files for submission to ENA via webin-cli.
"""

import argparse
import yaml
import pandas as pd
import os
import subprocess
import sys


def check_emblmygff3():
    """Check if EMBLmyGFF3 is available in the system PATH."""
    try:
        subprocess.run(['EMBLmyGFF3', '--help'], check=True, 
                      stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError:
        raise EnvironmentError(
            "EMBLmyGFF3 is not found in the system PATH. "
            "Please install it and ensure it's accessible."
        )


def setup_output_directory(output_path, overwrite):
    """Create or validate the output directory."""
    if os.path.exists(output_path) and not overwrite:
        raise FileExistsError(
            f"The output directory '{output_path}' already exists. "
            "Use --overwrite to overwrite or choose a different directory."
        )
    elif os.path.exists(output_path) and overwrite:
        print(f"Warning: The output directory '{output_path}' already exists and will be overwritten.")
    elif not os.path.exists(output_path):
        os.makedirs(output_path)


def process_sample(sample_id, scientific_name, sample_ena_id, config, 
                   output_path, project, platform, moleculetype):
    """Process a single sample and generate all required files."""
    print(f"Processing sample: {sample_id}")

    # Define expected paths
    fasta_path = f"results/assembled_sequence/{sample_id}.fasta"
    annotation_dir = f"results/annotations/{sample_id}/"

    # Check if assembly exists
    if not os.path.exists(fasta_path):
        print(f"   No assembly found\n")
        return
    else:
        print(f"   Assembly found")

    # Get the number of sequences in the fasta file
    with open(fasta_path, 'r') as fasta:
        num_sequences = sum(1 for line in fasta if line.startswith('>'))
        print(f"   Number of sequences: {num_sequences}")
        
        # If there is a single contig, check if it's circular
        topology = 'linear'
        if num_sequences == 1:
            fasta.seek(0)
            sequence_name = fasta.readline().strip()
            if 'circular' in sequence_name:
                topology = 'circular'
            print(f"   Topology: {topology}")

    # Search recursively for the gff annotation files
    list_gff_files = []
    for root, dirs, files in os.walk(annotation_dir):
        for file in files:
            if file.endswith('.gff'):
                annotation_path = os.path.join(root, file)
                list_gff_files.append(annotation_path)
    
    if len(list_gff_files) == 0:
        print(f"   No annotation found\n")
        return
    elif len(list_gff_files) != num_sequences:
        print(f"   Warning: Number of annotation files ({len(list_gff_files)}) "
              f"does not match number of sequences ({num_sequences})\n")
    else:
        print(f"   Found {len(list_gff_files)} annotation file(s)")

    # If multiple gff files found, write to a temporary file
    if len(list_gff_files) > 1:
        temp_gff_path = f"{output_path}/{sample_id}.gff"
        with open(temp_gff_path, 'w') as temp_gff:
            temp_gff.write("##gff-version 3\n")
            for gff_file in list_gff_files:
                with open(gff_file, 'r') as gff:
                    for line in gff:
                        if not line.startswith('#'):
                            temp_gff.write(line)
        gff_path = temp_gff_path
    else:
        gff_path = list_gff_files[0]

    # Run EMBLmyGFF3
    print(f"   Running EMBLmyGFF3 for sample {sample_id}")
    cmd = [
        'EMBLmyGFF3', gff_path, fasta_path,
        '--topology', topology,
        '--molecule_type', moleculetype,
        '--transl_table', str(config['mitos_code']),
        '--species', scientific_name,
        '--locus_tag', 'LOCUSTAG',
        '--project_id', project,
        '-o', f'{output_path}/{sample_id}.embl'
    ]

    log_file = f'{output_path}/{sample_id}.log'
    with open(log_file, 'w') as log:
        result = subprocess.run(cmd, stdout=log, stderr=log, text=True)
    print(f"   Log written to: {log_file}")

    # Remove temporary gff file if created
    if len(list_gff_files) > 1:
        os.remove(temp_gff_path)

    # Copy fasta to output directory
    subprocess.run(['cp', fasta_path, f'{output_path}/{sample_id}.fasta'])

    # Handle circular/chromosomal assemblies
    if topology == 'circular':
        print(f"   Treating upload as chromosome")
        
        print('   Writing a chromosome list')
        with open(f'{output_path}/{sample_id}_chromosome_list.txt', 'w') as chrom_file:
            chrom_file.write(f"{sample_id}\tMIT\tchromosome\tcircular\tMitochondrion\n")

        # Get sequence length
        with open(f'results/seqkit/{sample_id}.txt', 'r') as seqkit_file:
            seqkit_file.readline()  # skip header
            seq_len = int(seqkit_file.readline().strip('\n').split()[4].replace(',', ''))

        print('   Writing an AGP file')
        with open(f'{output_path}/{sample_id}_apg_file.txt', 'w') as apg_file:
            apg_file.write("##agp-version\t2.0\n")
            apg_file.write(f"MIT\t1\t{seq_len}\t1\tW\t{sample_id}\t1\t{seq_len}\t+\n")

    # Get sequence coverage
    with open(f'results/blobtools/{sample_id}/table.tsv', 'r') as blobtools_file:
        blobtools_file.readline()  # skip header
        seq_coverage = blobtools_file.readline().strip('\n').split()[4]

    # Write manifest file
    print('   Writing a manifest file')
    with open(f'{output_path}/{sample_id}_manifest.txt', 'w') as manifest:
        manifest.write(f'STUDY\t{project}\n')
        manifest.write(f'SAMPLE\t{sample_ena_id}\n')
        manifest.write(f'ASSEMBLYNAME\t{sample_id}\n')
        manifest.write(f'COVERAGE\t{seq_coverage}\n')
        manifest.write(f'PROGRAM\tGetOrganelle\n')
        manifest.write(f'PLATFORM\t{platform}\n')
        manifest.write(f'MINGAPLENGTH\t1\n')
        manifest.write(f'MOLECULETYPE\t{moleculetype}\n')
        manifest.write(f'FASTA\t{sample_id}.fasta\n')
    
    print()


def main():
    parser = argparse.ArgumentParser(
        description='Create ENA submission files for mitochondrial genome assemblies.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  %(prog)s -s ena_sample_list.txt -c config/config.yaml -o test_embl
  %(prog)s -s samples.txt -c config.yaml -o output --project PRJEB12345 --overwrite
        """
    )
    
    # Required arguments
    parser.add_argument(
        '-s', '--samples',
        required=True,
        help='Path to TSV file with sample information (ID, scientificName, SampleID)'
    )
    parser.add_argument(
        '-c', '--config',
        required=True,
        help='Path to snakemake config YAML file'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output directory for generated files'
    )
    
    # Optional arguments
    parser.add_argument(
        '--project',
        default='PRJXXX',
        help='ENA project accession (default: PRJXXX)'
    )
    parser.add_argument(
        '--platform',
        default='Illumina',
        help='Sequencing platform (default: Illumina)'
    )
    parser.add_argument(
        '--moleculetype',
        default='genomic DNA',
        help='Molecule type (default: genomic DNA)'
    )
    parser.add_argument(
        '--overwrite',
        action='store_true',
        help='Overwrite existing output directory'
    )
    
    args = parser.parse_args()
    
    # Check for EMBLmyGFF3
    check_emblmygff3()
    
    # Setup output directory
    setup_output_directory(args.output, args.overwrite)
    
    # Read samples file
    try:
        samples = pd.read_csv(args.samples, sep='\t')
        print(f"Loaded {len(samples)} samples from {args.samples}")
    except Exception as e:
        print(f"Error reading samples file: {e}")
        sys.exit(1)
    
    # Read config file
    try:
        with open(args.config, 'r') as f:
            config = yaml.safe_load(f)
    except Exception as e:
        print(f"Error reading config file: {e}")
        sys.exit(1)
    
    # Process each sample
    for _, row in samples.iterrows():
        sample_id = row['ID']
        scientific_name = row['scientificName']
        sample_ena_id = row['SampleID']
        
        process_sample(
            sample_id=sample_id,
            scientific_name=scientific_name,
            sample_ena_id=sample_ena_id,
            config=config,
            output_path=args.output,
            project=args.project,
            platform=args.platform,
            moleculetype=args.moleculetype
        )
    
    print("Processing complete!")


if __name__ == '__main__':
    main()
