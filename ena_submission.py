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


def split_circular_annotations(gff_path, fasta_path, output_path, sample_id):
    """
    Split annotations that span the junction in circular genomes.
    
    For circular genomes, some annotations may wrap around the origin.
    This function splits such features into two parts (part_a and part_b).
    
    Note: This assumes a single circular sequence, which should be validated
    before calling this function.
    """
    # Get sequence length from fasta (expect single sequence for circular genomes)
    seq_length = 0
    num_seqs = 0
    with open(fasta_path, 'r') as f:
        seq = ''
        for line in f:
            if line.startswith('>'):
                num_seqs += 1
                if num_seqs > 1:
                    raise ValueError(f"Circular genome should contain only one sequence, found {num_seqs}")
            else:
                seq += line.strip()
        seq_length = len(seq)
    
    if num_seqs != 1:
        raise ValueError(f"Expected 1 sequence for circular genome, found {num_seqs}")
    
    output_gff = f'{output_path}/{sample_id}_split.gff'
    
    with open(gff_path, 'r') as infile, open(output_gff, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                outfile.write(line)
                continue
            
            seqid, source, feature_type, start, end, score, strand, phase, attributes = fields
            start, end = int(start), int(end)
            
            # Check if annotation spans the junction (wraps around)
            if end < start:
                # print message about splitting
                print(f"   Splitting annotation with attribute {attributes} spanning junction: {start}-{end} on circular genome")
                # Split into two parts
                # Part A: from start to seq_length
                part_a = fields.copy()
                part_a[4] = str(seq_length)  # end at sequence end
                # Modify attributes to add -part_a suffix
                part_a[8] = attributes.replace('ID=', 'ID=').replace(';', '-part_a;', 1)
                if 'ID=' not in part_a[8]:
                    part_a[8] += ';Note=part_a'
                else:
                    part_a[8] = part_a[8].replace('ID=', 'ID=', 1).replace(';', '-part_a;', 1)
                outfile.write('\t'.join(part_a) + '\n')
                
                # Part B: from 1 to end
                part_b = fields.copy()
                part_b[3] = '1'  # start at 1
                # Modify attributes to add -part_b suffix
                part_b[8] = attributes.replace('ID=', 'ID=').replace(';', '-part_b;', 1)
                if 'ID=' not in part_b[8]:
                    part_b[8] += ';Note=part_b'
                else:
                    part_b[8] = part_b[8].replace('ID=', 'ID=', 1).replace(';', '-part_b;', 1)
                outfile.write('\t'.join(part_b) + '\n')
            else:
                # Normal annotation, write as is
                outfile.write(line)
    
    return output_gff


def process_sample(sample_id, scientific_name, sample_ena_id, config, 
                   output_path, project, platform, moleculetype, assemblytype, 
                   basedir, username, password):
    """Process a single sample and generate all required files."""
    print(f"Processing sample: {sample_id}")

    # Define expected paths
    fasta_path = os.path.join(basedir, f"results/assembled_sequence/{sample_id}.fasta")
    annotation_dir = os.path.join(basedir, f"results/annotations/{sample_id}/")

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

    # split annotations which span the break point for circular genomes
    if topology == 'circular':
        print(f"   Splitting annotations spanning the junction for circular genome")
        gff_path = split_circular_annotations(gff_path, fasta_path, output_path, sample_id)

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
        '--shame', '-o', f'{output_path}/{sample_id}.embl'
    ]

    log_file_emblmygff3 = f'{output_path}/{sample_id}_EMBLmyGFF3.log'
    with open(log_file_emblmygff3, 'w') as log:
        result = subprocess.run(cmd, stdout=log, stderr=log, text=True)
    print(f"   EMBLmyGFF3 log written to: {log_file_emblmygff3} and printed below:")
    print("   ----- EMBLmyGFF3 LOG START -----")
    with open(log_file_emblmygff3, 'r') as log:
        for line in log:
            print(f"      {line.strip()}")
    print("   ----- EMBLmyGFF3 LOG END -----")
    
    # Remove temporary gff file if created
    if len(list_gff_files) > 1:
        os.remove(temp_gff_path)

    # Copy fasta to output directory
    subprocess.run(['cp', fasta_path, f'{output_path}/{sample_id}.fasta'])

    # Gzip fasta 
    subprocess.run(['gzip', '-f', f'{output_path}/{sample_id}.fasta'])


    # Handle circular/chromosomal assemblies
    if topology == 'circular':
        print(f"   Treating upload as chromosome")
        
        print('   Writing a chromosome list')
        with open(f'{output_path}/{sample_id}_chromosome_list.txt', 'w') as chrom_file:
            chrom_file.write(f"{sample_id}\tMIT\tchromosome\tcircular\tMitochondrion\n")

        # Get sequence length
        seqkit_path = os.path.join(basedir, f'results/seqkit/{sample_id}.txt')
        with open(seqkit_path, 'r') as seqkit_file:
            seqkit_file.readline()  # skip header
            seq_len = int(seqkit_file.readline().strip('\n').split()[4].replace(',', ''))

        print('   Writing an AGP file')
        with open(f'{output_path}/{sample_id}_apg_file.txt', 'w') as apg_file:
            apg_file.write("##agp-version\t2.0\n")
            apg_file.write(f"MIT\t1\t{seq_len}\t1\tW\t{sample_id}\t1\t{seq_len}\t+\n")

    # Get sequence coverage
    blobtools_path = os.path.join(basedir, f'results/blobtools/{sample_id}/table.tsv')
    with open(blobtools_path, 'r') as blobtools_file:
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
        manifest.write(f'ASSEMBLYTYPE\t{assemblytype}\n')
        manifest.write(f'FASTA\t{sample_id}.fasta.gz\n')

    # Run ena-webin-cli to validate assembly
    print('   Validating the manifest file with ena-webin-cli')

    cmd2 = ['ena-webin-cli', '-context', 'genome', '-validate',
            '-manifest', f'{output_path}/{sample_id}_manifest.txt', 
            '-userName', f'{username}', '-password', f'{password}',
            '-inputDir', f'{output_path}']

    log_file_ena_webin = f'{output_path}/{sample_id}_ena_webin.log'
    with open(log_file_ena_webin, 'w') as log:
        result = subprocess.run(cmd2, stdout=log, stderr=log, text=True)
    print(f"   ena-webin-cli log written to: {log_file_ena_webin} and printed below:")
    print("   ----- ena-webin-cli LOG START -----")
    with open(log_file_ena_webin, 'r') as log:
        for line in log:
            print(f"      {line.strip()}")
    print("   ----- ena-webin-cli LOG END -----")   
    print()


def main():
    parser = argparse.ArgumentParser(
        description='Create ENA submission files for mitochondrial genome assemblies.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  %(prog)s -s ena_sample_list.txt -c config/config.yaml -o test_embl
  %(prog)s -s samples.txt -c config.yaml -o output --project PRJEB12345 --overwrite
  %(prog)s -s samples.txt -c config.yaml -o output -b /path/to/data
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
        '-b', '--basedir',
        default='.',
        help='Base directory for input files (default: current directory)'
    )
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
        '--assemblytype',
        default='isolate',
        help='Assembly type (default: isolate)'
    )
    parser.add_argument(
        '--ena_username',
        help='ENA username for webin-cli'
    )
    parser.add_argument(
        '--ena_password',
        help='ENA password for webin-cli'
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
            moleculetype=args.moleculetype,
            assemblytype=args.assemblytype,
            basedir=args.basedir,
            username=args.ena_username,
            password=args.ena_password
        )
    
    print("Processing complete!")

if __name__ == '__main__':
    main()
