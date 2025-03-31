#!/usr/bin/env python3
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
from datetime import datetime

if sys.stdout.isatty():
    RED = "\033[1;31m"
    GREEN = "\033[1;32m"
    YELLOW = "\033[1;33m"
    NC = "\033[0m"
else:
    RED = GREEN = YELLOW = NC = ""

def log_info(message):
    print(f"{GREEN}[{datetime.now():%Y-%m-%d %H:%M:%S}] [INFO]{NC} {message}")

def log_warn(message):
    print(f"{YELLOW}[{datetime.now():%Y-%m-%d %H:%M:%S}] [WARN]{NC} {message}")

def log_error(message):
    print(f"{RED}[{datetime.now():%Y-%m-%d %H:%M:%S}] [ERROR]{NC} {message}", file=sys.stderr)

def process_genbank(input_file: str, output_file: str) -> bool:
    """
    Process GenBank file to extract mat_peptide features and generate FASTA.
    
    Args:
        input_file: Path to input GenBank file
        output_file: Path for output FASTA file
        
    Returns:
        bool: True if successful, False if error occurs
    """
    try:
        # Read GenBank file
        record = SeqIO.read(input_file, "genbank")
        features = record.features
        accession = record.id

        # Extract mat_peptide features
        mat_peptides = [f for f in features if f.type == "mat_peptide"]
        
        if not mat_peptides:
            log_info(f"No mat_peptide features found in {input_file}")
            return True  # Valid outcome, not an error, simply go back to main bash function and retrieve cds
        cds_features = [f for f in features if f.type == "CDS"]
        fallback_gene, fallback_protein_id, fallback_locus_tag = "", "", ""

        # use CDS as fallback source if not features are found in any mat peptide
        if len(cds_features) == 1:
            # Use the single CDS as fallback for all mat_peptides
            fallback_cds = cds_features[0]
            fallback_gene = fallback_cds.qualifiers.get("gene", [""])[0]
            fallback_locus_tag = fallback_cds.qualifiers.get("locus_tag", [""])[0]
            fallback_protein_id = fallback_cds.qualifiers.get("protein_id", [""])[0]
            log_info(f'Using features from "{fallback_gene if fallback_gene else fallback_locus_tag}" CDS as fallback in case no corresponding features are found in any mat_peptide')

        # Process each mat_peptide
        subproteins = []
        for idx, feature in enumerate(mat_peptides, start=1):
            try:
                # Extract coordinates (0-based to 1-based conversion)
                start = int(feature.location.start)
                end = int(feature.location.end)
                location = f"{start + 1}..{end}"

                # Get metadata with fallback to polyprotein features
                gene = feature.qualifiers.get("gene", [fallback_gene])[0]
                protein_id = feature.qualifiers.get("protein_id", [fallback_protein_id])[0]
                product = feature.qualifiers.get("product", [None])[0]
                locus_tag = feature.qualifiers.get("locus_tag", [fallback_locus_tag])[0]

                #If werent found then get them for the corresponding cds that contains the mat_peptide
                if not (gene and protein_id and locus_tag):
                    for cds in cds_features:
                        cds_start = int(cds.location.start)
                        cds_end = int(cds.location.end)
                        if (start >= min(cds_start, cds_end)) and (end <= max(cds_start, cds_end)):
                            if not gene:
                                gene = cds.qualifiers.get("gene", [""])[0]
                            if not locus_tag:
                                locus_tag = cds.qualifiers.get("locus_tag", [""])[0]
                            if not protein_id:
                                protein_id = cds.qualifiers.get("protein_id", [""])[0]
                            break
                
                if not protein_id:
                    raise ValueError(f"Missing gene/protein_id at {location}")
                
                # Build FASTA header and description
                header = f"lcl|{accession}_cds_{protein_id}_{idx}"
                description_parts = [
                    f"[protein_id={protein_id}]",
                    f"[location={location}]"
                ]
                #Be sure that gene is first, or locus tag if there is not gene, or prodcut goes first if there is neither gene nor locustag
                if gene:
                    description_parts.insert(0, f"[gene={gene}]")
                if locus_tag:
                    if gene:
                        description_parts.append(f"[locus_tag={locus_tag}]")
                    else:
                        description_parts.insert(0, f"[locus_tag={locus_tag}]")
                if product:
                    if gene or locus_tag:
                        description_parts.insert(1, f"[protein={product}]")
                    else:    
                        description_parts.insert(0, f"[protein={product}]")

                # Create SeqRecord
                subproteins.append(SeqRecord(
                    seq=record.seq[start:end],
                    id=header,
                    description=" ".join(description_parts)
                ))
            except Exception as e:
                sys.stderr.write(f"[ERROR] Failed to process mat_peptide {idx}: {str(e)}\n")
                return False
        # Write output
        with open(output_file, "w") as f:
            SeqIO.write(subproteins, f, "fasta")
        log_info(f"Generated {len(subproteins)} subproteins in {output_file}")
        return True

    except Exception as e:
        log_error(f"{str(e)}", file=sys.stderr)
        return False

if __name__ == "__main__":
    if len(sys.argv) != 3:
        log_info("Usage: python genbank_to_fasta.py <input.gbk> <output.fasta>")
        sys.exit(1)
        
    input_path = sys.argv[1]
    output_path = sys.argv[2]
    
    if not process_genbank(input_path, output_path):
        sys.exit(1)
    sys.exit(0)
