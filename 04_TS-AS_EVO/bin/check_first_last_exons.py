import sys

def parse_gtf(gtf_file):
    """Parse GTF file and extract transcript information."""
    transcripts = {}
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if fields[2] == 'CDS':  # Consider only coding exons
                attributes = {}
                for item in fields[8].strip().split(';'):
                    if item.strip():
                        parts = item.strip().split(' ', 1)
                        if len(parts) == 2:
                            key, value = parts
                            attributes[key.strip()] = value.strip('" ')
                transcript_id = attributes.get('transcript_id')
                chromosome = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]
                if transcript_id:
                    if transcript_id not in transcripts:
                        transcripts[transcript_id] = {'chromosome': chromosome, 'exons': []}
                    transcripts[transcript_id]['exons'].append({'start': start, 'end': end, 'strand': strand})
    return transcripts

def label_coordinates(coord_str, transcripts):
    """Label coordinates as 'first', 'middle', 'last', or 'not_found'."""
    chromosome, coords = coord_str.split(':')
    start, end = map(int, coords.split('-'))

    labels = []
    for transcript in transcripts.values():
        first_exon = transcript['exons'][0]
        last_exon = transcript['exons'][-1]
        if first_exon['start'] == start and first_exon['end'] == end:
            labels.append('first')
        elif last_exon['start'] == start and last_exon['end'] == end:
            labels.append('last')
    
    if not labels:  # If no matching exon is found
        return 'not_found'

    if 'last' in labels:
        return 'last'
    elif 'first' in labels:
        return 'first'
    else:
        return 'middle'

def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py <GTF_file> <coordinates_file> <output_file>")
        sys.exit(1)

    gtf_file = sys.argv[1]
    coordinates_file = sys.argv[2]
    output_file = sys.argv[3]
    transcripts = parse_gtf(gtf_file)

    with open(coordinates_file, 'r') as f, open(output_file, 'w') as out_f:
        for line in f:
            coordinate = line.strip()
            label = label_coordinates(coordinate, transcripts)
            out_f.write(f"The coordinates {coordinate} are labeled as '{label}'.\n")

if __name__ == "__main__":
    main()
