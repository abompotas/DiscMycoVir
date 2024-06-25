from Bio.Blast import NCBIXML


def get_text_output(filepath, array=False):
    with open(filepath, 'r') as txt:
        if array:
            contents = txt.readlines()
        else:
            contents = txt.read()
    return contents


def parse_blast_xml(blast_file):
    records = []
    with open(blast_file, 'r') as xml:
        blast_records = NCBIXML.parse(xml)
        for r in blast_records:
            records.append(get_blast_record(r))
    return records


def get_blast_record(blast_record):
    return {
        'queryName': blast_record.query,
        'queryLetters': blast_record.query_letters,
        'descriptions': get_blast_descriptions(blast_record),
        'alignments': get_blast_alignments(blast_record)
    }


def get_blast_descriptions(blast_record):
    descriptions = []
    for d in blast_record.descriptions:
        descriptions.append({
            'title': d.title,
            'score': d.score,
            'bits': d.bits,
            'e': d.e,
            'numAlignments': d.num_alignments
        })
    return descriptions


def get_blast_alignments(blast_record):
    alignments = []
    for a in blast_record.alignments:
        alignments.append({
            'hitId': a.hit_id,
            'hitDef': a.hit_def,
            'length': a.length,
            'hsps': get_blast_hsps(a)
        })
    return alignments


def get_blast_hsps(blast_alignment):
    hsps = []
    for h in blast_alignment.hsps:
        hsps.append({
            'score': h.score,
            'bits': h.bits,
            'expect': h.expect,
            'numAlignments': h.num_alignments,
            'identities': h.identities,
            'positives': h.positives,
            'gaps': h.gaps,
            'alignLength': h.align_length,
            'strand': h.strand,
            'frame': h.frame,
            'query': h.query,
            'queryStart': h.query_start,
            'queryEnd': h.query_end,
            'match': h.match,
            'sbjct': h.sbjct,
            'sbjctStart': h.sbjct_start,
            'sbjctEnd': h.sbjct_end,
        })
    return hsps
