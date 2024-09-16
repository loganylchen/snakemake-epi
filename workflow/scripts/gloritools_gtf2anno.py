import sys
import re
import os


ATTRIBUTE_PATTERN=re.compile('([a-z_]+) "(.+?)"[;]?')

def _parse_attribute(attribute: str):
    attribute_dict = dict()
    for att,value in ATTRIBUTE_PATTERN.findall(attributs):
        attribute_dict[att]=value
    return attribute_dict


def _parse_gtf_line(line: str):
    seqname, source, feature, start, end, score, strand, frame, attribute_str = line.strip('\n').split('\t')
    if feature not in ['transcript', 'exon','CDS']:
        return None
    attributes = _parse_attribute(attribute_str)
    return {'seqname': seqname,
            'feature': feature,
            'start': int(start),
            'end': int(end),
            'strand': strand,
            'attributes': attributes}


def _sort_and_parse_exons(tid_dict):
    '''
    tid_dict = {
        'exons': [(start1,end1),(start2,end2),(start3,end3),...,(startn,endn)]
    }
    '''
    if 'exons' in tid_dict:
        tid_dict['exons'] = sorted(tid_dict['exons'],key=lambda x: x[0])
        
        related_start = 0
        transcript_length = 0
        tid_dict['related_exons'] = []
        
        
        for exon in tid_dict['exons']:
            exon_length = exon[1] - exon[0] + 1
            related_end = related_start + exon_length -1 
            tid_dict['related_exons'].append((related_start,related_end))
            related_start = related_end + 1
            transcript_length += exon_length    
        tid_dict['transcript_length'] = transcript_length 
        if tid_dict['strand'] == '-':
            tid_dict['related_exons'] =[(tid_dict['transcript_length']-exon[1]-1,tid_dict['transcript_length']-exon[0]-1,) for exon in tid_dict['related_exons'][::-1]]
    else:
        raise KeyError('exons was not in the dict: {tid_dict[attributes][transcript_id]}')
    if 'CDS' in tid_dict and len(tid_dict['CDS'])>0:
        tid_dict['CDS'] = sorted(tid_dict['CDS'],key=lambda x: x[0])
        related_start = 0
        first_cds = tid_dict['CDS'][0]
        exon_length = 0
        
        for exon in tid_dict['exons']:
            if exon[0]<=first_cds[0] and exon[1]>=first_cds[1]:
                related_start = first_cds[0]-exon[0]  + exon_length
                break
            else:
                exon_length += exon[1] - exon[0] + 1
                
        coding_length = 0
        tid_dict['related_CDS'] = []
        for cds in tid_dict['CDS']:
            cds_length = cds[1] - cds[0] + 1
            related_end = related_start + cds_length -1 
            tid_dict['related_CDS'].append((related_start,related_end))
            related_start = related_end + 1
            coding_length += cds_length    
        tid_dict['coding_length'] = coding_length 
        if tid_dict['strand'] == '-':
            tid_dict['related_CDS'] =[(tid_dict['transcript_length']-exon[1]-1,tid_dict['transcript_length']-exon[0]-1,) for exon in tid_dict['related_CDS'][::-1]]
    
def readGTF(gtf_file):
    '''
    GTF file is 1-based, and the end position is included.
    '''
    gtf_dict = dict()
    transcript_id = None
    
    with (gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                line_dict = _parse_gtf_line(line)
                if None is line_dict:
                    continue
                elif line_dict['feature'] == 'transcript':
                    if transcript_id is not None:
                        _sort_and_parse_exons(gtf_dict[transcript_id])
                    transcript_id = line_dict['attributes']['transcript_id']
                    gtf_dict[transcript_id] = line_dict
                    gtf_dict[transcript_id]['exons'] = []
                    gtf_dict[transcript_id]['CDS'] = []
                    gtf_dict[transcript_id]['coding_length']=0
                    # gtf_dict[transcript_id]['relate_exons'] = []
                    gtf_dict[transcript_id]['transcript_length'] = 0
                elif line_dict['feature'] == 'exon':
                    if line_dict['attributes']['transcript_id'] == transcript_id:
                        start = line_dict['start']
                        end = line_dict['end']
                        gtf_dict[transcript_id]['exons'].append((start, end))
                    else:
                        raise ValueError('Exon is not in the same transcript')
                elif line_dict['feature'] == 'CDS':
                    if line_dict['attributes']['transcript_id'] == transcript_id:
                        start = line_dict['start']
                        end = line_dict['end']
                        gtf_dict[transcript_id]['CDS'].append((start, end))
                    else:
                        raise ValueError('CDS is not in the same transcript')
                else:
                    continue
  
    if transcript_id is not None:
        _sort_and_parse_exons(gtf_dict[transcript_id])
    return gtf_dict


with open(snakemake.log[0], "w") as log_f:
    sys.stderr = sys.stdout = log_f
    print(f'reading {snakemake.input.gtf}')
    gtf_dict = readGTF(snakemake.input.gtf)
    print(f'Done loading')
    result_list = []
    for transcript in gtf_dict:
        gene_name = gtf_dict[transcript]['attributes'].get('gene_name','.')
    transcript_id = transcript
    chromosome = gtf_dict[transcript]['seqname']
    transcript_strand = gtf_dict[transcript]['strand']
    transcript_start = str(gtf_dict[transcript]['start'])
    transcript_end = str(gtf_dict[transcript]['end'])
    cds_start = ','.join([str(i[0]) for i in gtf_dict[transcript]['CDS']]) if len(gtf_dict[transcript]['CDS'])>0 else '.'
    cds_end = ','.join([str(i[1]) for i in gtf_dict[transcript]['CDS']]) if len(gtf_dict[transcript]['CDS'])>0 else '.'
    exon_count = str(len(gtf_dict[transcript]['exons']))
    exon_start = ','.join([str(i[0]) for i in gtf_dict[transcript]['exons']]) if len(gtf_dict[transcript]['exons'])>0 else '.'
    exon_end = ','.join([str(i[1]) for i in gtf_dict[transcript]['exons']]) if len(gtf_dict[transcript]['exons'])>0 else '.'
    transcript_length = str(gtf_dict[transcript]['transcript_length'])
    gene_id = gtf_dict[transcript]['attributes']['gene_id']
    coding_length = str(gtf_dict[transcript]['coding_length'])
    transcript_biotype = gtf_dict[transcript]['attributes'].get('transcript_biotype','.')
    transcript_version = str(gtf_dict[transcript]['attributes'].get('transcript_version','.'))
    print(gene_name,gene_id,transcript_id)
    result_list.append(
        '\t'.join([gene_name,transcript_id,chromosome,transcript_strand,transcript_start,
                    transcript_end,cds_start,cds_end,exon_count,exon_start,exon_end,transcript_length,
                    gene_id,transcript_biotype,transcript_version,coding_length
                    ])
        )
os.makedirs(os.path.dirname(output.anno),exists_ok=True)
with open(output.anno,'w') as f:
    f.write('\n'.join(result_list))
    





    



