import sys




with open(snakemake.log[0], "w") as log_f:
    sys.stderr = sys.stdout = log_f
    

def get_geneanno(list_gene):
    lw = open(options.output, 'w+')
    list_tt=[]
    with open(options.input, 'r') as gtf:
        line = gtf.readline()
        while (line):
            row = line.strip().split("\t")
            if line.startswith("#") or row[2] != 'gene':
                line = gtf.readline()
                continue
            chr = row[0]
            dir = row[6]
            # print(row)
            start, end = int(row[3]) - 1, int(row[4])
            # print(list_gene)
            try:
                row2 = [i.split(" ") for i in row if re.search('gene_id', i)][0]
                gene_id = row2[row2.index("gene_id") + 1].replace('\"', '').strip(";")
                if gene_id not in list_gene:
                    trans_id = gene_id
                    gene_biotype = row2[row2.index("gene_biotype") + 1].replace('\"', '').strip(";")
                    # print(gene_id,gene_biotype)
                    if not trans_id in dictInfo:
                        dictInfo[trans_id] = {'start': [],
                                              'end': [],
                                              'dir': "",
                                              'chr': "",
                                              'gene_id': ""}
                        dictInfo[trans_id]['dir'] = dir
                        dictInfo[trans_id]['chr'] = chr
                        dictInfo[trans_id]['gene_id'] = gene_id
                        dictInfo[trans_id]['start'] = dictInfo[trans_id]['start'] + [start]
                        dictInfo[trans_id]['end'] = dictInfo[trans_id]['end'] + [end]
                    line = gtf.readline()
                else:
                    line = gtf.readline()
            except ValueError:
                line = gtf.readline()
        for trans_id in dictInfo.keys():
            exonCount = str(len(dictInfo[trans_id]['start']))
            exonStarts = ','.join([str(i) for i in sorted(dictInfo[trans_id]['start'])]) + ","
            exonEnds = ','.join([str(i) for i in sorted(dictInfo[trans_id]['end'])]) + ","
            list_tt.append("\t".join(
                ['.', trans_id, dictInfo[trans_id]['chr'], dictInfo[trans_id]['dir'], 'txStart', 'txEnd', 'cdsStart',
                 'cdsEnd', exonCount, exonStarts, exonEnds, ".", dictInfo[trans_id]['gene_id'], '.', '.', '.']))
            # print(len(dictInfo),list_tt)
        lw.writelines("\n".join(list_tt)+"\n")


if __name__ == "__main__":
    # Parser
    usage = "Usage: python %prog -i <gtf> > output.anno"
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument("-i", dest="input", help="Input file")
    parser.add_argument("-o", dest="output", help="Output file")
    options = parser.parse_args()
    # format:
    # transcript_id gene_id gene_name chr dir 5_UTR_len CDS_len 3_UTR_len 5_UTR start_codon CDS stop_codon 3_UTR 5_len CDS_len 3_len Exon
    dictInfo = {}
    list_gene=[]
    with open(options.input, 'r') as gtf:
        line = gtf.readline()
        while (line):
            row = line.strip().split("\t")
            if line.startswith("#") or row[2] == 'gene':
                line = gtf.readline()
                continue
            chr = row[0]
            dir = row[6]
            # print(row)
            start, end = int(row[3]) - 1, int(row[4])
            try:
                row2 = [i.split(" ") for i in row if re.search('transcript_id', i)][0]
                trans_id = row2[row2.index("transcript_id") + 1].replace('\"', '').strip(";")
                gene_id = row2[row2.index("gene_id") + 1].replace('\"', '').strip(";")
                try:
                    gene_name = row2[row2.index("gene") + 1].replace('\"', '').strip(";")
                except ValueError:
                    gene_name = gene_id
                type = row[2]
                if not trans_id in dictInfo:
                    dictInfo[trans_id] = {'start': [],
                                          'end': [],
                                          'dir': "",
                                          'chr': "",
                                          'gene_id': ""}
                    dictInfo[trans_id]['dir'] = dir
                    dictInfo[trans_id]['chr'] = chr
                    dictInfo[trans_id]['gene_id'] = gene_id
                if type == "exon" or type == "tRNAscan":
                    dictInfo[trans_id]['start'] = dictInfo[trans_id]['start'] + [start]
                    dictInfo[trans_id]['end'] = dictInfo[trans_id]['end'] + [end]
                line = gtf.readline()
            except ValueError:
                line = gtf.readline()
    for trans_id in dictInfo.keys():
        gene_name2=dictInfo[trans_id]['gene_id']
        if gene_name2 not in list_gene:
            list_gene.append(gene_name2)
    get_geneanno(list_gene)
    
    
from rich import progress
import time


# %% ../../nbs/0002_utils_gtf.ipynb 3
import logging
logger = logging.getLogger('baleen')

# %% ../../nbs/0002_utils_gtf.ipynb 4
def _parse_attribute(attribute: str):
    attribute_dict = dict()
    for att in attribute.strip('\n').split(';'):
        res = att.split(' "')
        if len(res) == 2:
            key, value = res
            attribute_dict[key.strip()] = value.replace('"', '')
    return attribute_dict


# %% ../../nbs/0002_utils_gtf.ipynb 5



# %% ../../nbs/0002_utils_gtf.ipynb 6
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
    
    

# %% ../../nbs/0002_utils_gtf.ipynb 7
def readGTF(gtf_file):
    '''
    GTF file is 1-based, and the end position is included.
    '''
    gtf_dict = dict()
    transcript_id = None
    
    with progress.open(gtf_file, 'r',description='Loading GTF:') as f:
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
        time.sleep(1)
    if transcript_id is not None:
        _sort_and_parse_exons(gtf_dict[transcript_id])
    return gtf_dict