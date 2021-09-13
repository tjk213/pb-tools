#!/usr/bin/env python3

import sys

counts = {}

sample_files = [
    'S4wk_B5.csv',
    'S4wk_B14.csv',
    'S6wk_A6.csv',
    'S6wk_A15.csv',
    'S6wk_B7.csv',
    'S6wk_B16.csv',
    'S8wk_A8.csv',
    'S8wk_A17.csv',
    'S8wk_B18.csv'
]

num_samples = len(sample_files)

out_rnas = set([
    'TAAAAACGCTGGCGGCCTAG',  # Non Targeting
    'GGTCATGGGCGACATGGACG',  # Fam118b      
    'TCACAGACCTGACCCTGCTG',  # Nfyc         
    'TACTTGGCAAGAGTGTCCCG',  # Vmn1r197     
    'CCTCATCACCCTGCCCCGCG',  # Srgap3       
    'CGTGACCTACTCGAACGTGG',  # Hnrnpa0      
    'CGGGCGTGAAAATGTACCAA',  # Smim11       
    'TCAGCGCGAAGGGCACGCGG',  # Gdnf         
    'TCAGACCAACCCAGGACCAG',
    'AGGACTTGTACTTACCACAC',
    'CTACAGTCCCTTTGACGTGG',
    'GGATGAAAGCAGCGCGTCTG',
    'TCCTTGGCCTCTGATTCCCG',
    'CTTCCCTGAAGAAAGGACGT',
    'CATCTGCTGCGTAGTCAGCA',
    'TGACAGCAATCCACCTCCGA',
    'ACATGTTTGTTTATATGACA',
    'CGCACGTCCATGCTGCACAG',
    'CTACGTTTGTACCTTCAGTG',
    'GCAGTCACAAATCATTGTTG',
    'TATGCCCTACGACCTGTGTG',
    'CAGCTTGTGGCAAACTTACG',
    'ACTACATGAAATAGTTAGCA',
    'ACTATTCAAGCCAATGTCGT',
    'AGTACACGTTCTTCACGTTG',
    'TGACACGAACTCGCTTGTCA',
    'CCTCTGTGGAGGATCCACCA',
    'TAATTTAACATCCTTACCTG',
    'TTGCCATGCTATCAAAGTTG',
    'CCGCGCGACTAGGATCGCGT',
    'CTGGTCTCAAGCGCCAGTGG',
    'CAGATGCCCCCCAAGTCAAG',
    'TTCTGTACCAGGGGTCGTCG',
    'AAGTTGGCCGTTCAGTCGAT',
    'CATGACATCGTGCCACCCGT',
    'ACCAGACAGTCACTCTCTCG',
    'GAGTGCGCAGTAGTACTGAA',
    'GCTTGCTCTGCTGCTCACTG',
    'TGTGTAGCCAGATGAGACCA',
    'CCAGATGTATTCCATCCCCG',
    'TTCTGGGGTGGAACGCACGT',
    'AAGGCACGAAGACATAGACA',
    'GCCGTAAGCGGGCCGGTTGA',
    'GGATATTCGCGCGGTCTTCA',
    'GTACATGCGGCAAGTCGACT',
    'TCACAACCCCCGACTATCGC',
])

def add_count(idx, key, val):
    if not key in counts.keys():
        counts[key] = ['0' for i in range(num_samples)]
    counts[key][idx] = val
    return

def main():

    print('Parsing gRNA samples...',file=sys.stderr,end='',flush=True)

    for i,sample_file in enumerate(sample_files):
        samples = open(sample_file,'r')
        lines = [x for x in samples]

        # Verify header

        labels = lines[0].strip().split(',')
        assert labels[0] == 'gene',    f'Unexpected label [file={sample_file}]'
        assert labels[1] == 'sgRNA',   f'Unexpected label [file={sample_file}]'
        assert labels[2] == 'gRNASeq', f'Unexpected label [file={sample_file}]'
        assert labels[3] == 'Count',   f'Unexpected label [file={sample_file}]'

        # Parse samples

        for line in lines[1::]:
            id, gene, full_name, rna_seq, count = line.strip().split(',')
            add_count(i,(gene,rna_seq),count)

        samples.close()
    print('done.',file=sys.stderr)

    print('Filtering samples...',file=sys.stderr,end='',flush=True)
    final_samples = []
    for gene,rna_seq in counts.keys():
        if rna_seq in out_rnas:
            sample = [gene,rna_seq]
            for cnt in counts[(gene,rna_seq)]:
                sample.append(cnt)
            final_samples.append(sample)
    print('done.',file=sys.stderr,flush=True)

    print('Sorting samples...',file=sys.stderr,end='',flush=True)
    final_samples.sort(key=lambda x: '' if x[0].startswith('NonTargeting') else x[0])
    print('done.',file=sys.stderr)

    print('Saving results...',file=sys.stderr,end='',flush=True)
    columns = 'Gene, gRNA, ' + ','.join(sample_files)
    print(columns)

    for sample in final_samples:
        print(', '.join(sample))
    print('done.',file=sys.stderr)
    
    return False

if __name__ == '__main__': sys.exit(main())
