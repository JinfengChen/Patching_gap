#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def usage():
    test="name"
    message='''
python Patching_gap.py --fasta ../input/ACT/OS_Chr10_2616367_3045941/HEG4_chr10_2404092_2868093.fasta --gff ../input/ACT/OS_Chr10_2616367_3045941/HEG4_chr10_2404092_2868093.gene.gff --sequence ../input/sequence/mPing_inversion.gap.fa
Patching one gap in the sequence and update gff after the position of gap.
--fasta: fasta sequence to patch gaps
--gff: gene gff
--sequence: gap sequence by Sanger, need to on plus strand. seqid must be gapseq:141357-141503, 141357-141503 is a approximate position of gap in reference sequence.

    '''
    print message


'''

'''
def update_act(fasta, path, prefix):
    repeat = '/usr/local/bin/RepeatMasker -species rice -q -xsmall -nolow -no_is -norna %s\n' % (fasta)
    repeat += 'perl /rhome/cjinfeng/software/bin/repeat_to_gff.pl %s' % (fasta + '.out')
    os.system(repeat)
    
    embl  = 'perl ./act/GFF2embl.pl -gff %s -embl %s -fasta %s\n' % (path + '/' + prefix + '.gene.update.gff', path + '/' + prefix + '.gene.update.embl' , path + '/' + prefix + '.update.fasta')
    embl += 'perl ./act/gffrepeat2embl.pl -repeat %s -embl %s -title %s\n' % (fasta + '.out.gff', path + '/' + prefix + '.gene.update.embl', path + '/' + prefix + '.update')
    embl += 'mv %s %s' % (path + '/' + prefix + '.update.merge', path + '/' + prefix + '.update.embl')    
    #print embl
    os.system(embl)    

    act   = 'cp %s ./\n' % (path + '/' + '*.fasta')  
    act  += 'perl ./act/runblast2seq.pl\n'
    act  += 'perl ./act/run2act.pl\n'
    act  += 'rm *.fasta *.blast *.temp *.nhr *.nin *.nsq\n'
    act  += 'mv *4ACT %s' % (path)
    #print act
    os.system(act)
 
'''
chr10   maker   gene    9314    13560   .       +       .       ID=HEG4_Os10g12220;Name=queuine tRNA-ribosyltransferase,putative,expressed;
chr10   maker   mRNA    9314    13560   .       +       .       ID=HEG4_Os10g12220.1;Parent=HEG4_Os10g12220;
'''
def update_gff(gff, update, path, update_file):
    ofile = open (path + '/' + update_file, 'w')
    with open (gff, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if int(unit[3]) > update[0]:
                    unit[3] = str(int(unit[3]) + update[1])
                    unit[4] = str(int(unit[4]) + update[1])
                    newline = '\t'.join(unit)
                    print >> ofile, newline
                else:
                    print >> ofile, line   
    ofile.close()


def write_seq(seqid, seq, path):
    ofile = open (path + '/' + seqid + '.update.fasta', 'w')
    newrecord = SeqRecord(Seq(seq), id=seqid + '.update', description="")
    SeqIO.write(newrecord, ofile, "fasta")
    ofile.close()

'''
get fasta sequence in file, return dict with id->seq
'''
def fasta_seq(fastafile):
    data = defaultdict(str) 
    for record in SeqIO.parse(fastafile,"fasta"):
        data[str(record.id)] = str(record.seq)
    return data


'''
parse Sanger sequence and we know where the gap in the reference
'''
def gap_seq(fastafile):
    gap  = []
    s  = re.compile(r'\w+\:(\d+)\-(\d+)')
    for record in SeqIO.parse(fastafile,"fasta"):
        m = s.search(str(record.id))
        length = len(record.seq)
        print length
        if m:
            gap.append(int(m.groups(0)[0]))
            gap.append(int(m.groups(0)[1]))
            gap.append(int(length))
    return gap

'''
blast:
gapseq:141357-141503    HEG4_chr10_2404092_2868093      100.00  740     0       0       30      769     140623  141362  0.0     1467
gapseq:141357-141503    HEG4_chr10_2404092_2868093      99.86   691     0       1       338     1028    262508  261819  0.0     1354
gapseq:141357-141503    HEG4_chr10_2404092_2868093      100.00  435     0       0       335     769     261988  261554  0.0      862
gapseq:141357-141503    HEG4_chr10_2404092_2868093      99.60   252     0       1       862     1112    141490  141741  9e-138   484
gapseq:141357-141503    HEG4_chr10_2404092_2868093      100.00  170     0       0       859     1028    140928  141097  1e-93    337
gapseq:141357-141503	HEG4_chr10_2404092_2868093	100.00	167	0	0	338	504	141490	141656	8e-92	 331
gap:
[141357, 141503]
'''
def parse_blastm8(blast, gap, ref, qry, path):
    left  = []
    right = []
    
    with open (blast, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if float(unit[2]) >= 99:
                    start = int(unit[8]) if int(unit[8]) < int(unit[9]) else int(unit[9])
                    end   = int(unit[9]) if int(unit[8]) < int(unit[9]) else int(unit[8])
                    strand= '+' if int(unit[8]) < int(unit[9]) else '-'
                    '''blast hit near to within 20 bp of gap boundary'''
                    if abs(int(unit[6]) - 0) <= 50 and abs(end - gap[0]) < 20:
                        print 'left hit: %s' % ('\t'.join(unit))
                        left.append(unit)
                    elif abs (start - gap[1]) < 20 and abs(int(unit[7]) - gap[2]) <= 50:
                        print 'right hit: %s' % ('\t'.join(unit))
                        right.append(unit)
    
    filler_pos = []
    gap_pos    = []                
    if len(left) == 1:
        temp1 = '\t'.join(left[0])
        filler_pos.append(int(left[0][7])+1-1-5) ## str is 0-based and we replace 5 more base flanking the gap
        gap_pos.append(int(left[0][9])+1-1-5)
        print 'Left Anchor: %s\t%s\t%s\t%s' % (left[0][6], left[0][7], left[0][8], left[0][9])
        
    if len(right) == 1:
        temp1 = '\t'.join(right[0])
        filler_pos.append(int(right[0][6])-1-1+1+5) ## str is 0-based and str[1:2] will get the the second word, not including 2
        gap_pos.append(int(right[0][8])-1-1+1+5)
        print 'Right Anchor: %s\t%s\t%s\t%s' % (right[0][6], right[0][7], right[0][8], right[0][9])
  
        '''filler sequence'''
        filler = qry[right[0][0]][filler_pos[0]:filler_pos[1]]
        gapper  = ref[right[0][1]][gap_pos[0]:gap_pos[1]]
        print 'Filler pos and sequence: %s %s %s' % (filler_pos[0], filler_pos[1], filler)
        print 'Gaper pos and sequence: %s %s %s' % (gap_pos[0], gap_pos[1], gapper)
 
        '''patching'''
        before_gap = ref[right[0][1]][0:gap_pos[0]]
        after_gap  = ref[right[0][1]][gap_pos[1]:]
        filled     = before_gap + filler + after_gap
        boundary   = gap_pos[1] ### after this boundary, the gene gff are need to update
        update_len = int(len(filler) - len(gapper)) ### after boundary, the gene gff need to add update_len


        '''update files'''
        write_seq(right[0][1], filled, path)
        return [boundary, update_len] 

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta')
    parser.add_argument('-g', '--gff')
    parser.add_argument('-s', '--sequence')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.sequence) > 0
    except:
        usage()
        sys.exit(2)

    '''parse gap position'''
    gap = gap_seq(args.sequence)
    print '%s %s %s' % (gap[0], gap[1], gap[2])
    
    '''get sequence'''
    ref = fasta_seq(args.fasta)
    qry = fasta_seq(args.sequence) 

    '''blast gap sequence to reference'''
    directory  = os.path.abspath(os.path.dirname(args.fasta))
    ref_prefix = os.path.basename(os.path.splitext(args.fasta)[0])
    qry_prefix = os.path.basename(os.path.splitext(args.sequence)[0])
    blastout   = '%sVS%s.blastm8' % (qry_prefix, ref_prefix)
    print '%s %s %s' % (ref_prefix, qry_prefix, blastout)
    blast = 'formatdb -i %s -p F\n' % (args.fasta)
    blast += 'blastall -p blastn -i %s -d %s -e 1e-5 -o %s -m 8' % (args.sequence, args.fasta, blastout)
    #print blast
    os.system(blast)

    '''parse blast to patch gap'''
    update = parse_blastm8(blastout, gap, ref, qry, directory)
    os.system('rm %s %s' % (blastout, args.fasta + '.*'))
    update_gff(args.gff, update, directory, ref_prefix + '.gene.update.gff')
    update_act(directory + '/' + ref_prefix + '.update.fasta', directory, ref_prefix)

if __name__ == '__main__':
    main()

