#! /usr/bin/python

# Script written by Ivan Sovic.
# Copyright Ivan Sovic, 2016. All rights reserved.
#               www.sovic.org
#
# Script for left aligning indels in a SAM file. Provides a generic function for left alignemnt as well as a commandline usage implementation.

import os;
import sys;
import subprocess;

import fastqparser;

# VERBOSE_DEBUG = True;
VERBOSE_DEBUG = False;

from leftalign import *;

def split_cigar_string(cigar):
    ops = [];
    count_string = '';
    for i in xrange(0, len(cigar)):
        if (cigar[i] in 'M=XIDSH'):
            op = CigarOp(cigar[i], int(count_string));
            ops.append(op);
            count_string = '';
        elif (cigar[i].isdigit):
            count_string += cigar[i];
        else:
            sys.stderr.write('ERROR: Unknown CIGAR op! Skipping the op.\n');
            count_string = '';
    return ops;

def parse_sam_line(raw_line, rname_hash):
    split_line = raw_line.split('\t');

    sam_line = BamAlignment();

    sam_line.Filename = ''; # in_sam;
    sam_line.Name = split_line[0];
    sam_line.AlignmentFlag = int(split_line[1]);
    sam_line.RefName = split_line[2];
    sam_line.RefID = rname_hash[sam_line.RefName] if (sam_line.RefName in rname_hash) else -1;
    sam_line.Position = int(split_line[3]) - 1;         # BAM is 0-based, and SAM 1-based.
    sam_line.MapQuality = int(split_line[4]);
    sam_line.CigarString = split_line[5];
    sam_line.CigarData = split_cigar_string(sam_line.CigarString);
    sam_line.MateRefName = split_line[6];
    sam_line.MateRefID = rname_hash[sam_line.MateRefName] if (sam_line.MateRefName in rname_hash) else -1;
    sam_line.MatePosition = int(split_line[7]);
    sam_line.InsertSize = int(split_line[8]);
    sam_line.QueryBases = split_line[9];
    sam_line.Qualities = split_line[10];
    sam_line.Length = len(sam_line.QueryBases);

    sam_line.Bin = 0;
    sam_line.AlignedBases = '';
    sam_line.TagData = '';

    return sam_line;

def parse_sam(in_sam):
    try:
        fp_in = open(in_sam, 'r');
    except IOError, e:
        sys.stderr.write('ERROR: Could not open file %s for reading! Exiting.\n' % (in_sam));
        sys.stderr.write(str(e));
        exit(1);

    headers = [];
    sam_lines = [];
    rname_hash = {};
    for line in fp_in:
        line = line.strip();
        if (len(line) == 0 or (len(line) > 0 and line[0] == '@' and line.startswith('@SQ') == False)):
            headers.append(line);
            continue;
        
        if (line.startswith('@SQ') == True):
            split_line = line.split('\t');
            for val in split_line[1:]:
                if (val.startswith('SN:')): rname_hash[val.split('SN:')[-1]] = len(rname_hash.keys());
            headers.append(line);
            continue;

        sam_line = parse_sam_line(line, rname_hash);
        sam_lines.append(sam_line);

    return [headers, sam_lines];

def load_and_process_sam(in_sam, ref_file, fp_out):
    [sam_headers, sam_lines] = parse_sam(in_sam);
    [headers_ref, seqs_ref, quals_ref] = fastqparser.read_fastq(ref_file);

    seqs_ref_hash = {};
    for i in xrange(0, len(seqs_ref)):
        seqs_ref_hash[headers_ref[i]] = seqs_ref[i];
        seqs_ref_hash[headers_ref[i].split()[0]] = seqs_ref[i];

    for sam_line in sam_lines:
        if (sam_line.IsMapped() == False): continue;
        # print sam_line.verbose();
        seq_ref = seqs_ref_hash[sam_line.RefName];
        stablyLeftAlign(sam_line, seq_ref, 1, False);
        # print sam_line.verbose();
        # print '';
        # for c in sam_line.CigarData:
        #     print c.verbose();

def process_sam_on_the_fly(in_sam, ref_file, fp_out):
    if (not os.path.exists(in_sam)):
        sys.stderr.write('ERROR: File "%s" does not exist!\n' % (in_sam));
        exit(1);
    if (not os.path.exists(ref_file)):
        sys.stderr.write('ERROR: File "%s" does not exist!\n' % (ref_file));
        exit(1);

    [headers_ref, seqs_ref, quals_ref] = fastqparser.read_fastq(ref_file);
    seqs_ref_hash = {};
    for i in xrange(0, len(seqs_ref)):
        seqs_ref_hash[headers_ref[i]] = seqs_ref[i];
        seqs_ref_hash[headers_ref[i].split()[0]] = seqs_ref[i];

    try:
        fp_in = open(in_sam, 'r');
    except IOError, e:
        sys.stderr.write('ERROR: Could not open file %s for reading! Exiting.\n' % (in_sam));
        sys.stderr.write(str(e));
        exit(1);

    headers = [];
    sam_lines = [];
    rname_hash = {};
    num_alignments = 0;
    for line in fp_in:
        line = line.strip();
        if (len(line) == 0 or (len(line) > 0 and line[0] == '@' and line.startswith('@SQ') == False)):
            fp_out.write('%s\n' % (line));
            headers.append(line);
            continue;
        
        if (line.startswith('@SQ') == True):
            split_line = line.split('\t');
            for val in split_line[1:]:
                if (val.startswith('SN:')): rname_hash[val.split('SN:')[-1]] = len(rname_hash.keys());
            fp_out.write('%s\n' % (line));
            headers.append(line);
            continue;

        sam_line = parse_sam_line(line, rname_hash);
        num_alignments += 1;

        if (sam_line.QueryBases == '*'):
            fp_out.write('%s\n' % (line))
            continue

        if (sam_line.IsMapped() == False):
            fp_out.write('%s\n' % (line));
            continue;

        seq_ref = seqs_ref_hash[sam_line.RefName];
        stablyLeftAlign(sam_line, seq_ref, 50, False);

        split_line = line.split('\t');
        split_line[5] = sam_line.CigarString;
        new_line = '\t'.join(split_line);
        fp_out.write('%s\n' % (new_line));

        sys.stderr.write('\rProcessed %d alignments...' % (num_alignments));

    sys.stderr.write('done!\n\n');
    fp_in.close();

def run_bamleftalign(in_sam, ref_file):
    filename = os.path.splitext(in_sam)[0];
    command = '../../samscripts/src/convert_to_bam.sh %s' % (filename);
    sys.stderr.write('Running: "%s"\n' % (command));
    subprocess.call(command, shell=True);
    command = 'samtools faidx %s' % (ref_file);
    sys.stderr.write('Running: "%s"\n' % (command));
    subprocess.call(command, shell=True);
    command = 'cat %s-sorted.bam | test-leftalign/bamleftalign/bamleftalign -d -f %s > %s.leftaligned.bam' % (filename, ref_file, filename);
    sys.stderr.write('Running: "%s"\n' % (command));
    subprocess.call(command, shell=True);
    command = 'samtools view -h %s.leftaligned.bam > %s.leftaligned.sam' % (filename, filename);
    sys.stderr.write('Running: "%s"\n' % (command));
    subprocess.call(command, shell=True);



def main():
    if (len(sys.argv) < 3 or len(sys.argv) > 4):
        sys.stderr.write('Script for left-aligning indels in SAM files.\n');
        sys.stderr.write('Usage:\n');
        sys.stderr.write('  %s <ref.fasta> <in.sam> [<out.sam>]\n' % (sys.argv[0]));
        exit(1);

    in_ref = sys.argv[1];
    in_sam = sys.argv[2];
    fp_out = sys.stdout;
    if (len(sys.argv) >= 4):
        out_sam = sys.argv[3];
        try:
            fp_out = open(out_sam, 'w');
        except IOError, e:
            sys.stderr.write('ERROR: Could not open file %s for writing! Exiting.\n' % (out_sam));
            sys.stderr.write(str(e));
            exit(1);

    process_sam_on_the_fly(in_sam, in_ref, fp_out);
    # run_bamleftalign(in_sam, in_ref);

    if (fp_out != sys.stdout):
        fp_out.close();

if __name__ == "__main__":
    main();




# ../../samscripts/src/convert_to_bam.sh msa3.handmade
# samtools faidx msa3.handmade.ref.fa
# cat msa3.handmade-sorted.bam | test-leftalign/bamleftalign/bamleftalign -d -f msa3.handmade.ref.fa > msa3.handmade-sorted.leftaligned.bam
# samtools view -h msa3.handmade-sorted.leftaligned.bam > msa3.handmade-sorted.leftaligned.sam
