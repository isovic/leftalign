#! /usr/bin/python

# Script written by Ivan Sovic.
# Copyright Ivan Sovic, 2016. All rights reserved.
#               www.sovic.org
#
# Script for left aligning indels in a SAM file. Provides a generic function for left alignemnt as well as a commandline usage implementation.
#
# The script is based on the C++ code provided in the Freebayes package (https://github.com/ekg/freebayes),
# concretely, https://github.com/ekg/freebayes/blob/master/src/LeftAlign.cpp . The licence of the original implementation is listed below.
#
# This code is more-or-less a literal translation from the original C++ implementation to Python. Some things might be implemented more optimal.
#
# The original Freebayes licence:
# Copyright (c) 2010 Erik Garrison, Gabor Marth

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import os;
import sys;
import subprocess;

import fastqparser;

# VERBOSE_DEBUG = True;
VERBOSE_DEBUG = False;

class CigarOp:
    def __init__(self, op_type = '\0', op_len = 0):
        self.Type = op_type;			# char     Type;   //!< CIGAR operation type (MIDNSHPX=)
        self.Length = op_len;			# uint32_t Length; //!< CIGAR operation length (number of bases)

    def verbose(self):
        return '%d%s' % (self.Length, self.Type);

class BamAlignment:
    def __init__(self):
        self.Name = '';					# std::string Name;               # read name
        self.Length = 0;				# int32_t     Length;             # length of query sequence
        self.QueryBases = '';			# std::string QueryBases;         # 'original' sequence (contained in BAM file)
        self.AlignedBases = '';			# std::string AlignedBases;       # 'aligned' sequence (QueryBases plus deletion, padding, clipping chars)
        self.Qualities = '';			# std::string Qualities;          # FASTQ qualities (ASCII characters, not numeric values)
        self.TagData = '';				# std::string TagData;            # tag data (use provided methods to query/modify)
        self.RefID = 0;					# int32_t     RefID;              # ID number for reference sequence
        self.Position = 0;				# int32_t     Position;           # position (0-based) where alignment starts
        self.Bin = 0;					# uint16_t    Bin;                # BAM (standard) index bin number for this alignment
        self.MapQuality = 0;			# uint16_t    MapQuality;         # mapping quality score
        self.AlignmentFlag = 0;			# uint32_t    AlignmentFlag;      # alignment bit-flag (use provided methods to query/modify)
        self.CigarData = [];			# std::vector<CigarOp> CigarData; # CIGAR operations for this alignment
        self.MateRefID = 0;				# int32_t     MateRefID;          # ID number for reference sequence where alignment's mate was aligned
        self.MatePosition = 0;			# int32_t     MatePosition;       # position (0-based) where alignment's mate starts
        self.InsertSize = 0;			# int32_t     InsertSize;         # mate-pair insert size
        self.Filename = '';				# std::string Filename;           # name of BAM file which this alignment comes from

        self.RefName = '';
        self.MateRefName = '';
        self.CigarString = '';

    def IsMapped(self):
        if (self.AlignmentFlag & 0x0004): return False;
        return True;

    def IsReverse(self):
        if (self.AlignmentFlag & 0x0010): return True;
        return False;

    def verbose(self):
        return ('%s\t%d\t%s' % (self.Name, self.AlignmentFlag, self.CigarString));



class FBIndelAllele:
    def __init__(self, i, l, p, rp, s, n):
	    self.insertion = i;	# bool
	    self.length = l;		# int
	    self.position = p;	# int
	    self.readPosition = rp;	# int
	    self.sequence = s;		# string
	    self.splice = n;			# bool

    def verbose(self):
        t = 'i' if (self.insertion) else 'd';
        out = '%s:%d:%d:%s:%s' % (t, self.position, self.readPosition, self.sequence, ('splice' if (self.splice) else ''));
        return out;

    def homopolymer(self):
        if (len(self.sequence) == 0): return False;
        c = self.sequence[0];
        for s in self.sequence:
            if (s != c): return False;
        return True;

    def is_equal_to(self, op2):
        return (self.insertion == op2.insertion
            and self.length == op2.length
            and self.position == op2.position
            and self.sequence == op2.sequence
            and self.splice == op2.splice);

    def is_less_than(self, op2):
        return (self.verbose() < op2.verbose());

def FBhomopolymer(sequence):
    if (len(sequence) == 0): return False;
    c = sequence[0];
    for s in sequence:
        if (s != c): return False;
    return True;

def LEFTALIGN_DEBUG(msg):
    if (VERBOSE_DEBUG == True):
       sys.stderr.write(msg);

#include "LeftAlign.h"

#bool debug;

# Attempts to left-realign all the indels represented by the alignment cigar.
#
# This is done by shifting all indels as far left as they can go without
# mismatch, then merging neighboring indels of the same class.  leftAlign
# updates the alignment cigar with changes, and returns true if realignment
# changed the alignment cigar.
#
# To left-align, we move multi-base indels left by their own length as long as
# the preceding bases match the inserted or deleted sequence.  After this
# step, we handle multi-base homopolymer indels by shifting them one base to
# the left until they mismatch the reference.
#
# To merge neighboring indels, we iterate through the set of left-stabilized
# indels.  For each indel we add a new cigar element to the new cigar.  If a
# deletion follows a deletion, or an insertion occurs at the same place as
# another insertion, we merge the events by extending the previous cigar
# element.
#
# In practice, we must call this function until the alignment is stabilized.
#
# def leftAlign(BamAlignment& alignment, string& referenceSequence, bool debug):
def leftAlign(alignment, referenceSequence, debug):

    arsOffset = 0; # pointer to insertion point in aligned reference sequence
    alignedReferenceSequence = referenceSequence + '';
    aabOffset = 0;
    alignmentAlignedBases = alignment.QueryBases + '';

    # store information about the indels
    indels = [];			# vector<FBIndelAllele> indels;

    rp = 0;  # read position, 0-based relative to read
    sp = 0;  # sequence position

    softBegin = '';
    softEnd = '';

#    stringstream cigar_before, cigar_after;

    cigar_before = '';
    cigar_after = '';

#    for (vector<CigarOp>::const_iterator c = alignment.CigarData.begin(); c != alignment.CigarData.end(); ++c) {
    for c in alignment.CigarData:
        l = c.Length;
        t = c.Type;
        cigar_before += '%d%s' % (l, t);

        if (t == 'M' or t == 'X' or t == '='): # match or mismatch
            sp += l;
            rp += l;
        elif (t == 'D'): # deletion
            indels.append(FBIndelAllele(False, l, sp, rp, referenceSequence[sp:(sp+l)], False));
            alignmentAlignedBases = alignmentAlignedBases[0:(rp + aabOffset)] + '-'*l + alignmentAlignedBases[(rp + aabOffset):];	# alignmentAlignedBases.insert(rp + aabOffset, string(l, '-'));
            aabOffset += l;
            sp += l;  # update reference sequence position
        elif (t == 'N'):
            indels.append(FBIndelAllele(False, l, sp, rp, referenceSequence[sp:(sp+l)], True));
            alignmentAlignedBases = alignmentAlignedBases[0:(rp + aabOffset)] + '-'*l + alignmentAlignedBases[(rp + aabOffset):];	# alignmentAlignedBases.insert(rp + aabOffset, string(l, '-'));
            aabOffset += l;
            sp += l;  # update reference sequence position
        elif (t == 'I'): # insertion
            indels.append(FBIndelAllele(True, l, sp, rp, alignment.QueryBases[rp:(rp+l)], False));
            alignedReferenceSequence = alignedReferenceSequence[0:(sp + len(softBegin) + arsOffset)] + '-'*l + alignedReferenceSequence[(sp + len(softBegin) + arsOffset):];	# alignedReferenceSequence.insert(sp + softBegin.size() + arsOffset, string(l, '-'));
            arsOffset += l;
            rp += l;
        elif (t == 'S'): # soft clip, clipped sequence present in the read not matching the reference
            # remove these bases from the refseq and read seq, but don't modify the alignment sequence
            if (rp == 0):
                alignedReferenceSequence = '*'*l + alignedReferenceSequence;	# alignedReferenceSequence = string(l, '*') + alignedReferenceSequence;
                softBegin = alignmentAlignedBases[0:l];							# softBegin = alignmentAlignedBases.substr(0, l);
            else:
                alignedReferenceSequence = alignedReferenceSequence + '*'*l;		# alignedReferenceSequence = alignedReferenceSequence + string(l, '*');
                softEnd = alignmentAlignedBases[(len(alignmentAlignedBases)-l):];	# softEnd = alignmentAlignedBases.substr(alignmentAlignedBases.size() - l, l);
            rp += l;
        elif (t == 'H'): # hard clip on the read, clipped sequence is not present in the read
        	pass;
            #} else if (t == 'N') { # skipped region in the reference not present in read, aka splice
            #sp += l;


    alignedLength = sp + 0;
    LEFTALIGN_DEBUG('| ' + cigar_before + '\n' + '| ' + alignedReferenceSequence + '\n' + '| ' + alignmentAlignedBases + '\n');

    # if no indels, return the alignment
    if (len(indels) == 0):
    	return False;

    # for each indel, from left to right
    #     while the indel sequence repeated to the left and we're not matched up with the left-previous indel
    #         move the indel left

    previous = indels[0];	# vector<FBIndelAllele>::iterator previous = indels.begin();
#    for indel in indels:		# for (vector<FBIndelAllele>::iterator id = indels.begin(); id != indels.end(); ++id) {
    for indel_id in xrange(0, len(indels)):
        indel = indels[indel_id];

        # left shift by repeats
        #
        # from 1 base to the length of the indel, attempt to shift left
        # if the move would cause no change in alignment optimality (no
        # introduction of mismatches, and by definition no change in gap
        # length), move to the new position.
        # in practice this moves the indel left when we reach the size of
        # the repeat unit.
        #
        steppos = 0;		# int steppos;
        readsteppos = 0;	# int readsteppos;
#        FBIndelAllele& indel = *id;
        i = 1;
        while (i <= indel.length):
            steppos = indel.position - i;
            readsteppos = indel.readPosition - i;

            if (VERBOSE_DEBUG == True):
                if (debug):
                    if (steppos >= 0 and readsteppos >= 0):
                        sys.stderr.write('(1.0) %s\n' % referenceSequence[steppos:(steppos + indel.length)]);
                    sys.stderr.write('(1.1) %s\n' % alignment.QueryBases[readsteppos:(readsteppos+indel.length)]);
                    sys.stderr.write('(1.2) %s\n' % indel.sequence);
                    sys.stderr.write('i = %d, indel.length = %d\n' % (i, indel.length));

            while (steppos >= 0 and readsteppos >= 0
                   and indel.splice == False
                   and indel.sequence == referenceSequence[steppos:(steppos + indel.length)]
                   and indel.sequence == alignment.QueryBases[readsteppos:(readsteppos + indel.length)]
                   and (indel is indels[0]
                       or (previous.insertion and steppos >= previous.position)
                       or (previous.insertion == False and steppos >= (previous.position + previous.length)))):

            	LEFTALIGN_DEBUG('(1.3) %s%s shifting %dbp left\n' % (('insertion' if (indel.insertion) else 'deletion'), indel.verbose(), i));	# LEFTALIGN_DEBUG((indel.insertion ? "insertion " : "deletion ") << indel << " shifting " << i << "bp left" << endl);

                indel.position -= i;
                indel.readPosition -= i;
                steppos = indel.position - i;
                readsteppos = indel.readPosition - i;

            i += 1;
            while (i <= indel.length and (indel.length % i) != 0):
                i += 1;         # do { ++i; } while (i <= indel.length and (indel.length % i) != 0);



        # left shift indels with exchangeable flanking sequence
        #
        # for example:
        #
        #    GTTACGTT           GTTACGTT
        #    GT-----T   ---->   G-----TT
        #
        # GTGTGACGTGT           GTGTGACGTGT
        # GTGTG-----T   ---->   GTG-----TGT
        #
        # GTGTG-----T           GTG-----TGT
        # GTGTGACGTGT   ---->   GTGTGACGTGT
        #
        #
        steppos = indel.position - 1;
        readsteppos = indel.readPosition - 1;
        while (steppos >= 0 and readsteppos >= 0
               and alignment.QueryBases[readsteppos] == referenceSequence[steppos]
               and alignment.QueryBases[readsteppos] == indel.sequence[-1]
               and (indel is indels[0]
                   or (previous.insertion and (indel.position - 1) >= previous.position)
                   or (previous.insertion == False and (indel.position - 1) >= (previous.position + previous.length)))):
            LEFTALIGN_DEBUG('%s%s exchanging bases %dbp left\n' % (('insertion' if (indel.insertion) else 'deletion'), indel.verbose(), 1));	# LEFTALIGN_DEBUG((indel.insertion ? "insertion " : "deletion ") << indel << " exchanging bases " << 1 << "bp left" << endl);
            indel.sequence = indel.sequence[-1] + indel.sequence[0:-1];
            indel.position -= 1;
            indel.readPosition -= 1;
            steppos = indel.position - 1;
            readsteppos = indel.readPosition - 1;
        # tracks previous indel, so we don't run into it with the next shift
        previous = indels[indel_id]
#    }





    # bring together floating indels
    # from left to right
    # check if we could merge with the next indel
    # if so, adjust so that we will merge in the next step
    if (len(indels) > 1):
        previous = indels[0];
        for indel in indels[1:]:	# for (vector<FBIndelAllele>::iterator id = (indels.begin() + 1); id != indels.end(); ++id) { FBIndelAllele& indel = *id;
            # parsimony: could we shift right and merge with the previous indel?
            # if so, do it
            prev_end_ref = previous.position if (previous.insertion) else (previous.position + previous.length);
            prev_end_read = previous.readPosition if (not previous.insertion) else (previous.readPosition + previous.length);
            if (previous.splice == False and indel.splice == False and
                previous.insertion == indel.insertion
                and ((previous.insertion
                and (previous.position < indel.position
                and previous.readPosition + previous.readPosition < indel.readPosition))
                or
                (previous.insertion == False
                and (previous.position + previous.length < indel.position)
                and (previous.readPosition < indel.readPosition)
                ))):
                if (previous.homopolymer()):
                    seq = referenceSequence[prev_end_ref:(indel.position)];
                    readseq = alignment.QueryBases[prev_end_read:(prev_end_read + indel.position - prev_end_ref)];
                    LEFTALIGN_DEBUG('seq: ' + seq + '\nreadseq: ' + readseq + '\n');
                    if (previous.sequence[0] == seq[0]
                            and FBhomopolymer(seq)
                            and FBhomopolymer(readseq)):
                        LEFTALIGN_DEBUG('moving %s right to %d\n' % (previous.verbose(), (indel.position if (indel.insertion) else indel.position - previous.length)));
                        previous.position = indel.position if (indel.insertion) else indel.position - previous.length;
                else:
                    pos = previous.position;
                    while (pos < int(len(referenceSequence)) and
                            ((previous.insertion and pos + previous.length <= indel.position)
                            or
                            (previous.insertion == False and pos + previous.length < indel.position))
                            and previous.sequence 
                                == referenceSequence[(pos + previous.length):(pos + previous.length + previous.length)]):
                        pos += previous.length;
                    if (pos < previous.position and
                        ((previous.insertion and pos + previous.length == indel.position)
                        or
                        (previous.insertion == False and pos == indel.position - previous.length))
                       ):
                        LEFTALIGN_DEBUG('right-merging tandem repeat: moving %s right to  %d\n' (previous.verbose(), pos));
                        previous.position = pos;
            
            previous = indel;

    # for each indel
    #     if ( we're matched up to the previous insertion (or deletion) 
    #          and it's also an insertion or deletion )
    #         merge the indels
    #
    # and simultaneously reconstruct the cigar

    newCigar = [];  # vector<CigarOp> newCigar;

    if (len(softBegin) > 0):
        newCigar.append(CigarOp('S', len(softBegin)));

    indel_id = 0;
#    indel = indels[indel_id];  # vector<FBIndelAllele>::iterator id = indels.begin();
    last = indels[indel_id];       # FBIndelAllele last = *id++;
    indel_id += 1;
    if (last.position > 0):
        newCigar.append(CigarOp('M', last.position));
    if (last.insertion):
        newCigar.append(CigarOp('I', last.length));
    elif (last.splice):
        newCigar.append(CigarOp('N', last.length));
    else:
        newCigar.append(CigarOp('D', last.length));

    lastend = last.position if (last.insertion) else (last.position + last.length);
    LEFTALIGN_DEBUG('%s,' % (last.verbose()));

    while (indel_id < len(indels)):         #for (; id != indels.end(); ++id) {
        indel = indels[indel_id];            # FBIndelAllele& indel = *id;
        LEFTALIGN_DEBUG('%s,' % (indel.verbose()));
        if (indel.position < lastend):
            sys.stderr.write('impossibility?: indel realigned left of another indel\n');
            sys.stderr.write('%s %d\n%s\n' % (alignment.Name, alignment.Position, alignment.QueryBases));
            exit(1);
        elif (indel.position == lastend and indel.insertion == last.insertion):
            op = newCigar[-1];  # CigarOp& op = newCigar.back();
            op.Length += indel.length;
        elif (indel.position >= lastend): # also catches differential indels, but with the same position
            newCigar.append(CigarOp('M', indel.position - lastend));
            if (indel.insertion):
                newCigar.append(CigarOp('I', indel.length));
            elif (indel.splice):
                newCigar.append(CigarOp('N', indel.length));
            else: # deletion
                newCigar.append(CigarOp('D', indel.length));

        last = indels[indel_id];
        indel_id += 1;
        lastend = last.position if (last.insertion) else (last.position + last.length);
    
    if (lastend < alignedLength):
        newCigar.append(CigarOp('M', alignedLength - lastend));

    if (len(softEnd) > 0):
        newCigar.append(CigarOp('S', len(softEnd)));

    LEFTALIGN_DEBUG('\n');

    if (VERBOSE_DEBUG == True):
        if (debug):
            for c in alignment.CigarData:   # for (vector<CigarOp>::const_iterator c = alignment.CigarData.begin(); c != alignment.CigarData.end(); ++c) {
                sys.stderr.write('%d%s' % (c.Length, c.Type));
            sys.stderr.write('\n');

    alignment.CigarData = newCigar;

    # print 'newCigar:\n%s' % (newCigar);

    for c in alignment.CigarData:   # for (vector<CigarOp>::const_iterator c = alignment.CigarData.begin(); c != alignment.CigarData.end(); ++c) {
        cigar_after += '%d%s' % (c.Length, c.Type);

    alignment.CigarString = cigar_after;

    LEFTALIGN_DEBUG('%s\n' % (cigar_after));

    # check if we're realigned
    if (cigar_after == cigar_before):
        return False;
    else:
        return True;

# def countMismatches(BamAlignment& alignment, string referenceSequence):
def countMismatches(alignment, referenceSequence):
    mismatches = 0;
    sp = 0;
    rp = 0;
    for c in alignment.CigarData:
        l = c.Length;
        t = c.Type;
        if (t == 'M' or t == 'X' or t == '='): # match or mismatch
            for i in xrange(0, l):
                if (alignment.QueryBases[rp] != referenceSequence[sp]):
                    mismatches += 1;
                sp += 1;
                rp += 1;
        elif (t == 'D'): # deletion
            sp += l;  # update reference sequence position
        elif (t == 'I'): # insertion
            rp += l;  # update read position
        elif (t == 'S'): # soft clip, clipped sequence present in the read not matching the reference
            rp += l;
        elif (t == 'H'): # hard clip on the read, clipped sequence is not present in the read
        	pass;
        elif (t == 'N'): # skipped region in the reference not present in read, aka splice
            sp += l;

    return mismatches;

# # Iteratively left-aligns the indels in the alignment until we have a stable
# # realignment.  Returns true on realignment success or non-realignment.
# # Returns false if we exceed the maximum number of realignment iterations.
# #
# def stablyLeftAlign(BamAlignment& alignment, string referenceSequence, int maxiterations, bool debug):
def stablyLeftAlign(alignment, referenceSequence, maxiterations, debug):
    if (VERBOSE_DEBUG == True):
        mismatchesBefore = countMismatches(alignment, referenceSequence);

    if (not leftAlign(alignment, referenceSequence, debug)):
        LEFTALIGN_DEBUG('did not realign\n');
        return True;

    else:

        maxiterations -= 1;
        while (maxiterations > 0 and leftAlign(alignment, referenceSequence, debug)):
            maxiterations -= 1;
            LEFTALIGN_DEBUG('realigning ...\n');

        if (VERBOSE_DEBUG == True):
            mismatchesAfter = countMismatches(alignment, referenceSequence);

            if (mismatchesBefore != mismatchesAfter):
                sys.stderr.write(alignment.Name + '\n');
                sys.stderr.write('ERROR: found %d mismatches before, but %d after left realignment!\n' % (mismatchesBefore, mismatchesAfter));
                exit(1);

        if (maxiterations <= 0):
            return False;
        else:
            return True;

   	return False;

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

        if (sam_line.IsMapped() == False):
            headers.append(line);
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
        sys.stderr.write('  %s <ref.fasta> <in.sam> <out.sam>\n' % (sys.argv[0]));
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
    run_bamleftalign(in_sam, in_ref);

    if (fp_out != sys.stdout):
        fp_out.close();

if __name__ == "__main__":
    main();




# ../../samscripts/src/convert_to_bam.sh msa3.handmade
# samtools faidx msa3.handmade.ref.fa
# cat msa3.handmade-sorted.bam | test-leftalign/bamleftalign/bamleftalign -d -f msa3.handmade.ref.fa > msa3.handmade-sorted.leftaligned.bam
# samtools view -h msa3.handmade-sorted.leftaligned.bam > msa3.handmade-sorted.leftaligned.sam
