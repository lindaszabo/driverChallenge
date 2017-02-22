# The input to the script is at most 50 DNA sequences (i.e, the character
# set is limited to T/C/G/A) whose length does not exceed 1000 characters.
# The sequences are given in FASTA format. These sequences are all different
# fragments of one chromosome.
# The specific set of sequences you will get satisfy a very unique property:
# there exists a unique way to reconstruct the entire chromosome from these
# reads by gluing together pairs of reads that overlap by more than half their
# length.

# usage: parseFragments.py -f input.fa

# output: unique sequence that contains each input string as a substring

import sys
import os
import re
import logging
from Bio import SeqIO
import argparse
from collections import deque
from difflib import SequenceMatcher


class fragObj:
    def __init__(self, frag):
        self.frag = frag  # current extension of fragment string
        self.leftFragOL = len(frag)/2 + 1  # length overlap req to extend left
        self.rightFragOL = len(frag)/2 + 1  # len overlap req to extend right

    # Paste new seq onto left or right side of self.frag if enough overlap
    # param ext: fragObj containing sequence to attempt to glue to self.frag
    # return: True if able to paste ext on either side, else False
    def extendWithOverlap(self, ext):

        match = SequenceMatcher(None, self.frag, ext.frag).find_longest_match(
            0, len(self.frag), 0, len(ext.frag))

        # match.a: start pos of match in self.frag
        # match.b: start pos of match in ext.frag
        # match.size: length of match

        if (match.a == 0 and  # match covers start of my fragment
            match.size >= self.leftFragOL and  # match > half my left frag
            match.size >= ext.rightFragOL and  # match > half new right frag
            (match.b + match.size == len(ext.frag) or
             match.a + match.size == len(self.frag))):  # match end of a frag
            # update overlap needed for another seq to extend on the left
            self.leftFragOL = len(ext.frag)/2 + 1

            log.debug("will paste on LHS \nself: %s\next:%s" %
                      (self.frag, ext.frag))

            # new frag is part of extension before match, match, trailing seq
            # trailing seq may be from either sequence, <= 1 seq is non-empty
            self.frag = "".join([ext.frag[:match.b],
                                 self.frag[match.a: match.a + match.size],
                                 ext.frag[match.b+match.size:],
                                 self.frag[match.a+match.size:]])

            log.debug("new self.frag: %s" % (self.frag))

            return True

        elif (match.b == 0 and  # match covers start of extension
              match.size >= self.rightFragOL and  # match > half my right frag
              match.size >= ext.leftFragOL and  # match > half new left frag
              (match.b + match.size == len(ext.frag) or
               match.a + match.size == len(self.frag))):  # match end of a frag

            # update overlap needed for another seq to extend on the right
            self.rightFragOL = len(ext.frag)/2 + 1

            log.debug("will paste on RHS \nself: %s\next:%s" %
                      (self.frag, ext.frag))

            # new frag is part of extension before match, match, trailing seq
            # trailing seq may be from either sequence, <= 1 seq is non-empty
            self.frag = "".join([self.frag[:match.a],
                                 self.frag[match.a: match.a + match.size],
                                 ext.frag[match.b+match.size:],
                                 self.frag[match.a+match.size:]])

            log.debug("new self.frag: %s" % (self.frag))

            return True

        else:
            return False

    def __str__(self):
        return "; ".join([self.frag, str(self.leftFragOL),
                          str(self.rightFragOL)])


# Construct realignment from partial fragments.
# param cur_frags: deque containing partially constructed fragments
# param breakAfter: int max number of iterations before exiting while loop
# return: unique sequence that contains each fragment as a substring
def getChromosomeSeq(cur_frags, breakAfter):

    baseFrag = cur_frags.pop()

    cnt_loop = 1  # counter to bail if can't glue all pieces together

    while len(cur_frags) > 0 and cnt_loop < breakAfter:
        for frag in cur_frags:
            cnt_loop = cnt_loop + 1
            if baseFrag.extendWithOverlap(frag):
                # able to add it to the base fragment, move on
                cur_frags.remove(frag)
                break

    if cnt_loop == args.maxSeqCount:
        logging.error("Input sequences not reconstructable.")
        return "Input sequences not reconstructable."
    else:
        return baseFrag.frag


patt_DNA = re.compile("^[ACTG]+$")  # to check for unexpected input in fasta

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fastaFile',
                        help='fasta of overlapping chromosome fragments',
                        required=True)
    parser.add_argument('-c', '--maxSeqCount', type=int, default=50,
                        help='max number of sequences expected in input file')
    parser.add_argument('-l', '--maxSeqLen', type=int, default=1000,
                        help='max length of sequence allowed in input file')
    args = parser.parse_args()

    partial_frags = deque([])  # for efficient removal of fragments at end

    in_handle = None  # so can check close in finally block

    # set up logging
    logfile = "pasteFragments_%s.log" % (os.path.basename(args.fastaFile))
    logging.basicConfig(filename=logfile, level=logging.DEBUG, filemode='w')
    log = logging.getLogger()

    try:
        in_handle = open(args.fastaFile, "rU")
        fasta_sequences = SeqIO.parse(in_handle, "fasta")

        for fasta in fasta_sequences:
            # check a few basic assumptions, report errors to user
            if not patt_DNA.search(str(fasta.seq)):
                logging.error("INPUT ERROR. Non-DNA sequence: %s", (fasta.seq))
                sys.exit("INPUT ERROR. Non-DNA sequence: " + fasta.seq)
            if len(fasta.seq) > args.maxSeqLen:
                logging.error("INPUT ERROR. Sequence longer than: %d",
                              (args.maxSeqLen))
                sys.exit("INPUT ERROR. Sequence longer than: " +
                         str(args.maxSeqLen))

            newFrag = fragObj(str(fasta.seq))  # store seq info in obj

            # flag will be updated if able to extend existing frag with new
            didExtend = False

            # check each prev. fragment to see if we can extend with the new
            for frag in partial_frags:
                if frag.extendWithOverlap(newFrag):
                    # newFrag was added to an existing fragment
                    didExtend = True
                    break

            if not didExtend:
                # could not extend existing frag, add as a new partial fragment
                partial_frags.append(newFrag)
        in_handle.close()

        # construct realignment from partial fragments
        print getChromosomeSeq(partial_frags, args.maxSeqCount)

    except Exception as e:
        logging.error("%s %s", (e, sys.exc_info()[0]))
        print "Exception", e
        print "error:", sys.exc_info()[0]
        sys.exit(2)

    finally:
        if in_handle:
            in_handle.close()
