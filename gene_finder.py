# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Sarah Barden

"""

import random
import math
from load import load_seq
dna = load_seq("./data/X73525.fa")

from amino_acids import aa, codons, aa_table   # you may find these useful


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))


# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###
def get_complement(nucleotide):
    """ Returns the complementary nucleotide
        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return 'C'


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    rev_comp = ''
    for i in range(0, len(dna)):
        nucleo = dna[i]
        comp = get_complement(nucleo)
        rev_comp = comp + rev_comp
    return rev_comp


def split_into_codons(dna):
    """ Takes a DNA sequence (a string) and splits it into a list of codons
        as a list of strings.

        dna: a DNA sequence
        returns: DNA sequence split into codons
    >>> split_into_codons("ATGTGATAG")
    ['ATG', 'TGA', 'TAG']
    >>> split_into_codons("ATGTGATAGCC")
    ['ATG', 'TGA', 'TAG', 'CC']
    """
    dna_split = []
    length = math.ceil(len(dna)/3)
    for i in range(0, length):
        j = 3*i
        codon = dna[j:j+3]
        dna_split += [codon]
    return dna_split


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon (TGA, TAA, or TAG).  If there is
        no in frame stop codon, returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    index = -1
    dna_split = split_into_codons(dna)

    for i, j in enumerate(dna_split):
        if j == "TAG" or j == "TGA" or j == "TAA":
            index = i
    if index == -1:
        return ''.join(dna_split)  # if there is no stop codon
    return ''.join(dna_split[:index])


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("ATGCCCATGTTTTAG")
    ['ATGCCCATGTTT']
    """
    dna_split = split_into_codons(dna)
    i = 0
    all_ORFs_oneframe = []
    length_of_orf = 0
    while i < len(dna_split):
        if dna_split[i] == "ATG":
            orf = rest_of_ORF(dna[int(i)*3:])
            all_ORFs_oneframe += [orf]
            length_of_orf = math.ceil(len(orf)/3)
            i += length_of_orf
        else:
            i += 1
    return all_ORFs_oneframe


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    all_ORFs = []
    for i in range(0, 3):
        dna_new = dna[i:]
        all_ORFs += find_all_ORFs_oneframe(dna_new)
    return all_ORFs


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    orfs = []
    orfs = find_all_ORFs(dna) + find_all_ORFs(get_reverse_complement(dna))
    return orfs


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    orfs = find_all_ORFs_both_strands(dna)
    longest = max(orfs, key=len)
    return longest


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF
    """
    i = 0
    longest_each_trial = []
    while i < num_trials:
        shuffled_dna = shuffle_string(dna)
        longest_each_trial += [longest_ORF(shuffled_dna)]
        i += 1

    longest_longest = max(longest_each_trial, key=len)
    return len(longest_longest)


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    dna_codons = split_into_codons(dna)
    i = 0
    aa_string = ''
    while i < len(dna_codons):
        if len(dna_codons[i]) == 3:
            amino_acid = aa_table[dna_codons[i]]
            aa_string += amino_acid
        i += 1
    return aa_string


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    orfs = find_all_ORFs_both_strands(dna)
    threshold = longest_ORF_noncoding(dna, 1500)
    print('threshold is', threshold)
    # print(orfs)
    print(len(orfs))
    aa_sequences = []
    i = 0
    while i < len(orfs):
        print('i is', i)
        print(len(orfs[i]))
        if len(orfs[i]) > threshold:
            print('if')
            aa_sequences += [coding_strand_to_AA(orfs[i])]
        i += 1
    print(aa_sequences)


if __name__ == "__main__":
    import doctest
    gene_finder(dna)
    # doctest.testmod(verbose=True)
    # doctest.run_docstring_examples(coding_strand_to_AA, globals())
