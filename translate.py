#! /usr/bin/env python3

import sys

def translate_sequence(rna_sequence, genetic_code):
    """Translates a sequence of RNA into a sequence of amino acids.
    Translates `rna_sequence` into string of amino acids, according to the
    `genetic_code` given as a dict. Translation begins at the first position of
    the `rna_sequence` and continues until the first stop codon is encountered
    or the end of `rna_sequence` is reached.
    If `rna_sequence` is less than 3 bases long, or starts with a stop codon,
    an empty string is returned.
    """
    rna_seq = rna_sequence.upper()
    start = 0
    proteins = ''
    for i in range(start, len(rna_seq), 3):
        codon = rna_seq[i:i + 3]
        if codon in ['UAG', 'UAA', 'UGA'] or len(codon) != 3:
            break
        else: proteins += genetic_code[codon]
    return proteins
    pass

def get_all_translations(rna_sequence, genetic_code):
    """Get a list of all amino acid sequences encoded by an RNA sequence.

    All three reading frames of `rna_sequence` are scanned from 'left' to
    'right', and the generation of a sequence of amino acids is started
    whenever the start codon 'AUG' is found. The `rna_sequence` is assumed to
    be in the correct orientation (i.e., no reverse and/or complement of the
    sequence is explored).

    The function returns a list of all possible amino acid sequences that
    are encoded by `rna_sequence`.

    If no amino acids can be translated from `rna_sequence`, an empty list is
    returned.
    """
    rna_seq = rna_sequence.upper()
    amino_seq = []
    start = 0

    def translate(start,rna_seq,genetic_code):
        proteins = ''
        for i in range(start, len(rna_seq), 3):
            codon = rna_seq[i:i + 3]
            if codon in ['UAG', 'UAA', 'UGA'] or len(codon) != 3:
                break
            else: proteins += genetic_code[codon]
        return proteins

    while start < len(rna_seq):
        start_codon = rna_seq[start:start + 3]
        if start_codon == 'AUG':
            translation = translate(start, rna_seq, genetic_code)
            amino_seq.append(translation)
        start += 1
    return amino_seq
    pass

def get_reverse(sequence):
    """Reverse orientation of `sequence`.
    Returns a string with `sequence` in the reverse order.
    If `sequence` is empty, an empty string is returned.
    """
    return ''.join(reversed(sequence.upper()))
    pass

def get_complement(sequence):
    """Get the complement of `sequence`.
    Returns a string with the complementary sequence of `sequence`.
    If `sequence` is empty, an empty string is returned.
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'A': 'U', 'U': 'A'}
    return ''.join(complement.get(base, base) for base in sequence.upper())
    pass

def reverse_and_complement(sequence):
    """Get the reversed and complemented form of `sequence`.
    Returns a string that is the reversed and complemented sequence
    of `sequence`.
    If `sequence` is empty, an empty string is returned.
    """
    rev_comp = get_reverse(get_complement(sequence))
    return rev_comp
    pass

def get_longest_peptide(rna_sequence, genetic_code):
    """Get the longest peptide encoded by an RNA sequence.

    Explore six reading frames of `rna_sequence` (the three reading frames of
    `rna_sequence`, and the three reading frames of the reverse and complement
    of `rna_sequence`) and return (as a string) the longest sequence of amino
    acids that it encodes, according to the `genetic_code`.

    If no amino acids can be translated from `rna_sequence` nor its reverse and
    complement, an empty string is returned.
    """
    peptides = get_all_translations(rna_sequence = rna_sequence, genetic_code = genetic_code)
    reverse_comp = reverse_and_complement(rna_sequence)
    reverse_tran = get_all_translations(rna_sequence = reverse_comp, genetic_code = genetic_code)
    peptides += reverse_tran
    if not peptides:
        return ""
    if len(peptides) < 2:
        return peptides[0]
    mst_bses = -1
    lngst_pep = -1
    for peptide, aa_seq in enumerate(peptides):
        if len(aa_seq) > mst_bses:
            lngst_pep = peptide
            mst_bses = len(aa_seq)
    return peptides[lngst_pep]

if __name__ == '__main__':
    genetic_code = {'GUC': 'V', 'ACC': 'T', 'GUA': 'V', 'GUG': 'V', 'ACU': 'T', 'AAC': 'N', 'CCU': 'P', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'CAU': 'H', 'AAU': 'N', 'AGU': 'S', 'GUU': 'V', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'CCA': 'P', 'ACA': 'T', 'CCC': 'P', 'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'UGC': 'C', 'CAG': 'Q', 'GAU': 'D', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'AGG': 'R', 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': '*', 'GGA': 'G', 'UAC': 'Y', 'GAC': 'D', 'UAG': '*', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'GGC': 'G', 'AUG': 'M', 'CUG': 'L', 'GAG': 'E', 'CUC': 'L', 'AGA': 'R', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'GCU': 'A', 'GAA': 'E', 'AUU': 'I', 'UUG': 'L', 'UUA': 'L', 'UGA': '*', 'UUC': 'F'}
    rna_seq = ("AUG"
            "UAC"
            "UGG"
            "CAC"
            "GCU"
            "ACU"
            "GCU"
            "CCA"
            "UAU"
            "ACU"
            "CAC"
            "CAG"
            "AAU"
            "AUC"
            "AGU"
            "ACA"
            "GCG")
    longest_peptide = get_longest_peptide(rna_sequence = rna_seq,
            genetic_code = genetic_code)
    assert isinstance(longest_peptide, str), "Oops: the longest peptide is {0}, not a string".format(longest_peptide)
    message = "The longest peptide encoded by\n\t'{0}'\nis\n\t'{1}'\n".format(
            rna_seq,
            longest_peptide)
    sys.stdout.write(message)
    if longest_peptide == "MYWHATAPYTHQNISTA":
        sys.stdout.write("Indeed.\n")
