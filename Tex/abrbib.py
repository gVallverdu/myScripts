#!/usr/bin/env python

import re
import sys
import os

CASSI= {"physical review b": "Phys. Rev. B",
        "journal of applied crystallography": "J. Appl. Crystallogr.",
        "journal of the electrochemical society": "J. Electrochem. Soc.",
        "the journal of chemical physics": "J. Chem. Phys.",
        "journal of chemical physics": "J. Chem. Phys.",
        "solid state ionics": "Solid State Ionics",
        "computational materials science": "Comput. Mater. Sci.",
        "journal of materials chemistry": "J. Mater. Chem.",
        "vibrational spectroscopy": "Vib. Spectrosc.",
        "journal of power sources": "J. Power Sources",
        "physical review letters": "Phys. Rev. Lett.",
        "acta crystallographica section c": "Acta Crystallogr., Sect. C",
        "crystallography reports": "Crystallogr. Rep.",
        "philosophical magazine part b": "Philos. Mag. B",
        "journal of non-crystalline solids": "J. Non-Cryst. Solids",
        "journal of solid state chemistry": "J. Solid State Chem.",
        "solid state communications": "Solid State Commun.",
        "advanced energy materials": "Adv. Energy Mater.",
#        "Zeitschrift für Physik": "Z. Phys.",
        "Zeitschrift fuer Physik": "Z. Phys.",
        "surface science": "Surf. Sci.",
        "catalysis today": "Catal. Today"
        }

def cassi_abr(fname):
    """
    Modify the journal key of all bib entries in fname and replace it by its
    corresponding abreviation
    """
    j_patt = re.compile("^\s*journal\s*=\s*\{\s*(.+)\s*\}")

    with open(fname, "r") as f:
        newbib = ""
        for line in f:
            m = j_patt.match(line)
            if m:
                journal = m.group(1)
                newbib += "    journal = {{{0}}},\n".format(find_abbreviation(journal))
            else:
                newbib += line

    base, ext = os.path.splitext(fname)
    with open(base + "_abr" + ext, "w") as f:
        f.write(newbib)

def find_abbreviation(journal):
    """ 
    Look for abreviation of journal in CASSI conventions 
    
    Args:
        journal: full name of the journal
        
    Returns:
        abrreviation: abreviation of the journal or the journal if the 
                      abreviation was not found.
    """
    if journal.lower() in CASSI:
        return CASSI[journal.lower()]
    else:
        print("abreviation of {0} not found".format(journal))
        suggest = [j for j in CASSI if dameraulevenshtein(j, journal.lower()) > .5]
        if suggest:
            print("Similar Journals :")
            for j in suggest:
                print("    * {0}: {1}".format(j, CASSI[j]))
        return journal

def dameraulevenshtein(seq1, seq2):
    """Calculate the Damerau-Levenshtein distance between sequences.

    This distance is the number of additions, deletions, substitutions,
    and transpositions needed to transform the first sequence into the
    second. Although generally used with strings, any sequences of
    comparable objects will work.

    Transpositions are exchanges of *consecutive* characters; all other
    operations are self-explanatory.

    This implementation is O(N*M) time and O(M) space, for N and M the
    lengths of the two sequences.
    """
    # http://blog.developpez.com/philben/p11268/vba-access/similarite_entre_deux_chaines_de_caracte
    # http://mwh.geek.nz/2009/04/26/python-damerau-levenshtein-distance/ 
    # Conceptually, this is based on a len(seq1) + 1 * len(seq2) + 1 matrix.
    # However, only the current and two previous rows are needed at once,
    # so we only store those.
    oneago = None
    thisrow = list(range(1, len(seq2) + 1)) + [0]
    for x in range(len(seq1)):
        # Python lists wrap around for negative indices, so put the
        # leftmost column at the *end* of the list. This matches with
        # the zero-indexed strings and saves extra calculation.
        twoago, oneago, thisrow = oneago, thisrow, [0] * len(seq2) + [x + 1]
        for y in range(len(seq2)):
            delcost = oneago[y] + 1
            addcost = thisrow[y - 1] + 1
            subcost = oneago[y - 1] + (seq1[x] != seq2[y])
            thisrow[y] = min(delcost, addcost, subcost)
            # This block deals with transpositions
            if (x > 0 and y > 0 and seq1[x] == seq2[y - 1]
                and seq1[x-1] == seq2[y] and seq1[x] != seq2[y]):
                thisrow[y] = min(thisrow[y], twoago[y - 2] + 1)

    dls = 1 - thisrow[len(seq2) - 1] / max(len(seq1), len(seq2))
    return dls

if __name__ == "__main__":
    #print("barette ", dameraulevenshtein("tbalrette", "barette"))
    #print("tablette ", dameraulevenshtein("tbalrette", "tablette"))
    #print("turbulette ", dameraulevenshtein("tbalrette", "turbulette"))    
    cassi_abr(sys.argv[1])

