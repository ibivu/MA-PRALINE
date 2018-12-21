from __future__ import division, absolute_import, print_function

import codecs
import csv

from six.moves import range
from six.moves import zip

from praline import write_features_jalview
from praline.container import Alphabet, ScoreMatrix, PlainTrack
from praline.core import *


class MAPRALINEError(Exception):
    pass

class PrositeParserError(Exception):
    pass

def _generate_mappings(symbols):
    """Helper method to generate a mapping between symbols and indices.

    :param symbols: list containing the symbols to map
    :returns: a list of 2-tuples containing the generated mapping between
        symbols and indices
    """
    return list(zip(symbols, list(range(len(symbols)))))

_SYM_PROSITE = [u'*', u'M', u'S']
_MAP_PROSITE = _generate_mappings(_SYM_PROSITE)
ALPHABET_PROSITE = Alphabet('mapraline.alphabet.PrositePatternAnnotation', _MAP_PROSITE)

def write_motifs_jalview(f, sequences, trids, descriptions, match_color,
                         spacer_color):
    """Write motif annotations of a set of sequences in Jalview's sequence
    feature format.

    :param f: a filename or file object to write the sequences to
    :param sequences: the collection of sequences to write
    :param trids: the track id or a list of track ids containing the motif
        track(s) to write
    :param descriptions: the description of the sequence features or a list of
        descriptions of the sequence features, one per track
    :param match_color: the color of motif matches
    :param spacer_color: the color of motif spacers

    """
    colors = {u"M": match_color, u"S": spacer_color}
    if isinstance(trids, list):
        colors = [colors] * len(trids)

    write_features_jalview(f, sequences, trids, descriptions, colors)
