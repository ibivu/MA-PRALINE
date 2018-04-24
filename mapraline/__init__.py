import codecs
import csv

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
    return zip(symbols, range(len(symbols)))

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
    if not isinstance(trids, list):
        trids = [trids]
    if not isinstance(descriptions, list):
        descriptions = [descriptions]

    match_idx = ALPHABET_PROSITE.symbol_to_index(u"M")
    spacer_idx = ALPHABET_PROSITE.symbol_to_index(u"S")

    rows = []
    rows.append(["match", match_color])
    rows.append(["spacer", spacer_color])

    for i, seq in enumerate(sequences):
        name = seq.name.split()[0]
        for trid, description in zip(trids, descriptions):
            track = seq.get_track(trid)
            if track.tid != PlainTrack.tid:
                s = "can only write plain motif tracks to a Jalview " \
                    "formatted file"
                raise DataError(s)

            if track.alphabet.aid != ALPHABET_PROSITE.aid:
                s = "need a motif annotation alphabet to write to Jalview " \
                    "format, but got an alphabet of type '{0}'"
                raise DataError(s.format(track.alphabet.aid))

            for j, symbol_idx in enumerate(track.values):
                if symbol_idx == match_idx:
                    type_ = "match"
                elif symbol_idx == spacer_idx:
                    type_ = "spacer"
                else:
                    continue

                row = [description, name, "-1", str(j + 1), str(j + 1),
                       type_]
                rows.append(row)

    lines = ["\t".join(row) for row in rows]

    should_close = False
    if isinstance(f, str):
        f = codecs.open(f, 'w', 'ascii')
        should_close = True

    f.write("\n".join(lines))
    f.write("\n")

    if should_close:
        f.close()
