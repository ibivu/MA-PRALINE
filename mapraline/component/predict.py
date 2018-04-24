import numpy as np

from mapraline import *
from mapraline.prosite import pattern_to_re

from praline.core import *
from praline.container import ProfileTrack, Sequence, PlainTrack
from praline.container import Alphabet, ALPHABET_AA, ALPHABET_DNA

_ALLOWED_ALPHABETS = [ALPHABET_AA.aid, ALPHABET_DNA.aid]

class PrositePatternAnnotator(Component):
    tid = "mapraline.component.PrositePatternAnnotator"

    inputs = {'sequence': Port(Sequence.tid),
              'pattern': Port(str),
              'track_id': Port(str)}
    outputs = {'prediction_track': Port(PlainTrack.tid)}

    options = {}
    defaults = {}

    def execute(self, sequence, pattern, track_id):
        track = sequence.get_track(track_id)
        if track.alphabet.aid not in _ALLOWED_ALPHABETS:
            s = "can only predict prosite pattern " \
                "matches for AA and DNA alphabets"
            raise ComponentError(s)

        if track.tid != PlainTrack.tid:
            s = "can only predict prosite pattern matches for plain tracks"
            raise ComponentError(s)

        re_pat = pattern_to_re(pattern)

        indices = track.values
        annotation_syms_list = []
        sym_list = []
        for n in xrange(indices.shape[0]):
            sym_list.append(track.alphabet.index_to_symbol(indices[n]))
        aa_syms = "".join(sym_list)

        for match in re_pat.finditer(aa_syms, overlapped=True):
            annotation_syms = [u'*'] * indices.shape[0]
            # First annotate the optional spacer groups from the match.
            try:
                match_starts_spacer = match.starts('spacer')
                match_ends_spacer = match.ends('spacer')
                for spacer_start, spacer_end in zip(match_starts_spacer, match_ends_spacer):
                    for n in xrange(spacer_start, spacer_end):
                        if annotation_syms[n] != 'M':
                            annotation_syms[n] = 'S'
            except IndexError:
                pass


            # Then annotate anything that wasn't annotated as a spacer but did match
            # as a match.
            for n in xrange(match.start(), match.end()):
                if annotation_syms[n] == u'*':
                    annotation_syms[n] = u'M'

            annotation_syms_list.append(annotation_syms)

        # Now collapse the individual annotation sequences into a single consensus sequence.
        consensus_annotation_syms = []
        for n in xrange(indices.shape[0]):
            consensus_sym = u'*'
            for m in xrange(len(annotation_syms_list)):
                if annotation_syms_list[m][n] == u'M':
                    consensus_sym = u'M'
                elif annotation_syms_list[m][n] == u'S' and consensus_sym in [u'*', u'S']:
                    consensus_sym = u'S'

            consensus_annotation_syms.append(consensus_sym)

        annotation_track = PlainTrack(consensus_annotation_syms, ALPHABET_PROSITE)
        outputs = {'prediction_track': annotation_track}
        yield CompleteMessage(outputs=outputs)
