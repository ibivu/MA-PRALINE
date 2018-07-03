from __future__ import division, absolute_import, print_function

import sys
import re
import pickle
import itertools
import os
import os.path

import numpy as np
from six.moves import range

from praline import load_score_matrix, window
from praline.core import *
from praline.container import Alignment, ALPHABET_AA, TRACK_ID_INPUT
from praline.container import PlainTrack

from mapraline import *
from mapraline.cmd import get_prosite_score_matrix, open_resource

_TRACK_ID_BASE = "mapraline.track.PrositePatternAnnotation"

def main():
    path = sys.argv[1]
    aa_score_matrix_path = sys.argv[2]
    gap_penalties = [-float(penalty) for penalty in sys.argv[3].split(",")]

    # Load inputs and other data.
    with open_resource(aa_score_matrix_path, "matrices") as f:
        aa_score_matrix = load_score_matrix(f, alphabet=ALPHABET_AA)

    aa_scores, motif_scores = get_aa_motif_scores(path, aa_score_matrix,
                                                  gap_penalties, 1.0, 0.0)

    #print aa_scores[0]
    print(sum(aa_scores[0]))
    print()
    #print motif_scores
    print(sum(sum(track_scores) for track_scores in motif_scores))

def get_aa_motif_scores(alignment_path, aa_score_matrix, gap_penalties,
                        motif_match_score, motif_mismatch_score):
    motif_score_matrix = get_prosite_score_matrix(motif_match_score,
                                                  motif_mismatch_score)

    with open(alignment_path, 'rb') as fi:
        alignment = pickle.load(fi)

    track_ids = []

    seq = alignment.items[0]
    for trid, tracks in seq.tracks:
        if trid.startswith(_TRACK_ID_BASE):
            track_ids.append(trid)

    aa_scores = calc_sum_of_pairs(alignment, [TRACK_ID_INPUT], aa_score_matrix, gap_penalties)
    motif_scores = calc_sum_of_pairs(alignment, track_ids, motif_score_matrix, [0.0])

    return (aa_scores, motif_scores)


def calc_sum_of_pairs(alignment, track_ids, score_matrix, gap_penalties):
    track_scores = []
    for track_idx, trid in enumerate(track_ids):
        print(track_idx + 1, len(track_ids))

        alphabet = None
        tracks = []
        for sequence in alignment.items:
            track = sequence.get_track(trid)
            if track.tid != PlainTrack.tid:
                s = "can only calculate sum of pairs score for plain tracks"
                raise DataError(s)
            alphabet = track.alphabet
            tracks.append(track)

        path = alignment.path
        positions = []
        for i, i_next in window(list(range(path.shape[0]))):
            inc_cols = (path[i_next, :]-path[i, :]) > 0
            position = []
            for j, inc_col in enumerate(inc_cols):
                if path[i_next, j] == (-1):
                    s = "the sum of pairs calculation does not currently " \
                        "support local alignments"
                    raise DataError(s)

                if inc_col:
                    seq_idx = path[i_next, j]
                    symbol = alphabet.index_to_symbol(tracks[j].values[seq_idx-1])
                    position.append(symbol)
                else:
                    position.append(None)

            positions.append(position)

        gap_lengths = precalc_gap_lengths(positions)

        scores = []
        for pos_idx, position in enumerate(positions):
            score = 0.0
            itr = itertools.combinations_with_replacement(enumerate(position), 2)
            for (idx_a, sym_a), (idx_b, sym_b) in itr:
                if idx_a == idx_b:
                    continue

                # Need to make this a bit faster, so cache the gap penalty
                # for a position instead of constantly looking it up.
                cached_gap_penalties = {}

                if sym_a is None:
                    if not idx_a in cached_gap_penalties:
                        gap_len = gap_lengths[pos_idx, idx_a]
                        cached_gap_penalties[idx_a] = calc_gap_penalty(gap_penalties, gap_len)
                    score += cached_gap_penalties[idx_a]
                elif sym_b is None:
                    if not idx_b in cached_gap_penalties:
                        gap_len = gap_lengths[pos_idx, idx_b]
                        cached_gap_penalties[idx_b] = calc_gap_penalty(gap_penalties, gap_len)
                    score += cached_gap_penalties[idx_b]
                else:
                    score += score_matrix.score((sym_a, sym_b))
            scores.append(score)

        track_scores.append(scores)

    return track_scores

def precalc_gap_lengths(positions):
    gap_lengths = np.empty((len(positions), len(positions[0])), dtype=np.uint32)

    for n in range(gap_lengths.shape[0]):
        for m in range(gap_lengths.shape[1]):
            if n == 0:
                if positions[n][m] is None:
                    gap_lengths[n, m] = 1
                else:
                    gap_lengths[n, m] = 0
            else:
                if positions[n][m] is None:
                    gap_lengths[n, m] = gap_lengths[n-1, m] + 1
                else:
                    gap_lengths[n, m] = 0

    return gap_lengths

def calc_gap_penalty(penalties, length):
    try:
        return penalties[length - 1]
    except IndexError:
        return penalties[-1]

if __name__ == '__main__':
    main()
