from __future__ import division, absolute_import, print_function

import sys
import argparse
import os
import uuid
import shutil
import tarfile
import os.path
import pickle

import six
import six.moves.urllib.parse

from mapraline import *
from mapraline.component import PrositePatternAnnotator

from praline import load_score_matrix, load_sequence_fasta, open_builtin
from praline import write_alignment_clustal, write_alignment_fasta
from praline.core import *
from praline.container import ALPHABET_AA, TRACK_ID_INPUT, TRACK_ID_PREPROFILE
from praline.container import ALPHABET_DNA, PlainTrack
from praline.component import GlobalMasterSlaveAligner, PairwiseAligner
from praline.component import LocalMasterSlaveAligner, ProfileBuilder
from praline.component import GuideTreeBuilder, TreeMultipleSequenceAligner
from praline.component import AdHocMultipleSequenceAligner
from praline.component import DummyMasterSlaveAligner
from praline.util import run, write_log_structure

ROOT_TAG = "__ROOT__"
_TRACK_ID_BASE_PATTERN = "mapraline.track.MotifAnnotationPattern"
_TRACK_ID_BASE_FILE = "mapraline.track.MotifAnnotationFile"

def main():
    # Parse arguments.
    args = parse_args()
    verbose = not args.quiet or args.verbose

    # See if we're doing DNA or protein alignments.
    # TODO: if unspecified, autodetect this based on input file contents?
    alphabet = None
    if args.input_dna:
        alphabet = ALPHABET_DNA
    elif args.input_protein:
        alphabet = ALPHABET_AA

    # Setup the execution manager.
    index = TypeIndex()
    index.autoregister()
    if args.num_threads > 1:
        manager = ParallelExecutionManager(index, args.num_threads)
    else:
        manager = Manager(index)

    seqs = load_sequence_fasta(args.input, alphabet)

    # Load inputs and other data.
    if args.score_matrix is not None:
        score_matrix_file = args.score_matrix
    else:
        if alphabet == ALPHABET_AA:
            score_matrix_file = 'blosum62'
        elif alphabet == ALPHABET_DNA:
            score_matrix_file = 'nucleotide'

    # Read score parameters.
    with open_resource(score_matrix_file, "matrices") as f:
        score_matrices = [load_score_matrix(f, alphabet=alphabet)]
    gap_series = [-float(x) for x in args.gap_penalties.split(",")]

    # Setup environment.
    keys = {}
    keys['gap_series'] = gap_series
    keys['score_threshold'] = args.preprofile_score
    keys['linkage_method'] = args.tree_linkage
    keys['waterman_eggert_iterations'] = args.num_preprofile_alignments
    keys['debug'] = args.debug
    if args.merge_semiglobal:
        keys['merge_mode'] = 'semiglobal'
        keys['dist_mode'] = 'semiglobal'
    else:
        keys['merge_mode'] = 'global'
        keys['dist_mode'] = 'global'

    if args.no_accelerate:
        keys['accelerate'] = False
    else:
        keys['accelerate'] = True
    env = Environment(keys=keys)

    # Initialize root node for output
    root_node = TaskNode(ROOT_TAG)

    # Annotate the motifs from the files and patterns.
    track_scores = do_motif_annotation(args, env, manager, seqs,
                                             verbose, root_node)
    # Build score matrices.
    motif_score_matrices = {}
    for trid, score in six.iteritems(track_scores):
        if score is None:
            score = args.motif_match_score

        motif_score_matrices[trid] = get_motif_score_matrix(score,
                                                            args.score_spacers)

    # Add all the new annotation tracks to the list of tracks to use
    # in the alignment.
    track_id_sets = [[TRACK_ID_INPUT]]
    for trid, track in seqs[0].tracks:
        if trid in motif_score_matrices:
            track_id_sets.append([trid])
            score_matrices.append(motif_score_matrices[trid])

    # Build initial sets of which sequences to align against every master
    # sequence. By default, we want to align every input sequence against
    # every other input sequence.
    master_slave_seqs = []
    all_seqs = list(seqs)
    for master_seq in seqs:
        slave_seqs = []
        for slave_seq in seqs:
            if slave_seq is not master_seq:
                slave_seqs.append(slave_seq)
        master_slave_seqs.append((master_seq, slave_seqs))

    master_slave_alignments = do_master_slave_alignments(args, env,
                                                         manager,
                                                         master_slave_seqs,
                                                         track_id_sets,
                                                         score_matrices,
                                                         verbose,
                                                         root_node)

    # Build preprofiles from master-slave alignments.
    do_preprofiles(args, env, manager, master_slave_alignments, seqs,
                   verbose, root_node)
    msa_track_id_sets = _replace_input_track_id(track_id_sets)

    # Do multiple sequence alignment from preprofile-annotated sequences.
    alignment = do_multiple_sequence_alignment(args, env, manager, seqs,
                                               msa_track_id_sets, score_matrices,
                                               verbose, root_node)

    # Write alignment to output file.
    outfmt = args.output_format
    if outfmt == 'fasta':
        write_alignment_fasta(args.output, alignment, TRACK_ID_INPUT)
    elif outfmt == "clustal":
        write_alignment_clustal(args.output, alignment, TRACK_ID_INPUT,
                                score_matrix)
    else:
        raise DataError("unknown output format: '{0}'".format(outfmt))

    # Dump pickled alignment object if user asked for it.
    if args.dump_alignment is not None:
        with open(args.dump_alignment, 'wb') as fo:
            pickle.dump(alignment, fo)

    if args.dump_all_tracks is not None:
        try:
            os.mkdir(args.dump_all_tracks)
        except OSError:
            pass

        all_trids = []
        for trid, track in alignment.items[0].tracks:
            if track.tid == PlainTrack.tid:
                all_trids.append(trid)

        for trid in all_trids:
            filename = "dump-{0}.aln".format(trid)
            path = os.path.join(args.dump_all_tracks, filename)

            if outfmt == "fasta":
                write_alignment_fasta(path, alignment, trid)
            elif outfmt == "clustal":
                write_alignment_clustal(path, alignment, trid, None)
            else:
                raise DataError("unknown output format: '{0}'".format(outfmt))

    if verbose:
        sys.stdout.write('\n')

    # Collect log bundles
    if args.debug > 0:
        write_log_structure(root_node)

def _replace_input_track_id(track_id_sets):
    new_track_id_sets = []
    for s in track_id_sets:
        new_s = []
        for tid in s:
            if tid == TRACK_ID_INPUT:
                new_s.append(TRACK_ID_PREPROFILE)
            else:
                new_s.append(tid)
        new_track_id_sets.append(new_s)
    return new_track_id_sets

def get_motif_score_matrix(match_score, score_spacers):
    all_symbols = {u"*", u"M", u"S"}

    d = {}
    for sym_one in all_symbols:
        for sym_two in all_symbols:
            if sym_one == u"M" and sym_two == u"M":
                score = match_score
            elif score_spacers and (sym_one == u"S" and sym_two == u"S"):
                score = match_score
            else:
                score = 0.0
            d[(sym_one, sym_two)] = score

    return ScoreMatrix(d, [ALPHABET_PROSITE, ALPHABET_PROSITE])


def do_motif_annotation(args, env, manager, seqs, verbose, root_node):
    FMT_TRACK_ID = "{0}_{1}"

    track_scores = {}

    execution = Execution(manager, ROOT_TAG)
    seq_patterns = []
    for seq in seqs:
        for pair in args.patterns:
            pattern = pair[0]
            if len(pair) > 1:
                score = float(pair[1])
            else:
                score = None
            seq_patterns.append((seq, pattern, score))

            component = PrositePatternAnnotator
            task = execution.add_task(component)
            task.environment(env)
            task.inputs(sequence=seq, pattern=pattern,
                        track_id=TRACK_ID_INPUT)

    outputs = run(execution, verbose=verbose, root_node=root_node)
    for n, output in enumerate(outputs):
        seq, pattern, score = seq_patterns[n]

        track = output['prediction_track']

        trid = FMT_TRACK_ID.format(_TRACK_ID_BASE_PATTERN, pattern)
        seq.add_track(trid, track)
        track_scores[trid] = score

    for pair in args.annotation_files:
        annotation_file = pair[0]
        if len(pair) > 1:
            score = float(pair[1])
        else:
            score = None

        annotation_seqs = load_sequence_fasta(annotation_file,
                                              ALPHABET_PROSITE)
        name_tracks = {}
        for annotation_seq in annotation_seqs:
            track = annotation_seq.get_track(TRACK_ID_INPUT)
            name_tracks[annotation_seq.name] = track

        for seq in seqs:
            track = name_tracks[seq.name]

            trid = FMT_TRACK_ID.format(_TRACK_ID_BASE_FILE, annotation_file)
            seq.add_track(trid, track)
            track_scores[trid] = score

    return track_scores

def do_master_slave_alignments(args, env, manager, seqs,
                               track_id_sets, score_matrices, verbose,
                               root_node):
    execution = Execution(manager, ROOT_TAG)

    master_slave_alignments = [None for seq in seqs]
    for master_seq, slave_seqs in seqs:
        if args.preprofile_global:
            component = GlobalMasterSlaveAligner
        elif args.preprofile_local:
            component = LocalMasterSlaveAligner
        else:
            component = DummyMasterSlaveAligner

        task = execution.add_task(component)
        task.environment(env)
        task.inputs(master_sequence=master_seq, slave_sequences=slave_seqs,
                    track_id_sets=track_id_sets, score_matrices=score_matrices)

    outputs = run(execution, verbose=verbose, root_node=root_node)
    for n, output in enumerate(outputs):
        master_slave_alignments[n] = output['alignment']

    return master_slave_alignments


def do_multiple_sequence_alignment(args, env, manager, seqs,
                                   track_id_sets, score_matrices,
                                   verbose, root_node):
    if args.pregen_tree:
        # Dummy preprofiles, so we can safely align by sequence.
        sub_env = Environment(parent=env)

        if not args.preprofile_local and not args.preprofile_global:
            sub_env.keys['squash_profiles'] = True

        # Build guide tree
        component = GuideTreeBuilder
        execution = Execution(manager, ROOT_TAG)
        task = execution.add_task(component)
        task.environment(sub_env)
        task.inputs(sequences=seqs, track_id_sets=track_id_sets,
                    score_matrices=score_matrices)

        outputs = run(execution, verbose=verbose, root_node=root_node)[0]

        # Build MSA
        component = TreeMultipleSequenceAligner
        execution = Execution(manager, ROOT_TAG)
        task = execution.add_task(component)
        task.environment(env)
        task.inputs(sequences=seqs, guide_tree=outputs['guide_tree'],
                    track_id_sets=track_id_sets, score_matrices=score_matrices)

        outputs = run(execution, verbose=verbose, root_node=root_node)[0]
    else:
        component = AdHocMultipleSequenceAligner
        execution = Execution(manager, ROOT_TAG)
        task = execution.add_task(component)
        task.environment(env)
        task.inputs(sequences=seqs, track_id_sets=track_id_sets,
                    score_matrices=score_matrices)

        outputs = run(execution, verbose=verbose, root_node=root_node)[0]

    return outputs['alignment']


def do_preprofiles(args, env, manager, alignments, seqs, verbose, root_node):
    for i, alignment in enumerate(alignments):
        component = ProfileBuilder
        execution = Execution(manager, ROOT_TAG)
        task = execution.add_task(component)
        task.environment(env)
        task.inputs(alignment=alignment, track_id=TRACK_ID_INPUT)

        outputs = run(execution, verbose=verbose, root_node=root_node)[0]
        track = outputs['profile_track']
        seqs[i].add_track(TRACK_ID_PREPROFILE, track)

def pair(arg):
    return arg.split(":")

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="input file in FASTA format")
    parser.add_argument("output", help="output alignment")
    parser.add_argument("-g", "--gap-penalties",
                        help="comma separated list of positive gap penaties",
                        default="11,1", dest="gap_penalties")
    parser.add_argument("--motif-match-score",
                        default=1.0, dest="motif_match_score", type=float,
                        help="on matching motif status, boost score by " \
                             "this amount")
    parser.add_argument("--spacers-as-match",
                        default=False, action="store_true",  dest="score_spacers",
                        help="treat spacers as real matches for the scoring")
    parser.add_argument("-m", "--score-matrix",
                        help="score matrix to use for alignment",
                        default=None, dest="score_matrix")
    parser.add_argument("-t", "--threads", help="number of threads to use",
                        default=1,
                        dest="num_threads", type=int)
    parser.add_argument("-s", "--preprofile-score", default=None,
                        dest="preprofile_score", type=float,
                        help="exclude preprofile alignments by score")
    parser.add_argument("-f", "--output-format", default="fasta",
                        dest="output_format",
                        help="write the alignment in the specified format")
    parser.add_argument("--tree-linkage", default="average",
                        dest="tree_linkage",
                        help="use this linkage method when building tree")
    parser.add_argument("--no-accelerate", default=False,
                        dest="no_accelerate", action="store_true",
                        help="disable the fast pairwise aligner")
    parser.add_argument("--preprofile-alignments", default=2,
                        dest="num_preprofile_alignments", type=int,
                        help="local preprofile alignments per sequence")
    parser.add_argument("--debug", "-d", action="count", dest="debug",
                        default=0, help="enable debugging output")
    parser.add_argument("--dump-alignment-obj", type=str, default=None,
                        dest="dump_alignment",
                        help="dump final alignment object to filename")
    parser.add_argument("--dump-all-tracks", type=str, default=None,
                        dest="dump_all_tracks",
                        help="write alignment files for all the tracks")

    parser.add_argument('-p', '--pattern', action='append', dest="patterns",
                    type=pair, default=[],
                    help="annotate this prosite pattern in the sequence " \
                         "(specify a number after a colon to override " \
                         "the global match boost score)")
    parser.add_argument('-a', '--annotation-file', action='append',
                    type=pair, dest="annotation_files", default=[],
                    help="read motif annotation tracks from a FASTA file " \
                         "(specify a number after a colon to override " \
                         "the global match boost score)")

    group = parser.add_mutually_exclusive_group()
    group.add_argument("--input-dna",
                       help="input contains DNA sequences",
                       action="store_true", dest="input_dna",
                       default=False)
    group.add_argument("--input-protein",
                       help="input contains protein sequences",
                       action="store_true", dest="input_protein",
                       default=True)

    group = parser.add_mutually_exclusive_group()
    group.add_argument("--preprofile-none", help="build no preprofiles",
                       action="store_true", dest="preprofile_none",
                       default=False)
    group.add_argument("--preprofile-local", help="build local preprofiles",
                       action="store_true", dest="preprofile_local",
                       default=False)
    group.add_argument("--preprofile-global", help="build global preprofiles",
                       action="store_true", dest="preprofile_global",
                       default=False)

    group = parser.add_mutually_exclusive_group()
    group.add_argument('--msa-tree', dest="pregen_tree",
                        action="store_true", default=True,
                        help="do pre-generated tree mutliple alignment")
    group.add_argument('--msa-adhoc',  dest="pregen_tree",
                        action="store_false", default=True,
                        help="do ad-hoc tree multiple alignment")

    group = parser.add_mutually_exclusive_group()
    group.add_argument('--merge-global', dest="merge_global",
                        action="store_true", default=True,
                        help="merge with a global alignment")
    group.add_argument('--merge-semiglobal',  dest="merge_semiglobal",
                        action="store_true", default=False,
                        help="merge with a semiglobal alignment")

    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v", "--verbose", dest="verbose", help="be verbose",
                        action="store_true", default=False)
    group.add_argument("-q", "--quiet", dest="quiet", help="be quiet",
                       action="store_true", default=True)

    return parser.parse_args()

def open_resource(filename, prefix):
    try:
        return open(filename)
    except IOError:
        return open_builtin(os.path.join(prefix, filename))
