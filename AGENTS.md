# AGENTS.md

This file provides guidance to LLM agents when working with code in
this repository.

## Project Overview

JAABA (the Janelia Automatic Animal Behavior Annotator) is a
machine-learning system for training classifiers that detect discrete
animal behaviors (e.g. walking, grooming, following) from video.  A
user labels some frames in an interactive GUI, and JAABA trains a
boosting-based classifier that then annotates behavior across large
video sets.

It is a MATLAB codebase with a GUI frontend.

## Inputs and data model

A JAABA classifier operates on per-frame features derived from
tracking data for one or more animals:

- **`.trx` / `trx.mat`**: a centroid-and-heading track (x, y, theta,
  a, b, ...) for each animal, the classic JAABA input.
- **`.trk`**: keypoint ("landmark") tracks for each animal, produced by
  APT.  JAABA can generate a `trx` on the fly from a `.trk` (see APT
  integration below).
- **Movie**: read via `filehandling/get_readframe_fcn`, which supports
  many formats (`.ufmf`, `.fmf`, `.sbfmf`, `.avi`, `.mov`, `.mp4`,
  `.m4v`, `.seq`, `.mjpg`, sequential `.tif`/`.png`, `.h5`, `.klb`).

An **experiment directory** holds one experiment's files under names
fixed by the project config: the movie (e.g. `movie.ufmf`), the trx
(e.g. `trx.mat`), a `perframe/` directory of per-frame feature files, a
`clips/` directory, optional `scores_<behavior>.mat`, and (for APT
projects) the `.trk` file(s).

A **project** is stored in a `.jab` file: a `.mat` file holding one
variable `x`, which is a `Macguffin` object
(`perframe/Macguffin.m`).  It records the feature lexicon, behavior
name(s), per-experiment file names, window-feature params, labels, and
the trained classifier.

## Architecture

- **Model**: `JLabelData` (`perframe/JLabelData.m`) holds the
  domain state (experiments, labels, features, classifier) and is
  usable without a GUI for batch work.
- **GUI**: `JLabel` (`perframe/JLabel.m`) plus `JLabelGUIData`.  The
  code predates a strict model-view-presenter split.  When practical,
  keep model logic in `JLabelData` and presentation logic in the GUI
  layer, but do not attempt large architectural rewrites as part of
  unrelated changes.
- **Project object**: `Macguffin` (`perframe/Macguffin.m`) is what a
  `.jab` serializes to; `saveJabFile` / `openJabFile` on `JLabelData`
  read and write it.
- **Feature computation**: per-frame features live in
  `perframe/compute_perframe_features/` (and
  `perframe/larva_compute_perframe_features/`); window features and the
  boosting classifier are in `perframe/`.

## Key directories

- `perframe/`: the main JAABA application, model, GUI, features, and
  classifier.  Almost everything of interest lives here.
- `filehandling/`: movie readers (`get_readframe_fcn` and friends),
  trx/trk I/O, and container-format parsing.
- `misc/`: general-purpose helper functions used throughout.
- `spaceTime/`: space-time (optical-flow-based) feature code and
  third-party toolboxes.
- `tests/`: the test suite (see below).
- `compiled/`, `docs/`, `demo/`, `dev/`, `figurecode/`: build
  artifacts, documentation, and assorted scripts.

## Important files

- `perframe/StartJAABA.m`: main GUI entry point.
- `perframe/SetUpJAABAPath.m`: sets up the MATLAB path for JAABA
  (adds `perframe`, `misc`, `filehandling`, `compute_perframe_features`,
  `params`, `tests`, `spaceTime`, ...).  `StartJAABA` runs this first.
- `perframe/JLabelData.m`: the model.
- `perframe/JLabel.m`: the GUI.
- `perframe/Macguffin.m`: the `.jab` project object.
- `filehandling/get_readframe_fcn.m`: the movie-reading entry point.
- `tests/testAll.m`: the test runner.

## Setup and launch

```matlab
SetUpJAABAPath;   % configure the MATLAB path for JAABA only
StartJAABA;       % launch the GUI
```

`SetUpJAABAPath` is a script in `perframe/`; run it from that
directory (or ensure `perframe/` is on the path first).  It sets up
*only* JAABA -- it does not add APT or other sibling projects to the
path.

## Testing

### Structure

- Tests live in `tests/`.  Each test is a MATLAB function
  `function success = testFoo()` that returns a logical `success`
  (true on pass) and/or throws an error on failure.
- `tests/testAll.m` is the runner: it globs every `.m` file in
  `tests/` (except itself and hidden/backup files), calls each by its
  base name, treats a thrown error or a false return as a failure, and
  prints a summary.  Because it globs by filename, a test's **filename
  must equal its function name**.
- Tests do not need to print their own "PASSED" message; reporting is
  `testAll`'s job.

### Test data

Test fixtures live on the Janelia filesystem under
`/groups/branson/bransonlab/projects/JAABA/test_data/` (a few tests
also reference `/groups/branson/home/kabram/...`).  These absolute
paths are hard-coded in the tests, so the suite requires that share to
be mounted.  Movie-reader test movies, including the MP4 fixtures, are
under `test_data/mp4/`.

### Running tests

```matlab
SetUpJAABAPath; testAll()          % run the whole suite
SetUpJAABAPath; testFoo()          % run a single test for debugging
```

Run headless (so GUI tests do not pop windows or steal focus) with:

```
xvfb-run -a matlab -batch "SetUpJAABAPath; testAll()"
```

Do this from the `perframe/` directory so `SetUpJAABAPath` is found.
Note that MATLAB startup plus a parpool spin-up can take a couple of
minutes before any test runs.

## Coding conventions

- Indent with two spaces; never use tabs.  Top-level `function`s are
  not indented.
- Match the style of the surrounding code.  This is a large, diverse 
  codebase; prefer minimal, local changes over sweeping restyling.
- New tests must follow the `function success = testFoo()` convention
  above.

## Git conventions

- Use the 50/72 rule for commit messages; give non-trivial commits a
  body explaining the why.
- Do not add a `Co-Authored-By: Claude` line to commit messages.

## Miscellaneous notes

- `filehandling/get_readframe_fcn` returns
  `[readframe, nframes, fid, headerinfo]`.  For `.mp4/.mov/.m4v` the
  frame count comes from container metadata (`mp4_read_nframes`) and
  random access uses a fast keyframe seek (`mp4_read_index` /
  `mp4_seek_read_frame`) rather than decoding from the start.
- Frame reading must be order-independent: reading a frame out of order
  must return the same pixels as reading it in order (the
  `testUFMFReading` / `testAVIReading` / `testMOVReading` tests guard
  this).
