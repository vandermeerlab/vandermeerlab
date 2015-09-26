## HOWTO_ExpKeys_Metadata

### Purpose

This document describes how to construct **ExpKeys** and **Metadata**
files, essential components of
[good data and analysis project management](http://ctnsrv.uwaterloo.ca/vandermeerlab/doku.php?id=analysis:nsb2015:week0).

### Background

Freshly recorded raw data generally needs to go through a number of
preprocessing and annotation steps before it can be analyzed and
shared. These steps should be such that a competent colleague will be
able to use your data.

When preprocessing and annotation are complete for a given data set
(typically a single recording session) that session can be
**promoted**. [This](http://ctnsrv.uwaterloo.ca/vandermeerlab/doku.php?id=analysis:nsb2015:week2#data_files_overview)
page describes the files that make up a promoted session. Two of these
are **ExpKeys** and **Metadata**.

Both files contain descriptive information about a single data set,
including basic properties such as the experimenter, subject, and
task (ExpKeys) as well as more session- and experiment-specific data
such as trial intervals, event times, and linearized maze paths. 

The rationale for including this latter information as part of a
promoted data set (rather than as, say, intermediate analysis files
that are not part of the data) is that this information (1) changes
infrequently, and (2) tends to be used in many different analyses.

### ExpKeys

#### Overview

The **ExpKeys** file (named `Rxxx_YYYY_MM_DD_keys.m`, note
underscores) is a MATLAB script that defines a [struct](http://www.mathworks.com/help/matlab/examples/create-a-structure-array.html), with fields
that can be divided into _three_ sections:

1. **Required** fields (all ExpKeys files must have these)
2. **Optional** fields (some ExpKeys files have these; when they do,
   the field names are standardized)
3. **Wildcard** fields (made up by the experimenter; if used, must be
   commented in-line, and explained in the
   [experiment description](http://ctnsrv.uwaterloo.ca/vandermeerlab/doku.php?id=analysis:dataanalysis#task_descriptions_and_metadata))

Because `*keys.m` files are
[scripts](http://www.mathworks.com/help/matlab/learn_matlab/scripts.html)
(this forces the use of underscores rather than dashes in the file
name), they are human-readable. This is typically not true for
information that goes in metadata.

#### Loading ExpKeys

To load a `*keys.m` file, use
[LoadExpKeys.m](https://github.com/mvdm/vandermeerlab/blob/master/code-matlab/shared/io/LoadExpKeys.m).

Notes:

* `.notes`: [string] Located at the top of the keys file when present. 
Any comments noting anything specific to this
session that is not clear from the overall experiment description,
such as 'Multiple headstage detachments, excluded from data --
MvdM'. The contents of this field are displayed by default when
ExpKeys are loaded (this is one reason why it is important and not
simply convenient to use `LoadExpKeys`). 

#### Required ExpKeys fields

* `.species`: [string] 'Rat', 'Mouse', 'Chinchilla', etc.
* `.behavior`: [string] identifier for task used, e.g. 'LinearTrack'
  (each task is described
  [here](http://ctnsrv.uwaterloo.ca/vandermeerlab/doku.php?id=analysis:dataanalysis#task_descriptions_and_metadata)
* `.target`: [1 x nTargets cell array] recording targets, e.g. {'dCA1'}
* `.experimenter`: [string] callsign of experimenter, e.g. 'MvdM'


* `.prerecord`: [2 x 1 double] start and end times of prerecord
* `.postrecord`: [2 x 1 double] start and end times of postrecord
* `.task`: [2 x nBlocks double] start and end times of task; note that
  if your task has multiple blocks that you want to describe here
  (e.g. standard and reversal), use the optional `.taskBlocks` field.


#### Optional Expkeys fields

If video data is included in this data set:

* `.VTConvFactor`: [2 x 1 double] conversion factor for obtaining
  centimeters from raw camera pixels (in the `.nvt` or `.tsp` file) in
  the [x,y] dimension

If more than one target:

* `.electrodeTarget`: [1 x nElectrodes integer] mapping from electrode number
  to target label in `.target`, e.g. `[1 1 2 NaN 2]` maps tetrode 1 and
  2 to the first target, and tetrode 3 and 5 to the second target

If more than one task block

* `.taskBlocks`: [1 x nBlocks cell array] labels of task blocks,
  e.g. `{'Standard','Reversal'}`; note the order of these labels
  corresponds to the order of the start and end times in `.task`.

Other common fields:

* `.day`: [integer] day of recording in a sequence
* `.weight`: [integer] weight of subject that day
* `.age`: [integer] age of subject in days postnatal


* `.goodSWR`: [string] cell array of filenames with good SWR quality
* `.goodTheta`: [string] cell array of filenames with good theta
* `.goodGamma`: [string] cell array of filenames with good gamma


* `.tetrodeDepths`: [1 x nElectrodes integer] estimated electrode
  depths in um from surface of cortex

#### Wildcard Expkeys fields

These can be invented at the discretion of the experimenter to include
data not covered by the above fields.

### Metadata

#### Overview

Unlike ExpKeys files, which are `.m` files readable using a text
editor, Metadata files are binary `.mat` files, and can therefore only
be read after a few more steps. Therefore, Metadata files are more
suitable for larger data (e.g. coordinates of linearized
trajectories), or data that you don't need to access outside of a
MATLAB analysis anyway (e.g. event times).

The metadata file should be named `Rxxx-YYYY-MM-DD-metadata.mat`, and
can be loaded with [LoadMetadata.m]().

There are no required metadata fields, but common fields are described
below.

#### Optional metadata fields

* `.Coord`: linearized paths, generated by [LinearizePos.m]().
* `.evt`: event times, loaded by [LoadEvents.m]().

#### Wildcard metadata fields

### ExpKeys vs. Metadata

Should a given descriptor, when not covered by the required or optional
categories above, go into ExpKeys or into Metadata? This is up to the
experimenter, but some considerations are:

- Would it be helpful if the information is _human-readable_ (from a
  text editor) without having to open up MATLAB? If yes, consider
  ExpKeys. An example is the optional `.goodSWR` field: this contains
  filenames which the user might want to view in
  [Neuraview](http://neuralynx.com/research_software/data_analysis_software/)
  without having to open up MATLAB to find out what the good files
  are.
- Information that is impossible or impractical for humans to read
  (e.g. long lists of trial interval times), or annoying to convert to
  text, Metadata is best.
