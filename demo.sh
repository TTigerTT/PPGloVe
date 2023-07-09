

#!/bin/bash
#set -e

# Makes programs, downloads sample data, trains a GloVe model, and then evaluates it.
# One optional argument can specify the language used for eval script: matlab, octave or [default] python

make
#if [ ! -e text8 ]; then
 # if hash wget 2>/dev/null; then
 #   wget http://mattmahoney.net/dc/text8.zip
 # else
   # curl -O http://mattmahoney.net/dc/text8.zip
 # fi
 # unzip text8.zip
  #rm text8.zip
#fi

CORPUS=enwiki-10.txt
VOCAB_FILE=mergedvocab_20mb.txt
COOCCURRENCE_FILE=pmerge_20mb.bin
COOCCURRENCE_SHUF_FILE=pmergeshuff_20mb.bin

CIPCOOCCURRENCE_FILE=merge_20mb.bin
CIPCOOCCURRENCE_SHUF_FILE=secmergeshuf_20series.bin

BUILDDIR=build
SAVE_FILE=vectors2
VERBOSE=2
MEMORY=4.0
VOCAB_MIN_COUNT=5
VECTOR_SIZE=100
MAX_ITER=15
WINDOW_SIZE=15
BINARY=2
NUM_THREADS=50
X_MAX=10
if hash python 2>/dev/null; then
    PYTHON=python
else
    PYTHON=python3
fi
# echo "$ gcc generateTriple.c -o generateTriple "
#       gcc generateTriple.c -o generateTriple
# echo "$ ./ generateTriple "
#       ./generateTriple

# echo "$ $PYTHON splitdata.py"
#        $PYTHON splitdata.py

# for i in {0..9}; do
#     infile="user_${i}_20mb.txt"
#     outfile="user${i}_vocab_20mb.txt"
#     echo "$ $BUILDDIR/vocab_count -min-count $VOCAB_MIN_COUNT -verbose $VERBOSE < $infile > $outfile"
#     $BUILDDIR/vocab_count -min-count $VOCAB_MIN_COUNT -verbose $VERBOSE < $infile > $outfile
# done


# echo "$ $PYTHON mergevocab.py"
#        $PYTHON mergevocab.py

# for i in {0..9}; do
#     uservocab="user${i}_vocab_20mb.txt"
#     infile="user_${i}_20mb.txt"
#     outfile="pusercoocur${i}_20mb.bin"
#     echo "$ $BUILDDIR/cooccur_persuser -memory $MEMORY -vocab-file mergedvocab_20mb.txt -uservocab $uservocab -verbose $VERBOSE -window-size $WINDOW_SIZE < $infile > $outfile"
#     $BUILDDIR/cooccur_persuser -memory $MEMORY -vocab-file mergedvocab_20mb.txt -uservocab $uservocab -verbose $VERBOSE -window-size $WINDOW_SIZE < $infile > $outfile
# done

# for i in {0..9}; do
#     uservocab="user${i}_vocab_20mb.txt"
#     infile="user_${i}_20mb.txt"
#     outfile="usercoocur${i}_20mb.bin"
#     echo "$ $BUILDDIR/per_file_corr -memory $MEMORY -vocab-file mergedvocab_20mb.txt -uservocab $uservocab -verbose $VERBOSE -window-size $WINDOW_SIZE < $infile > $outfile"
#     $BUILDDIR/per_file_corr -memory $MEMORY -vocab-file mergedvocab_20mb.txt -uservocab $uservocab -verbose $VERBOSE -window-size $WINDOW_SIZE < $infile > $outfile
# done


# echo "$ $BUILDDIR/realmerge   -verbose $VERBOSE   > merge_20mb.bin"
# $BUILDDIR/realmerge -verbose $VERBOSE  > merge_20mb.bin

# echo "$ $BUILDDIR/mergeuserdata   -verbose $VERBOSE   > pmerge_20mb.bin"
# $BUILDDIR/mergeuserdata -verbose $VERBOSE  > pmerge_20mb.bin



# echo "$ $BUILDDIR/shuffle -memory $MEMORY -verbose $VERBOSE < pmerge_20mb.bin > pmergeshuff_20mb.bin"
# $BUILDDIR/shuffle -memory $MEMORY -verbose $VERBOSE < pmerge_20mb.bin > pmergeshuff_20mb.bin


# echo "$ $BUILDDIR/shuffcipdata -memory $MEMORY -verbose $VERBOSE < merge_20mb.bin > mergeshuff_20mb.bin"
# $BUILDDIR/shuffcipdata -memory $MEMORY -verbose $VERBOSE < merge_20mb.bin > mergeshuff_20mb.bin
# 



echo "$ $BUILDDIR/glove -save-file $SAVE_FILE -threads $NUM_THREADS -input-file $COOCCURRENCE_SHUF_FILE -x-max $X_MAX -iter $MAX_ITER -vector-size $VECTOR_SIZE -binary $BINARY -vocab-file $VOCAB_FILE -verbose $VERBOSE -eta 0.09"
$BUILDDIR/glove -save-file $SAVE_FILE -threads $NUM_THREADS -input-file $COOCCURRENCE_SHUF_FILE -x-max $X_MAX -iter $MAX_ITER -vector-size $VECTOR_SIZE -binary $BINARY -vocab-file $VOCAB_FILE -verbose $VERBOSE -eta 0.09


echo "$ $BUILDDIR/ciplogandf"
$BUILDDIR/ciplogandf 

echo "$ $BUILDDIR/traincip_triple -save-file $SAVE_FILE -threads $NUM_THREADS -input-file $CIPCOOCCURRENCE_SHUF_FILE -x-max $X_MAX -iter $MAX_ITER -vector-size $VECTOR_SIZE -binary $BINARY -vocab-file $VOCAB_FILE -verbose $VERBOSE -eta 0.13"

$BUILDDIR/traincip_triple -save-file $SAVE_FILE -threads $NUM_THREADS -input-file $CIPCOOCCURRENCE_SHUF_FILE -x-max $X_MAX -iter $MAX_ITER -vector-size $VECTOR_SIZE -binary $BINARY -vocab-file $VOCAB_FILE -verbose $VERBOSE -eta 0.13 



if [ "$CORPUS" = 'enwiki-10.txt' ]; then
   if [ "$1" = 'matlab' ]; then
       matlab -nodisplay -nodesktop -nojvm -nosplash < ./eval/matlab/read_and_evaluate.m 1>&2 
   elif [ "$1" = 'octave' ]; then
       octave < ./eval/octave/read_and_evaluate_octave.m 1>&2
   else
       echo "$ $PYTHON eval/python/evaluate.py"
       $PYTHON eval/python/evaluate.py
   fi
fi

