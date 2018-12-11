
#!/bin/bash

set -e

odtmroot=$(pwd)

debug=""
npes=1
while getopts 'dj:' flag; do
    case "${flag}" in
    d) debug=".debug" ;;
    j) npes=$OPTARG ;;
    esac
done

shift $(($OPTIND - 1))

opts=$@

EXE="odtm.exe"

execdir="$odtmroot/exec"
mkmf="$odtmroot/bin/mkmf"

mkmftemplate="$odtmroot/bin/mkmf.template$debug"

FMS_UTILS=$odtmroot/src/fms_shared

FMS_UTILITIES="$FMS_UTILS/include \
			   $FMS_UTILS/platform \
               $FMS_UTILS/constants \
			   $FMS_UTILS/fms \
			   $FMS_UTILS/time_manager \
				$FMS_UTILS/mpp \
				$FMS_UTILS/diag_manager  \
				$FMS_UTILS/memutils \
				$FMS_UTILS/constants \
				$FMS_UTILS/mpp/include \
				$FMS_UTILS/data_override \
				$FMS_UTILS/horiz_interp \
				$FMS_UTILS/time_interp \
				$FMS_UTILS/axis_utils \
				$FMS_UTILS/mosaic"

paths="$odtmroot/src/odtm $odtmroot/src/nemuro"



echo '...............Compiling mppnccombine.....................'

cd $odtmroot/src/postproc/mppnccombine
make

echo '...............Done Compiling mppnccombine.....................'



echo '...............Compiling datefunc.....................'

cd $odtmroot/src/postproc/datefunc
make

echo '...............Done Compiling datefunc.....................'


mkdir -p $execdir/lib_fms

echo '...............Compiling lib_fms.....................'
cd $execdir/lib_fms

$mkmf -f -p lib_fms.a -t $mkmftemplate $FMS_UTILITIES

make -j 16
echo '...............Done compiling lib_fms.....................'


echo '...............Compiling ODTM.....................'

mkdir -p $execdir/odtm

cd $execdir/odtm

$mkmf -f -p $EXE -t $mkmftemplate -o "-I$execdir/lib_fms" -l "$execdir/lib_fms/lib_fms.a" $paths

make -j $npes $opts

echo '...............Done Compiling ODTM.....................'


