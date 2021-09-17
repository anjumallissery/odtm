
#!/bin/bash

echo '...............Setting up environment.....................'


if [ ! -f .env ]; then
	echo ".env file does not exist. Run init.sh first."
	exit
fi

. .env

. $rootdir/bin/env.$MACH

set -e

debug=""
npes=1
workdir='none'
while getopts 'dj:w:' flag; do
    case "${flag}" in
    d) debug=".debug" ;;
    j) npes=$OPTARG ;;
    w) workdir=$OPTARG ;;
    esac
done

shift $(($OPTIND - 1))

opts=$@

EXE="odtm.exe"

execdir="$rootdir/exec"
mkmf="$rootdir/bin/mkmf"

mkmftemplate="$rootdir/bin/mkmf.template$debug"

FMS_UTILS=$rootdir/src/fms_shared

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

paths="$rootdir/src/odtm"

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


echo "#-------------------------MAKE RUN_NCCOMBINEP2R--------------------------------------"
cppDef="-Dlib_mppnccp2r -Duse_libMPI"
exe=run_mppnccp2r
paths="$rootdir/src/postproc/mppnccombinep2r"
export LD=$FC
mkdir -p $execdir/$exe
cd $execdir/$exe

OPTS="-I$execdir/lib_fms"

LIBS="$execdir/lib_fms/lib_fms.a"

$mkmf -c "$cppDef" -f -p ${exe} -t $mkmftemplate -o "$OPTS" -l "$LIBS"  $paths
make -j $npes
echo "#--------------------------------------------------------------------------------"

filestocopy="data_table diag_table input.nml run_mppnccombine.sh odtm_submit.pbs"

if [ "$workdir" != "none" ]; then
	wrkdir="$rootdir/work/$workdir"
	if [ -d "$wrkdir" ]; then
		echo "Work directory $wrkdir already exist!!"
		exit
	fi
	mkdir -p $wrkdir
  mkdir -p $wrkdir/INPUT 
  mkdir -p $wrkdir/RESTART
  mkdir -p $wrkdir/OUTPUT

	cd $wrkdir/INPUT 
	wget https://www.dropbox.com/s/zwg6839nq0d5sxb/input.tar.gz
	tar -zxvf input.tar.gz
	cp -f input.nml $wrkdir/
	cp -f data_table $wrkdir/
	cp -f diag_table $wrkdir/

	for f in $filestocopy; do
		cp $rootdir/scripts/$f $wrkdir/
		sed -i "s|_ROOTDIR_|$rootdir|g" $wrkdir/$f
		sed -i "s|_EXPNAME_|$workdir|g" $wrkdir/$f
  done
	echo 
	echo "Experiment directory is created: $wrkdir"
	echo
fi
