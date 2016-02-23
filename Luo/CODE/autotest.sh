#! /bin/bash -u
 
# Watches a Fortran source file (*.[fF]90) and its corresponding
# unit test (*.fun), and runs unit test with Funit if either changes.

host=`hostname`

#Load modules containing funit and ifort on K and cypher
if [[ $host = "K"* ]] ; then
  . /usr/local/pkgs/modules/init/bash
  module load ruby_2.0.0_p247
  module load intel_2013.4.183
elif [[ $host = "cypher"* ]] ; then
  . /usr/local/pkgs/modules/init/bash
  module load ruby_2.0.0_p247
  module load intel_2013.4.183
fi

funit_exec=funit

if [ ! "$(${funit_exec} --version 2>/dev/null)" ]; then
  echo ${funit_exec} rubygem required
  echo \ \ sudo gem install ${funit_exec}
  exit 1
fi

fortran_dirs="-s ../libcore -s ../libdefs -s ../libinit -s ../PHYSICS_MODULES -s ../FuncLib90 -s ../libturb -s ../libddfb -s ../PHYSICS_DEPS -s ../LibF90"
c_dirs="-l ../libbc/enginesim -l ../libddfb"

run_funit(){
  #env FC=ifort \
      #FCFLAGS="-O0 -warn -check -traceback -g -DDATADIR=\'path\'" \
  env FC=gfortran \
      FCFLAGS="-O0 -Wall -fbounds-check -fbacktrace -g -DDATADIR=\'path\'" \
      LDFLAGS="-lc++" \
      FSFLAG=-I \
      CC=gcc \
      CXX=g++ \
      MAKE_OPTS='-j' \
    ${funit_exec} ${fortran_dirs} ${c_dirs} $1
}

trap \
  "${funit_exec} ${c_dirs} --clean; \
   rm -f *.{o,mod,MOD} ../{libcore,libdefs,libinit,libturb,PHYSICS_MODULES,FuncLib90,PHYSICS_DEPS}/*.{o,mod,MOD}; \
   exit" \
  INT

if test $# -ne 1; then
  echo Usage: `basename $0` unit_test_name
  exit 1
else
  test_name=$1
  run_funit $test_name 
  while true; do
    if [ $(find $test_name.{fun,[fF]90} -maxdepth 0 -newer TestRunner.f90 | wc -l) -ne 0 ]
    then
      run_funit $test_name 
    fi
    sleep 1
  done
fi
