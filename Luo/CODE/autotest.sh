#! /bin/bash -u
 
# Watches a Fortran source file (*.[fF]90) and its corresponding
# unit test (*.fun), and runs unit test with Funit if either changes.

export CC=gcc \
       CXX=g++

funit_exec=funit

if [ ! "$(${funit_exec} --version 2>/dev/null)" ]; then
  echo ${funit_exec} rubygem required
  echo \ \ sudo gem install ${funit_exec}
  exit 1
fi

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
      ${funit_exec} -s ./unit_tests $1
}

trap \
  "${funit_exec} --clean; \
   rm -f *.{o,mod,MOD} 
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
