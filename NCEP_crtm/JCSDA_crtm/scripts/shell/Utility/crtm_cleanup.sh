#!/bin/bash

# Script to clean up a CRTM working copy, removing the library build products
# as well as the regression builds and results.
#
# $Id$

script_id()
{
  REVISION='$Revision$'
  LAST_CHANGED_DATE='$LastChangedDate$'
  echo
  echo "${SCRIPT_NAME} ${REVISION} ${LAST_CHANGED_DATE}"
  echo " "`date`
  echo " Support email: NCEP.List.EMC.JCSDA_CRTM.Support@noaa.gov"
}

usage()
{
  echo
  echo " Usage: crtm_cleanup.sh [-hx]"
  echo
  echo "   Script to clean up a CRTM working copy, removing the library build"
  echo "   products as well as the regression builds and results."
  echo
  echo " Options:"
  echo "   -h"
  echo "         Print this message and exit"
  echo
  echo "   -x"
  echo "         Turn on execution tracing"
  echo
  echo
}

error_message()
{
  MESSAGE=$1
  echo >&2
  echo "  *********************" >&2
  echo "  ${SCRIPT_NAME}(ERROR): ${MESSAGE}" >&2
  echo "  *********************" >&2
}


info_message()
{
  MESSAGE=$1
  echo "  ${SCRIPT_NAME}(INFORMATION): ${MESSAGE}"
}



########################################################################
#                           MAIN SCRIPT BEGINS                         #
########################################################################

# Setup
SCRIPT_NAME=`basename $0`
# ...Definitions
SUCCESS=0; TRUE=0
FAILURE=1; FALSE=1
REMOVE="rm -fr"


# Parse the command line options
while getopts :hx OPTVAL; do
  # Exit if option argument looks like another option
  case ${OPTARG} in
    -*) break ;;
  esac
  # Parse the valid options
  case ${OPTVAL} in
    h)  usage | more; exit ${SUCCESS} ;;
    x)  script_id; set -x ;;
    \?) OPTVAL=${OPTARG}; break ;;
  esac
done
# ...Remove the options processed
shift $((OPTIND - 1))
# ...Output invalidities based on OPTVAL
case ${OPTVAL} in
  # If OPTVAL contains nothing, then all options
  # have been successfully parsed
  \?) : ;;
  # Invalid option
  ?) usage
     error_message "Invalid option '-${OPTARG}'"
     exit ${FAILURE} ;;
esac


# =================
# Start the cleanup
# =================

CURRENT_DIR=${PWD}

# Determine the CRTM version
CRTM_VERSION=`tr -d "'" < ${CRTM_SOURCE_ROOT}/CRTM_Version.inc`
CRTM_LIB_DIR="crtm_${CRTM_VERSION}"


# Clean up the library
cd ${CRTM_SOURCE_ROOT}/Build
if [ $? -ne 0 ]; then
  error_message "Error changing to CRTM build directory, ${CRTM_SOURCE_ROOT}/Build"
  cd ${CURRENT_DIR}; exit ${FAILURE}
fi
# ...Use makefile if it exists
if [ -f Makefile ]; then
  # Uninstall library
  make uninstall
  if [ $? -ne 0 ]; then
    error_message "Error uninstalling CRTM library"
    cd ${CURRENT_DIR}; exit ${FAILURE}
  fi
  # Clean up build products
  make distclean
  if [ $? -ne 0 ]; then
    error_message "Error cleaning up CRTM library build products"
    cd ${CURRENT_DIR}; exit ${FAILURE}
  fi
# ...Otherwise just delete the directory
elif [ -d ${CRTM_LIB_DIR} ]; then
  # Explicitly delete the library directory
  ${REMOVE} ${CRTM_LIB_DIR}
  if [ $? -ne 0 ]; then
    error_message "Error deleting the CRTM library directory, ${CRTM_LIB_DIR}"
    cd ${CURRENT_DIR}; exit ${FAILURE}
  fi
else
  info_message "No library build makefile or ${CRTM_LIB_DIR} directory. Nothing to cleanup!"
fi


# Clean up the source links
cd ${CRTM_SOURCE_ROOT}
if [ $? -ne 0 ]; then
  error_message "Error changing to CRTM source directory, ${CRTM_SOURCE_ROOT}"
  cd ${CURRENT_DIR}; exit ${FAILURE}
fi
make realclean
if [ $? -ne 0 ]; then
  error_message "Error cleaning up CRTM library source links"
  cd ${CURRENT_DIR}; exit ${FAILURE}
fi


# Clean up the regression build
CRTM_REGTEST_DIR="${CRTM_TEST_ROOT}/Main/regression_tests"
cd ${CRTM_REGTEST_DIR}
if [ $? -ne 0 ]; then
  error_message "Error changing to CRTM regression test directory, ${CRTM_REGTEST_DIR}"
  cd ${CURRENT_DIR}; exit ${FAILURE}
fi

# ...Use makefile if it exists
if [ -f Makefile ]; then
  make realclean
  if [ $? -ne 0 ]; then
    error_message "Error cleaning up CRTM regression test build products"
    cd ${CURRENT_DIR}; exit ${FAILURE}
  fi
else
  info_message "No regression test makefile. Nothing to cleanup!"
fi


# Return to original directory
cd ${CURRENT_DIR}
