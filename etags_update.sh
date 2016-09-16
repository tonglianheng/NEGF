#!/bin/bash

# --------------------------------------------------------------------
# Script for recrusively updating TAGS file used by emacs for a source
# tree.
# --------------------------------------------------------------------
# Lianheng Tong (2014), tonglianheng@gmail.com
# --------------------------------------------------------------------

######################################################################
# You should modify the following to suit your source code
######################################################################

ETAGS=etags
SRC_ROOT=.
INCLUDE_REGEXP='.*\.(F|f90|h)'
IGNORE_REGEXP='[^a-zA-Z]*svn[^a-zA-Z]*'

# If $ETAGS_REGEXP_FILE is found, then -r option of etags is used.
ETAGS_REGEXP_FILE=./etags_regexp

######################################################################
# You do not need to modify the code below
######################################################################

if [ "x$ETAGS_REGEXP_FILE" != "x" ] && \
   [ -f "$ETAGS_REGEXP_FILE" ] ; then
    cat >&2 <<EOF
[32m$ETAGS_REGEXP_FILE file found[0m, using -r option of etags
EOF
    etags_command="$ETAGS -r @etags_regexp"
else
    etags_command="$ETAGS"
fi

# Linux find and Mac OS X find has different options, so need to call
# them differently
if (man find | grep "\-E" >& /dev/null) ; then
    find -L -E . -regex $INCLUDE_REGEXP -and \
         -not -regex $IGNORE_REGEXP \
         -print | xargs $etags_command
else
    find -L . -regextype posix-extended -regex $INCLUDE_REGEXP -and \
         -not -regex $IGNORE_REGEXP \
         -print | xargs $etags_command
fi

echo "done"
