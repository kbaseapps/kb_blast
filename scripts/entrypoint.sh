#!/bin/bash

Make_BLAST_DB=/kb/module/blast/bin/makeblastdb

. /kb/deployment/user-env.sh

python ./scripts/prepare_deploy_cfg.py ./deploy.cfg ./work/config.properties

if [ -f ./work/token ] ; then
  export KB_AUTH_TOKEN=$(<./work/token)
fi

if [ $# -eq 0 ] ; then
  sh ./scripts/start_server.sh
elif [ "${1}" = "test" ] ; then
  echo "Run Tests"
  make test
elif [ "${1}" = "async" ] ; then
  sh ./scripts/run_async.sh
elif [ "${1}" = "init" ] ; then
  echo "Initialize module"

  mkdir -p /data/blast_dbs
  cd /data/blast_dbs

  DB=Archaea-RS
  echo "downloading: https://portal.nersc.gov/project/kbase/GTDB/r207.0/blast_dbs/v1.0.0/$DB.faa.gz"
  curl -o $DB.faa.gz https://portal.nersc.gov/project/kbase/GTDB/r207.0/blast_dbs/v1.0.0/$DB.faa.gz
  gunzip $DB.faa.gz
  $Make_BLAST_DB -in $DB.faa -parse_seqids -dbtype prot -out $DB.faa
  gzip $DB.faa
  
  DB=Archaea-GB
  echo "downloading: https://portal.nersc.gov/project/kbase/GTDB/r207.0/blast_dbs/v1.0.0/$DB.faa.gz"
  curl -o $DB.faa.gz https://portal.nersc.gov/project/kbase/GTDB/r207.0/blast_dbs/v1.0.0/$DB.faa.gz
  gunzip $DB.faa.gz
  $Make_BLAST_DB -in $DB.faa -parse_seqids -dbtype prot -out $DB.faa
  gzip $DB.faa

  DB=Bacteria-RS
  echo "downloading: https://portal.nersc.gov/project/kbase/GTDB/r207.0/blast_dbs/v1.0.0/$DB.faa.gz"
  curl -o $DB.faa.gz https://portal.nersc.gov/project/kbase/GTDB/r207.0/blast_dbs/v1.0.0/$DB.faa.gz
  gunzip $DB.faa.gz
  $Make_BLAST_DB -in $DB.faa -parse_seqids -dbtype prot -out $DB.faa
  gzip $DB.faa

  DB=Bacteria-GB
  echo "downloading: https://portal.nersc.gov/project/kbase/GTDB/r207.0/blast_dbs/v1.0.0/$DB.faa.gz"
  curl -o $DB.faa.gz https://portal.nersc.gov/project/kbase/GTDB/r207.0/blast_dbs/v1.0.0/$DB.faa.gz
  gunzip $DB.faa.gz
  $Make_BLAST_DB -in $DB.faa -parse_seqids -dbtype prot -out $DB.faa
  gzip $DB.faa

  if [ -s "/data/blast_dbs/Archaea-RS.psq" -a -s "/data/blast_dbs/Archaea-GB.psq" -a -s "/data/blast_dbs/Bacteria-RS.10.psq" -a -s "/data/blast_dbs/Bacteria-GB.10.psq" ]; then
    echo "DATA DOWNLOADED SUCCESSFULLY"
    touch /data/__READY__
  else
    echo "Init failed"
  fi
elif [ "${1}" = "bash" ] ; then
  bash
elif [ "${1}" = "report" ] ; then
  export KB_SDK_COMPILE_REPORT_FILE=./work/compile_report.json
  make compile
else
  echo Unknown
fi
