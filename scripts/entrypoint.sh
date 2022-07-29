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

  for DB in Archaea-RS Archaea-GB Bacteria-RS Bacteria-GB; do
    echo "downloading: https://portal.nersc.gov/project/kbase/GTDB/r207.0/blast_dbs/v1.0.0/$DB-blast_2.13.0_dbs.tgz"
    curl -o $DB-blast_2.13.0_dbs.tgz https://portal.nersc.gov/project/kbase/GTDB/r207.0/blast_dbs/v1.0.0/$DB-blast_2.13.0_dbs.tgz
    echo "tar xfz $DB-blast_2.13.0_dbs.tgz"
    tar xfz $DB-blast_2.13.0_dbs.tgz
    echo "rm $DB-blast_2.13.0_dbs.tgz"
    rm $DB-blast_2.13.0_dbs.tgz
  done

  # validate complete
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
