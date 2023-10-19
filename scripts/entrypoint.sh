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

  export GTDB_VER_INT=214
  export GTDB_VER_FLT=214.1
  export DATA_VER=1.1.0
  export BLAST_VER=2.14.1
  export DL_URL=https://portal.nersc.gov/project/kbase/GTDB/r${GTDB_VER_FLT}/blast_dbs/v${DATA_VER}
  for DB in Archaea-RS Archaea-GB Bacteria-RS Bacteria-GB; do
      export TARGET_FILE=GTDB_r${GTDB_VER_INT}-sprep-${DB}-blast_${BLAST_VER}_dbs.tgz
      export FULL_URL=${DL_URL}/${TARGET_FILE}
      echo "downloading: ${FULL_URL}"
    curl -o ${TARGET_FILE} ${FULL_URL}
    echo "tar xfz ${TARGET_FILE}"
    tar xfz ${TARGET_FILE}
    echo "rm ${TARGET_FILE}"
    rm ${TARGET_FILE}
  done

  # validate complete
  if [ -s "/data/blast_dbs/Archaea-RS.psq" -a -s "/data/blast_dbs/Archaea-GB.psq" -a -s "/data/blast_dbs/Bacteria-RS.10.psq" -a -s "/data/blast_dbs/Bacteria-GB.15.psq" ]; then
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
