#/bin/bash
SCRIPT_DIR_NO_CD=$(dirname $0)
SCRIPT_DIR=$(cd $(dirname $0) && pwd)

echo "SCRIPT_DIR_NO_CD: $SCRIPT_DIR_NO_CD"
echo "SCRIPT_DIR: $SCRIPT_DIR"
echo "PWD: $PWD"
exit 0
