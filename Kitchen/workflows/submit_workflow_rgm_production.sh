MODEL="rec"
TAGDIR="calib1.0"
TAG="${TAGDIR}_t0"
RUNS="15359"
TARGET="LAr"
#INPUTS="/mss/clas12/rg-m/data/"
INPUTS="/volatile/clas12/rg-m/48Ca/prod1.4/decoded/015359/"
OUTDIR="/volatile/clas12/rg-m/${TARGET}/${TAGDIR}"
CLARA="/group/clas12/packages/clara/5.0.2_7.0.1"
SCHEMA="dst"
YAML="/home/clas12-5/users/tkutz/kitchen/yaml/rgm_data_recon_${SCHEMA}schema_prod.yaml"
#TRAIN="calib"
#clas12-workflow.py --model $MODEL --runGroup rgm --tag $TAG --runs $RUNS --inputs $INPUTS --outDir $OUTDIR --clara $CLARA --reconYaml $YAML --trainYaml $TRAIN --submit
clas12-workflow.py --model $MODEL --runGroup rgm --tag $TAG --threads 32 --reconSize 2 --runs $RUNS --inputs $INPUTS --outDir $OUTDIR --clara $CLARA --reconYaml $YAML --submit
