for((i=$1;i<$2;i++))
do
    echo "Clearing model_output..."
    rm -r model_output/*
    mr=` printf "%02d/" $i `
    echo "Copying data/model_run${mr} into model_output..."
    cp -r data/model_run${mr}* model_output/
    echo "Running PARAMETERS.sh"
    cd scripts/
    ./PARAMETERS.sh
    cd ..
    echo "Running model..."
    ./run.sh 0 1000
    echo "Copying model_output into data/model_run${mr} ..."
    rm -r data/model_run${mr}*
    cp -r model_output/* data/model_run${mr}
    echo "Running ./COLLECT.sh for ../model_run${mr} ..."
    cd scripts/
    ./COLLECT.sh ../data/model_run${mr} ../data/model_run${mr}
    cd ../
    echo "Copying /data/model_run${mr}I* into scripts/calculations/balancefunction/model_run${mr}/"
    cp -v data/model_run${mr}I* scripts/calculations/balancefunction/model_run${mr}
done