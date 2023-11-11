#!/bin/bash

# Helper script to run the bench tool on all models and print the results in share-able format

printf "Usage: ./bench.sh [n_threads] [encoder-only]\n"

if [ -z "$1" ]; then
    n_threads=4
else
    n_threads=$1
fi

encoder_only=0
if [ -z "$2" ]; then
    encoder_only=0
else
    encoder_only=$2
fi

models=(                                               \
      "tiny"   "tiny-q5_0"   "tiny-q5_1"   "tiny-q8_0" \
      "base"   "base-q5_0"   "base-q5_1"   "base-q8_0" \
     "small"  "small-q5_0"  "small-q5_1"  "small-q8_0" \
    "medium" "medium-q5_0" "medium-q5_1" "medium-q8_0" \
     "large"  "large-q5_0"  "large-q5_1"  "large-q8_0" \
)

if [ "$encoder_only" -eq 0 ]; then
    printf "\n"
    printf "Running memcpy benchmark\n"
    printf "\n"

    ./bench -w 1 -t $n_threads 2>&1

    printf "\n"
    printf "Running ggml_mul_mat benchmark with $n_threads threads\n"
    printf "\n"

    ./bench -w 2 -t $n_threads 2>&1

    printf "\n"
    printf "Running benchmark for all models\n"
    printf "This can take a while!\n"
    printf "\n"
fi

printf "| %6s | %6s | %16s | %11s | %3s | %7s | %7s | %7s | %7s |\n" "CPU" "OS" "Config" "Model" "Th" "Enc." "Dec." "PP" "Commit"
printf "| %6s | %6s | %16s | %11s | %3s | %7s | %7s | %7s | %7s |\n" "---" "---" "---" "---" "---" "---" "---" "---" "---"

for model in "${models[@]}"; do
    # actual run
    # store stderr output in a variable in order to parse it later
    output=$(./bench -m ./models/ggml-$model.bin -t $n_threads 2>&1)
    ret=$?

    # parse the output:
    encode_time=$(echo "$output" | grep "encode time" | awk '{print $11}')
    decode_time=$(echo "$output" | grep "decode time" | awk '{print $11}')
    prompt_time=$(echo "$output" | grep "prompt time" | awk '{print $11}')
    system_info=$(echo "$output" | grep "system_info")
    n_threads=$(echo "$output" | grep "system_info" | awk '{print $4}')

    # floor to milliseconds
    #encode_time=${encode_time%.*}
    #decode_time=${decode_time%.*}
    #prompt_time=${prompt_time%.*}

    config=""

    if [[ $system_info == *"AVX2 = 1"* ]]; then
        config="$config AVX2"
    fi

    if [[ $system_info == *"NEON = 1"* ]]; then
        config="$config NEON"
    fi

    if [[ $system_info == *"BLAS = 1"* ]]; then
        config="$config BLAS"
    fi

    if [[ $system_info == *"COREML = 1"* ]]; then
        config="$config COREML"
    fi

    if [[ $system_info == *"METAL = 1"* ]]; then
        config="$config METAL"
    fi

    commit=$(git rev-parse --short HEAD)

    if [ $ret -eq 0 ]; then
        printf "| <todo> | <todo> | %16s | %11s | %3s | %7s | %7s | %7s | %7s |\n" "$config" "$model" "$n_threads" "$encode_time" "$decode_time" "$prompt_time" "$commit"
    fi
done
