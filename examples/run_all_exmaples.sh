# for each bin in example_bins, run it

# run all examples
run_all() {
    for bin in $(ls $examples_dir); do
        echo "Running $bin"
        $examples_dir/$bin
    done
}

main() {
    # read args
    if [ $# -ne 1 ]; then
        echo "Usage: $0 <examples_dir>"
        exit 1
    fi

    examples_dir=$1

    echo "Running all examples in $examples_dir"
    run_all
}

main $@