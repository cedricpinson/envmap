dir=$1
for file in $(ls ${dir}/sample*.data)
do
    numSample=$(head -1 "${file}" | cut -d' ' -f3)
    title=$(head -1 "${file}")
    gnuplot -e "filename='${file}';output_file='${file}.png';numSample='${numSample}';titleGraph='${title}" distrib.gplot
done
