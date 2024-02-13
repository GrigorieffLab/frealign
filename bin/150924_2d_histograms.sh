#!/bin/bash
#
# Alexis Rohou, September 2015
#

occ_threshold="80.0"


rm -f tmp_all.txt

for input_par_fn in atp_166_r{1,2,4,5,7,9,10}.par; do
  output_fig_fn="${input_par_fn%%.par}_angular_distribution.pdf"
  # Get theta & phi out of the parameter file
  awk '$1 !~ /C/{if ($12 > '$occ_threshold') {print $3" "$4}}' $input_par_fn > tmp.txt
  cat tmp.txt >> tmp_all.txt
  # Compute the plot
  150923_2d_histogram_03.py tmp.txt $output_fig_fn
  echo "Output figure: $output_fig_fn"
done

output_fig_fn="atp_166_angular_distribution.pdf"
150923_2d_histogram_03.py tmp_all.txt $output_fig_fn
echo "Output figure: $output_fig_fn"
