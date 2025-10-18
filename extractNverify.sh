#!/bin/bash


# This program extracts datasets downloaded from ncbi or other databases in zip format and verifies them using the provided md5sum.txt


if [ -z "$1" ]; then
  echo "Usage: $0 <genomedataset.zip>"
  exit 1
fi

input_file="$1"
input_dir="$(dirname "$input_file")"
filename_with_ext=$(basename "$input_file")
filename="${filename_with_ext%.*}"
output_dir="$input_dir/$filename"

if unzip "$input_file" -d "$output_dir"; then
  cd "$output_dir" || exit
  ls
  md5sum -c md5sum.txt
else
  echo "Failed to unzip $input_file"
  exit 1
fi
