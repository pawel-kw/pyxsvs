#! /bin/bash
rm -rf *.png
rm -rf *.edf
rm -rf *.p
sed --in-place '/mask\.edf/d' xsvs_input.txt 
